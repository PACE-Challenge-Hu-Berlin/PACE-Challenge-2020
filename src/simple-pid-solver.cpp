#include <cstring>
#include <iostream>

#include "connected-components.hpp"
#include "simple-pid-solver.hpp"

namespace {
	constexpr bool debug_trees = false;

	// Helper class to print object during debugging.
	struct print_composition {
		friend std::ostream &operator<< (std::ostream &os, const print_composition &self) {
			os << "({";
			for(size_t i = 0; i < self.p_->vertices.size(); i++) {
				if(i > 0)
					os << ", ";
				os << self.p_->vertices[i];
			}
			os << "}, h=" << self.p_->h;
			if(self.p_->trivial)
				os << ", trivial";
			os << ")";
			return os;
		}

		print_composition(const simple_pid_solver::feasible_composition &p)
		: p_{&p} { }

	private:
		const simple_pid_solver::feasible_composition *p_;
	};

	struct print_memory {
		friend std::ostream &operator<< (std::ostream &os, const print_memory &self) {
			if(self.n_ > 10 * (size_t(1) << 30)) {
				os << self.n_ << " GiB";
			}else if(self.n_ >= 10 * (size_t(1) << 20)) {
				os << (self.n_ >> 20) << " MiB";
			}else if(self.n_ > 10 * (size_t(1) << 10)) {
				os << (self.n_ >> 10) << " KiB";
			}else{
				os << self.n_ << " B";
			}
			return os;
		}

		print_memory(size_t n)
		: n_{n} { }

	private:
		size_t n_;
	};

	template<typename C>
	vertex_span copy_to_arena(const C &vec, memory_arena &memory) {
		arena_allocator<vertex> alloc{memory};
		auto p = alloc.allocate(vec.size());
		memcpy(p, vec.data(), sizeof(vertex) * vec.size());
		return vertex_span{p, p + vec.size()};
	}

	template<typename C>
	vertex_span copy_to_queue(const C &vec, queue_memory &memory) {
		queue_allocator<vertex> alloc{memory};
		auto p = alloc.allocate(vec.size());
		memcpy(p, vec.data(), sizeof(vertex) * vec.size());
		return vertex_span{p, p + vec.size()};
	}

	void free_in_arena(vertex_span vs, memory_arena &memory) {
		arena_allocator<vertex> alloc{memory};
		alloc.deallocate(const_cast<vertex *>(vs.data()), vs.size());
	}

	void free_in_queue(vertex_span vs, queue_memory &memory) {
		queue_allocator<vertex> alloc{memory};
		alloc.deallocate(const_cast<vertex *>(vs.data()), vs.size());
	}
}

simple_pid_solver::simple_pid_solver(graph &g)
: g_{&g} { }

int simple_pid_solver::compute_treedepth() {
	connected_components cc;
	cc.compute(*g_);

	// Eliminate all but the first component.
	for(size_t j = 1; j < cc.num_components(); j++)
		for(vertex v : cc.component(j))
			g_->eliminate(v);

	int td = 0;

	// Invariant: Only the i-th component is uneliminated.
	size_t i;
	for(i = 0; /* exit through break below */ ; i++) {
		int k = 0;
		std::cerr << "solving CC " << (i + 1) << " of " << cc.num_components() << std::endl;
		while(true) {
			if(decide_treedepth_(k))
				break;
			k++;
		}
		td = std::max(td, k);

		// Eliminate the i-th component, uneliminate the (i+1)-th.
		if(i + 1 >= cc.num_components())
			break;
		for(vertex v : cc.component(i))
			g_->eliminate(v);
		for(vertex v : cc.component(i + 1))
			g_->uneliminate(v);
	}

	// Uneliminate all but the i-th component.
	for(size_t j = 0; j < cc.num_components(); j++) {
		if(j == i)
			continue;
		for(vertex v : cc.component(j))
			g_->uneliminate(v);
	}

	return td;
}

bool simple_pid_solver::decide_treedepth_(int k) {
	for(const feasible_tree &tree : active_trees)
		free_in_arena(tree.vertices, eternal_arena_);
	for(const feasible_tree &tree : inactive_trees)
		free_in_arena(tree.vertices, eternal_arena_);
	active_trees.clear();
	inactive_trees.clear();
	staged_trees.clear();
	eternal_arena_.reset();

	for(vertex v : g_->vertices()) {
		if(g_->degree(v) + 1 > k)
			continue;
		workset_.clear();
		workset_.push_back(v);
		staged_tree staged{copy_to_arena(workset_, eternal_arena_), 1, true};
		staged_trees.emplace(staged.vertices, staged);
	}

	for(int h = 1; h < k; h++) {
		if(staged_trees.empty())
			break;

		std::cerr << "constructing k = " << k << ", h = " << h << std::endl;

		// Move trees for the current h from the staging buffer to the set of feasible trees.
		auto sit = staged_trees.begin();
		while(sit != staged_trees.end()) {
			const auto &staged = sit->second;
			if(staged.h != h) {
				++sit;
				continue;
			}

			// Find the neighbor set of the tree to turn it into an atomic forest.
			pivot_marker_.reset(g_->id_limit());
			pivot_neighbor_marker_.reset(g_->id_limit());
			separator_.clear();
			for(vertex v : staged.vertices)
				pivot_marker_.mark(v);
			for(vertex v : staged.vertices)
				for(vertex w : g_->neighbors(v)) {
					if(pivot_marker_.is_marked(w))
						continue;
					if(pivot_neighbor_marker_.is_marked(w))
						continue;
					pivot_neighbor_marker_.mark(w);
					separator_.push_back(w);
				}

			auto rv = staged.vertices.front();
			join_q_.push(feasible_forest{rv, rv,
					copy_to_queue(staged.vertices, join_memory_),
					copy_to_queue(separator_, join_memory_),
					true, staged.trivial});
			join_memory_.seal();
			active_trees.insert(staged.vertices, feasible_tree{staged.vertices, staged.h});
			sit = staged_trees.erase(sit);
		}

		num_join_ = 0;
		num_compose_ = 0;
		num_stage_ = 0;
		num_unimproved_ = 0;
		std::cerr << "    k = " << k << ", h = " << h << ": there are "
				<< inactive_trees.size() << " inactive trees" << std::endl;
		std::cerr << "    k = " << k << ", h = " << h << ": there are "
				<< join_q_.size() << " initial forests" << std::endl;

		while(!join_q_.empty()) {
			feasible_forest forest = join_q_.front();
			join_q_.pop();
			join_(k, h, forest);
			num_join_++;

			while(!compose_q_.empty()) {
				feasible_composition comp = compose_q_.front();
				compose_q_.pop();
				compose_(k, h, comp);
				num_compose_++;

				free_in_queue(comp.vertices, compose_memory_);
				free_in_queue(comp.prefix, compose_memory_);
				compose_memory_.reclaim();
			}

			free_in_queue(forest.vertices, join_memory_);
			free_in_queue(forest.separator, join_memory_);
			join_memory_.reclaim();
		}

		std::cerr << "    staged " << num_stage_
				<< " trees, discarded " << num_unimproved_
				<< " (" << num_join_ << " joins, " << num_compose_ << " compositions)"
				<< std::endl;
		std::cerr << "    empty separator: " << stats_.num_empty_separator << std::endl;
		std::cerr << "    pruned forests: " << stats_.num_pruned_forests
				<< ", compositions: " << stats_.num_pruned_compositions << std::endl;
		std::cerr << "    eternal memory: " << print_memory{eternal_arena_.used_space()}
				<< ", allocations: " << eternal_arena_.num_allocations()
				<< std::endl;
		std::cerr << "    join memory: " << print_memory{join_memory_.max_used_space()}
				<< std::endl;
		std::cerr << "    compose memory: " << print_memory{compose_memory_.max_used_space()}
				<< std::endl;

		for(const feasible_tree &tree : active_trees)
			inactive_trees.insert(tree.vertices, tree);
		active_trees.clear();
	}

	// Move trees for the final h from the staging buffer to the set of feasible trees.
	auto sit = staged_trees.begin();
	while(sit != staged_trees.end()) {
		const auto &staged = sit->second;
		if(staged.h != k) {
			++sit;
			continue;
		}

		active_trees.insert(staged.vertices, feasible_tree{staged.vertices, staged.h});
		sit = staged_trees.erase(sit);
	}

	for(feasible_tree &tree : active_trees) {
		if(tree.vertices.size() == g_->num_vertices())
			return true;
	}
	for(feasible_tree &tree : inactive_trees) {
		if(tree.vertices.size() == g_->num_vertices())
			return true;
	}
	return false;
}

void simple_pid_solver::join_(int k, int h, feasible_forest &forest) {
	// Mark and count the neighbors of the feasible forest.
	pivot_marker_.reset(g_->id_limit());
	pivot_neighbor_marker_.reset(g_->id_limit());
	unsigned int num_forest_neighbors = 0;
	for(vertex v : forest.vertices)
		pivot_marker_.mark(v);
	for(vertex v : forest.vertices) {
		for(vertex w : g_->neighbors(v)) {
			if(pivot_marker_.is_marked(w))
				continue;
			if(pivot_neighbor_marker_.is_marked(w))
				continue;
			pivot_neighbor_marker_.mark(w);
			++num_forest_neighbors;
		}
	}

	auto join_with = [&] (const feasible_tree &tree, vertex min_rv, vertex max_rv) {
		int num_total_neighbors = num_forest_neighbors;

		// Determine the neighbors of the tree that we are trying to join.
		associate_marker_.reset(g_->id_limit());
		associate_neighbor_marker_.reset(g_->id_limit());
		for(vertex v : tree.vertices) {
			assert(!pivot_marker_.is_marked(v));
			associate_marker_.mark(v);
		}

		for(vertex v : tree.vertices) {
			for(vertex w : g_->neighbors(v)) {
				assert(!pivot_marker_.is_marked(w));
				if(associate_marker_.is_marked(w))
					continue;
				if(associate_neighbor_marker_.is_marked(w))
					continue;
				associate_neighbor_marker_.mark(w);
				if(!pivot_neighbor_marker_.is_marked(w))
					++num_total_neighbors;
			}
		}

		if(h + num_total_neighbors > k) {
			++stats_.num_pruned_forests;
			return;
		}

		separator_.clear();
		for(vertex v : forest.separator)
			if(associate_neighbor_marker_.is_marked(v))
				separator_.push_back(v);

		if(!separator_.size()) {
			++stats_.num_empty_separator;
			return;
		}

		workset_.clear();
		workset_.insert(workset_.end(), forest.vertices.begin(), forest.vertices.end());
		workset_.insert(workset_.end(), tree.vertices.begin(), tree.vertices.end());
		join_q_.push(feasible_forest{min_rv, max_rv,
				copy_to_queue(workset_, join_memory_),
				copy_to_queue(separator_, join_memory_),
				false, false});
		join_memory_.seal();
	};

	auto non_disjoint_predicate = [&] (vertex v) {
		return pivot_marker_.is_marked(v) || pivot_neighbor_marker_.is_marked(v);
	};

	// To avoid duplicate joins, we only join trees into this forest that either
	//    - are ordered before all components of this forest,
	// or - are ordered after all components of this forest
	// w.r.t. the order induced by representative vertices.

	// Join only lesser active trees into this forest.
	for(feasible_tree &tree : active_trees.disjoint_less(sieve_query_,
				forest.min_rv, non_disjoint_predicate)) {
		auto rv = tree.vertices.front();
		join_with(tree, rv, forest.max_rv);
	}
	// Join both lesser and greater inactive trees into this forest.
	for(feasible_tree &tree : inactive_trees.disjoint_less(sieve_query_,
				forest.min_rv, non_disjoint_predicate)) {
		auto rv = tree.vertices.front();
		join_with(tree, rv, forest.max_rv);
	}
	for(feasible_tree &tree : inactive_trees.disjoint_greater(sieve_query_,
				forest.max_rv, non_disjoint_predicate)) {
		auto rv = tree.vertices.front();
		join_with(tree, forest.min_rv, rv);
	}

	if(!forest.atomic || forest.trivial) {
		compose_q_.push(feasible_composition{
				copy_to_queue(forest.vertices, compose_memory_),
				copy_to_queue(forest.separator, compose_memory_),
				h, forest.trivial});
		compose_memory_.seal();
	}
}

void simple_pid_solver::compose_(int k, int h, feasible_composition &comp) {
	pivot_marker_.reset(g_->id_limit());
	pivot_neighbor_marker_.reset(g_->id_limit());
	int num_forest_neighbors = 0;
	for(vertex v : comp.vertices)
		pivot_marker_.mark(v);
	for(vertex v : comp.vertices) {
		for(vertex w : g_->neighbors(v)) {
			if(pivot_marker_.is_marked(w))
				continue;
			if(pivot_neighbor_marker_.is_marked(w))
				continue;
			pivot_neighbor_marker_.mark(w);
			++num_forest_neighbors;
		}
	}

	for(vertex u : comp.prefix) {
		// Count number of neighbors of the composition.
		// Since u is a neighbor of the vertex set, subtract one to begin with.
		int num_composed_neighbors = num_forest_neighbors - 1;
		for(vertex w : g_->neighbors(u)) {
			if(pivot_marker_.is_marked(w) || pivot_neighbor_marker_.is_marked(w))
				continue;
			++num_composed_neighbors;
		}

		if(comp.h + 1 + num_composed_neighbors > k) {
			++stats_.num_pruned_compositions;
			continue;
		}

		separator_.clear();
		for(vertex v : comp.prefix)
			if(v > u)
				separator_.push_back(v);

		workset_.clear();
		workset_.insert(workset_.end(), comp.vertices.begin(), comp.vertices.end());
		workset_.push_back(u);
		std::sort(workset_.begin(), workset_.end());

		auto improves_staged = [&] (const staged_tree &staged) {
			if(comp.h + 1 < staged.h)
				return true;
			if(comp.h + 1 > staged.h)
				return false;
			return !comp.trivial && staged.trivial;
		};

		auto it = staged_trees.find(vertex_span{workset_});
		if(it == staged_trees.end()) {
			++num_stage_;
			staged_tree staged{copy_to_arena(workset_, eternal_arena_),
					comp.h + 1, comp.trivial};
			staged_trees.emplace(staged.vertices, staged);
		}else if(improves_staged(it->second)) {
			++num_stage_;
			it->second.h = comp.h + 1;
			it->second.trivial = comp.trivial;
		}else{
			++num_unimproved_;
		}

		compose_q_.push(feasible_composition{
				copy_to_queue(workset_, compose_memory_),
				copy_to_queue(separator_, compose_memory_),
				comp.h + 1, comp.trivial});
		compose_memory_.seal();
	}
}
