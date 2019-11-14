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

	template<typename T, typename A>
	std::vector<T> copy_from(const std::vector<T, A> &vec) {
		std::vector<T> res;
		res.reserve(vec.size());
		res.insert(res.end(), vec.begin(), vec.end());
		return res;
	}

	template<typename T>
	std::vector<T, arena_allocator<T>> copy_to_arena(const std::vector<T> &vec,
			memory_arena &arena) {
		std::vector<T, arena_allocator<T>> res{arena_allocator<T>{arena}};
		res.reserve(vec.size());
		res.insert(res.end(), vec.begin(), vec.end());
		return res;
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
	feasible_trees.clear();
	staged_trees.clear();
	eternal_arena_.reset();

	for(vertex v : g_->vertices()) {
		if(g_->degree(v) + 1 > k)
			continue;
		workset_.clear();
		workset_.push_back(v);
		staged_tree staged{copy_to_arena(workset_, eternal_arena_), 1, true};
		vertex_span vs{staged.vertices}; // Careful with the move below.
		staged_trees.emplace(vs, std::move(staged));
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

			size_t idx = feasible_trees.size();
			join_q_.push(feasible_forest{idx,
					copy_from(staged.vertices), // Careful with the move below.
					std::move(separator_), true, staged.trivial});
			feasible_trees.push_back(feasible_tree{std::move(staged.vertices), staged.h});
			sit = staged_trees.erase(sit);
		}

		num_join_ = 0;
		num_compose_ = 0;
		std::cerr << "    k = " << k << ", h = " << h << ": there are "
				<< join_q_.size() << " initial forests" << std::endl;

		while(!join_q_.empty()) {
			feasible_forest forest = std::move(join_q_.front());
			join_q_.pop();
			join_(k, h, forest);
			num_join_++;

			if(!forest.atomic || forest.trivial)
				compose_q_.push({std::move(forest.vertices), std::move(forest.separator),
						h, forest.trivial});
		}

		std::cerr << "    k = " << k << ", h = " << h << ": there are "
				<< compose_q_.size() << " initial compositions" << std::endl;

		while(!compose_q_.empty()) {
			feasible_composition comp = std::move(compose_q_.front());
			compose_q_.pop();
			compose_(k, h, comp);
			num_compose_++;

			if(comp.h > h) {
				auto improves_staged = [&] (const staged_tree &staged) {
					if(comp.h < staged.h)
						return true;
					if(comp.h > staged.h)
						return false;
					return !comp.trivial && staged.trivial;
				};

				std::sort(comp.vertices.begin(), comp.vertices.end());

				auto it = staged_trees.find(vertex_span{comp.vertices});
				if(it == staged_trees.end()) {
					staged_tree staged{copy_to_arena(comp.vertices, eternal_arena_),
							comp.h, comp.trivial};
					vertex_span vs{staged.vertices}; // Careful with the move below.
					staged_trees.emplace(vs, std::move(staged));
				}else if(improves_staged(it->second)) {
					it->second.h = comp.h;
					it->second.trivial = comp.trivial;
				}
			}
		}

		std::cerr << "    k = " << k << ", h = " << h
				<< ": staged " << staged_trees.size() << " trees"
				<< " (" << num_join_ << " joins, " << num_compose_ << " compositions)"
				<< std::endl;
		std::cerr << "    non-disjoint objects: " << stats_.num_non_disjoint
				<< ", non-separated: " << stats_.num_non_separated
				<< ", empty separator: " << stats_.num_empty_separator << std::endl;
		std::cerr << "    eternal memory: " << print_memory{eternal_arena_.used_space()}
				<< ", allocations: " << eternal_arena_.num_allocations()
				<< std::endl;
	}

	// Move trees for the final h from the staging buffer to the set of feasible trees.
	auto sit = staged_trees.begin();
	while(sit != staged_trees.end()) {
		const auto &staged = sit->second;
		if(staged.h != k) {
			++sit;
			continue;
		}
		feasible_trees.push_back(feasible_tree{std::move(staged.vertices), staged.h});
		sit = staged_trees.erase(sit);
	}

	for(const auto &tree : feasible_trees)
		if(tree.vertices.size() == g_->num_vertices())
			return true;
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

	for(size_t i = 0; i < forest.idx; i++) {
		// We only need to join forests with smaller indices into this forest.
		const auto &tree = feasible_trees[i];

		int num_total_neighbors = num_forest_neighbors;

		// Helper function. Return true iff the forest and the tree are disjoint.
		// Marks the vertices of the tree.
		auto validate_tree_vertices = [&] () -> bool {
			associate_marker_.reset(g_->id_limit());
			associate_neighbor_marker_.reset(g_->id_limit());
			for(vertex v : tree.vertices) {
				associate_marker_.mark(v);
				if(pivot_marker_.is_marked(v))
					return false;
			}
			return true;
		};

		// Helper function. Return true iff the forest and the tree are separated.
		// Marks the neighbors of the tree.
		auto validate_tree_neighbors = [&] () -> bool {
			for(vertex v : tree.vertices) {
				for(vertex w : g_->neighbors(v)) {
					if(pivot_marker_.is_marked(w))
						return false;
					if(associate_marker_.is_marked(w))
						continue;
					if(associate_neighbor_marker_.is_marked(w))
						continue;
					associate_neighbor_marker_.mark(w);
					if(!pivot_neighbor_marker_.is_marked(w))
						++num_total_neighbors;
				}
			}
			return true;
		};

		// Determine the neighbors of the tree that we are trying to join.
		// Do not join if the objects are not disjoint or not separated.
		if(!validate_tree_vertices()) {
			++stats_.num_non_disjoint;
			continue;
		}
		if(!validate_tree_neighbors()) {
			++stats_.num_non_separated;
			continue;
		}

		if(h + num_total_neighbors > k)
			continue;

		separator_.clear();
		for(vertex v : forest.separator)
			if(associate_neighbor_marker_.is_marked(v))
				separator_.push_back(v);

		if(!separator_.size()) {
			++stats_.num_empty_separator;
			continue;
		}

		workset_.clear();
		workset_.insert(workset_.end(), forest.vertices.begin(), forest.vertices.end());
		workset_.insert(workset_.end(), tree.vertices.begin(), tree.vertices.end());
		join_q_.push(feasible_forest{forest.idx, std::move(workset_), std::move(separator_),
				false, false});
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

		if(comp.h + 1 + num_composed_neighbors > k)
			continue;

		separator_.clear();
		for(vertex v : comp.prefix)
			if(v > u)
				separator_.push_back(v);

		workset_.clear();
		workset_.insert(workset_.end(), comp.vertices.begin(), comp.vertices.end());
		workset_.push_back(u);
		compose_q_.push({std::move(workset_), std::move(separator_), comp.h + 1, comp.trivial});
	}
}
