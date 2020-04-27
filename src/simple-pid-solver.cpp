#include <cstring>
#include <iostream>

#include "connected-components.hpp"
#include "kernelization.hpp"
#include "precedence-by-inclusion.hpp"
#include "simple-pid-solver.hpp"
#include "subgraph.hpp"

namespace {
	constexpr bool debug_trees = false;

	// Helper class to print object during debugging.
	struct print_staged {
		friend std::ostream &operator<< (std::ostream &os, const print_staged &self) {
			os << "({";
			for(size_t i = 0; i < self.p_->vertices.size(); i++) {
				if(i > 0)
					os << ", ";
				os << self.p_->vertices[i];
			}
			os << "}, h=" << self.p_->h << ")";
			return os;
		}

		print_staged(const simple_pid_solver::staged_tree &p)
		: p_{&p} { }

	private:
		const simple_pid_solver::staged_tree *p_;
	};

	struct print_forest {
		friend std::ostream &operator<< (std::ostream &os, const print_forest &self) {
			os << "({";
			for(size_t i = 0; i < self.p_->vertices.size(); i++) {
				if(i > 0)
					os << ", ";
				os << self.p_->vertices[i];
			}
			os << "}, R=" << self.p_->rv << ", S=" << self.p_->sweep_rv << ")";
			return os;
		}

		print_forest(const simple_pid_solver::feasible_forest &p)
		: p_{&p} { }

	private:
		const simple_pid_solver::feasible_forest *p_;
	};

	struct print_tree {
		friend std::ostream &operator<< (std::ostream &os, const print_tree &self) {
			os << "({";
			for(size_t i = 0; i < self.p_->vertices.size(); i++) {
				if(i > 0)
					os << ", ";
				os << self.p_->vertices[i];
			}
			os << "})";
			return os;
		}

		print_tree(const simple_pid_solver::feasible_tree &p)
		: p_{&p} { }

	private:
		const simple_pid_solver::feasible_tree *p_;
	};

	struct print_composition {
		friend std::ostream &operator<< (std::ostream &os, const print_composition &self) {
			os << "({";
			for(size_t i = 0; i < self.p_->vertices.size(); i++) {
				if(i > 0)
					os << ", ";
				os << self.p_->vertices[i];
			}
			os << "}, {";
			for(size_t i = 0; i < self.p_->separator.size(); i++) {
				if(i > 0)
					os << ", ";
				os << self.p_->separator[i];
			}
			os << "}, candidates={";
			for(size_t i = 0; i < self.p_->candidates.size(); i++) {
				if(i > 0)
					os << ", ";
				os << self.p_->candidates[i];
			}
			os << "})";
			return os;
		}

		print_composition(const simple_pid_solver::feasible_composition &p)
		: p_{&p} { }

	private:
		const simple_pid_solver::feasible_composition *p_;
	};

	template<typename C>
	span<typename C::value_type> copy_to_arena(const C &vec, memory_arena &memory) {
		arena_allocator<typename C::value_type> alloc{memory};
		auto p = alloc.allocate(vec.size());
		memcpy(p, vec.data(), sizeof(typename C::value_type) * vec.size());
		return span<typename C::value_type>{p, p + vec.size()};
	}

	template<typename C>
	span<typename C::value_type> copy_to_queue(const C &vec, queue_memory &memory) {
		queue_allocator<typename C::value_type> alloc{memory};
		auto p = alloc.allocate(vec.size());
		memcpy(p, vec.data(), sizeof(typename C::value_type) * vec.size());
		return span<typename C::value_type>{p, p + vec.size()};
	}

	template<typename T>
	void free_in_arena(span<T> vs, memory_arena &memory) {
		arena_allocator<T> alloc{memory};
		alloc.deallocate(const_cast<T *>(vs.data()), vs.size());
	}

	template<typename T>
	void free_in_queue(span<T> vs, queue_memory &memory) {
		queue_allocator<T> alloc{memory};
		alloc.deallocate(const_cast<T *>(vs.data()), vs.size());
	}
}

simple_pid_solver::simple_pid_solver(graph &g)
: g_{&g} { }

int simple_pid_solver::compute_treedepth() {
	decomp_.resize(g_->id_limit(), nil_vertex());

	kernelization kern;
	if(no_kernelization)
		kern.compute_trivial(*g_);
	else
		kern.compute(*g_);

	int td = 0;
	for(size_t i = 0; i < kern.num_components(); ++i) {
		sg_ = kern.component(i);

		int k = 0;
		std::cerr << "solving component " << (i + 1)
				<< " of " << kern.num_components() << std::endl;
		while(true) {
			try {
				if(decide_treedepth_(k))
					break;
			} catch(const std::bad_alloc &) {
				std::cerr << "allocation failure in simple-pid solver\n"
						<< "    eternal memory is at: "
						<< print_memory(eternal_arena_.used_space()) << "\n"
						<< "    join memory is at: "
						<< print_memory(join_memory_.used_space()) << "\n"
						<< "    compose memory is at: "
						<< print_memory(compose_memory_.used_space()) << std::endl;
				throw;
			}
			k++;
		}
		td = std::max(td, k);
	}

	return td;
}

bool simple_pid_solver::decide_treedepth_(int k) {
	active_trees.clear();
	inactive_trees.clear();
	for(const auto &entry : staged_trees) {
		const auto &staged = entry.second;
		free_in_arena(staged.vertices.as_span(), eternal_arena_);
		free_in_arena(staged.separator, eternal_arena_);
	}
	staged_trees.clear();
	eternal_arena_.reset();

	if(no_inclusion_precedence)
	{
		inclusion_precedence_.compute_trivial(sg_);
	}
	else
	{
		inclusion_precedence_.compute(sg_);
	}

	pivot_marker_.reset(inclusion_precedence_.num_components());
	for(vertex  v : sg_.vertices()) {
		auto cmp_id = inclusion_precedence_.component_id(v);

		if(pivot_marker_.is_marked(cmp_id))
			continue;
		pivot_marker_.mark(cmp_id);

		if(!inclusion_precedence_.is_maximal(cmp_id))
		{
			++stats_.num_pruned_by_precedence;
			continue;
		}

		pivot_neighbor_marker_.reset(g_->id_limit());
		int neighbors_counter = 0;
		for(auto u : inclusion_precedence_.component(cmp_id))
		{
			if(!pivot_neighbor_marker_.is_marked(u))
			{
				pivot_neighbor_marker_.mark(u);
				neighbors_counter++;
			}

			for(auto v : g_->neighbors(u))
			{
				if(!pivot_neighbor_marker_.is_marked(v))
				{
					pivot_neighbor_marker_.mark(v);
					neighbors_counter++;
				}
			}
		}
		if(neighbors_counter > k)
		{
			continue;
		}

		workset_.clear();
		for(auto u : inclusion_precedence_.component(cmp_id))
		{
			workset_.push_back(u);
		}
		std::sort(workset_.begin(), workset_.end());

		staged_tree staged{copy_to_arena(workset_, eternal_arena_),
				static_cast<int>(inclusion_precedence_.component_size(cmp_id)),
				copy_to_arena(workset_, eternal_arena_)};
		staged_trees.emplace(staged.vertices, staged);
	}
	
	std::cerr <<  "Pruned components from dynamic programming base " << stats_.num_pruned_by_precedence << std::endl;

	for(int h = 1; h < k; h++) {
		std::cerr << "constructing k = " << k << ", h = " << h << std::endl;

		// Move trees for the current h from the staging buffer to the set of feasible trees.
		bool all_trees_inactive = true;
		for(const auto &entry : staged_trees) {
			const auto &staged = entry.second;
			if(staged.h >= h)
				all_trees_inactive = false;
			if(staged.h != h)
				continue;
			active_trees.insert(sg_, staged.vertices, feasible_tree{staged.vertices, staged.h});
		}

		if(all_trees_inactive)
			break;

		num_join_ = 0;
		num_compose_ = 0;
		num_stage_ = 0;
		num_unimproved_ = 0;
		std::cerr << "    there are "
				<< active_trees.size() << " active trees, "
				<< inactive_trees.size() << " inactive trees" << std::endl;

		for(const auto &tree : active_trees) {
			if(tree.h != h)
				continue;

			// Find the neighbor set of the tree to turn it into an atomic forest.
			pivot_marker_.reset(sg_.id_limit());
			pivot_neighbor_marker_.reset(sg_.id_limit());
			candidates_.clear();
			for(vertex v : tree.vertices)
				pivot_marker_.mark(v);
			for(vertex v : tree.vertices)
				for(vertex w : sg_.neighbors(v)) {
					if(pivot_marker_.is_marked(w))
						continue;
					if(pivot_neighbor_marker_.is_marked(w))
						continue;
					pivot_neighbor_marker_.mark(w);
					candidates_.push_back(w);
				}

			auto trivial = (static_cast<int>(tree.vertices.size()) == h);
			if(trivial) {
				for(vertex u : candidates_) {
					associate_neighbor_marker_.reset(sg_.id_limit());
					int tree_neighbors = static_cast<int>(candidates_.size()) - 1;
					for(vertex v : sg_.neighbors(u)) {
						if(pivot_marker_.is_marked(v) || pivot_neighbor_marker_.is_marked(v))
							continue;
						assert(v != u);
						if(associate_neighbor_marker_.is_marked(v))
							continue;
						associate_neighbor_marker_.mark(v);
						++tree_neighbors;
					}

					if(h + 1 + tree_neighbors > k)
						continue;

					workset_.clear();
					workset_.insert(workset_.end(), tree.vertices.begin(), tree.vertices.end());
					workset_.push_back(u);
					std::sort(workset_.begin(), workset_.end());

					separator_.clear();
					separator_.push_back(u);

					do_stage_(h + 1, workset_, separator_);
				}
			}

			auto rv = tree.vertices.front();
			join_q_.push({feasible_forest{rv, sg_.id_limit(),
					tree.vertices.as_span(),
					copy_to_queue(candidates_, join_memory_),
					true}, ownership::borrowed});
			join_memory_.seal();
		}

		while(!join_q_.empty()) {
			feasible_forest forest = join_q_.front().first;
			ownership forest_ownership = compose_q_.front().second;
			join_q_.pop();
			join_(k, h, forest);
			num_join_++;

			while(!compose_q_.empty()) {
				feasible_composition comp = compose_q_.front().first;
				ownership comp_ownership = compose_q_.front().second;
				compose_q_.pop();
				compose_(k, h, comp);
				num_compose_++;

				if(comp_ownership == ownership::owned) {
					free_in_queue(comp.separator, compose_memory_);
					free_in_queue(comp.candidates, compose_memory_);
				}
				compose_memory_.reclaim();
			}

			if(forest_ownership == ownership::owned)
				free_in_queue(forest.vertices, join_memory_);
			free_in_queue(forest.candidates, join_memory_);
			join_memory_.reclaim();
		}

		std::cerr << "    staged " << num_stage_
				<< " trees, discarded " << num_unimproved_
				<< " (" << num_join_ << " joins, " << num_compose_ << " compositions)"
				<< std::endl;
		std::cerr << "    empty separator: " << stats_.num_empty_separator << std::endl;
		std::cerr << "    pruned forests: " << stats_.num_pruned_forests
				<< ", compositions: " << stats_.num_pruned_compositions << std::endl;
		std::cerr << "    protected separators: " << stats_.num_protected_separators
				<< std::endl;
		std::cerr << "    eternal memory: " << print_memory{eternal_arena_.used_space()}
				<< ", allocations: " << eternal_arena_.num_allocations()
				<< std::endl;
		std::cerr << "    join memory: " << print_memory{join_memory_.max_used_space()}
				<< ", compose memory: " << print_memory{compose_memory_.max_used_space()}
				<< std::endl;
		std::cerr << "    join time: " << print_time{stats_.time_join}
				<< ", assembling: " << print_time{stats_.time_join_assemble}
				<< ", compose: " << print_time{stats_.time_compose} << std::endl;

		for(const feasible_tree &tree : active_trees)
			inactive_trees.insert(sg_, tree.vertices, tree);
		active_trees.clear();
	}

	return recover_decomposition_();
}

bool simple_pid_solver::recover_decomposition_() {
	workset_.clear();
	for(vertex v : sg_.vertices())
		workset_.push_back(v);
	auto root_it = staged_trees.find(vertex_key{workset_});
	if(root_it == staged_trees.end())
		return false;

	for(vertex v : sg_.vertices())
		decomp_[v] = nil_vertex();

	std::stack<std::tuple<vertex, span<vertex>, span<vertex>>> stack;
	stack.push({0, root_it->second.vertices.as_span(), root_it->second.separator});

	connected_components cc;
	while(!stack.empty()) {
		auto e = stack.top();
		stack.pop();

		// u is the root of the current subtree.
		// Arrange the separator in an arbitrary order.
		vertex u = std::get<0>(e);
		for(vertex v : std::get<2>(e)) {
			decomp_[v] = u;
			u = v;
		}

		pivot_marker_.reset(sg_.id_limit());
		for(vertex v : std::get<1>(e))
			pivot_marker_.mark(v);
		for(vertex v : std::get<2>(e)) {
			assert(pivot_marker_.is_marked(v));
			pivot_marker_.unmark(v);
		}

		cc.compute(induced_subgraph{sg_, pivot_marker_});
		if(!cc.num_components())
			continue;

		for(size_t i = 0; i < cc.num_components(); ++i) {
			workset_.clear();
			for(vertex v : cc.component(i))
				workset_.push_back(v);
			std::sort(workset_.begin(), workset_.end());

			auto it = staged_trees.find(vertex_key{workset_});
			assert(it != staged_trees.end());
			stack.push({u, it->second.vertices.as_span(), it->second.separator});
		}
	}

	return true;
}

void simple_pid_solver::join_(int k, int h, feasible_forest &forest) {
	profiling_timer join_timer;

	// Mark and count the neighbors of the feasible forest.
	pivot_marker_.reset(sg_.id_limit());
	pivot_neighbor_marker_.reset(sg_.id_limit());
	unsigned int num_forest_neighbors = 0;
	for(vertex v : forest.vertices)
		pivot_marker_.mark(v);
	for(vertex v : forest.vertices) {
		for(vertex w : sg_.neighbors(v)) {
			if(pivot_marker_.is_marked(w))
				continue;
			if(pivot_neighbor_marker_.is_marked(w))
				continue;
			pivot_neighbor_marker_.mark(w);
			++num_forest_neighbors;
		}
	}

	protected_marker_.reset(sg_.id_limit());
	bool has_protected_separators = false;
	if(!no_protected_separators)
		for(vertex v : forest.candidates) {
			bool protect = true;
			for(vertex w : sg_.neighbors(v))
				if(!pivot_marker_.is_marked(w) && !pivot_neighbor_marker_.is_marked(w)) {
					protect = false;
					break;
				}
			if(protect) {
				protected_marker_.mark(v);
				has_protected_separators = true;
			}
		}

	auto join_with = [&] (const feasible_tree &tree) {
		int num_total_neighbors = num_forest_neighbors;

		// Determine the neighbors of the tree that we are trying to join.
		associate_marker_.reset(sg_.id_limit());
		associate_neighbor_marker_.reset(sg_.id_limit());
		for(vertex v : tree.vertices) {
			assert(!pivot_marker_.is_marked(v));
			associate_marker_.mark(v);
		}

		for(vertex v : tree.vertices) {
			for(vertex w : sg_.neighbors(v)) {
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
		for(vertex v : forest.candidates)
			if(associate_neighbor_marker_.is_marked(v)) {
				separator_.push_back(v);
			}else if(protected_marker_.is_marked(v)) {
				++stats_.num_protected_separators;
				return;
			}

		if(!separator_.size()) {
			++stats_.num_empty_separator;
			return;
		}

		workset_.clear();
		workset_.insert(workset_.end(), forest.vertices.begin(), forest.vertices.end());
		workset_.insert(workset_.end(), tree.vertices.begin(), tree.vertices.end());

		// For active trees, it holds that rv < forest.rv; however, for inactive trees,
		// it it only guaranteed that rv < forest.sweep_rv.
		// It always holds that forest.rv < forest.sweep_rv.
		auto rv = tree.vertices.front();
		assert(rv < forest.sweep_rv);
		join_q_.push({feasible_forest{std::min(rv, forest.rv), rv,
				copy_to_queue(workset_, join_memory_),
				copy_to_queue(separator_, join_memory_),
				false}, ownership::owned});
		join_memory_.seal();
	};
	auto join_with_other = [&] (const feasible_tree &tree) {
		profiling_timer assemble_timer;
		join_with(tree);
		stats_.time_join_assemble += assemble_timer.elapsed();
	};

	auto member_predicate = [&] (vertex v) {
		return pivot_marker_.is_marked(v);
	};
	auto neighbor_predicate = [&] (vertex v) {
		return pivot_neighbor_marker_.is_marked(v);
	};

	// To avoid duplicate joins, we only join trees into this forest that
	// are ordered before all components of this forest.
	// Note that we have to distinguish active and inactive forests.

	// Join only preceeding active trees into this forest.
	active_trees.list_predecessors(sg_, sieve_query_,
				k - h - num_forest_neighbors, forest.rv,
				member_predicate, neighbor_predicate, join_with_other);
	// Join both preceeding and succeeding inactive trees into this forest.
	inactive_trees.list_predecessors(sg_, sieve_query_,
				k - h - num_forest_neighbors, forest.sweep_rv,
				member_predicate, neighbor_predicate, join_with_other);
	stats_.time_join += join_timer.elapsed();

	if(!forest.atomic) {
		if(has_protected_separators) {
			separator_.clear();
			candidates_.clear();
			for(vertex v : forest.candidates) {
				if(protected_marker_.is_marked(v)) {
					separator_.push_back(v);
				}else{
					candidates_.push_back(v);
				}
			}

			int ch = h + static_cast<int>(separator_.size());
			assert(ch + static_cast<int>(candidates_.size()) <= k);

			workset_.clear();
			workset_.insert(workset_.end(), forest.vertices.begin(), forest.vertices.end());
			workset_.insert(workset_.end(), separator_.begin(), separator_.end());
			std::sort(workset_.begin(), workset_.end());

			auto st = do_stage_(ch, workset_, separator_);

			compose_q_.push({feasible_composition{st.vertices.as_span(),
					copy_to_queue(separator_, compose_memory_),
					copy_to_queue(candidates_, compose_memory_)}, ownership::owned});
			compose_memory_.seal();
		}else{
			compose_q_.push({feasible_composition{forest.vertices,
					span<vertex>{},
					forest.candidates}, ownership::borrowed});
			compose_memory_.seal();
		}
	}
}

void simple_pid_solver::compose_(int k, int h, feasible_composition &comp) {
	profiling_timer compose_timer;

	pivot_marker_.reset(sg_.id_limit());
	pivot_neighbor_marker_.reset(sg_.id_limit());
	int num_forest_neighbors = 0;
	for(vertex v : comp.vertices)
		pivot_marker_.mark(v);
	for(vertex v : comp.vertices) {
		for(vertex w : sg_.neighbors(v)) {
			if(pivot_marker_.is_marked(w))
				continue;
			if(pivot_neighbor_marker_.is_marked(w))
				continue;
			pivot_neighbor_marker_.mark(w);
			++num_forest_neighbors;
		}
	}

	for(vertex u : comp.candidates) {
		// Count number of neighbors of the composition.
		// Since u is a neighbor of the vertex set, subtract one to begin with.
		int num_composed_neighbors = num_forest_neighbors - 1;
		for(vertex w : sg_.neighbors(u)) {
			if(pivot_marker_.is_marked(w) || pivot_neighbor_marker_.is_marked(w))
				continue;
			++num_composed_neighbors;
		}

		int ch = h + static_cast<int>(comp.separator.size()) + 1;

		if(ch + num_composed_neighbors > k) {
			++stats_.num_pruned_compositions;
			continue;
		}

		bool precedence_okay = true;
		for(auto v : inclusion_precedence_.vertex_successors(u))
		{
			if(!pivot_marker_.is_marked(v) && !pivot_neighbor_marker_.is_marked(v))
			{
				precedence_okay = false;
				break;
			}
		}
		if(!precedence_okay)
		{
			continue;
		}

		workset_.clear();
		workset_.insert(workset_.end(), comp.vertices.begin(), comp.vertices.end());
		workset_.push_back(u);
		std::sort(workset_.begin(), workset_.end());

		separator_.clear();
		separator_.insert(separator_.end(), comp.separator.begin(), comp.separator.end());
		separator_.push_back(u);

		auto st = do_stage_(ch, workset_, separator_);

		candidates_.clear();
		for(vertex v : comp.candidates)
			if(v > u)
				candidates_.push_back(v);

		compose_q_.push({feasible_composition{st.vertices.as_span(),
				copy_to_queue(separator_, compose_memory_),
				copy_to_queue(candidates_, compose_memory_)}, ownership::owned});
		compose_memory_.seal();
	}
	stats_.time_compose += compose_timer.elapsed();
}

// Helper function to stage trees.
auto simple_pid_solver::do_stage_(int h,
		const std::vector<vertex> &vertices,
		const std::vector<vertex> &separator) -> staged_tree {
	auto it = staged_trees.find(vertex_key{vertices});
	if(it == staged_trees.end()) {
		bool success;
		auto key = vertex_key{copy_to_arena(vertices, eternal_arena_)};
		std::tie(it, success) = staged_trees.emplace(key, staged_tree{});
		assert(success);
		it->second.vertices = key;
		it->second.h = h;
		it->second.separator = copy_to_arena(separator, eternal_arena_);
		++num_stage_;
	}else if(it->second.h > h) {
		free_in_arena(it->second.separator, eternal_arena_);
		it->second.h = h;
		it->second.separator = copy_to_arena(separator, eternal_arena_);
		++num_stage_;
	}else{
		++num_unimproved_;
	}
	return it->second;
}
