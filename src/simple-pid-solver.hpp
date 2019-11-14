#pragma once

#include <queue>
#include <unordered_map>

#include "data-structures.hpp"
#include "graph.hpp"
#include "mm.hpp"

struct simple_pid_solver {
	struct statistics {
		int64_t num_non_disjoint = 0;
		int64_t num_non_separated = 0;
		int64_t num_empty_separator = 0;
	};

	// Lightweight reference to a set of vertices.
	struct vertex_span {
		template<typename A>
		vertex_span(const std::vector<vertex, A> &vec)
		: begin{vec.data()}, end{vec.data() + vec.size()} { }

		const vertex *begin;
		const vertex *end;
	};

	struct feasible_tree {
		std::vector<vertex, arena_allocator<vertex>> vertices;
		int h;
	};

	struct feasible_forest {
		size_t idx;
		std::vector<vertex> vertices;
		std::vector<vertex> separator;
		// Forest only consists of a single tree.
		bool atomic;
		// Forest is atomic and the only tree is a path.
		bool trivial;
	};

	struct feasible_composition {
		std::vector<vertex> vertices;
		std::vector<vertex> prefix;
		int h;
		bool trivial;
	};

	struct staged_tree {
		std::vector<vertex, arena_allocator<vertex>> vertices;
		int h;
		// Vertices can only be arranged in a path.
		bool trivial;
	};

	struct staging_hash {
		// Possibly the worst hash function ever imagined.
		// TODO: Pick a reasonable hash function.
		size_t operator() (const vertex_span &vs) const {
			size_t hash = 0x12345678;
			for(auto p = vs.begin; p != vs.end; ++p)
				hash = 13 * hash + *p;
			return hash;
		}
	};

	struct staging_equals {
		bool operator() (const vertex_span &vs1, const vertex_span &vs2) const {
			size_t sz1 = vs1.end - vs1.begin;
			size_t sz2 = vs1.end - vs1.begin;
			if(sz1 != sz2)
				return false;
			for(size_t i = 0; i < sz1; i++)
				if(vs1.begin[i] != vs2.begin[i])
					return false;
			return true;
		}
	};

	simple_pid_solver(graph &g);

	int compute_treedepth();

	const statistics &stats() {
		return stats_;
	}

private:
	bool decide_treedepth_(int k);
	void join_(int k, int h, feasible_forest &forest);
	void compose_(int k, int h, feasible_composition &comp);

	graph *g_;
	memory_arena eternal_arena_;
	statistics stats_;

	// Set of all feasible trees of height <= current h.
	std::vector<feasible_tree> feasible_trees;

	// (Incomplete) set of feasible trees of height > current h.
	std::unordered_map<
		vertex_span,
		staged_tree,
		staging_hash,
		staging_equals
	> staged_trees;

	// Queues to process newly found feasible forests and feasible compositions.
	std::queue<feasible_forest> join_q_;
	std::queue<feasible_composition> compose_q_;

	// Auxiliary data structures to compute neighbor sets and separators etc.
	boolean_marker pivot_marker_;
	boolean_marker pivot_neighbor_marker_;
	boolean_marker associate_marker_;
	boolean_marker associate_neighbor_marker_;
	std::vector<vertex> workset_;
	std::vector<vertex> separator_;

	int64_t num_join_;
	int64_t num_compose_;
};
