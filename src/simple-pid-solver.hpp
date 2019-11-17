#pragma once

#include <queue>
#include <unordered_map>

#include "block-sieve.hpp"
#include "data-structures.hpp"
#include "graph.hpp"
#include "mm.hpp"

struct simple_pid_solver {
	struct statistics {
		int64_t num_empty_separator = 0;
		int64_t num_pruned_forests = 0;
		int64_t num_pruned_compositions = 0;
	};

	struct feasible_tree {
		vertex_span vertices;
		int h;
	};

	struct feasible_forest {
		// A representative of a component is the vertex with smallest vertex ID in that component.
		// As the representatives are ordered (by their vertex ID), they induce an order
		// on the component.
		vertex min_rv; // Minimal representative of all components.
		vertex max_rv; // Maximal representative of all components.
		vertex_span vertices;
		vertex_span separator;
		// Forest only consists of a single tree.
		bool atomic;
		// Forest is atomic and the only tree is a path.
		bool trivial;
	};

	struct feasible_composition {
		vertex_span vertices;
		vertex_span prefix;
		int h;
		bool trivial;
	};

	struct staged_tree {
		vertex_span vertices;
		int h;
		// Vertices can only be arranged in a path.
		bool trivial;
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
	queue_memory join_memory_;
	queue_memory compose_memory_;
	statistics stats_;

	// The following data structure store all feasible trees of height <= current h.
	// Set of trees that are currently being expanded. Those have height == current h.
	block_sieve<feasible_tree> active_trees;
	// Set of trees that were already expanded. Those have height < current h.
	block_sieve<feasible_tree> inactive_trees;

	// (Incomplete) set of feasible trees of height > current h.
	std::unordered_map<vertex_span, staged_tree> staged_trees;

	// Queues to process newly found feasible forests and feasible compositions.
	std::queue<feasible_forest> join_q_;
	std::queue<feasible_composition> compose_q_;

	// Auxiliary data structures to compute neighbor sets and separators etc.
	block_sieve_query<feasible_tree> sieve_query_;
	boolean_marker pivot_marker_;
	boolean_marker pivot_neighbor_marker_;
	boolean_marker associate_marker_;
	boolean_marker associate_neighbor_marker_;
	std::vector<vertex> workset_;
	std::vector<vertex> separator_;

	int64_t num_join_;
	int64_t num_compose_;
	int64_t num_stage_;
	int64_t num_unimproved_;
};
