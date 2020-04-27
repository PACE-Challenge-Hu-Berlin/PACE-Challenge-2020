#pragma once

#include <queue>
#include <unordered_map>

#include "block-sieve.hpp"
#include "data-structures.hpp"
#include "graph.hpp"
#include "mm.hpp"
#include "utilities.hpp"
#include "precedence-by-inclusion.hpp"

struct simple_pid_solver {
	struct statistics {
		int64_t num_empty_separator = 0;
		int64_t num_pruned_forests = 0;
		int64_t num_pruned_compositions = 0;
		int64_t num_pruned_by_precedence = 0;
		int64_t num_protected_separators = 0;
		profiling_duration time_join{};
		profiling_duration time_join_assemble{};
		profiling_duration time_compose{};
	};

	struct feasible_tree {
		vertex_key vertices;
		int h;
	};

	struct feasible_forest {
		// A representative of a component is the vertex with smallest vertex ID in that component.
		// As the representatives are ordered (by their vertex ID), they induce an order
		// on the components.
		vertex rv; // Minimal representative of all trees.
		vertex sweep_rv; // Representative used to sweep over inactive trees.
		span<vertex> vertices;
		span<vertex> separator;
		// Forest only consists of a single tree.
		bool atomic;
		// Forest is atomic and the only tree is a path.
		bool trivial;
	};

	struct feasible_composition {
		span<vertex> vertices;
		span<vertex> prefix;
		span<vertex> separator;
	};

	struct staged_tree {
		vertex_key vertices;
		int h;
		span<vertex> separator;
	};


	simple_pid_solver(graph &g);

	int compute_treedepth();

	const std::vector<vertex> &decomposition() {
		return decomp_;
	}

	const statistics &stats() {
		return stats_;
	}

	bool no_kernelization = false;
	bool no_inclusion_precedence = false;
	bool no_protected_separators = false;

private:
	enum class ownership {
		none,
		owned,
		borrowed
	};

	bool decide_treedepth_(int k);
	bool recover_decomposition_();
	void join_(int k, int h, feasible_forest &forest);
	void compose_(int k, int h, feasible_composition &comp);

	staged_tree do_stage_(int h, const std::vector<vertex> &vertices,
			const std::vector<vertex> &separator);

	graph *g_;
	graph sg_;
	memory_arena eternal_arena_;
	queue_memory join_memory_;
	queue_memory compose_memory_;
	std::vector<vertex> decomp_;
	statistics stats_;

	precedence_by_inclusion inclusion_precedence_;

	// The following data structure store all feasible trees of height <= current h.
	// Set of trees that are currently being expanded. Those have height == current h.
	separated_sieve<feasible_tree> active_trees;
	// Set of trees that were already expanded. Those have height < current h.
	separated_sieve<feasible_tree> inactive_trees;

	// (Incomplete) set of feasible trees of height > current h.
	std::unordered_map<vertex_key, staged_tree> staged_trees;

	// Queues to process newly found feasible forests and feasible compositions.
	std::queue<std::pair<feasible_forest, ownership>> join_q_;
	std::queue<std::pair<feasible_composition, ownership>> compose_q_;

	// Auxiliary data structures to compute neighbor sets and separators etc.
	separated_sieve_query<feasible_tree> sieve_query_;
	boolean_marker pivot_marker_;
	boolean_marker pivot_neighbor_marker_;
	boolean_marker protected_marker_;
	boolean_marker associate_marker_;
	boolean_marker associate_neighbor_marker_;
	std::vector<vertex> workset_;
	std::vector<vertex> separator_;
	std::vector<vertex> candidates_;

	int64_t num_join_;
	int64_t num_compose_;
	int64_t num_stage_;
	int64_t num_unimproved_;
};
