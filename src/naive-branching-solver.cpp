#include <climits>
#include <iostream>
#include <tuple>

#include "connected-components.hpp"
#include "naive-branching-solver.hpp"

// Note: this algorithm is intended to be a simple baseline
//       (e.g., to verify that other algorithms yield correct results).
//       It is *not* meant to be as sophisticated as possible.

naive_branching_solver::naive_branching_solver(graph &g)
: g_{&g} { }

int naive_branching_solver::compute_treedepth() {
	if(!g_->num_vertices())
		return 0;

	// Compute a non-trivial elimination tree if it exists.
	bool feasible;
	int s;
	std::tie(feasible, s) = recurse_(g_->num_vertices() - 1);
	std::cerr << "evaluated " << stats_.num_branches << " branches" << std::endl;
	if(!feasible)
		return g_->num_vertices();
	return s;
}

// Return true and compute the treedepth of the graph if it is <= k.
// Otherwise, return false.
std::pair<bool, int> naive_branching_solver::recurse_(int k) {

	// We should should never recursively call this with k < 0.
	// (In particular, the k = 0 case is handled without recursion.)
	assert(k >= 0);

	// Handle k = 0 (graph must be empty).
	if(!g_->num_vertices())
		return {true, 0};
	if(!k)
		return {false, k + 1};

	// TODO: for k = 1, we could explicitly check that the graph is an independent set.
	// TODO: for k = 2, we could explicitly check that the graph is a star.

	connected_components cc;
	cc.compute(*g_);

	// TODO: during the DFS that computes the CC, we can extract the maximal DFS depth.
	//       log_2(DFS depth) is a lower bound for treedepth since we find a long path
	//       as a subgraph in the input graph.

	bool feasible;
	int s;

	if(cc.num_components() > 1) {
		int worst = 0;
		for(size_t i = 0; i < cc.num_components(); i++) {
			for(vertex v : cc.component(i))
				g_->eliminate(v);
			stats_.num_branches++;
			std::tie(feasible, s) = recurse_(k);
			for(vertex v : cc.component(i))
				g_->uneliminate(v);
			if(!feasible)
				return {false, k + 1};
			worst = std::max(worst, s);
		}
		return {true, worst};
	}else{
		int best = k;
		for(vertex v : g_->vertices()) {
			g_->eliminate(v);
			stats_.num_branches++;
			std::tie(feasible, s) = recurse_(best - 1);
			g_->uneliminate(v);
			if(feasible)
				best = std::min(best, s);
		}
		// If all vertices v yield TD >= k, the TD of the current graph is >= k + 1.
		if(best == k)
			return {false, k + 1};
		return {true, 1 + best};
	}
}
