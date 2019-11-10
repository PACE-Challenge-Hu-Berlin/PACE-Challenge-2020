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
// Otherwise, return false and (k + 1) as treedepth.
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

		// Eliminate all but the first component.
		for(size_t j = 1; j < cc.num_components(); j++)
			for(vertex v : cc.component(j))
				g_->eliminate(v);

		// Invariant: Only the i-th component is uneliminated.
		size_t i;
		for(i = 0; /* exit through break below */ ; i++) {
			stats_.num_branches++;
			std::tie(feasible, s) = recurse_(k);
			worst = std::max(worst, s);

			if(!feasible)
				break;

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

		if(worst == k + 1)
			return {false, k + 1};
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
