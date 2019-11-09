#pragma once

#include <utility>

#include "graph.hpp"

struct naive_branching_solver {
	struct statistics {
		int num_branches = 0;
	};

	naive_branching_solver(graph &g);

	int compute_treedepth();

	const statistics &stats() {
		return stats_;
	}

private:
	std::pair<bool, int> recurse_(int k);

	graph *g_;
	statistics stats_;
};
