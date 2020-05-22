#pragma once

#include <queue>

#include "graph.hpp"

struct kernelization {

	struct statistics {
		int64_t num_universal = 0;
	};

	void compute(const graph &ig);
	void compute_trivial(const graph &ig);

	size_t num_components() {
		return components_.size();
	}

	const std::tuple<graph, int, vertex> &component(size_t i) {
		return components_[i];
	}

	int get_depth()
	{
		return maximal_depth_;
	}

	const std::vector<vertex>& get_kernel_decomposition()
	{
		return kernel_tree_;
	}

private:
	bool decompose_(graph &g, int depth, vertex root);
	bool reduce_universal_(graph &g, int depth, vertex root);

	std::vector<std::tuple<graph, int, vertex>> components_;
	statistics stats_;

	std::queue<std::tuple<graph, int, vertex>> q_;
	boolean_marker marker_;

	std::vector<vertex> kernel_tree_;
	
	int maximal_depth_;
};
