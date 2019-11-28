#pragma once

#include <queue>

#include "graph.hpp"

struct kernelization {
	struct statistics {
		int64_t num_universal = 0;
	};

	void compute(const graph &ig);

	size_t num_components() {
		return components_.size();
	}

	const graph &component(size_t i) {
		return components_[i];
	}

private:
	bool decompose_(graph &g);
	bool reduce_universal_(graph &g);

	std::vector<graph> components_;
	statistics stats_;

	std::queue<graph> q_;
	boolean_marker marker_;
};
