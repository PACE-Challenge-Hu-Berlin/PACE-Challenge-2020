#pragma once

#include <vector>

// Data structure to mark integers (e.g., vertices of a graph).
// Use a timestamping technique to accelerate resets.
struct boolean_marker {
	void reset(size_t n) {
		vec_.resize(n);
		if(ts_ == UINT8_MAX) {
			std::fill(vec_.begin(), vec_.end(), 0);
			ts_ = 0;
		}
		++ts_;
	}

	void mark(size_t x) {
		vec_[x] = ts_;
	}

	bool is_marked(size_t x) {
		return vec_[x] == ts_;
	}

private:
	std::vector<uint8_t> vec_;
	uint8_t ts_ = 0;
};
