#pragma once

#include "data-structures.hpp"
#include "graph.hpp"

struct connected_components {
	struct vertex_range {
		vertex_range(std::vector<vertex>::const_iterator begin,
				std::vector<vertex>::const_iterator end)
		: begin_{begin}, end_{end} { }

		auto begin() const {
			return begin_;
		}

		auto end() const {
			return end_;
		}

	private:
		std::vector<vertex>::const_iterator begin_;
		std::vector<vertex>::const_iterator end_;
	};

	size_t num_components() const {
		assert(cc_indices_.size());
		return cc_indices_.size() - 1;
	}

	vertex_range component(size_t i) const {
		return vertex_range{cc_lists_.begin() + cc_indices_[i],
				cc_lists_.begin() + cc_indices_[i + 1]};
	}

	void compute(const graph &g) {
		visited_.reset(g.id_limit());
		cc_lists_.clear();
		cc_indices_.clear();

		for(vertex r : g.vertices()) {
			if(visited_.is_marked(r))
				continue;

			// This vertex starts a new CC.
			cc_indices_.push_back(cc_lists_.size());

			stack_.push_back({r, g.neighbors(r).begin()});
			visited_.mark(r);
			cc_lists_.push_back(r);

			while(!stack_.empty()) {
				auto v = stack_.back().first;
				auto &it = stack_.back().second;

				// Skip neighbors that are already visited.
				while(it != g.neighbors(v).end() && visited_.is_marked(*it))
					++it;

				if(it == g.neighbors(v).end()) {
					stack_.pop_back();
					continue;
				}

				// Note: the reference to it is invalidated by the push_back!
				auto w = *it;
				++it;
				stack_.push_back({w, g.neighbors(w).begin()});
				visited_.mark(w);
				cc_lists_.push_back(w);
			}
		}
		// Store number of CCs in last element to avoid branch in component().
		cc_indices_.push_back(cc_lists_.size());
	}

private:
	// Store all CCs in a single array to avoid memory allocations.
	// Index into this array with a second array.
	std::vector<vertex> cc_lists_;
	std::vector<size_t> cc_indices_;

	std::vector<std::pair<vertex, graph::neighbor_range::iterator>> stack_;
	boolean_marker visited_;
};
