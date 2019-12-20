#pragma once

#include "data-structures.hpp"
#include "graph.hpp"
#include "directed-graph.hpp"

#include <stack>

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

	template<typename G>
	void compute(const G &g) {
		// We would like to avoid the allocation here (by storing the stack in the class)
		// but that is not possible due to the template parameter G.
		// TODO: add a caching memory allocator that allows us to reuse the allocation.
		//       Right now, that would not really help since connected_components is not
		//       on a hot path in the PID code and the naive branching code duplicates the
		//       connected_components object anyway.
		std::vector<std::pair<vertex, typename G::neighbor_range::iterator>> stack;

		visited_.reset(g.id_limit());
		cc_lists_.clear();
		cc_indices_.clear();

		for(vertex r : g.vertices()) {
			if(visited_.is_marked(r))
				continue;

			// This vertex starts a new CC.
			cc_indices_.push_back(cc_lists_.size());

			stack.push_back({r, g.neighbors(r).begin()});
			visited_.mark(r);
			cc_lists_.push_back(r);

			while(!stack.empty()) {
				auto v = stack.back().first;
				auto &it = stack.back().second;

				// Skip neighbors that are already visited.
				while(it != g.neighbors(v).end() && visited_.is_marked(*it))
					++it;

				if(it == g.neighbors(v).end()) {
					stack.pop_back();
					continue;
				}

				// Note: the reference to it is invalidated by the push_back!
				auto w = *it;
				++it;
				stack.push_back({w, g.neighbors(w).begin()});
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

	boolean_marker visited_;
};


struct strongly_connected_components
{
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

	vertex_range component(size_t i) const {
		return vertex_range{scc_lists_.begin() + scc_indices_[i],
				scc_lists_.begin() + scc_indices_[i + 1]};
	}

	size_t component_id(vertex v) const {
		return scc_id_[v];
	}

	size_t num_components() const {
		assert(scc_indices_.size());
		return scc_indices_.size() - 1;
	}

	size_t component_size(size_t i) const {
		return scc_indices_[i+1] - scc_indices_[i];
	}

	void compute(const directed_graph& g_)
	{
		n = g_.num_vertices();
	       	done.assign(n, 0);
	       	in_stack.assign(n, 0),
	       	depth.resize(n);
	       	back_depth.resize(n);
	       	scc_id_.resize(n);
		scc_lists_.clear();
		scc_indices_.clear();
		for(unsigned int v = 0; v < n; v++)
		{
			assert(!in_stack[v]);
			if(done[v])
			{
				continue;
			}
			dfs(v, 0, g_);
		}
		assert(scc_lists_.size() == n);
		scc_indices_.push_back(n);
	}


private:
	unsigned int n;

	std::stack<unsigned int> stack_;

	std::vector<bool> done;
	std::vector<bool> in_stack;
	std::vector<unsigned int> depth;
	std::vector<unsigned int> back_depth;

	std::vector<vertex> scc_lists_;
	std::vector<size_t> scc_indices_;
	std::vector<size_t> scc_id_;

	void dfs(unsigned int v, unsigned int depth_, const directed_graph& g_){
		stack_.push(v);
		in_stack[v] = 1;
		depth[v] = back_depth[v] = depth_;

		for(auto u : g_.neighbors(v))
		{
			if(!in_stack[u] && !done[u])
			{
				dfs(u, depth_ + 1, g_); 
				back_depth[v] = std::min(back_depth[v], back_depth[u]);
			}
			else if(in_stack[u])
			{
				
				back_depth[v] = std::min(back_depth[v], depth[u]);
			}
		}
		if(depth[v] == back_depth[v])
		{
			scc_indices_.push_back(scc_lists_.size());
			bool found_v = false;
			while(!found_v)
			{
				unsigned int u = stack_.top();
				if(u == v)
				{
					found_v = true;
				}
				stack_.pop();
				scc_lists_.push_back(u);
				scc_id_[u] = scc_indices_.size() - 1;

				done[u] = 1;
				in_stack[u] = 0;

			}
		}
	}
};
