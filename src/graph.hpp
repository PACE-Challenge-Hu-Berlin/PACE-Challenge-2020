#pragma once

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <vector>

using vertex = unsigned int;

// TODO: vertex nil_vertex = static_cast<vertex>(-1);

struct graph {
	// Range that iterates over all present vertices of the graph.
	struct vertex_range {
		struct sentinel { };

		struct iterator {
			explicit iterator(const graph *g)
			: g_{g} {
				skip_();
			}

			bool operator!= (sentinel) const {
				return v_ != g_->id_limit_;
			}

			iterator &operator++ () {
				assert(v_ < g_->id_limit_);
				++v_;
				skip_();
				return *this;
			}

			vertex operator* () const {
				return v_;
			}

		private:
			void skip_() {
				while(v_ != g_->id_limit_ && !g_->present_[v_])
					++v_;
			}

			const graph *g_;
			vertex v_ = 0;
		};

		explicit vertex_range(const graph *g)
		: g_{g} { }

		iterator begin() const {
			return iterator{g_};
		}

		sentinel end() const {
			return sentinel{};
		}

	private:
		const graph *g_;
	};

	// Range that iterates over all edges incident to present vertices.
	struct neighbor_range {
		struct sentinel { };

		struct iterator {
			explicit iterator(const graph *g, const vertex *c, const vertex *end)
			: g_{g}, c_{c}, end_{end} {
				skip_();
			}

			bool operator== (sentinel) const {
				return c_ == end_;
			}
			bool operator!= (sentinel) const {
				return c_ != end_;
			}
			bool operator!= (const iterator &other) const {
				return c_ != other.c_;
			}

			iterator &operator++ () {
				assert(c_);
				assert(c_ < end_);
				++c_;
				skip_();
				return *this;
			}

			const vertex &operator* () const {
				return *c_;
			}

		private:
			void skip_() {
				while(c_ != end_ && !g_->present_[*c_])
					++c_;
			}

			const graph *g_;
			const vertex *c_;
			const vertex *end_;
		};

		explicit neighbor_range(const graph *g, vertex v)
		: g_{g}, v_{v} { }

		iterator begin() const {
			const auto &adj = g_->adj_lists_[v_];
			return iterator{g_, adj.data(), adj.data() + adj.size()};
		}

		sentinel end() const {
			return sentinel{};
		}

	private:
		const graph *g_;
		vertex v_;
	};

	// Largest vertex ID + 1.
	// Useful if vertex IDs are used to index into arrays.
	unsigned int id_limit() const {
		return id_limit_;
	}

	// Number of present vertices.
	size_t num_vertices() const {
		return num_present_;
	}

	vertex add_vertex() {
		vertex v = id_limit_++;
		assert(present_.size() == v);
		assert(adj_lists_.size() == v);
		++num_present_;
		present_.push_back(1);
		adj_lists_.emplace_back();
		return v;
	}

	void add_edge(vertex u, vertex v) {
		assert(u < id_limit_);
		assert(v < id_limit_);
		assert(present_[u]);
		assert(present_[v]);
		auto &u_adj = adj_lists_[u];
		auto &v_adj = adj_lists_[v];
		assert(std::find(u_adj.begin(), u_adj.end(), v) == u_adj.end());
		assert(std::find(v_adj.begin(), v_adj.end(), u) == v_adj.end());
		adj_lists_[u].push_back(v);
		adj_lists_[v].push_back(u);
	}

	// Eliminate a vertex, i.e., temporarily remove it and all incident edges from the graph.
	void eliminate(vertex v) {
		assert(v < id_limit_);
		assert(present_[v]);
		assert(num_present_ > 0);
		present_[v] = 0;
		--num_present_;
	}

	// Undo a vertex elimination.
	void uneliminate(vertex v) {
		assert(v < id_limit_);
		assert(!present_[v]);
		present_[v] = 1;
		++num_present_;
	}

	vertex_range vertices() const {
		return vertex_range{this};
	}

	neighbor_range neighbors(vertex v) const {
		assert(v < id_limit_);
		assert(present_[v]);
		return neighbor_range{this, v};
	}

	// Compute the degree of a vertex.
	// Note that this is not O(1) because we need to check for eliminated vertices.
	int degree(vertex v) const {
		int d = 0;
		for(vertex w : neighbors(v)) {
			(void)w;
			++d;
		}
		return d;
	}

private:
	unsigned int id_limit_ = 0;
	unsigned int num_present_ = 0;
	std::vector<uint8_t> present_;
	std::vector<std::vector<vertex>> adj_lists_;
};
