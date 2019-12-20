#pragma once

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <vector>

#include "data-structures.hpp"

struct directed_graph {
	// Range that iterates over all vertices of the directed_graph.
	struct vertex_range {
		struct sentinel { };

		struct iterator {
			explicit iterator(const directed_graph *g)
			: g_{g} {}

			bool operator!= (sentinel) const {
				return v_ != g_->num_vertices_;
			}

			iterator &operator++ () {
				assert(v_ < g_->num_vertices_);
				++v_;
				return *this;
			}

			vertex operator* () const {
				return v_;
			}

		private:

			const directed_graph *g_;
			vertex v_ = 0;
		};

		explicit vertex_range(const directed_graph *g)
		: g_{g} { }

		iterator begin() const {
			return iterator{g_};
		}

		sentinel end() const {
			return sentinel{};
		}

	private:
		const directed_graph *g_;
	};

	// Range that iterates over all edges neighbors of the parameter-vertex. 
	struct neighbor_range {
		struct sentinel { };

		struct iterator {
			explicit iterator(const directed_graph *g, const vertex *c, const vertex *end)
			: g_{g}, c_{c}, end_{end} {}

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
				return *this;
			}

			const vertex &operator* () const {
				return *c_;
			}

		private:

			const directed_graph *g_;
			const vertex *c_;
			const vertex *end_;
		};

		explicit neighbor_range(const directed_graph *g, vertex v)
		: g_{g}, v_{v} { }

		iterator begin() const {
			const auto &adj = g_->adj_lists_[v_];
			return iterator{g_, adj.data(), adj.data() + adj.size()};
		}

		sentinel end() const {
			return sentinel{};
		}

	private:
		const directed_graph *g_;
		vertex v_;
	};

	directed_graph(size_t n):num_vertices_(n), adj_lists_(n){}
	directed_graph() = default;

	void add_edge(vertex u, vertex v) {
		assert(u < num_vertices_);
		assert(v < num_vertices_);
		assert(u != v);
		adj_lists_[u].push_back(v);
	}

	size_t num_vertices() const {
		return num_vertices_;
	}

	vertex_range vertices() const {
		return vertex_range{this};
	}

	neighbor_range neighbors(vertex v) const {
		assert(v < num_vertices_);
		return neighbor_range{this, v};
	}

	int degree(vertex v) const
	{
		return adj_lists_[v].size();
	}

	void remove_multiedges()
	{
		for(auto &v : adj_lists_)
		{
			std::sort(v.begin(), v.end());
			v.erase(std::unique(v.begin(), v.end()),
					v.end());
		}
	}

	friend std::ostream& operator<<(std::ostream& out, const directed_graph& g)
	{
		out << "c Graph data\n";
		unsigned long edge_counter = 0;
		for(const auto x : g.vertices())
		{
			edge_counter += g.degree(x);
		}
		out << g.num_vertices() << " " << edge_counter << "\n";
		out << "c Adjacencies:\n";
		for(auto x : g.vertices())
		{
			for(auto y : g.neighbors(x))
			{
				out << x << " " << y << "\n";
			}
		}
		out << std::endl;//to flush output
		return out;
	}

private:
	unsigned int num_vertices_;
	std::vector<std::vector<vertex>> adj_lists_;
};
