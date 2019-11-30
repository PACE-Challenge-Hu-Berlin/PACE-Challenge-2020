#pragma once

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <functional> // For std::hash.
#include <vector>

#include "data-structures.hpp"


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
				while(v_ != g_->id_limit_ && !g_->exists_[v_])
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
			: g_{g}, c_{c}, end_{end} { }

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
				assert(c_ != end_);
				++c_;
				return *this;
			}

			const vertex &operator* () const {
				return *c_;
			}

		private:
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

	friend void swap(graph &g1, graph &g2) {
		using std::swap;
		swap(g1.id_limit_, g2.id_limit_);
		swap(g1.num_vertices_, g2.num_vertices_);
		swap(g1.exists_, g2.exists_);
		swap(g1.adj_lists_, g2.adj_lists_);
	}

	graph() = default;

	graph(const graph &other) = default;

	graph(graph &&other)
	: graph{} {
		swap(*this, other);
	}

	graph &operator= (graph other) {
		swap(*this, other);
		return *this;
	}

	void resize_ids(unsigned int n) {
		// TODO: support shrinking, but only if there are no vertices.
		assert(id_limit_ <= n);
		exists_.resize(n, false);
		adj_lists_.resize(n);
		id_limit_ = n;
	}

	// Largest vertex ID + 1.
	// Useful if vertex IDs are used to index into arrays.
	unsigned int id_limit() const {
		return id_limit_;
	}

	// Number of present vertices.
	size_t num_vertices() const {
		return num_vertices_;
	}

	vertex add_vertex(vertex v) {
		assert(v < id_limit_);
		assert(!exists_[v]);
		++num_vertices_;
		exists_[v] = 1;
		return v;
	}

	void add_edge(vertex u, vertex v) {
		assert(u < id_limit_);
		assert(v < id_limit_);
		assert(exists_[u]);
		assert(exists_[v]);
		auto &u_adj = adj_lists_[u];
		auto &v_adj = adj_lists_[v];
		assert(std::find(u_adj.begin(), u_adj.end(), v) == u_adj.end());
		assert(std::find(v_adj.begin(), v_adj.end(), u) == v_adj.end());
		adj_lists_[u].push_back(v);
		adj_lists_[v].push_back(u);
	}

	void clear() {
		std::fill(exists_.begin(), exists_.end(), 0);
		for(auto &adj : adj_lists_)
			adj.clear();
		num_vertices_ = 0;
	}

	vertex_range vertices() const {
		return vertex_range{this};
	}

	neighbor_range neighbors(vertex v) const {
		assert(v < id_limit_);
		assert(exists_[v]);
		return neighbor_range{this, v};
	}

	int degree(vertex v) const {
		return static_cast<int>(adj_lists_[v].size());
	}

	friend std::ostream& operator<<(std::ostream& out, const graph& g)
	{
		out << "c Graph data\n";
		unsigned long edge_counter = 0;
		for(const auto x : g.vertices())
		{
			edge_counter += g.degree(x);
		}
		assert(!(edge_counter & 2));
		out << g.num_vertices() << " " << edge_counter/2 << "\n";
		out << "c Adjacencies:\n";
		for(auto x : g.vertices())
		{
			for(auto y : g.neighbors(x))
			{
				if(x > y)continue;
				out << x << " " << y << "\n";
			}
		}
		out << std::endl;//to flush output
		return out;
	}
private:
	unsigned int id_limit_ = 0;
	unsigned int num_vertices_ = 0;
	std::vector<uint8_t> exists_;
	std::vector<std::vector<vertex>> adj_lists_;
};

// Lightweight reference to a set of vertices.
struct vertex_span {
	friend bool operator== (vertex_span vs1, vertex_span vs2) {
		if(vs1.size() != vs2.size())
			return false;
		for(size_t i = 0; i < vs1.size(); i++)
			if(vs1[i] != vs2[i])
				return false;
		return true;
	}

	friend bool operator!= (vertex_span vs1, vertex_span vs2) {
		return !(vs1 == vs2);
	}

	vertex_span()
	: begin_{nullptr}, end_{nullptr} { }

	vertex_span(const vertex *begin, const vertex *end)
	: begin_{begin}, end_{end} { }

	vertex_span(span<vertex> s)
	: begin_{s.begin()}, end_{s.end()} { }

	template<typename A>
	explicit vertex_span(const std::vector<vertex, A> &vec)
	: begin_{vec.data()}, end_{vec.data() + vec.size()} { }

	span<vertex> as_span() const {
		return span<vertex>{begin_, end_};
	}

	vertex operator[] (size_t i) const {
		return begin_[i];
	}

	const vertex *begin() const {
		return begin_;
	}

	const vertex *end() const {
		return end_;
	}

	vertex front() const {
		assert(begin_ != end_);
		return *begin_;
	}

	const vertex *data() const {
		return begin_;
	}

	bool empty() const {
		return !size();
	}

	size_t size() const {
		return end_ - begin_;
	}

private:
	const vertex *begin_;
	const vertex *end_;
};

namespace std {
	template<>
	struct hash<vertex_span> {
		// Possibly the worst hash function ever imagined.
		// TODO: Pick a reasonable hash function.
		size_t operator() (vertex_span vs) const {
			size_t hash = 0x12345678;
			for(auto p = vs.begin(); p != vs.end(); ++p)
				hash = 13 * hash + *p;
			return hash;
		}
	};
}
