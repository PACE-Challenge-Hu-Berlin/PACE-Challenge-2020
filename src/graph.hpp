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

struct induced_subgraph {
	// Range that iterates over all present vertices of the subgraph.
	struct vertex_range {
		struct sentinel { };

		struct iterator {
			explicit iterator(const boolean_marker *marker, graph::vertex_range::iterator c)
			: marker_{marker}, c_{c} {
				skip_();
			}

			bool operator!= (sentinel) const {
				return c_ != graph::vertex_range::sentinel{};
			}

			iterator &operator++ () {
				++c_;
				skip_();
				return *this;
			}

			vertex operator* () const {
				return *c_;
			}

		private:
			void skip_() {
				while(c_ != graph::vertex_range::sentinel{} && !marker_->is_marked(*c_))
					++c_;
			}

			const boolean_marker *marker_;
			graph::vertex_range::iterator c_;
		};

		explicit vertex_range(const graph *g, const boolean_marker *marker)
		: g_{g}, marker_{marker} { }

		iterator begin() const {
			return iterator{marker_, g_->vertices().begin()};
		}

		sentinel end() const {
			return sentinel{};
		}

	private:
		const graph *g_;
		const boolean_marker *marker_;
	};

	// Range that iterates over all edges incident to present vertices.
	struct neighbor_range {
		struct sentinel { };

		struct iterator {
			explicit iterator(const boolean_marker *marker, graph::neighbor_range::iterator c)
			: marker_{marker}, c_{c} {
				skip_();
			}

			bool operator== (sentinel) const {
				return c_ == graph::neighbor_range::sentinel{};
			}
			bool operator!= (sentinel) const {
				return c_ != graph::neighbor_range::sentinel{};
			}
			bool operator!= (const iterator &other) const {
				return c_ != other.c_;
			}

			iterator &operator++ () {
				++c_;
				skip_();
				return *this;
			}

			const vertex &operator* () const {
				return *c_;
			}

		private:
			void skip_() {
				while(c_ != graph::neighbor_range::sentinel{} && !marker_->is_marked(*c_))
					++c_;
			}

			const boolean_marker *marker_;
			graph::neighbor_range::iterator c_;
		};

		explicit neighbor_range(const graph *g, const boolean_marker *marker, vertex v)
		: g_{g}, marker_{marker}, v_{v} { }

		iterator begin() const {
			return iterator{marker_, g_->neighbors(v_).begin()};
		}

		sentinel end() const {
			return sentinel{};
		}

	private:
		const graph *g_;
		const boolean_marker *marker_;
		vertex v_;
	};

	induced_subgraph(const graph &g, const boolean_marker &marker)
	: g_{&g}, marker_{&marker} { }

	size_t id_limit() const {
		return g_->id_limit();
	}

	// Number of present vertices.
	// Note that this is not O(1) because we need to check for eliminated vertices.
	size_t num_vertices() const {
		assert(g_->id_limit() <= marker_->size());
		int n = 0;
		for(vertex v : vertices()) {
			(void)v;
			++n;
		}
		return n;
	}

	vertex_range vertices() const {
		assert(g_->id_limit() <= marker_->size());
		return vertex_range{g_, marker_};
	}

	neighbor_range neighbors(vertex v) const {
		assert(g_->id_limit() <= marker_->size());
		assert(v < marker_->size());
		assert(marker_->is_marked(v));
		return neighbor_range{g_, marker_, v};
	}

	// Compute the degree of a vertex.
	// Note that this is not O(1) because we need to check for eliminated vertices.
	int degree(vertex v) const {
		assert(g_->id_limit() <= marker_->size());
		int d = 0;
		for(vertex w : neighbors(v)) {
			(void)w;
			++d;
		}
		return d;
	}

	void materialize(graph &mg) const {
		assert(g_->id_limit() <= marker_->size());
		mg.clear();
		mg.resize_ids(g_->id_limit());
		for(vertex v : vertices())
			mg.add_vertex(v);
		for(vertex v : vertices())
			for(vertex w : neighbors(v))
				if(v <= w)
					mg.add_edge(v, w);
	}

private:
	const graph *g_;
	const boolean_marker *marker_;
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

	template<typename A>
	vertex_span(const std::vector<vertex, A> &vec)
	: begin_{vec.data()}, end_{vec.data() + vec.size()} { }

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
