#pragma once

#include <cassert>
#include <cstdint>

#include "data-structures.hpp"
#include "graph.hpp"

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

struct complement_induced_subgraph {
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
				while(c_ != graph::vertex_range::sentinel{} && marker_->is_marked(*c_))
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
				while(c_ != graph::neighbor_range::sentinel{} && marker_->is_marked(*c_))
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

	complement_induced_subgraph(const graph &g, const boolean_marker &marker)
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
		assert(!marker_->is_marked(v));
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
