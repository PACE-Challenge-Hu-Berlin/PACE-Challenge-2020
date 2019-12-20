#pragma once
#include "graph.hpp"
#include "connected-components.hpp"
#include "data-structures.hpp"
#include "directed-graph.hpp"

struct precedence_by_inclusion
{

	//a mask over graph_neighbors that returns vertices that are not in the same scc
	struct dependent_vertex_range
	{
		struct sentinel{};
		struct iterator{

			using graph_sentinel = directed_graph::neighbor_range::sentinel;
			explicit iterator(const directed_graph& g_, const strongly_connected_components& scc_, vertex v):
				g_(g_), scc_(scc_), v(v), git_(g_.neighbors(v).begin())
			{
				skip();
			}

			bool operator!=(sentinel)
			{
				return git_ != graph_sentinel{};
			}

			iterator& operator++()
			{
				++git_;
				skip();
				return *this;
			}

			const vertex& operator* ()
			{
				return *git_;
			}

		private:
			const directed_graph& g_;
			const strongly_connected_components& scc_;
			vertex v;
			directed_graph::neighbor_range::iterator git_;

			void skip()
			{
				while(git_ != graph_sentinel{} &&
					       	scc_.component_id(*git_) == scc_.component_id(v))
				{
					++git_;
				}
			}

		};

		explicit dependent_vertex_range(const directed_graph& g_,
				const strongly_connected_components& scc_, vertex v):
				g_(g_), scc_(scc_), v(v){}

		iterator begin()
		{
			return iterator{g_, scc_, v};
		}
		sentinel end()
		{
			return sentinel{};
		}

		const directed_graph& g_;
		const strongly_connected_components& scc_;
		vertex v;

		iterator begin() const {
			return iterator{g_, scc_, v};
		}
		sentinel end() const
		{
			return sentinel{};
		}
	};

	precedence_by_inclusion() = default;

	void compute(const graph& g_) 
	{
		n = g_.id_limit();
		adj_lists_ = std::move(directed_graph(n));// [todo] avoid reallocation? possibly a reset function
		for(auto u : g_.vertices())
		{
			for(auto v : g_.neighbors(u))
			{
				if(covers(g_, u, v))
				{
					adj_lists_.add_edge(u, v);
				}
			}
		}
		adj_lists_.remove_multiedges();

		run_inclusion_scc();
	}
	void compute_trivial(const graph& g_) 
	{
		n = g_.id_limit();
		adj_lists_ = std::move(directed_graph(n));// [todo] avoid reallocation? possibly a reset function

		run_inclusion_scc();
	}

	unsigned int num_components() const
	{
		return scc_.num_components();
	}

	strongly_connected_components::vertex_range component(size_t i) const 
	{
		return scc_.component(i);
	}

	unsigned int component_id(vertex u) const
	{
		return scc_.component_id(u);
	}

	size_t component_size(size_t i) const
	{
		return scc_.component_size(i);
	}

	bool is_minimal(size_t i) const
	{
		return scc_minimal_[i];
	}

	bool is_maximal(size_t i) const
	{
		return scc_dag_.degree(i) == 0;
	}

	directed_graph::neighbor_range component_successors(vertex v) const
       	{
		return scc_dag_.neighbors(v);
	}

	dependent_vertex_range vertex_successors(vertex v) const
	{
		return dependent_vertex_range(adj_lists_, scc_, v);
	}

	friend std::ostream& operator << (std::ostream& out, const precedence_by_inclusion& prec)
	{
		out << "num of components: " << prec.num_components() << std::endl;
		for(unsigned int i = 0; i < prec.num_components(); i++)
		{
			for(auto v : prec.component(i))
			{
				out << v << " ";
			}
			out << std::endl;
		}
		out << "component ids" << std::endl;
		for(unsigned int i = 0; i < prec.n; i++)
		{
			out << prec.component_id(i) << std::endl;
		}
		out << "covering graph" << std::endl;
		out << prec.adj_lists_ << std::endl;
		out << "scc dag" << std::endl;
		out << prec.scc_dag_ << std::endl;
		return out;
	}
		
private:
	unsigned int n;
	directed_graph adj_lists_;

	strongly_connected_components scc_;
	directed_graph scc_dag_;
	std::vector<bool> scc_minimal_;

	bool covers(const graph& g_, vertex u, vertex v)
	{
		static boolean_marker u_neighbors;
		u_neighbors.reset(n);
		for(auto w : g_.neighbors(u))
		{
			assert(w < n);
			u_neighbors.mark(w);
		}

		for(auto w : g_.neighbors(v))
		{
			assert(w < n);
			if(w == u) 
			{
				continue;
			}
			if(!u_neighbors.is_marked(w))
				return false;
		}
		return true;
	}

	void run_inclusion_scc()
	{
		scc_.compute(adj_lists_);
		scc_dag_ = std::move(directed_graph{scc_.num_components()});// [todo] also avoid reallocatoin
		scc_minimal_.assign(scc_.num_components(), true);

		for(unsigned int u = 0; u  < n; u++)
		{
			unsigned int idu = scc_.component_id(u);
			for(const auto& v : adj_lists_.neighbors(u))
			{
				unsigned int idv = scc_.component_id(v);
				if(idu == idv)continue;
				scc_dag_.add_edge(idu, idv);
				scc_minimal_[idv] = false;
			}

		}
		scc_dag_.remove_multiedges();
	}

};
