#pragma once
#include "graph.hpp"
#include "connected-components.hpp"
#include "data-structures.hpp"
#include "directed-graph.hpp"

struct precedence_by_inclusion
{
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

	unsigned int num_components()
	{
		return scc_.num_components();
	}

	strongly_connected_components::vertex_range component(size_t i) const 
	{
		return scc_.component(i);
	}

	unsigned int component_id(vertex u)
	{
		return scc_.component_id(u);
	}

	size_t component_size(size_t i) const
	{
		return scc_.component_size(i);
	}

	directed_graph::neighbor_range component_neighbors(vertex v) const
       	{
		return scc_dag_.neighbors(v);
		}
		
private:
	unsigned int n;
	directed_graph adj_lists_;

	strongly_connected_components scc_;
	directed_graph scc_dag_;

	bool covers(const graph& g_, vertex u, vertex v)
	{
		static boolean_marker u_neighbors;
		u_neighbors.reset(n);
		for(auto w : g_.neighbors(u))
		{
			u_neighbors.mark(w);
		}

		for(auto w : g_.neighbors(v))
		{
			if(!u_neighbors.is_marked(w))
				return false;
		}
		return true;
	}

	void run_inclusion_scc()
	{
		scc_.compute(adj_lists_);
		scc_dag_ = std::move(directed_graph{scc_.num_components()});// [todo] also avoid reallocatoin

		for(unsigned int u = 0; u  < n; u++)
		{
			unsigned int idu = scc_.component_id(u);
			for(const auto& v : adj_lists_.neighbors(u))
			{
				unsigned int idv = scc_.component_id(v);
				if(idu == idv)continue;
				scc_dag_.add_edge(idu, idv);
			}

		}
		scc_dag_.remove_multiedges();
	}

};
