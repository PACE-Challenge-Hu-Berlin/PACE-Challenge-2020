#include <iostream>
#include <tuple>

#include "connected-components.hpp"
#include "kernelization.hpp"
#include "utilities.hpp"
#include "subgraph.hpp"

void kernelization::compute(const graph &ig) {
	std::cerr << "computing kernel" << std::endl;
	coarse_profiling_timer timer;
	kernel_tree_.assign(ig.id_limit(), nil_vertex());
	maximal_depth_ = 0;
	q_.push(std::make_tuple(ig, 0, 0));

	while(!q_.empty()) {
		graph g = std::move(std::get<0>(q_.front()));
		int depth = std::get<1>(q_.front());
		vertex root = std::get<2>(q_.front());
		q_.pop();
		std::cerr << "considering subgraph of size " << g.num_vertices() << std::endl;

		if(decompose_(g, depth, root))
			continue;
		if(reduce_universal_(g, depth, root))
			continue;

		components_.push_back(std::make_tuple(std::move(g), depth, root));
	}
	auto elapsed = timer.elapsed();

	std::cerr << "kernalization finished" << std::endl;
	std::cerr << "    components: " << components_.size() << std::endl;
	std::cerr << "    universal vertices: " << stats_.num_universal << std::endl;
	std::cerr << "    kernelization time: " << print_time(elapsed) << std::endl;
}

void kernelization::compute_trivial(const graph &ig) {
	std::cerr << "computing kernel" << std::endl;
	coarse_profiling_timer timer;
	kernel_tree_.assign(ig.id_limit(), nil_vertex());
	maximal_depth_ = 0;
	q_.push(std::make_tuple(ig, 0, 0));

	while(!q_.empty()) {
		graph g = std::move(std::get<0>(q_.front()));
		int depth = std::get<1>(q_.front());
		vertex root = std::get<2>(q_.front());
		q_.pop();
		std::cerr << "considering subgraph of size " << g.num_vertices() << std::endl;

		if(decompose_(g, depth, root))
			continue;

		components_.push_back(std::make_tuple(std::move(g), depth, root));
	}
	auto elapsed = timer.elapsed();

	std::cerr << "kernalization finished" << std::endl;
	std::cerr << "    components: " << components_.size() << std::endl;
	std::cerr << "    kernelization time: " << print_time(elapsed) << std::endl;
}

bool kernelization::decompose_(graph &g, int depth, vertex root) {
	connected_components cc;
	cc.compute(g);

	if(cc.num_components() == 1)
		return false;

	std::cerr << "    decomposition yields " << cc.num_components()
			<< " components" << std::endl;

	for(size_t i = 0; i < cc.num_components(); ++i) {
		marker_.reset(g.id_limit());
		for(vertex v : cc.component(i))
			marker_.mark(v);

		induced_subgraph isg{g, marker_};
		graph sg;
		isg.materialize(sg);
		q_.push(std::make_tuple(std::move(sg), depth, root));
	}
	return true;
}

bool kernelization::reduce_universal_(graph &g, int depth, vertex root) {
	int count = 0;
	marker_.reset(g.id_limit());
	for(vertex v : g.vertices()) {
		if(static_cast<size_t>(g.degree(v)) + 1 == g.num_vertices()) {
			marker_.mark(v);
			kernel_tree_[v] = root;
			root = v;
			++depth;
			++count;
		}
	}

	maximal_depth_ = std::max(depth, maximal_depth_);

	if(!count)
		return false;

	std::cerr << "    elimination of " << count << " universal vertices" << std::endl;

	complement_induced_subgraph isg{g, marker_};
	graph sg;
	isg.materialize(sg); // TODO: add materialize_inline.
	q_.push(std::make_tuple(std::move(sg), depth, root));
	stats_.num_universal += count;
	return true;
}
