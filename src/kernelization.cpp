#include <iostream>

#include "connected-components.hpp"
#include "kernelization.hpp"
#include "utilities.hpp"
#include "subgraph.hpp"

void kernelization::compute(const graph &ig) {
	std::cerr << "computing kernel" << std::endl;
	coarse_profiling_timer timer;
	q_.push(ig);

	while(!q_.empty()) {
		graph &g = q_.back();

		if(decompose_(g)) {
			q_.pop();
			continue;
		}

		if(reduce_universal_(g))
			continue;

		components_.push_back(std::move(g));
		q_.pop();
	}
	auto elapsed = timer.elapsed();

	std::cerr << "    components: " << components_.size() << std::endl;
	std::cerr << "    universal vertices: " << stats_.num_universal << std::endl;
	std::cerr << "    kernelization time: " << print_time(elapsed) << std::endl;
}

bool kernelization::decompose_(graph &g) {
	connected_components cc;
	cc.compute(g);

	if(cc.num_components() == 1)
		return false;

	for(size_t i = 0; i < cc.num_components(); ++i) {
		marker_.reset(g.id_limit());
		for(vertex v : cc.component(i))
			marker_.mark(v);

		induced_subgraph isg{g, marker_};
		graph sg;
		isg.materialize(sg);
		q_.push(std::move(sg));
	}
	return true;
}

bool kernelization::reduce_universal_(graph &g) {
	int count = 0;
	marker_.reset(g.id_limit());
	for(vertex v : g.vertices()) {
		if(static_cast<size_t>(g.degree(v)) + 1 == g.num_vertices()) {
			marker_.mark(v);
			++count;
		}
	}

	if(!count)
		return false;

	complement_induced_subgraph isg{g, marker_};
	graph sg;
	isg.materialize(sg); // TODO: add materialize_inline.
	g = sg;
	stats_.num_universal += count;
	return true;
}
