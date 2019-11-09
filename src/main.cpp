#include <fstream>

#include "io.hpp"
#include "naive-branching-solver.hpp"

static const char *usage_text =
	"Usage: td [OPTIONS] INSTANCE\n"
	"Possible OPTIONS are:\n"
	"    None\n";

struct options {
	const char *instance = nullptr;
};

void parse_options(char **p, options &opts) {
	auto error = [] (const char *text) {
		std::cerr << "Usage error: " << text << "\n" << usage_text << std::flush;
		exit(2);
	};

	// TODO: Here, we can parse command line arguments.

	if(!(*p))
		error("not enough positional arguments");
	opts.instance = *(p++);
	if(*p)
		error("too many positional arguments");
}

int main(int argc, char *argv[]) {
	options opts;
	parse_options(argv + 1, opts);

	graph g;
	{
		std::ifstream ifs{opts.instance};
		if(!ifs)
			throw std::runtime_error("could not open instance file");
		read_pace_graph(ifs, g);
	}
	std::cerr << "graph has " << g.num_vertices() << " vertices" << std::endl;

	naive_branching_solver solver{g};
	int solution = solver.compute_treedepth();
	std::cout << "s " << solution << std::endl;
}
