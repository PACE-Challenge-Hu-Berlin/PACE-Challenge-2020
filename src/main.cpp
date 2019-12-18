#include <cstring>
#include <fstream>

#include "io.hpp"
#include "naive-branching-solver.hpp"
#include "simple-pid-solver.hpp"
#include "utilities.hpp"

static const char *usage_text =
	"Usage: td [OPTIONS] INSTANCE\n"
	"Possible OPTIONS are:\n"
	"    --solver SOLVER\n"
	"        Select a solver algorithm {naive-branching, simple-pid}.\n"
	"    --no-precedence\n"
	"        Prune simple-pid to choose separators more carefully.\n";

enum class solver_algorithm {
	none,
	naive_branching,
	simple_pid,
};

bool no_precedence_ = false;

struct options {
	const char *instance = nullptr;
	solver_algorithm solver = solver_algorithm::naive_branching;
};

void parse_options(char **p, options &opts) {
	auto error = [] (const char *text) {
		std::cerr << usage_text << "Usage error: " << text << std::endl;
		exit(2);
	};

	// Argument for unary options.
	const char *arg;

	auto handle_nullary_option = [&] (const char *name) -> bool {
		assert(*p);
		if(std::strcmp(*p, name))
			return false;
		++p;
		return true;
	};

	auto handle_unary_option = [&] (const char *name) -> bool {
		assert(*p);
		if(std::strcmp(*p, name))
			return false;
		++p;
		if(!(*p))
			error("expected argument for unary option");
		arg = *p;
		++p;
		return true;
	};

	// TODO: Here, we can parse command line arguments.
	while(*p && !std::strncmp(*p, "--", 2)) {
		if(handle_nullary_option("--profile")) {
			global_profiling_flag = true;
		}else if(handle_unary_option("--solver")) {
			if(!std::strcmp(arg, "naive-branching")) {
				opts.solver = solver_algorithm::naive_branching;
			}else if(!std::strcmp(arg, "simple-pid")) {
				opts.solver = solver_algorithm::simple_pid;
			}else{
				error("unknown solver algorithm");
			}
		} else if(handle_nullary_option("--no-precedence")){
			no_precedence_ = true;
		} else{
			error("unknown command line option");
		}
	}

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

	int solution;
	profiling_duration solve_duration;
	if(opts.solver == solver_algorithm::naive_branching) {
		naive_branching_solver solver{g};
		coarse_profiling_timer solve_timer;
		solution = solver.compute_treedepth();
		solve_duration = solve_timer.elapsed();
	}else{
		assert(opts.solver == solver_algorithm::simple_pid);
		simple_pid_solver solver{g, no_precedence_};
		coarse_profiling_timer solve_timer;
		solution = solver.compute_treedepth();
		solve_duration = solve_timer.elapsed();
	}
	std::cout << "s " << solution << std::endl;
	std::cerr << "time: " << print_time(solve_duration) << std::endl;
}
