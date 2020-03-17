#pragma once

#include <iostream>
#include <set>
#include <sstream>
#include <string>

#include "graph.hpp"

// AvdG: Parsing is always tedious.
// The following parser works but it is not very robust against syntax errors.
bool read_pace_graph(std::istream &ifs, graph &g) {
	std::string l; // Current line.
	std::string w; // Current word (for prefixed lines).

	// Checks if the line consists entirely of whitespace.
	auto is_whitespace = [&] () -> bool {
		for(auto c : l)
			if(!isspace(c))
				return false;
		return true;
	};

	// Checks if the line has a whitespace-word-whitespace prefix.
	auto is_prefixed = [&] () -> bool {
		const char *s = l.data();
		while(isspace(*s))
			++s;
		if(!isalpha(*s))
			return false;
		while(isalpha(*s))
			++s;
		return isspace(*s);
	};

	// Parse the initial "p"-prefixed line.
	bool found_p = false;
	int n, m;
	while(std::getline(ifs, l)) {
		if(is_whitespace())
			continue;
		std::istringstream iss{l};
		if(!is_prefixed()) {
			std::cerr << "expected only prefixed lines before 'p'" << std::endl;
			return false;
		}
		iss >> w;
		if(w == "c") {
			continue;
		}else if(w == "p") {
			iss >> w;
			iss >> n;
			iss >> m;
			if(!iss) {
				std::cerr << "i/o error" << std::endl;
				return false;
			}
			// TODO: Make sure that values are parsed correctly.
			// TODO: Make sure that line ends here (except for whitespace).

			g.resize_ids(n + 1);
			for(int i = 1; i <= n; i++)
				g.add_vertex(i);
			found_p = true;
			break;
		}else{
			std::cerr << "expected only 'c' prefixes before 'p'" << std::endl;
			return false;
		}
	}

	if(!found_p) {
		std::cerr << "found no 'p' line" << std::endl;
		return false;
	}

	// Parse all edges.
	std::set<std::pair<int, int>> edge_set;
	while(std::getline(ifs, l)) {
		if(is_whitespace())
			continue;
		std::istringstream iss{l};
		if(!is_prefixed()) {
			int u, v;
			iss >> u;
			iss >> v;
			if(!iss) {
				std::cerr << "i/o error" << std::endl;
				return false;
			}
			// TODO: Make sure that values are parsed correctly.
			// TODO: Make sure that line ends here (except for whitespace).

			if(u <= 0 || u > n) {
				std::cerr << "vertex " << u << " out of bounds" << std::endl;
				return false;
			}
			if(v <= 0 || v > n) {
				std::cerr << "vertex " << v << " out of bounds" << std::endl;
				return false;
			}

			if(u == v) {
				std::cerr << "ignoring self loop at " << u << std::endl;
				continue;
			}
			if(u > v)
				std::swap(u, v);

			bool unique;
			std::tie(std::ignore, unique) = edge_set.insert({u, v});
			if(!unique) {
				std::cerr << "ignoring duplicate edge {" << u << ", " << v << "}" << std::endl;
				continue;
			}

			g.add_edge(u, v);
		}else{
			iss >> w;
			if(w == "c") {
				continue;
			}else{
				std::cerr << "expected only 'c' prefixes after 'p'" << std::endl;
				return false;
			}
		}
	}

	return true;
}
