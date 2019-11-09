#pragma once

#include <iostream>
#include <sstream>
#include <string>

#include "graph.hpp"

// AvdG: Parsing is always tedious.
// The following parser works but it is not very robust against syntax errors.
void read_pace_graph(std::istream &ifs, graph &g) {
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
	while(std::getline(ifs, l)) {
		if(is_whitespace())
			continue;
		std::istringstream iss{l};
		if(!is_prefixed())
			throw std::runtime_error("expected only prefixed lines before 'p'");
		iss >> w;
		if(w == "c") {
			continue;
		}else if(w == "p") {
			int n, m;
			iss >> w;
			iss >> n;
			iss >> m;
			// TODO: Make sure that values are parsed correctly.
			// TODO: Make sure that line ends here (except for whitespace).

			for(int i = 0; i <= n; i++)
				g.add_vertex();
			g.eliminate(0); // Do not use vertex ID zero.
			found_p = true;
			break;
		}else{
			throw std::runtime_error("expected only 'c' prefixes before 'p'");
		}
	}

	if(!found_p)
		throw std::runtime_error("found no 'p' line");

	// Parse all edges.
	while(std::getline(ifs, l)) {
		if(is_whitespace())
			continue;
		std::istringstream iss{l};
		if(!is_prefixed()) {
			int u, v;
			iss >> u;
			iss >> v;
			// TODO: Make sure that values are parsed correctly.
			// TODO: Make sure that line ends here (except for whitespace).

			g.add_edge(u, v);
		}else{
			iss >> w;
			if(w == "c") {
				continue;
			}else{
				throw std::runtime_error("expected only 'c' prefixes after 'p'");
			}
		}
	}
}
