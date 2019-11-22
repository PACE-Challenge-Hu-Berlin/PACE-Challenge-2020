#pragma once

#include <cstddef>
#include <vector>

#include "data-structures.hpp"
#include "graph.hpp"

template<typename E>
struct separated_sieve {
private:
	// TODO: improve cache efficiency by storing all siblings in a contiguous buffer.
	// TODO: store the first vertex of the prefix inside this struct to improve locality.
	struct trie {
		vertex_span vs;
		E element;
		trie *parent = nullptr;
		trie *next = nullptr;
		trie *child = nullptr;
		// Length of prefix that distinguishes this element from siblings.
		unsigned int plen = 0;
		std::vector<vertex> new_neighbors;
		unsigned int self_overlap;
	};

	enum class comparator {
		none,
		less,
		greater
	};

public:
	// Query state, required for list_less() / list_greater().
	// Only one operation can use an object of this class at a time.
	struct query {
		std::vector<std::tuple<trie *, unsigned int, int>> stack;
	};

private:
	struct sieve_sentinel { };

	struct sieve_iterator {
		explicit sieve_iterator(separated_sieve *sieve)
		: sieve_{sieve} {
			skip_roots_();
		}

		bool operator== (sieve_sentinel) const {
			return !cur_;
		}
		bool operator!= (sieve_sentinel) const {
			return cur_;
		}

		sieve_iterator &operator++ () {
			assert(cur_);

			// Go downwards if possible.
			if(cur_->child) {
				cur_ = cur_->child;
				return *this;
			}

			// Go sidewards if we are at the lowest level already.
			if(cur_->next) {
				cur_ = cur_->next;
				return *this;
			}

			// Go upwards (multiple times) and at least once sidewards.
			while(true) {
				cur_ = cur_->parent;
				if(!cur_) {
					++r_;
					skip_roots_();
					return *this;
				}
				if(cur_->next) {
					cur_ = cur_->next;
					return *this;
				}
			}
		}

		E &operator* () const {
			assert(cur_);
			return cur_->element;
		}

	private:
		void skip_roots_() {
			assert(!cur_);
			while(r_ < sieve_->roots_.size()) {
				if(sieve_->roots_[r_]) {
					cur_ = sieve_->roots_[r_];
					break;
				}
				r_++;
			}
		}

		separated_sieve *sieve_;
		trie *cur_ = nullptr;
		size_t r_ = 0;
	};

public:
	size_t size() {
		return size_;
	}

	void insert(const graph &g, vertex_span vs, E element) {
		assert(vs.size());

		auto nd = new trie{vs, element};
		unsigned int nlen = nd->vs.size();

		auto r = log_ceil2int(vs.size());
		if(r >= roots_.size())
			roots_.resize(r + 1, nullptr);

		if(!roots_[r]) {
			nd->plen = nlen;
			recompute_neighbors_(g, nd, 0);
			roots_[r] = nd;
			++size_;
			return;
		}

		unsigned int pfx = 0;
		auto cur = roots_[r];
		trie *previous = nullptr;
		while(true) {
			// Order sets by their representative (this is only important in the first level).
			if(!pfx && vs[0] < cur->vs[0]) {
				trie *parent = cur->parent;
				nd->plen = nlen;
				nd->parent = parent;
				nd->next = cur;
				recompute_neighbors_(g, nd, pfx);

				if(previous) {
					assert(previous->next == cur);
					previous->next = nd;
				}else if(parent) {
					assert(parent->child == cur);
					parent->child = nd;
				}else{
					assert(roots_[r] == cur);
					roots_[r] = nd;
				}
				++size_;
				return;
			}

			unsigned int i;
			for(i = pfx; i < std::min(cur->plen, nlen); i++)
				if(cur->vs[i] != vs[i])
					break;

			if(i == pfx) {
				// No match. Proceed sideways.
				if(!cur->next) {
					nd->plen = nlen;
					nd->parent = cur->parent;
					recompute_neighbors_(g, nd, pfx);
					cur->next = nd;
					++size_;
					return;
				}

				previous = cur;
				cur = cur->next;
			}else if(i == cur->plen) {
				// Full match. Proceed downwards.
				if(!cur->child) {
					nd->plen = nlen;
					nd->parent = cur;
					recompute_neighbors_(g, nd, cur->plen);
					cur->child = nd;
					++size_;
					return;
				}

				previous = nullptr;
				pfx = cur->plen;
				cur = cur->child;
			}else{
				assert(i < cur->plen);

				// Partial match. Exchange the current node.
				assert(vs[0] == cur->vs[0]);
				trie *parent = cur->parent;
				nd->plen = i;
				nd->parent = parent;
				nd->child = cur;
				nd->next = cur->next;
				recompute_neighbors_(g, nd, pfx);
				cur->parent = nd;
				cur->next = nullptr;
				recompute_neighbors_(g, cur, i);

				if(previous) {
					assert(previous->next == cur);
					previous->next = nd;
				}else if(parent) {
					assert(parent->child == cur);
					parent->child = nd;
				}else{
					assert(roots_[r] == cur);
					roots_[r] = nd;
				}
				++size_;
				return;
			}
		}
	}

	void clear() {
		for(size_t r = 0; r < roots_.size(); ++r) {
			trie *cur = roots_[r];
			roots_[r] = nullptr;
			while(cur) {
				if(cur->child) {
					auto child = cur->child;
					cur->child = nullptr;
					cur = child;
				}else{
					auto parent = cur->parent;
					auto next = cur->next;
					delete cur;

					if(next) {
						cur = next;
					}else{
						cur = parent;
					}
				}
			}
		}
		size_ = 0;
	}

	sieve_iterator begin() {
		return sieve_iterator{this};
	}

	sieve_sentinel end() {
		return sieve_sentinel{};
	}

	template<typename C, typename N, typename F>
	void list_predecessors(const graph &g, query &query,
			int m, vertex rv, C contains, N neighbor, F functor) {
		return do_list_<comparator::less>(g, query, m, rv, contains, neighbor, functor);
	}

	template<typename C, typename N, typename F>
	void list_successors(const graph &g, query &query,
			int m, vertex rv, C contains, N neighbor, F functor) {
		return do_list_<comparator::greater>(g, query, m, rv, contains, neighbor, functor);
	}

private:
	void recompute_neighbors_(const graph &g, trie *nd, unsigned int pfx) {
		assert(!nd->parent || nd->parent->plen == pfx);

		marker_.reset(g.id_limit());
		for(unsigned int k = 0; k < pfx; k++)
			marker_.mark(nd->vs[k]);
		for(unsigned int k = 0; k < pfx; k++)
			for(vertex w : g.neighbors(nd->vs[k]))
				marker_.mark(w);

		nd->new_neighbors.clear();
		nd->self_overlap = 0;
		for(unsigned int k = pfx; k < nd->plen; k++) {
			if(marker_.is_marked(nd->vs[k]))
				continue;
			marker_.mark(nd->vs[k]);
			nd->self_overlap++; // We assume that vertices are not duplicated.
		}
		for(unsigned int k = pfx; k < nd->plen; k++)
			for(vertex w : g.neighbors(nd->vs[k])) {
				if(marker_.is_marked(w))
					continue;
				marker_.mark(w);
				nd->new_neighbors.push_back(w);
			}
	}

	// Lists all vertex sets that are
	//     (1) disjoint w.r.t. the given membership predicate,
	//     (2) separated w.r.t. the given neighborhood predicate,
	// and (3) have a representative (= smallest vertex) satisfying the given comparison.
	// Space complexity: O(n^2) where n is the number of vertices.
	template<comparator Cmp, typename C, typename N, typename F>
	void do_list_(const graph &g, query &query,
			int m, vertex rv, C contains, N neighbor, F functor) {
		auto is_disjoint_separated = [&] (trie *nd, unsigned int i, unsigned int j) -> bool {
			for(unsigned int k = i; k < j; k++)
				if(contains(nd->vs[k]) || neighbor(nd->vs[k]))
					return false;
			return true;
		};

		query.stack.clear();
		for(size_t r = 0; r < roots_.size(); ++r) {
			trie *cur = roots_[r];
			unsigned int pfx = 0;
			int pb = 0;
			while(true) {
				while(cur) {
					assert(!cur->parent || cur->parent->plen == pfx);

					if(Cmp == comparator::less) {
						if(cur->vs[0] >= rv) {
							cur = nullptr;
							break;
						}
					}else{
						assert(Cmp == comparator::greater);
						if(cur->vs[0] <= rv) {
							cur = cur->next;
							continue;
						}
					}

					if(is_disjoint_separated(cur, pfx, cur->plen)) {
						assert(cur->self_overlap <= cur->plen - pfx);
						int b = pb - cur->self_overlap + (cur->plen - pfx);
						for(vertex w : cur->new_neighbors) {
							if(!neighbor(w))
								b++;
						}
						if(b <= m + (1 << r)) {
							if(is_disjoint_separated(cur, cur->plen, cur->vs.size()))
								functor(cur->element);

							if(cur->child)
								query.stack.push_back({cur->child, cur->plen, b});
						}
					}
					cur = cur->next;
				}

				if(query.stack.empty())
					break;
				std::tie(cur, pfx, pb) = query.stack.back();
				query.stack.pop_back();
			}
		}
	}

	std::vector<trie *> roots_;
	size_t size_ = 0;
	boolean_marker marker_;
};

template<typename E>
using separated_sieve_query = typename separated_sieve<E>::query;
