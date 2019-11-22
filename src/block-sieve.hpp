#pragma once

#include <cstddef>
#include <vector>

#include "data-structures.hpp"
#include "graph.hpp"

template<typename E>
struct block_sieve {
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
	// Query state, required for disjoint_less().
	// Only one operation can use an object of this class at a time.
	struct query {
		std::vector<std::tuple<size_t, trie *, unsigned int, int>> stack;
	};

private:
	struct sieve_sentinel { };

	struct sieve_iterator {
		explicit sieve_iterator(block_sieve *sieve)
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

		block_sieve *sieve_;
		trie *cur_ = nullptr;
		size_t r_ = 0;
	};

	// Range that contains all vertex sets that are
	//     (1) disjoint w.r.t. the given membership predicate
	// and (2) have a representative (= smallest vertex) satisfying the given comparison.
	// Note: the iterators of this range are *single pass* only.
	// Also note that the range requires a (mutable) query state.
	// Space complexity: O(n^2) where n is the number of vertices.
	template<comparator Cmp, typename C, typename N>
	struct disjoint_range {
		struct sentinel { };

		struct iterator {
			explicit iterator(disjoint_range *range)
			: range_{range} { }

			bool operator== (sentinel) const {
				return range_->done_();
			}
			bool operator!= (sentinel) const {
				return !range_->done_();
			}

			iterator &operator++ () {
				range_->advance_();
				return *this;
			}

			E &operator* () const {
				return range_->get_();
			}

		private:
			disjoint_range *range_;
		};

		disjoint_range(block_sieve &sieve, const graph &g, query &query,
				int m, vertex rv, C contains, N neighbor)
		: sieve_{&sieve}, g_{&g}, query_{&query}, m_{m}, rv_{rv},
				contains_{contains}, neighbor_{neighbor} {
			query_->stack.clear();
			for(size_t r = 0; r < sieve_->roots_.size(); ++r) {
				if(!sieve_->roots_[r])
					continue;
				query_->stack.push_back({r, sieve_->roots_[r], 0, 0});
			}
			skip_();
		}

		disjoint_range(const disjoint_range &) = delete;

		disjoint_range &operator= (const disjoint_range &) = delete;

		iterator begin() {
			return iterator{this};
		}

		sentinel end() {
			return sentinel{};
		}

	private:
		// The following functions implement the iterators.

		bool done_() {
			return !cur_ && query_->stack.empty();
		}

		void advance_() {
			assert(cur_);
			assert(!cur_->parent || cur_->parent->plen == pfx_);
			assert(cur_->self_overlap <= cur_->plen - pfx_);
			int b = pb_ - cur_->self_overlap + (cur_->plen - pfx_);
			for(vertex w : cur_->new_neighbors) {
				if(!neighbor_(w))
					b++;
			}
			if(cur_->child)
				query_->stack.push_back({r_, cur_->child, cur_->plen, b});
			cur_ = cur_->next;
			skip_();
		}

		E &get_() {
			assert(cur_);
			return cur_->element;
		}

	private:
		// Skips to the next element that satisfies the query.
		void skip_() {
			auto is_disjoint_separated = [&] (unsigned int i, unsigned int j) -> bool {
				for(unsigned int k = i; k < j; k++)
					if(contains_(cur_->vs[k]) || neighbor_(cur_->vs[k]))
						return false;
				return true;
			};

			while(true) {
				while(cur_) {
					assert(!cur_->parent || cur_->parent->plen == pfx_);
					if(Cmp == comparator::less) {
						if(cur_->vs[0] >= rv_) {
							cur_ = nullptr;
							break;
						}
					}else{
						assert(Cmp == comparator::greater);
						if(cur_->vs[0] <= rv_) {
							cur_ = cur_->next;
							continue;
						}
					}

					if(is_disjoint_separated(pfx_, cur_->plen)) {
						assert(cur_->self_overlap <= cur_->plen - pfx_);
						int b = pb_ - cur_->self_overlap + (cur_->plen - pfx_);
						for(vertex w : cur_->new_neighbors) {
							if(!neighbor_(w))
								b++;
						}
						if(b <= m_ + (1 << r_)) {
							if(is_disjoint_separated(cur_->plen, cur_->vs.size()))
								return;

							if(cur_->child)
								query_->stack.push_back({r_, cur_->child, cur_->plen, b});
						}
					}
					cur_ = cur_->next;
				}

				if(query_->stack.empty())
					break;
				std::tie(r_, cur_, pfx_, pb_) = query_->stack.back();
				query_->stack.pop_back();
			}
		}

		block_sieve *sieve_;
		const graph *g_;
		query *query_;
		int m_;
		vertex rv_;
		C contains_;
		N neighbor_;
		size_t r_;
		trie *cur_ = nullptr;
		unsigned int pfx_;
		int pb_;
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

	template<typename C, typename N>
	disjoint_range<comparator::less, C, N> disjoint_less(const graph &g, query &query,
			int m, vertex rv, C contains, N neighbor) {
		return {*this, g, query, m, rv, contains, neighbor};
	}

	template<typename C, typename N>
	disjoint_range<comparator::greater, C, N> disjoint_greater(const graph &g, query &query,
			int m, vertex rv, C contains, N neighbor) {
		return {*this, g, query, m, rv, contains, neighbor};
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

	std::vector<trie *> roots_;
	size_t size_ = 0;
	boolean_marker marker_;
};

template<typename E>
using block_sieve_query = typename block_sieve<E>::query;
