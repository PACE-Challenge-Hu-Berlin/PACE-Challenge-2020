#pragma once

#include <cstddef>
#include <vector>

#include "graph.hpp"

// TODO: implement Tamaki's bound on the "margin" of vertex sets.

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
		std::vector<std::pair<trie *, unsigned int>> stack;
	};

private:
	struct sieve_sentinel { };

	struct sieve_iterator {
		explicit sieve_iterator(trie *root)
		: cur_{root} { }

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
				if(!cur_)
					return *this;
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
		trie *cur_;
	};

	// Range that contains all vertex sets that are
	//     (1) disjoint w.r.t. the given membership predicate
	// and (2) have a representative (= smallest vertex) satisfying the given comparison.
	// Note: the iterators of this range are *single pass* only.
	// Also note that the range requires a (mutable) query state.
	// Space complexity: O(n^2) where n is the number of vertices.
	template<comparator Cmp, typename C>
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

		disjoint_range(block_sieve &sieve, query &query, vertex rv, C contains)
		: sieve_{&sieve}, query_{&query}, rv_{rv}, contains_{contains} {
			query_->stack.clear();
			query_->stack.push_back({sieve_->root_, 0});
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
			if(cur_->child)
				query_->stack.push_back({cur_->child, cur_->plen});
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
			auto range_disjoint = [&] (unsigned int i, unsigned int j) -> bool {
				for(unsigned int k = i; k < j; k++)
					if(contains_(cur_->vs[k]))
						return false;
				return true;
			};

			while(true) {
				while(cur_) {
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

					if(range_disjoint(pfx_, cur_->plen)) {
						if(range_disjoint(cur_->plen, cur_->vs.size()))
							return;

						if(cur_->child)
							query_->stack.push_back({cur_->child, cur_->plen});
					}
					cur_ = cur_->next;
				}

				if(query_->stack.empty())
					break;
				std::tie(cur_, pfx_) = query_->stack.back();
				query_->stack.pop_back();
			}
		}

		block_sieve *sieve_;
		query *query_;
		vertex rv_;
		C contains_;
		trie *cur_ = nullptr;
		unsigned int pfx_;
	};

public:
	size_t size() {
		return size_;
	}

	void insert(vertex_span vs, E element) {
		assert(vs.size());

		auto nd = new trie{vs, element};
		unsigned int nlen = nd->vs.size();

		if(!root_) {
			nd->plen = nlen;
			root_ = nd;
			++size_;
			return;
		}

		unsigned int pfx = 0;
		auto cur = root_;
		trie *previous = nullptr;
		while(true) {
			// Order sets by their representative (this is only important in the first level).
			if(!pfx && vs[0] < cur->vs[0]) {
				trie *parent = cur->parent;
				nd->plen = nlen;
				nd->parent = parent;
				nd->next = cur;

				if(previous) {
					assert(previous->next == cur);
					previous->next = nd;
				}else if(parent) {
					assert(parent->child == cur);
					parent->child = nd;
				}else{
					assert(root_ == cur);
					root_ = nd;
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
				cur->parent = nd;
				cur->next = nullptr;

				if(previous) {
					assert(previous->next == cur);
					previous->next = nd;
				}else if(parent) {
					assert(parent->child == cur);
					parent->child = nd;
				}else{
					assert(root_ == cur);
					root_ = nd;
				}
				++size_;
				return;
			}
		}
	}

	void clear() {
		trie *cur = root_;
		root_ = nullptr;
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

	sieve_iterator begin() {
		return sieve_iterator{root_};
	}

	sieve_sentinel end() {
		return sieve_sentinel{};
	}

	template<typename C>
	disjoint_range<comparator::less, C> disjoint_less(query &query, vertex rv, C contains) {
		return {*this, query, rv, contains};
	}

	template<typename C>
	disjoint_range<comparator::greater, C> disjoint_greater(query &query, vertex rv, C contains) {
		return {*this, query, rv, contains};
	}

private:
	trie *root_ = nullptr;
	size_t size_ = 0;
};

template<typename E>
using block_sieve_query = typename block_sieve<E>::query;
