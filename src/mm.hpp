#pragma once

#include <cstddef>
#include <cstdint>
#include <deque>

// Memory area that is always freed entirely at once.
struct memory_arena {
	memory_arena();

	memory_arena(const memory_arena &) = delete;

	memory_arena &operator= (const memory_arena &) = delete;

	~memory_arena();

	void *allocate(size_t align, size_t size) {
		size_t offset = (frontier_ + align - 1) & ~(align - 1);
		size_t n_frontier = offset + size;
		if(limit_ < n_frontier)
			extend(n_frontier);
		frontier_ = n_frontier;
		++live_;
		++num_allocations_;
		return base_ + offset;
	}

	void deallocate() {
		--live_;
	}

	void extend(size_t threshold);

	void reset();

	size_t used_space() {
		return frontier_;
	}

	int64_t num_allocations() {
		return num_allocations_;
	}

private:
	size_t frontier_ = 0;
	size_t limit_ = 0;
	char *base_ = nullptr;
	// For debugging: Make sure there are no live objects when we free the arena.
	int64_t live_ = 0;
	// For debugging: Allow users to make sure the number of allocations is as expected.
	int64_t num_allocations_ = 0;
};

template <class T>
struct arena_allocator {
	using value_type = T;

	using propagate_on_container_copy_assignment = std::true_type;
	using propagate_on_container_move_assignment = std::true_type;
	using propagate_on_container_swap = std::true_type;

	arena_allocator(memory_arena &arena) noexcept
	: arena_{&arena} { }

	template <class U>
	arena_allocator(arena_allocator<U> const &other) noexcept
	: arena_{other.arena_} { }

	T *allocate(size_t n) {
		return reinterpret_cast<T *>(arena_->allocate(alignof(T), sizeof(T) * n));
	}

	void deallocate(T *, size_t) noexcept {
		arena_->deallocate();
	}

private:
	memory_arena *arena_;
};

template <class T, class U>
bool operator== (arena_allocator<T> const &lhs, arena_allocator<U> const &rhs) noexcept {
	return lhs.arena_ == rhs.arena_;
}

template <class T, class U>
bool operator!= (arena_allocator<T> const& lhs, arena_allocator<U> const& rhs) noexcept {
	return !(lhs == rhs);
}

// --------------------------------------------------------------------------------------

// Memory that is freed in the same order in which it is allocated (FIFO).
// Usage: all allocations done before a call to seal() are freed by the corresponding
// call to reclaim().
struct queue_memory {
private:
	struct slice {
		char *base;
		size_t frontier_limit = SIZE_MAX;
	};

	struct chunk {
		size_t limit;
	};

	// Minimal amount of memory that must be available in a slice
	// (if this is not available, the next chunk is put into a new slice).
	static constexpr size_t slice_threshold = 1024 * 1024 * 1024;

	// Amount of memory before a commit/uncommit is triggered.
	static constexpr size_t commit_threshold = 64 * 1024 * 1024;

public:
	queue_memory();

	queue_memory(const queue_memory &) = delete;

	queue_memory &operator= (const queue_memory &) = delete;

	~queue_memory();

	void *allocate(size_t align, size_t size) {
		if(!seal_ptr_) {
			if(slices_.empty() || slices_.back().frontier_limit < SIZE_MAX)
				open_chunk_();

			size_t head_offset = (frontier_ + alignof(chunk) - 1) & ~(alignof(chunk) - 1);
			size_t object_offset = (head_offset + sizeof(chunk) + align - 1) & ~(align - 1);
			size_t n_frontier = object_offset + size;
			if(commit_limit_ < n_frontier)
				commit_(n_frontier);
			seal_ptr_ = slices_.back().base + head_offset;
			frontier_ = n_frontier;
			return slices_.back().base + object_offset;
		}else{
			size_t offset = (frontier_ + align - 1) & ~(align - 1);
			size_t n_frontier = offset + size;
			if(commit_limit_ < n_frontier)
				commit_(n_frontier);
			frontier_ = n_frontier;
			return slices_.back().base + offset;
		}
	}

	void deallocate() {
		// Do nothing for now.
	}

	void seal() {
		// TODO: Handle zero allocations per seal().
		assert(!slices_.empty());

		new (seal_ptr_) chunk{frontier_};
		seal_ptr_ = nullptr;

		if(frontier_ + slice_threshold > slice_space_)
			slices_.back().frontier_limit = frontier_;
	}

	void reclaim() {
		assert(!slices_.empty());

		chunk *rc = reclaim_chunk_;
		if(!reclaim_chunk_)
			rc = reinterpret_cast<chunk *>(slices_.front().base);
		auto rc_limit = rc->limit;
		rc->~chunk();

		if(rc_limit == slices_.front().frontier_limit) {
			close_chunk_();
			reclaim_chunk_ = nullptr;
		}else{
			assert(rc_limit < slices_.front().frontier_limit);
			size_t head_offset = (rc_limit + alignof(chunk) - 1) & ~(alignof(chunk) - 1);
			if(rc_limit > uncommit_limit_ + commit_threshold)
				uncommit_(rc_limit);
			reclaim_chunk_ = reinterpret_cast<chunk *>(slices_.front().base + head_offset);
		}
	}

private:
	void open_chunk_();
	void close_chunk_();
	void commit_(size_t threshold);
	void uncommit_(size_t threshold);

	std::deque<slice> slices_;
	void *seal_ptr_ = nullptr;
	chunk *reclaim_chunk_ = nullptr;
	size_t slice_space_;
	size_t frontier_ = 0;
	size_t commit_limit_ = 0;
	size_t uncommit_limit_ = 0;
};

template <class T>
struct queue_allocator {
	using value_type = T;

	queue_allocator(queue_memory &memory) noexcept
	: memory_{&memory} { }

	template <class U>
	queue_allocator(queue_allocator<U> const &other) noexcept
	: memory_{other.memory_} { }

	T *allocate(size_t n) {
		return reinterpret_cast<T *>(memory_->allocate(alignof(T), sizeof(T) * n));
	}

	void deallocate(T *, size_t) noexcept {
		memory_->deallocate();
	}

private:
	queue_memory *memory_;
};

template <class T, class U>
bool operator== (queue_allocator<T> const &lhs, queue_allocator<U> const &rhs) noexcept {
	return lhs.memory_ == rhs.memory_;
}

template <class T, class U>
bool operator!= (queue_allocator<T> const& lhs, queue_allocator<U> const& rhs) noexcept {
	return !(lhs == rhs);
}
