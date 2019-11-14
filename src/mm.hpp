#pragma once

#include <cstddef>
#include <cstdint>

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
