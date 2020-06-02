#include <sys/mman.h>

#include <cassert>
#include <iostream>
#include <new>
#include <type_traits>

#include "config.hpp"
#include "data-structures.hpp"
#include "mm.hpp"

namespace {
	// Ridiculously large virtual memory area.
	static constexpr size_t vspace = 1024 * size_t(1024 * 1024 * 1024);

	// Prefer 2 MiB pages.
	static constexpr size_t page_size = 0x200000;

	bool reserve_vm(size_t size, void *&window) {
		// PROT_NONE areas do not count towards overcommitting.
		// There is also MAP_NORESERVE to disable overcomitting heuristics.
		auto res = mmap(nullptr, size, PROT_NONE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
		if(res == MAP_FAILED)
			return false;
		window = res;
		return true;
	}

	void unreserve_vm(void *window, size_t size) {
		if(munmap(window, size)) {
			std::cerr << "munmap() failure: cannot unreserve virtual memory" << std::endl;
			abort();
		}
	}

	bool commit_vm(void *window, size_t size) {
		if(mprotect(window, size, PROT_READ | PROT_WRITE))
			return false;
		return true;
	}

	void uncommit_vm(void *window, size_t size) {
		// TODO: should we also call mprotect to PROT_NONE the area again?
#if defined(MADV_FREE) && !defined(TD_AVOID_MADVFREE)
		if(madvise(window, size, MADV_FREE)) {
			std::cerr << "madvise(MADV_FREE) failure:"
					" cannot uncommit virtual memory" << std::endl;
			abort();
		}
#else
#warning MADV_FREE is not available
		if(madvise(window, size, MADV_DONTNEED)) {
			std::cerr << "madvise(MADV_DONTNEED) failure:"
					" cannot uncommit virtual memory" << std::endl;
			abort();
		}
#endif
	}
};

// --------------------------------------------------------------------------------------

memory_arena::memory_arena() {
	void *base_ptr;
	if(!reserve_vm(vspace, base_ptr)) {
		std::cerr << "failed to reserve virtual memory for memory arena" << std::endl;
		throw std::bad_alloc{};
	}
	base_ = reinterpret_cast<char *>(base_ptr);
}

memory_arena::~memory_arena() {
	unreserve_vm(base_, vspace);
}

void memory_arena::extend(size_t threshold) {
	if(threshold <= limit_)
		return;
	threshold = ceil2pow(threshold);
	if(!commit_vm(base_, threshold)) {
		std::cerr << "failed to extend memory arena to " << print_memory(threshold) << std::endl;
		throw std::bad_alloc{};
	}
	limit_ = threshold;
}

void memory_arena::reset() {
	assert(!live_);
	uncommit_vm(base_, limit_);
	frontier_ = 0;
	limit_ = 0;
	num_allocations_ = 0;
}

// --------------------------------------------------------------------------------------

queue_memory::queue_memory()
: slice_space_{vspace} { }

queue_memory::~queue_memory() {
	while(!slices_.empty()) {
		unreserve_vm(slices_.front().base, vspace);
		slices_.pop_front();
	}
}

void queue_memory::open_chunk_() {
	assert(!seal_ptr_);

	void *base_ptr;
	if(!reserve_vm(vspace, base_ptr)) {
		std::cerr << "failed to reserve virtual memory for queue memory" << std::endl;
		throw std::bad_alloc{};
	}
	slices_.push_back(slice{reinterpret_cast<char *>(base_ptr)});
	frontier_ = 0;
	commit_limit_ = 0;
}

void queue_memory::close_chunk_() {
	assert(!slices_.empty());

	unreserve_vm(slices_.front().base, vspace);
	slices_.pop_front();
	uncommit_limit_ = 0;
}

void queue_memory::commit_(size_t threshold) {
	assert(!slices_.empty());
	assert(threshold > commit_limit_);

	size_t lim = (threshold + commit_threshold + page_size - 1) & ~(page_size - 1);
	assert(lim > commit_limit_);
	if(!commit_vm(slices_.back().base + commit_limit_, lim - commit_limit_)) {
		std::cerr << "failed to extend queue memory to " << print_memory(lim) << std::endl;
		throw std::bad_alloc{};
	}
	commit_limit_ = lim;
}

void queue_memory::uncommit_(size_t threshold) {
	assert(!slices_.empty());
	assert(threshold > uncommit_limit_);

	size_t lim = threshold & ~(page_size - 1);
	assert(lim > uncommit_limit_);
	uncommit_vm(slices_.front().base + uncommit_limit_, lim - uncommit_limit_);
	uncommit_limit_ = lim;
}
