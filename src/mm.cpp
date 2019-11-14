#include <sys/mman.h>

#include <cassert>
#include <iostream>
#include <new>
#include <type_traits>

#include "mm.hpp"

namespace {
	template<typename T>
	struct dependent_false : std::false_type { };

	template<typename T>
	T ceil2pow(T x) {
		static_assert(dependent_false<T>::value, "ceil2pow() is not implemented for this type");
	}

	unsigned long ceil2pow(unsigned long x) {
		if(x <= 1)
			return 1;
		assert(x <= ((UINT64_MAX / 2) + 1));
		return (1 << (64 - __builtin_clzl(x - 1)));
	}

	// Ridiculously large virtual memory area.
	static size_t vspace = 1024 * size_t(1024 * 1024 * 1024);

	void *reserve_vm(size_t size) {
		// PROT_NONE areas do not count towards overcommitting.
		// There is also MAP_NORESERVE to disable overcomitting heuristics.
		auto window = mmap(nullptr, size, PROT_NONE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
		if(window == MAP_FAILED)
			throw std::bad_alloc{};
		return window;
	}

	void unreserve_vm(void *window, size_t size) {
		if(munmap(window, size)) {
			std::cerr << "munmap() failure: cannot unreserve virtual memory" << std::endl;
			abort();
		}
	}

	void commit_vm(void *window, size_t size) {
		if(mprotect(window, size, PROT_READ | PROT_WRITE))
			throw std::bad_alloc{};
	}

	void uncommit_vm(void *window, size_t size) {
		// TODO: should we also call mprotect to PROT_NONE the area again?
		if(madvise(window, size, MADV_FREE)) {
			std::cerr << "madvise(MADV_FREE) failure: cannot uncommit virtual memory" << std::endl;
			abort();
		}
	}
};

memory_arena::memory_arena() {
	base_ = reinterpret_cast<char *>(reserve_vm(vspace));
}

memory_arena::~memory_arena() {
	unreserve_vm(base_, vspace);
}

void memory_arena::extend(size_t threshold) {
	if(threshold <= limit_)
		return;
	threshold = ceil2pow(threshold);
	commit_vm(base_, threshold);
	limit_ = threshold;
}

void memory_arena::reset() {
	assert(!live_);
	uncommit_vm(base_, limit_);
	limit_ = 0;
	num_allocations_ = 0;
}
