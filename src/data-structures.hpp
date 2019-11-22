#pragma once

#include <cassert>
#include <vector>

// Data structure to mark integers (e.g., vertices of a graph).
// Use a timestamping technique to accelerate resets.
struct boolean_marker {
	void reset(size_t n) {
		vec_.resize(n);
		if(ts_ == UINT8_MAX) {
			std::fill(vec_.begin(), vec_.end(), 0);
			ts_ = 0;
		}
		++ts_;
	}

	void mark(size_t x) {
		vec_[x] = ts_;
	}

	bool is_marked(size_t x) {
		return vec_[x] == ts_;
	}

private:
	std::vector<uint8_t> vec_;
	uint8_t ts_ = 0;
};

template<typename T>
struct dependent_false : std::false_type { };

template<typename T>
inline T bit_clz(T x) {
	static_assert(dependent_false<T>::value, "bit_clz() is not implemented for this type");
}

inline unsigned int bit_clz(unsigned int x) {
	return __builtin_clz(x);
}

inline unsigned long bit_clz(unsigned long x) {
	return __builtin_clzl(x);
}

inline unsigned long ceil2pow(unsigned long x) {
	static_assert(sizeof(unsigned long) == 8, "unexpected sizeof(long)");
	if(x <= 1)
		return 1;
	assert(x <= ((UINT64_MAX / 2) + 1));
	return (1UL << (64 - bit_clz(x - 1)));
}

inline unsigned int log2int(unsigned int x) {
	static_assert(sizeof(unsigned int) == 4, "unexpected sizeof(int)");
	assert(x >= 1);
	return (32 - bit_clz(x)) - 1;
}

inline unsigned int log_ceil2int(unsigned int x) {
	static_assert(sizeof(unsigned int) == 4, "unexpected sizeof(int)");
	if(x <= 1)
		return 0;
	return (32 - bit_clz(x - 1));
}
