#pragma once

#include <chrono>
#include <iomanip>
#include <iostream>

struct print_time {
	friend std::ostream &operator<< (std::ostream &os, const print_time &self) {
		std::ios ios{nullptr};
		ios.copyfmt(os);
		os << std::fixed << std::setprecision(3) << self.dur_.count() << " s";
		os.copyfmt(ios);
		return os;
	}

	template<typename R, typename P>
	print_time(std::chrono::duration<R, P> d)
	: dur_{std::chrono::duration_cast<std::chrono::duration<double>>(d)} { }

private:
	std::chrono::duration<double> dur_;
};

using profiling_duration = std::chrono::high_resolution_clock::duration;

extern bool global_profiling_flag;

struct profiling_timer {
	profiling_timer() {
		if(global_profiling_flag)
			t0_ = std::chrono::high_resolution_clock::now();
	}

	profiling_timer(const profiling_timer &) = delete;

	profiling_timer &operator= (const profiling_timer &) = delete;

	profiling_duration elapsed() {
		if(global_profiling_flag)
			return std::chrono::high_resolution_clock::now() - t0_;
		return profiling_duration{};
	}

private:
	std::chrono::high_resolution_clock::time_point t0_{};
};

// Like profiling timer. This class is not disabled by global_profiling_flag.
struct coarse_profiling_timer {
	coarse_profiling_timer() {
		t0_ = std::chrono::high_resolution_clock::now();
	}

	coarse_profiling_timer(const coarse_profiling_timer &) = delete;

	coarse_profiling_timer &operator= (const coarse_profiling_timer &) = delete;

	profiling_duration elapsed() {
		return std::chrono::high_resolution_clock::now() - t0_;
	}

private:
	std::chrono::high_resolution_clock::time_point t0_{};
};
