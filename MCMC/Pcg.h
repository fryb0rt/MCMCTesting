#pragma once
#include "Config.h"
// *Really* minimal PCG32 code / (c) 2014 M.E. O'Neill / pcg-random.org
// Licensed under Apache License 2.0 (NO WARRANTY, etc. see website)

class Pcg {
	uint64_t mState;  
	uint64_t mInc;
public:
	INLINE Pcg(const uint64_t initstate = 1337, const uint64_t initseq = 1) {
		reset(initstate, initseq);
	}

	INLINE void reset(const uint64_t initstate = 1337, const uint64_t initseq = 1) {
		mState = 0U;
		mInc = (initseq << 1u) | 1u;
		uint();
		mState += initstate;
		uint();
	}

	INLINE Float operator()() {
		return uint() / Float(UINT32_MAX);
	}

	INLINE uint32_t uint() {
		uint64_t oldstate = mState;
		// Advance internal state
		mState = oldstate * 6364136223846793005ULL + (mInc | 1);
		// Calculate output function (XSH RR), uses old state for max ILP
		uint32_t xorshifted = uint32_t(((oldstate >> 18u) ^ oldstate) >> 27u);
		uint32_t rot = uint32_t(oldstate >> 59u);
		return (xorshifted >> rot) | (xorshifted << (uint32_t(-int32_t(rot)) & 31));
	}

	INLINE uint32_t uint(const uint32_t limit) {
		const uint32_t end = (UINT32_MAX / limit) * limit;
		uint32_t value = uint();
		while (value >= end) {
			value = uint();
		}
		return value % limit;
	}

};