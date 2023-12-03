#pragma once
#include <iostream>
#include <vector>
#include <array>
#include <algorithm>

class PermutationSampler {
	bool mSkipPartialIdentity;
	uint32_t mSize;
	std::vector<Float> mCache;
	std::vector<uint32_t> mPermut;
	uint32_t mOpCount, mSampleOp;
	Float mPositive;
public:

	PermutationSampler(const uint32_t size, const bool skipPartialIdentity):mSkipPartialIdentity(skipPartialIdentity), mSize(size) {
		mCache.resize(1i64 << size, -1.f);
		mPermut.resize(size);
		mPositive = Float(1);
	}

	template <typename TFunctor>
	const std::vector<uint32_t> & sample(const TFunctor& mFunctor, Float rnd, Float &count) {
		for (uint32_t i = 0; i < mSize; ++i) {
			mPermut[i] = i;
		}
		mOpCount = 0;
		mSampleOp = 0;
		count = countPermutation(mFunctor, 0, 0);
		for (uint32_t i = 0; i < mSize; ++i) {
			mPermut[i] = i;
		}
		samplePermutation(mFunctor, rnd * count, 0, 0);
		//std::cout << mOpCount << " " << mSampleOp << std::endl;
		mPositive *= Float(-1);
		return mPermut;
	}

	template <typename TFunctor>
	Float normalization(const TFunctor& mFunctor) {
		for (uint32_t i = 0; i < mSize; ++i) {
			mPermut[i] = i;
		}
		mOpCount = 0;
		Float result = countPermutation(mFunctor, 0, 0);
		mPositive *= Float(-1);
		return result;
	}

private:
	template <typename TFunctor>
	void samplePermutation(const TFunctor& mFunctor, Float rnd, const uint32_t level, const uint32_t cache) {
		for (uint32_t i = level; i < mSize; ++i) {
			++mSampleOp;
			std::swap(mPermut[i], mPermut[level]);
			if (!mSkipPartialIdentity || level != mPermut[level]) {
				const uint32_t newCache = cache | (1 << mPermut[level]);
				const Float funcValue = mFunctor(level, mPermut[level]);
				Float value = funcValue;
				if (level + 1 < mSize - 1) {
					value *= mCache[newCache] * mPositive;
				}
				else {
					value *= (!mSkipPartialIdentity || (mSize - 1) != mPermut[mSize - 1]) ? mFunctor(mSize - 1, mPermut[mSize - 1]) : 0.f;
				}
				if (value < rnd && i + 1 < mSize) {
					rnd -= value;
				}
				else {
					if (level + 1 < mSize - 1) {
						samplePermutation(mFunctor, rnd / funcValue, level + 1, newCache);
					}
					assert(mPermut[mSize - 1] != mSize - 1);
					return;
				}
			}
			std::swap(mPermut[i], mPermut[level]);
		}
		assert(false);
	}

	template <typename TFunctor>
	Float countPermutation(const TFunctor& mFunctor, const uint32_t level, const uint32_t cache) {
		if (level == mSize - 1) {
			return (!mSkipPartialIdentity || level != mPermut[level]) ? mFunctor(level, mPermut[level]) : 0.f;
		}
		Float & result = mCache[cache];
		if (result * mPositive < 0.f) {
			result = 0.f;
			for (uint32_t i = level; i < mSize; ++i) {
				++mOpCount;
				std::swap(mPermut[i], mPermut[level]);
				if (!mSkipPartialIdentity || level != mPermut[level]) {
					const uint32_t newCache = cache | (1 << mPermut[level]);
					result += mFunctor(level, mPermut[level]) * countPermutation(mFunctor, level + 1, newCache);
				}
				std::swap(mPermut[i], mPermut[level]);
			}
			result *= mPositive;
		}
		return abs(result);
	}

};