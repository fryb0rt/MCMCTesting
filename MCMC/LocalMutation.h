#pragma once
#include "Utils.h"

template<uint32_t TDimension>
class LocalMutation {
	Pcg& mRandom;
	bool mAdaptive;
	Float mGoalAcceptance;
	uint32_t mAcceptedAll, mRejectedAll, mAcceptedLocal, mRejectedLocal, mUpdates;
	Float mS1,mS2, mLogRatio;
	Float mLocalMutationSize;
	Vector<TDimension> mState;
public:
	INLINE LocalMutation(Pcg& random) :mRandom(random), mAdaptive(false) {
		mS1 = 1.0f / 1024.0f;
		mS2 = 1.0f / 64.0f;
		mLogRatio = -log(mS2 / mS1);
		mAcceptedAll = 0;
		mRejectedAll = 0;
	}

	INLINE Vector<TDimension> mutateState() {
		Vector<TDimension> temp;
		for (uint32_t d = 0; d < TDimension; ++d) {
			temp[d] = mutate1D(mState[d]);
		}
		return temp;
	}

	INLINE const Vector<TDimension>& getState() const {
		return mState;
	}

	INLINE void setState(const Vector<TDimension>& state) {
		mState = state;
	}

	INLINE void startAdaptation(const Float goalAcceptance) {
		mGoalAcceptance = goalAcceptance;
		mAdaptive = true;
		mAcceptedLocal = 0;
		mRejectedLocal = 0;
		mUpdates = 0;
		mLocalMutationSize = mS2;
	}

	INLINE void mutationWasAccepted(const Vector<TDimension>& proposed, const bool localMutation) {
		mState = proposed;
		++mAcceptedAll;
		++mAcceptedLocal;
		if (!mAdaptive || !localMutation) {
			return;
		}
		adaptMutation();
	}

	INLINE void mutationWasRejected(const bool localMutation) {
		++mRejectedAll;
		++mRejectedLocal;
		if (!mAdaptive || !localMutation) {
			return;
		}
		adaptMutation();
	}

	INLINE Float acceptanceRate() const {
		return mAcceptedAll / Float(mAcceptedAll + mRejectedAll);
	}
private:
	INLINE Float mutate1D(Float value) {
		Float sample = mRandom();
		bool add;

		if (sample < 0.5f) {
			add = true;
			sample *= 2.0f;
		}
		else {
			add = false;
			sample = 2.0f * (sample - 0.5f);
		}
		const Float dv = mAdaptive ? pow(sample, (1 / mLocalMutationSize) + Float(1)) :
			mS2 * exp(sample * mLogRatio);
		if (add) {
			value += dv;
			if (value >= 1)
				value -= 1;
		}
		else {
			value -= dv;
			if (value < 0)
				value += 1;
		}
		return value;
	}

	INLINE void adaptMutation() {
		const uint32_t count = mAcceptedLocal + mRejectedLocal;
		assert(count > 0);
		const Float ratio = mAcceptedLocal / Float(count);
		const Float newSize = mLocalMutationSize + (ratio - mGoalAcceptance) / (mUpdates + 1);
		if (newSize > 0 && newSize < 1) {
			++mUpdates;
			mLocalMutationSize = newSize;
		}
	}
};