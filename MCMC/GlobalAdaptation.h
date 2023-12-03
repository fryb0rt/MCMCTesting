#pragma once
#include "Algorithm.h"
#include "LocalMutation.h"

template<uint32_t TDimension>
class GlobalAdaptationAlgorithm : public Algorithm<TDimension> {
	LocalMutation<TDimension> mMutator;
	struct Sample {
		INLINE Sample() = default;
		INLINE Sample(Float value) :value(value) {} // For CDF search
		Vector<TDimension> sample;
		Float value;
	};
	std::vector<Sample> mSamples;
	Float mLargeStepProb;
public:
	INLINE GlobalAdaptationAlgorithm(const Integrand<TDimension>& integrand) : Algorithm(integrand), mMutator(this->mRandom) {
		mLargeStepProb = 0.3f;
	}

	virtual void run(std::vector<SampleAndPdf<TDimension>>& samples) override {
		do {
			mMutator.setState(randomVector<TDimension>(mRandom));
		} while ((mIntegrand.value(mMutator.getState())) == Float(0));
		// Init algorithm 
		/*mSamples.clear();
		mSamples.reserve(100);
		for (uint32_t i = 0; i < uint32_t(mSamples.capacity()); ++i) {
			addSample(haltonVector<TDimension>(i));
		}*/
		//mMutator.startAdaptation(0.6f);
		int adaptiveStep = 0;
		for (uint32_t i = 0; i < uint32_t(samples.size()) + MCMC_BURN_PERIOD; ++i) {
			const bool localStep = i % 3 != 0;
			bool alwaysAccept = false;
			Float probRatio = 1;
			Vector<TDimension> proposed;
			if (localStep) {
				proposed = mMutator.mutateState();
			}
			else {
				if (++adaptiveStep % 2 == 0 && mSamples.size() > 1000) {
					Float temp = mRandom(), pdf;
					typename std::vector<Sample>::iterator it;
					std::tie(it, pdf) = sampleDiscrete(mSamples.begin(), mSamples.end(), [](const Sample& a) { return a.value; }, temp);
					/*proposed = it->sample + randomVector<TDimension>(mRandom) * 0.1f;
					probRatio = pdf / mIntegrand.value(mMutator.getState());*/
					//alwaysAccept = true;
				}
				else {
					proposed = randomVector<TDimension>(mRandom);
				}
			}
			if (acceptRatio(mMutator.getState(), proposed, probRatio) > mRandom() || alwaysAccept) {
				mMutator.mutationWasAccepted(proposed, localStep);
			}
			else {
				mMutator.mutationWasRejected(localStep);
			}
			if (i >= MCMC_BURN_PERIOD) {
				const uint32_t index = i - MCMC_BURN_PERIOD;
				samples[index].sample = mMutator.getState();
				samples[index].pdf = mIntegrand.value(samples[index].sample);
				addSample(mMutator.getState());
			}
		}
	}

	static std::string sName() {
		return "GlobalAdaptation";
	}

	virtual std::string name() const override {
		return sName();
	}

	virtual bool hasNormalizedPdf() const override { return false; }

	virtual void printStats() const override {
		std::cout << "Acceptance rate: " << Float(100) * mMutator.acceptanceRate() << " %" << std::endl;
	}
private:
	INLINE Float acceptRatio(const Vector<TDimension>& current, const Vector<TDimension>& proposed, const Float probRatio) const {
		const Float valCurr = mIntegrand.value(current);
		const Float valProp = mIntegrand.value(proposed);
		return valCurr > Float(0) ? std::min(Float(1), valProp * probRatio / valCurr) : Float(1);
	}

	INLINE void addSample(const Vector<TDimension>&v) {
		Sample s;
		s.sample = v;
		s.value = mIntegrand.value(v);
		if (s.value > 0) {
			Float prev = mSamples.size() > 0 ? mSamples.back().value : 0;
			s.value = 1 + prev;
			mSamples.push_back(s);
		}
	}
};