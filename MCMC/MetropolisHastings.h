#pragma once
#include "Algorithm.h"
#include "LocalMutation.h"

template<uint32_t TDimension>
class MetropolisHastingsAlgorithm : public Algorithm<TDimension> {
	LocalMutation<TDimension> mMutator;
	Float mLargeStepProb;
public:
	INLINE MetropolisHastingsAlgorithm(const Integrand<TDimension>& integrand) : Algorithm<TDimension>(integrand), mMutator(this->mRandom) {
		mLargeStepProb = 0.3f;
	}

	virtual void run(std::vector<SampleAndPdf<TDimension>>& samples) override {
		do {
			mMutator.setState(randomVector<TDimension>(this->mRandom));
		} while ((this->mIntegrand.value(mMutator.getState(), Float(1))) == Float(0));
		mMutator.startAdaptation(0.3f);
		for (uint32_t i = 0; i < uint32_t(samples.size()) + MCMC_BURN_PERIOD; ++i) {
			const bool localStep = i % 3 != 0;
			Vector<TDimension> proposed;
			if (localStep) {
				proposed = mMutator.mutateState();
			}
			else {
				proposed = randomVector<TDimension>(this->mRandom);
			}
			if (acceptRatio(mMutator.getState(), proposed) > this->mRandom()) {
				mMutator.mutationWasAccepted(proposed, localStep);
			}
			else {
				mMutator.mutationWasRejected(localStep);
			}
			if (i >= MCMC_BURN_PERIOD) {
				const uint32_t index = i - MCMC_BURN_PERIOD;
				samples[index].sample = mMutator.getState();
				samples[index].pdf = this->mIntegrand.value(samples[index].sample, Float(1));
			}
		}
	}

	static std::string sName() {
		return "MetropolisHastings";
	}

	virtual std::string name() const override {
		return sName();
	}

	virtual bool hasNormalizedPdf() const override { return false; }

	virtual void printStats() const override {
		std::cout << "Acceptance rate: " << Float(100) * mMutator.acceptanceRate() << " %" << std::endl;
	}
private:
	INLINE Float acceptRatio(const Vector<TDimension>& current, const Vector<TDimension>& proposed) const {
		const Float valCurr = this->mIntegrand.value(current, Float(1));
		const Float valProp = this->mIntegrand.value(proposed, Float(1));
		return valCurr > Float(0) ? std::min(Float(1),valProp / valCurr) : Float(1);
	}
};