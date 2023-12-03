#pragma once
#include "Algorithm.h"

template<uint32_t TDimension>
class ReferenceAlgorithm : public Algorithm<TDimension> {
public:
	INLINE ReferenceAlgorithm(const Integrand<TDimension>& integrand) : Algorithm<TDimension>(integrand) {}

	virtual void run(std::vector<SampleAndPdf<TDimension>>& samples) override {
		const auto& dist = this->mIntegrand.getDistribution();
		for (uint32_t i = 0; i < uint32_t(samples.size()); ++i) {
			samples[i].sample = dist.sample(randomVector<TDimension>(this->mRandom));
			samples[i].pdf = dist.pdf(samples[i].sample);
		}
	}

	static std::string sName() {
		return "Reference";
	}

	virtual std::string name() const override {
		return sName();
	}
};