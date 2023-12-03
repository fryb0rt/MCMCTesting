#pragma once
#include "Algorithm.h"

template<uint32_t TDimension>
class UniformAlgorithm : public Algorithm<TDimension> {
public:
	INLINE UniformAlgorithm(const Integrand<TDimension>& integrand) : Algorithm<TDimension>(integrand) {}

	virtual void run(std::vector<SampleAndPdf<TDimension>>& samples) override {
		for (uint32_t i = 0; i < uint32_t(samples.size()); ++i) {
			samples[i].sample = randomVector<TDimension>(this->mRandom);
			samples[i].pdf = Float(1);
		}
	}

	static std::string sName() {
		return "Uniform";
	}

	virtual std::string name() const override {
		return sName();
	}
};