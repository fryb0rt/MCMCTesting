#pragma once
#include "Algorithm.h"

template<uint32_t TDimension>
class HaltonAlgorithm : public Algorithm<TDimension> {
public:
	INLINE HaltonAlgorithm(const Integrand<TDimension>& integrand) : Algorithm<TDimension>(integrand) {}

	virtual void run(std::vector<SampleAndPdf<TDimension>>& samples) override {
		for (uint32_t i = 0; i < uint32_t(samples.size()); ++i) {
			samples[i].sample = haltonVector<TDimension>(i);
			samples[i].pdf = Float(1);
		}
	}

	static std::string sName() {
		return "Halton";
	}

	virtual std::string name() const override {
		return sName();
	}
};