#pragma once
#include "Integrand.h"
#include "SampleAndPdf.h"
#include <string>

constexpr uint32_t MCMC_BURN_PERIOD = 1000;

template<uint32_t TDimension>
class Algorithm {
protected:
	const Integrand<TDimension>& mIntegrand;
	Pcg mRandom;
public:
	INLINE Algorithm(const Integrand<TDimension>& integrand):
		mIntegrand(integrand) {
		mRandom = Pcg(0xDEAD, 1);
	}
	virtual void run(std::vector<SampleAndPdf<TDimension>>& samples) = 0;

	virtual std::string name() const = 0;

	virtual bool hasNormalizedPdf() const { return true; }

	virtual void printStats() const {};

	virtual ~Algorithm() {}
};