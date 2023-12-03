#pragma once
#include "Algorithm.h"
#include "ReferenceAlgorithm.h"
#include "UniformAlgorithm.h"
#include "HaltonAlgorithm.h"
#include "MetropolisHastings.h"
//#include "GlobalAdaptation.h"
#include "ParallelTempering.h"
#include "AdaptiveEES.h"
#include "Permutations.h"
#include "EEM.h"
#include "Statistics.h"
#include "SampledSwaps.h"
#include <sstream>
#include <iostream>

template<uint32_t TDimension>
class TestSuite {
	static_assert(TDimension >= 2, "Dimension must be divisible by 2");
	const Integrand<TDimension> mIntegrand;
	Statistics<TDimension> mStats;
	std::vector<Algorithm<TDimension> *> mAlgorithms;
public:
	TestSuite(const uint32_t countDistributions,
		const Float minWeight,
		const Float maxWeight,
		const Float avgScale,
		const Float diffScale,
		const uint64_t seed):
		mIntegrand(randomMixture<TDimension>(countDistributions, minWeight, maxWeight, avgScale, diffScale, seed)), 
		mStats(mIntegrand, 1024) {
	}

	template<typename TAlgorithm, typename ... Types, typename = std::enable_if_t<std::is_base_of_v<Algorithm<TDimension>, TAlgorithm>>>
	void addAlgorithm(Types&& ... params) {
		mAlgorithms.push_back(new TAlgorithm(mIntegrand, std::forward<Types>(params)...));
	}

	void runAll(const uint32_t runCount, const uint32_t samplesPerRun) {
		std::vector<SampleAndPdf<TDimension>> samples(samplesPerRun);
		for (auto * alg : mAlgorithms) {
			mStats.clear();
			std::cout << "Executing " << runCount << " runs of " << alg->name();
			for (uint32_t r = 0; r < runCount; ++r) {
				alg->run(samples);
				mStats.addRunResult(samples, alg->hasNormalizedPdf());
				std::cout << ".";
			}
			std::cout << std::endl;
			mStats.histogram(alg->name());
			mStats.computeAndPrintStats(alg->name());
			alg->printStats();
		}
		std::cout << std::endl;
		mStats.printOverallStats();
	}

	~TestSuite() {
		for (auto * alg : mAlgorithms) {
			delete alg;
		}
	}
private:
	template<typename TAlgorithm, typename ... Types>
	bool tryCreateAlgorithm(const std::string &name, Types&& ... params) {
		if (TAlgorithm::sName() == name) {
			mAlgorithms.push_back(new TAlgorithm(mIntegrand, std::forward<Types>(params)...));
			return true;
		}
		return false;
	}
};