#pragma once
#include "Vector.h"
#include "Utils.h"
#include "Bitmap.h"
#include "Integrand.h"
#include "SampleAndPdf.h"
#include <iostream>
#include <string>
#include <iomanip>
#include <sstream>

struct AlgStats {
	Float avgSecMomentDiff, worstSecMomentDiff, bestSecMomentDiff;
	Float avgMiss, worstMiss, bestMiss;
	Float avgModesDiff, worstModesDiff, bestModesDiff;
	std::string algName;
};

template<uint32_t TDimension>
class Statistics {
	class RunStats {
		Float mSquares;
		uint32_t mSamples;
		std::vector<uint32_t> mModes;
	public:
		INLINE RunStats() = default;

		INLINE RunStats(const uint32_t modeCount): mSquares(0), mSamples(0), mModes(modeCount, 0) {
		}
		
		INLINE uint32_t sampleCount() const {
			return mSamples;
		}

		INLINE uint32_t modesHit() const {
			uint32_t miss = 0;
			for (uint32_t hits : mModes) {
				miss += uint32_t(hits > 0);
			}
			return miss;
		}

		INLINE uint32_t mostModeVisits() const {
			uint32_t hitsMost = 0;
			for (uint32_t hits : mModes) {
				hitsMost = std::max(hits, hitsMost);
			}
			return hitsMost;
		}

		INLINE uint32_t leastModeVisits() const {
			uint32_t hitsleast = mModes[0];
			for (uint32_t hits : mModes) {
				hitsleast = std::min(hits, hitsleast);
			}
			return hitsleast;
		}

		INLINE Float secondMoment() const {
			return mSquares / mSamples;
		}

		INLINE void addSample(const Float integrandValue, const Float pdf, const uint32_t mode) {
			if (integrandValue > 0) {
				++mModes[mode];
			}
			mSquares += integrandValue * integrandValue / pdf;
			++mSamples;
		}
	};
	std::vector<Vector<TDimension>> mSamples;
	const Integrand<TDimension> mIntegrand;
	std::vector<RunStats> mRuns;
	Float mRefSecondMoment, mRefFirstMoment;
	const uint32_t mResolution;
	std::vector<AlgStats> mOverAllStats;
public:
	INLINE Statistics(const Integrand<TDimension> &integrand, const uint32_t resolution):mIntegrand(integrand), mResolution(resolution) {
		mRefSecondMoment = Float(0);
		mRefFirstMoment = Float(0);
#ifdef _DEBUG
		mSamples.reserve(1000);
#else
		mSamples.reserve(1000000);
#endif
		
		for (uint32_t i = 0; i < uint32_t(mSamples.capacity()); ++i) {
			const auto s = integrand.getDistribution().sample(haltonVector<TDimension>(i));
			const Float v = integrand.value(s, Float(1));
			mRefFirstMoment += ((v / integrand.getDistribution().pdf(s)) - mRefFirstMoment) / (i + 1);
			mRefSecondMoment += ((v * v / integrand.getDistribution().pdf(s)) - mRefSecondMoment) / (i + 1);
			if (mIntegrand.inside(s)) {
				mSamples.push_back(s);
			}
		}
		histogram("Integrand");
		mSamples.clear();
	}

	INLINE void addRunResult(const std::vector<SampleAndPdf<TDimension>>& samples, const bool hasNormalizedPdf) {
		mRuns.push_back(RunStats(mIntegrand.modeCount()));
		mSamples.clear();
		const Float pdfNormalization = hasNormalizedPdf ? Float(1) : mRefFirstMoment;
		for (const auto& s : samples) {
			if (mIntegrand.inside(s.sample)) {
				mSamples.push_back(s.sample);
			}
			const Float value = mIntegrand.value(s.sample, Float(1));
			mRuns.back().addSample(value, s.pdf / pdfNormalization, mIntegrand.mode(s.sample));
		}
	}

	INLINE void clear() {
		mSamples.clear();
		mRuns.clear();
	}

	INLINE uint32_t runCount() const {
		return uint32_t(mRuns.size());
	}

	INLINE void computeAndPrintStats(const std::string& algName) {
		std::cout << std::right;
		std::cout << "============================================================" << std::endl;
		std::cout << algName << std::endl << std::endl;
		uint32_t i = 0;
		std::cout << std::fixed << std::setfill(' ');
		auto printLine = [this](const Float secondMomentDiff, const Float miss, const Float modesDiff, const bool endLine = true) {
			std::cout << std::setprecision(5) << " Second moment diff: " << std::setw(11) << secondMomentDiff << " %";
			std::cout << std::setprecision(5) << " Modes miss: " << std::setw(11) << miss << " %";
			std::cout << std::setprecision(5) << " Modes diff: " << std::setw(11) << modesDiff << " %";
			if (endLine) {
				std::cout << std::endl;
			}
		};

		AlgStats result;
		result.algName = algName;
		result.avgSecMomentDiff = Float(0);
		result.worstSecMomentDiff = Float(0);
		result.bestSecMomentDiff = Float(100);
		result.avgMiss = Float(0);
		result.bestMiss = Float(100);
		result.worstMiss = Float(0);
		result.avgModesDiff = Float(0); 
		result.worstModesDiff = Float(0);
		result.bestModesDiff = 10e7;
		for (const auto& run : mRuns) {
			const Float diff = Float(100) * abs(run.secondMoment()/ mRefSecondMoment - 1.f);
			const Float miss = Float(100) * Float(mIntegrand.modeCount() - run.modesHit()) / mIntegrand.modeCount();
			const Float modesDiff = Float(100) * Float(run.mostModeVisits() - run.leastModeVisits()) / run.sampleCount();
			std::cout << "Run #" << std::setw(3) << (++i);
			printLine(diff, miss, modesDiff, false);
			std::cout << "(M: " << std::setw(5) << run.mostModeVisits() << " L: " << std::setw(5) << run.leastModeVisits() << ")" << std::endl;
			
			result.avgSecMomentDiff += diff;
			result.avgMiss += miss;
			result.avgModesDiff += modesDiff;
			result.worstSecMomentDiff = std::max(diff, result.worstSecMomentDiff);
			result.worstMiss = std::max(result.worstMiss, miss);
			result.worstModesDiff = std::max(result.worstModesDiff, modesDiff);
			result.bestSecMomentDiff = std::min(diff, result.bestSecMomentDiff);
			result.bestMiss = std::min(result.bestMiss, miss);
			result.bestModesDiff = std::min(result.bestModesDiff, modesDiff);
		}
		result.avgSecMomentDiff /= mRuns.size();
		result.avgMiss /= mRuns.size();
		result.avgModesDiff /= mRuns.size();
		std::cout << std::endl;
		std::cout << "BEST    ";
		printLine(result.bestSecMomentDiff, result.bestMiss, result.bestModesDiff);
		std::cout << "AVERAGE ";
		printLine(result.avgSecMomentDiff, result.avgMiss, result.avgModesDiff);
		std::cout << "WORST   ";
		printLine(result.worstSecMomentDiff, result.worstMiss, result.worstModesDiff);
		std::cout << std::endl;
		mOverAllStats.push_back(result);
	}

	INLINE void printOverallStats() const {
		// Gather overall stats
		AlgStats bestAvgDiff = mOverAllStats[0], bestAvgMiss = mOverAllStats[0], bestBestDiff = mOverAllStats[0], bestBestMiss = mOverAllStats[0], bestWorstDiff = mOverAllStats[0], bestWorstMiss = mOverAllStats[0],
			bestAvgModesDiff = mOverAllStats[0], bestBestModesDiff = mOverAllStats[0], bestWorstModesDiff = mOverAllStats[0];
		AlgStats worstAvgDiff = mOverAllStats[0], worstAvgMiss = mOverAllStats[0], worstBestDiff = mOverAllStats[0], worstBestMiss = mOverAllStats[0], worstWorstDiff = mOverAllStats[0], worstWorstMiss = mOverAllStats[0],
			worstAvgModesDiff = mOverAllStats[0], worstBestModesDiff = mOverAllStats[0], worstWorstModesDiff = mOverAllStats[0];
		for (const auto& algStats : mOverAllStats) {
			//BEST
			if (bestAvgDiff.avgSecMomentDiff > algStats.avgSecMomentDiff) {
				bestAvgDiff = algStats;
			}
			if (bestAvgMiss.avgMiss > algStats.avgMiss) {
				bestAvgMiss = algStats;
			}
			if (bestAvgModesDiff.avgModesDiff > algStats.avgModesDiff) {
				bestAvgModesDiff = algStats;
			}
			if (bestBestDiff.bestSecMomentDiff > algStats.bestSecMomentDiff) {
				bestBestDiff = algStats;
			}
			if (bestBestMiss.bestMiss > algStats.bestMiss) {
				bestBestMiss = algStats;
			}
			if (bestBestModesDiff.bestModesDiff > algStats.bestModesDiff) {
				bestBestModesDiff = algStats;
			}
			if (bestWorstDiff.worstSecMomentDiff > algStats.worstSecMomentDiff) {
				bestWorstDiff = algStats;
			}
			if (bestWorstMiss.worstMiss > algStats.worstMiss) {
				bestWorstMiss = algStats;
			}
			if (bestWorstModesDiff.worstModesDiff > algStats.worstModesDiff) {
				bestWorstModesDiff = algStats;
			}
			//WORST
			if (worstAvgDiff.avgSecMomentDiff < algStats.avgSecMomentDiff) {
				worstAvgDiff = algStats;
			}
			if (worstAvgMiss.avgMiss < algStats.avgMiss) {
				worstAvgMiss = algStats;
			}
			if (worstAvgModesDiff.avgModesDiff < algStats.avgModesDiff) {
				worstAvgModesDiff = algStats;
			}
			if (worstBestDiff.bestSecMomentDiff < algStats.bestSecMomentDiff) {
				worstBestDiff = algStats;
			}
			if (worstBestMiss.bestMiss < algStats.bestMiss) {
				worstBestMiss = algStats;
			}
			if (worstBestModesDiff.bestModesDiff < algStats.bestModesDiff) {
				worstBestModesDiff = algStats;
			}
			if (worstWorstDiff.worstSecMomentDiff < algStats.worstSecMomentDiff) {
				worstWorstDiff = algStats;
			}
			if (worstWorstMiss.worstMiss < algStats.worstMiss) {
				worstWorstMiss = algStats;
			}
			if (worstWorstModesDiff.worstModesDiff < algStats.worstModesDiff) {
				worstWorstModesDiff = algStats;
			}
		}

		auto printLine = [this](const std::string& category, const AlgStats& stats, const Float value) {
			std::cout << std::left << std::setw(20) << std::setfill(' ') << category << std::setw(20) << stats.algName
				      << " " << std::right << std::setw(11) << std::setprecision(5) << value << " %" << std::endl;
		};


		std::cout << std::fixed << std::left;
		std::cout << "============================================================" << std::endl;
		std::cout << "OVERALL COMPARISON" << std::endl << std::endl;
		std::cout << "BEST" << std::endl;
		printLine("Average diff: ", bestAvgDiff, bestAvgDiff.avgSecMomentDiff);
		printLine("Average miss: ", bestAvgMiss, bestAvgMiss.avgMiss);
		printLine("Average modes diff: ", bestAvgModesDiff, bestAvgModesDiff.avgModesDiff);
		printLine("Best diff: ", bestBestDiff, bestBestDiff.bestSecMomentDiff);
		printLine("Best miss: ", bestBestMiss, bestBestMiss.bestMiss);
		printLine("Best modes diff: ", bestBestModesDiff, bestBestModesDiff.bestModesDiff);
		printLine("Worst diff: ", bestWorstDiff, bestWorstDiff.worstSecMomentDiff);
		printLine("Worst miss: ", bestWorstMiss, bestWorstMiss.worstMiss);
		printLine("Worst modes diff: ", bestWorstModesDiff, bestWorstModesDiff.worstModesDiff);
		std::cout << std::endl;

		std::cout << "WORST" << std::endl;
		printLine("Average diff: ", worstAvgDiff, worstAvgDiff.avgSecMomentDiff);
		printLine("Average miss: ", worstAvgMiss, worstAvgMiss.avgMiss);
		printLine("Average modes diff: ", worstAvgModesDiff, worstAvgModesDiff.avgModesDiff);
		printLine("Best diff: ", worstBestDiff, worstBestDiff.bestSecMomentDiff);
		printLine("Best miss: ", worstBestMiss, worstBestMiss.bestMiss);
		printLine("Best modes diff: ", worstBestModesDiff, worstBestModesDiff.bestModesDiff);
		printLine("Worst diff: ", worstWorstDiff, worstWorstDiff.worstSecMomentDiff);
		printLine("Worst miss: ", worstWorstMiss, worstWorstMiss.worstMiss);
		printLine("Worst modes diff: ", worstWorstModesDiff, worstWorstModesDiff.worstModesDiff);
		std::cout << std::endl;
		std::cout << "============================================================" << std::endl;
	}

	INLINE void histogram(const std::string& name) const {
		for (int d = 0; d < TDimension / 2; ++d) {
			std::ostringstream oss;
			oss << (2 * d + 1) << "," << (2*d + 2) << "-" << name << ".bmp";
			histogram(2 * d, 2 * d + 1).save(oss.str().c_str());
		}
	}

	INLINE Bitmap integrand(const uint32_t dim1, const uint32_t dim2, const uint32_t mResolution, const Float invTemperature) const {
		Bitmap out(mResolution, mResolution);
		Float highestValue(0);
		auto getPos = [this, dim1, dim2, mResolution](const uint32_t i1, const uint32_t i2) {
			Vector<TDimension> v(0.5f);
			v[dim1] = (i1 + Float(0.5)) / Float(mResolution);
			v[dim2] = (i2 + Float(0.5)) / Float(mResolution);
			return v;
		};
		for (uint32_t i1 = 0; i1 < mResolution; ++i1) {
			for (uint32_t i2 = 0; i2 < mResolution; ++i2) {
				highestValue = std::max(mIntegrand.value(getPos(i1, i2), invTemperature), highestValue);
			}
		}
		for (uint32_t i1 = 0; i1 < mResolution; ++i1) {
			for (uint32_t i2 = 0; i2 < mResolution; ++i2) {
				out.set(i1, i2, mapSampleCount(mIntegrand.value(getPos(i1, i2), invTemperature), highestValue));
			}
		}
		return out;
	}

private:

	INLINE Bitmap histogram(const uint32_t dim1, const uint32_t dim2) const {
		Bitmap out(mResolution, mResolution);
		std::vector<uint32_t> bins(mResolution * mResolution, 0);
		const auto index = [this](const Float value) {
			return std::min(uint32_t(value * mResolution), mResolution - 1);
		};
		for (const auto &v : mSamples) {
			++bins[index(v[dim1]) + index(v[dim2]) * mResolution];
		}
		uint32_t highestValue = 0;
		for (auto c : bins) {
			highestValue = std::max(c, highestValue);
		}
		for (uint32_t i1 = 0; i1 < mResolution; ++i1) {
			for (uint32_t i2 = 0; i2 < mResolution; ++i2) {
				out.set(i1, i2, mapSampleCount(Float(bins[i1 + i2 * mResolution]), Float(highestValue)));
			}
		}
		return out;
	}

	INLINE Rgb mapSampleCount(const Float value, const Float highestValue) const {
		if (value == 0.f) {
			return Rgb(Float(1));
		}
		const Float fraction = value / highestValue;
		Vector<3> hsv;
		hsv[0] = 240 - fraction * Float(240);
		hsv[1] = Float(1);
		hsv[2] = Float(1);
		return hsvToRgb(hsv);
	}
};