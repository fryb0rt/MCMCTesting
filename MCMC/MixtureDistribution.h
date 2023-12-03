#pragma once
#include "Utils.h"
#include "NormalDistribution.h"
#include <vector>

template<typename TDistribution>
class MixtureDistribution {
	std::vector<TDistribution> mDistributions;
	std::vector<Float> mCdf;
public:
	INLINE MixtureDistribution() = default;

	INLINE MixtureDistribution(const std::vector<TDistribution>& distributions, const std::vector<Float>& weights): mDistributions(distributions) {
		assert(weights.size() == distributions.size());
		mCdf.resize(weights.size());
		Float accum(0);
		for (int i = 0; i < weights.size(); ++i) {
			accum += weights[i];
			mCdf[i] = accum;
		}
	}

	template<typename TVector>
	INLINE Float pdf(const TVector& x) const {
		Float total(0), prevCdf(0);
		for (int i = 0; i < mDistributions.size(); ++i) {
			total += mDistributions[i].pdf(x) * (mCdf[i] - prevCdf);
			prevCdf = mCdf[i];
		}
		return total / prevCdf;
	}

	template<typename TVector>
	INLINE Float pdfTempered(const TVector& x, const Float invTemperature) const {
		Float total(0), prevCdf(0);
		for (int i = 0; i < mDistributions.size(); ++i) {
			total += mDistributions[i].pdfTempered(x, invTemperature) * (mCdf[i] - prevCdf);
			prevCdf = mCdf[i];
		}
		return total / prevCdf;
	}

	INLINE Float highestPdf(const Float invTemperature) const {
		Float highest(0), prevCdf(0);;
		for (int i = 0; i < mDistributions.size(); ++i) {
			highest = std::max(highest, mDistributions[i].highestPdf(invTemperature) * (mCdf[i] - prevCdf));
			prevCdf = mCdf[i];
		}
		return highest;
	}

	template<typename TVector>
	INLINE TVector sample(TVector rnd) const {
		std::vector<Float>::const_iterator it;
		Float pdf;
		std::tie(it, pdf) = sampleDiscrete(mCdf.begin(), mCdf.end(), rnd[0]);
		return mDistributions[it - mCdf.begin()].sample(rnd);
	}

	INLINE uint32_t modeCount() const {
		return uint32_t(mDistributions.size());
	}

	template<typename TVector>
	INLINE uint32_t mode(const TVector& x) const {
		uint32_t m = modeCount();
		Float prevCdf(0), maxP(0);
		for (int i = 0; i < mDistributions.size(); ++i) {
			const Float p = mDistributions[i].pdf(x) * (mCdf[i] - prevCdf);
			if (p > maxP) {
				m = i;
				maxP = p;
			}
		}
		return m;
	}
};

template<uint32_t TDimension>
INLINE MixtureDistribution<NormalDistribution<TDimension>> randomMixture(
	const uint32_t count, 
	const Float minWeight, 
	const Float maxWeight, 
	const Float avgScale, 
	const Float diffScale, 
	const uint64_t seed) {
	assert(count >= 1);
	static_assert(TDimension >= 2, "Dimension must be divisible by 2");
	std::vector<Float> weights;
	weights.resize(count);
	using Distribution = NormalDistribution<TDimension>;
	using Vec = Vector<TDimension>;
	using VecHalf = Vector<TDimension/2>;
	std::vector<Distribution> distributions;
	distributions.resize(count);
	Pcg random(seed, 1337);
	const Float dist = sqrt(avgScale);
	std::vector<Vec> means(count);
	nRooks(random, means);
	for (uint32_t i = 0; i < count; ++i) {
		weights[i] = minWeight + random() * (maxWeight - minWeight);
		Vec scale;
		for (uint32_t d = 0; d < TDimension / 2; ++d) {
			const Float r = 1 + random() * diffScale;
			scale[2 * d] = avgScale / r;
			scale[2 * d + 1] = avgScale * r;
		}
		distributions[i] = Distribution(
			Vec(dist * 3) + means[i] * Vec(1 - dist * 6),
			scale,
			VecHalf(0) + randomVector<TDimension/2>(random) * Float(2 * M_PI));
	}
	return MixtureDistribution<Distribution>(distributions, weights);
;}
