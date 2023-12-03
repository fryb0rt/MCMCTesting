#pragma once
#include "MixtureDistribution.h"
template<uint32_t TDimension>
class Integrand {
public:
	using TDist = MixtureDistribution<NormalDistribution<TDimension>>;
private:
	 TDist mDist;
public:
	INLINE Integrand(const TDist &dist) : mDist(dist) {}

	INLINE Float value(const Vector<TDimension>& v, const Float invTemperature) const {
		if (inside(v)) {
			return mDist.pdfTempered(v, invTemperature);
		} else {
			return 0.f;
		}
	}

	INLINE uint32_t modeCount() const {
		return mDist.modeCount();
	}

	INLINE uint32_t mode(const Vector<TDimension>& v) const {
		return mDist.mode(v);
	}

	INLINE const TDist& getDistribution() const {
		return mDist;
	}

	INLINE bool inside(const Vector<TDimension>& v) const {
		bool result = true;
		for (int i = 0; i < TDimension; ++i) {
			result &= v[i] >= 0.f && v[i] <= 1.f;
		}
		return result;
	}

	INLINE Float maxValue(const Float invTemperature) const {
		return mDist.highestPdf(invTemperature);
	}
};