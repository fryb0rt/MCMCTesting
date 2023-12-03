#pragma once
#include "Pcg.h"
#include "Vector.h"
#include <tuple>
#include <vector>

template<uint32_t TDim>
INLINE Vector<TDim> randomVector(Pcg& rnd) {
	Vector<TDim> temp;
	for (uint32_t i = 0; i < TDim; ++i) {
		temp[i] = rnd();
	}
	return temp;
}

template<uint32_t TDim>
INLINE Vector<TDim> randomVectorExponential(Pcg& rnd, const Float b) {
	Vector<TDim> temp;
	for (uint32_t i = 0; i < TDim; ++i) {
		temp[i] = -log(1 - rnd() * (1 - exp(-b * 5))) / 5;
	}
	return temp;
}

template<typename TIterator, typename TCdf>
INLINE std::tuple<TIterator, Float> sampleDiscrete(const TIterator begin, const TIterator end, const TCdf cdf, Float& rnd) {
	using ValueType = typename std::iterator_traits<TIterator>::value_type;
	const ValueType searched = ValueType(rnd * cdf(*(end - 1)));
	TIterator found = std::upper_bound(begin, end, searched, [&cdf](const ValueType& a, const ValueType& b) { return cdf(a) < cdf(b); });
	if (found == end) {
		--found;
	}
	Float prev = found == begin ? Float(0) : cdf(*(found - 1));
	Float pdf = cdf(*found) - prev;
	rnd = (cdf(searched) - prev) / pdf;
	return std::make_tuple(found, pdf);
}

template<typename TIterator>
INLINE std::tuple<TIterator, Float> sampleDiscrete(const TIterator begin, const TIterator end, Float& rnd) {
	return sampleDiscrete(begin, end, [](const Float a) { return a; }, rnd);
}

template<uint32_t TPower>
INLINE Float intpow(const Float value) {
	const Float t = intpow< TPower / 2>(value);
	if ((TPower & 1) == 0) {
		return t * t;
	}
	else {
		return t * t * value;
	}
}

template<>
INLINE Float intpow<0>(const Float value) {
	return Float(1);
}

template<>
INLINE Float intpow<1>(const Float value) {
	return value;
}

INLINE Vector<3> hsvToRgb(const Vector<3>& hsv)
{
	auto rgb = [](const Float r, const Float g, const Float b) {
		Vector<3> out;
		out[0] = r;
		out[1] = g;
		out[2] = b;
		return out;
	};
	Float H = hsv[0] / Float(60), S = hsv[1], V = hsv[2];
	while (H >= Float(6)) {
		H -= Float(6);
	}
	while (H < Float(0)) {
		H += Float(6);
	}

	const Float fract = H - floor(H);

	const Float P = V*(Float(1) - S);
	const Float Q = V*(Float(1) - S*fract);
	const Float T = V*(Float(1) - S*(Float(1) - fract));

	if (Float(0) <= H && H < Float(1))
		return rgb(V, T, P);
	else if (Float(1) <= H && H < Float(2))
		return rgb(Q, V, P);
	else if (Float(2) <= H && H < Float(3))
		return rgb(P, V, T);
	else if (Float(3) <= H && H < Float(4))
		return rgb(P, Q, V);
	else if (Float(4) <= H && H < Float(5))
		return rgb(T, P, V);
	else if (Float(5) <= H && H < Float(6))
		return rgb(V, P, Q);
	else
		return Vector<3>(0);
}


template<uint32_t TDim>
INLINE Vector<TDim> log(const Vector<TDim>& v) {
	Vector<TDim> temp;
	for (uint32_t i = 0; i < TDim; ++i) {
		temp[i] = log(v[i]);
	}
	return temp;
}

template<uint32_t TDim>
INLINE Vector<TDim> exp(const Vector<TDim>& v) {
	Vector<TDim> temp;
	for (uint32_t i = 0; i < TDim; ++i) {
		temp[i] = exp(v[i]);
	}
	return temp;
}

INLINE Float halton(const uint32_t base, uint32_t i) {
	Float f = 1, r = 0;
	while (i > 0) {
		f /= base;
		r = r + f * (i % base);
		i = uint32_t(i / base);
	}
	return r;
}

constexpr static uint32_t PRIMES[] = { 2, 3, 5, 11, 17, 23, 31, 43, 59, 71, 89, 107, 131, 149, 163, 191, 211, 227};
template<uint32_t TDim>
INLINE Vector<TDim> haltonVector(uint32_t index) {
	Vector<TDim> temp;
	for (uint32_t d = 0; d < TDim; ++d) {
		temp[d] = halton(PRIMES[d], index);
	}
	return temp;
}

template<uint32_t TDim>
INLINE void nRooks(Pcg& rnd, std::vector<Vector<TDim>>& out) {
	const Vector<TDim> strata(Float(1) / out.size());
	// Generate initial positions
	for (uint32_t i = 0; i < out.size(); ++i) {
		out[i] = randomVector<TDim>(rnd) * strata + strata * Float(i);
	}
	// Shuffle each dimension (except the first one)
	for (uint32_t d = 1; d < TDim; ++d) {
		for (uint32_t i = 0; i < (out.size()-1); ++i) {
			std::swap(out[i][d], out[i + rnd.uint(uint32_t(out.size()) - i)][d]);
		}
	}
}

template<typename T>
T clamp(const T& value, const T& minValue, const T& maxValue) {
	return std::max(std::min(value, maxValue), minValue);
}