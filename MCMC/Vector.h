#pragma once
#include "Config.h"

template<uint32_t TDim>
class Vector {
	static_assert(TDim >= 1, "Dimension must be positive.");
	Float mData[TDim];
public:

	INLINE Vector() = default;

	INLINE Vector(const Float x) {
		for (int i = 0; i < TDim; ++i) {
			mData[i] = x;
		}
	}

	INLINE Vector operator+(const Vector& v) const {
		Vector temp;
		for (int i = 0; i < TDim; ++i) {
			temp[i] = mData[i] + v.mData[i];
		}
		return temp;
	}

	INLINE Vector operator-(const Vector& v) const {
		Vector temp;
		for (int i = 0; i < TDim; ++i) {
			temp[i] = mData[i] - v.mData[i];
		}
		return temp;
	}

	INLINE Vector operator*(const Vector& v) const {
		Vector temp;
		for (int i = 0; i < TDim; ++i) {
			temp[i] = mData[i] * v.mData[i];
		}
		return temp;
	}

	INLINE Vector operator/(const Vector& v) const {
		Vector temp;
		for (int i = 0; i < TDim; ++i) {
			temp[i] = mData[i] / v.mData[i];
		}
		return temp;
	}

	INLINE Vector operator+=(const Vector& v) {
		for (int i = 0; i < TDim; ++i) {
			mData[i] += v.mData[i];
		}
		return *this;
	}

	INLINE Vector operator-=(const Vector& v) {
		Vector temp;
		for (int i = 0; i < TDim; ++i) {
			mData[i] -= v.mData[i];
		}
		return *this;
	}

	INLINE Vector operator*=(const Vector& v) {
		for (int i = 0; i < TDim; ++i) {
			mData[i] *= v.mData[i];
		}
		return *this;
	}

	INLINE Vector operator/=(const Vector& v) {
		for (int i = 0; i < TDim; ++i) {
			mData[i] /= v.mData[i];
		}
		return *this;
	}

	INLINE Vector operator-() const {
		Vector temp;
		for (int i = 0; i < TDim; ++i) {
			temp = -mData[i];
		}
		return temp;
	}

	INLINE Float operator[](const uint32_t index) const {
		assert(index >= 0 && index < TDim);
		return mData[index];
	}

	INLINE Float& operator[](const uint32_t index) {
		assert(index >= 0 && index < TDim);
		return mData[index];
	}

	INLINE Vector operator*=(const Float f) {
		for (int i = 0; i < TDim; ++i) {
			mData[i] *= f;
		}
		return *this;
	}

	INLINE Vector operator/=(const Float f) {
		for (int i = 0; i < TDim; ++i) {
			mData[i] /= f;
		}
		return *this;
	}

	INLINE Vector operator*(const Float f) const {
		Vector temp;
		for (int i = 0; i < TDim; ++i) {
			temp[i] = mData[i] * f;
		}
		return temp;
	}

	INLINE Vector operator/(const Float f) const {
		Vector temp;
		for (int i = 0; i < TDim; ++i) {
			temp[i] = mData[i] / f;
		}
		return temp;
	}
};

template<uint32_t TDim>
INLINE Vector<TDim> magnitudeSqr(const Vector<TDim>& v) {
	Float sqrSum = Float(0);
	for (int i = 0; i < TDim; ++i) {
		sqrSum += v[i] * v[i];
	}
	return sqrSum;
}

template<uint32_t TDim>
INLINE Vector<TDim> magnitude(const Vector<TDim>& v) {
	return sqrt(magnitudeSqr(v));
}

template<uint32_t TDim>
INLINE Vector<TDim> normalize(const Vector<TDim>& v) {
	return v / magnitude(v);
}

template<uint32_t TDim>
INLINE Float dot(const Vector<TDim>& v1, const Vector<TDim>& v2) {
	Float product = Float(0);
	for (int i = 0; i < TDim; ++i) {
		product += v1[i] * v2[i];
	}
	return product;
}

using Vector2 = Vector<2>;

INLINE Vector2 makeVector2(const Float v0, const Float v1) {
	Vector2 temp;
	temp[0] = v0;
	temp[1] = v1;
	return temp;
}

template<uint32_t TDim>
INLINE Vector2 pickVector2(const Vector<TDim> x, const uint32_t startIndex) {
	assert(startIndex >= 0 && startIndex + 1 < TDim);
	return makeVector2(x[startIndex], x[startIndex + 1]);
}