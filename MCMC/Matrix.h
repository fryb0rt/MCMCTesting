#pragma once
#include "Config.h"

template<uint32_t TRows, uint32_t TColumns>
class Matrix {
	static_assert(TRows >= 1, "Rows must be positive.");
	static_assert(TColumns >= 1, "Colums must be positive.");
	Vector<TColumns> mData[TRows];
public:

	INLINE Matrix() = default;

	INLINE Matrix(const Float x) {
		for (int i = 0; i < TRows; ++i) {
			mData[i] = Vector<TColumns>(x);
		}
	}
	INLINE Matrix operator+(const Matrix& v) const {
		Matrix temp;
		for (int i = 0; i < TRows; ++i) {
			temp[i] = mData[i] + v.mData[i];
		}
		return temp;
	}

	INLINE Matrix operator-(const Matrix& v) const {
		Matrix temp;
		for (int i = 0; i < TRows; ++i) {
			temp[i] = mData[i] - v.mData[i];
		}
		return temp;
	}

	template<uint32_t TColumns2>
	INLINE Matrix<TRows, TColumns2> operator*(const Matrix<TColumns, TColumns2>& v) const {
		Matrix temp;
		for (int r = 0; r < TRows; ++r) {
			for (int c = 0; c < TColumns2; ++c) {
				temp[r][c] = Float(0);
				for (int i = 0; i < TColumns; ++i) {
					temp[r][c] += mData[r][i] * v.mData[i][c];
				}
			}
		}
		return temp;
	}

	INLINE Matrix operator+=(const Matrix& v) {
		for (int i = 0; i < TRows; ++i) {
			mData[i] += v.mData[i];
		}
		return *this;
	}

	INLINE Matrix operator-=(const Matrix& v) {
		Matrix temp;
		for (int i = 0; i < TRows; ++i) {
			mData[i] -= v.mData[i];
		}
		return *this;
	}

	INLINE Matrix operator-() const {
		Matrix temp;
		for (int i = 0; i < TRows; ++i) {
			temp = -mData[i];
		}
		return temp;
	}

	INLINE Vector<TColumns> operator[](const uint32_t index) const {
		assert(index >= 0 && index < TRows);
		return mData[index];
	}

	INLINE Vector<TColumns>& operator[](const uint32_t index) {
		assert(index >= 0 && index < TRows);
		return mData[index];
	}

	INLINE Vector<TRows> operator*(const Vector<TColumns>& v) const {
		Vector<TRows> temp;
		for (int r = 0; r < TRows; ++r) {
			temp[r] = Float(0);
			for (int i = 0; i < TColumns; ++i) {
				temp[r] += mData[r][i] * v[i];
			}
		}
		return temp;
	}

	INLINE Matrix operator*=(const Float f) {
		for (int i = 0; i < TRows; ++i) {
			mData[i] *= f;
		}
		return *this;
	}

	INLINE Matrix operator/=(const Float f) {
		for (int i = 0; i < TRows; ++i) {
			mData[i] /= f;
		}
		return *this;
	}

	INLINE Matrix operator*(const Float f) const {
		Matrix temp;
		for (int i = 0; i < TRows; ++i) {
			temp[i] = mData[i] * f;
		}
		return temp;
	}

	INLINE Matrix operator/(const Float f) const {
		Matrix temp;
		for (int i = 0; i < TRows; ++i) {
			temp[i] = mData[i] / f;
		}
		return temp;
	}
};

using Matrix2x2 = Matrix<2, 2>;

INLINE Float determinant(const Matrix2x2 &m) {
	return m[0][0] * m[1][1] - m[0][1] * m[1][0];
}

INLINE Matrix2x2 invert(const Matrix2x2 &m) {
	const Float det = determinant(m);
	assert(det != Float(0));
	Matrix2x2 temp = m;
	std::swap(temp[0][0], temp[1][1]);
	temp[1][0] = -temp[1][0];
	temp[0][1] = -temp[0][1];
	return temp / det;
}

INLINE Matrix2x2 makeMatrix2x2(const Float m00, const Float m01, const Float m10, const Float m11) {
	Matrix2x2 temp;
	temp[0][0] = m00;
	temp[0][1] = m01;
	temp[1][0] = m10;
	temp[1][1] = m11;
	return temp;
}

INLINE Matrix2x2 cholesky(const Matrix2x2 &m) {
	const Float a = sqrt(m[0][0]);
	const Float b = m[0][1] / a;
	const Float c = sqrt(m[1][1] - b * b);
	return makeMatrix2x2(a, 0, b, c);
}