#pragma once
#include "Config.h"
#include "Vector.h"
#include "Matrix.h"

class NormalDistribution2D {
	Vector2 mMean;
	Matrix2x2  mSigmaCholesky, mInvertedSigma;
	Float mNormalization;
public:
	INLINE NormalDistribution2D() = default;

	INLINE NormalDistribution2D(const Vector2 mean, const Vector2 scale, const Float rotationRadians) :mMean(mean) {
		assert(scale[0] != Float(0) || scale[1] != Float(0));
		const Float sinRot = sin(rotationRadians);
		const Float cosRot = cos(rotationRadians);
		const Matrix2x2 rotMatrix = makeMatrix2x2(cosRot, -sinRot, sinRot, cosRot);
		const Matrix2x2 scaleMatrix = makeMatrix2x2(scale[0], Float(0), Float(0), scale[1]);
		const Matrix2x2 sigma = invert(rotMatrix) * scaleMatrix * rotMatrix;
		const Float det = determinant(sigma);
		assert(det > Float(0));
		mInvertedSigma = invert(sigma);
		mNormalization = Float(2 * M_PI) * sqrt(det);
		mSigmaCholesky = cholesky(sigma);
	}

	INLINE Float pdf(const Vector2 x) const {
		const Vector2 rel = x - mMean;
		const Float expArg = -0.5f * dot(rel, mInvertedSigma * rel);
		return exp(expArg) / mNormalization;
	}

	INLINE Float pdfTempered(const Vector2 x, const Float invTemperature) const {
		const Vector2 rel = x - mMean;
		const Float expArg = -0.5f * dot(rel, mInvertedSigma * rel) * invTemperature;
		return exp(expArg) / mNormalization * invTemperature;
	}
	
	INLINE Float highestPdf(const Float invTemperature) const {
		return invTemperature / mNormalization;
	}

	INLINE Vector2 sample(const Vector2 rnd) const {
		assert(1.f - rnd[0] > 0.f);
		const Float scale = sqrt(-2.f * log(1.f - rnd[0]));
		const Float angle = Float(2 * M_PI) * rnd[1];
		return mSigmaCholesky * (makeVector2(cos(angle), sin(angle)) * scale) + mMean;
	}
};

template<uint32_t TDim>
class NormalDistribution {
	static_assert(TDim >= 2 && TDim % 2 ==0, "Dimension must be positive and divisible by two");
	NormalDistribution2D mSubDistributions[TDim / 2];
public:
	INLINE NormalDistribution() = default;

	INLINE NormalDistribution(const Vector<TDim>& mean, const Vector<TDim>& scale, const Vector<TDim / 2>& rotationRadians) {
		for (int i = 0; i < TDim / 2; ++i) {
			mSubDistributions[i] = NormalDistribution2D(pickVector2(mean, 2 * i), pickVector2(scale, 2 * i), rotationRadians[i]);
		}
	}

	INLINE Float pdf(const Vector<TDim>& x) const {
		Float product = Float(1);
		for (int i = 0; i < TDim / 2; ++i) {
			product *= mSubDistributions[i].pdf(pickVector2(x, 2 * i));
		}
		return product;
	}

	INLINE Float pdfTempered(const Vector<TDim>& x, const Float invTemperature) const {
		Float product = Float(1);
		for (int i = 0; i < TDim / 2; ++i) {
			product *= mSubDistributions[i].pdfTempered(pickVector2(x, 2 * i), invTemperature);
		}
		return product;
	}

	INLINE Float highestPdf(const Float invTemperature) const {
		Float product = Float(1);
		for (int i = 0; i < TDim / 2; ++i) {
			product *= mSubDistributions[i].highestPdf(invTemperature);
		}
		return product;
	}

	INLINE Vector<TDim> sample(const Vector<TDim>& rnd) const {
		Vector<TDim> result;
		for (int i = 0; i < TDim / 2; ++i) {
			Vector2 temp = mSubDistributions[i].sample(pickVector2(rnd, 2 * i));
			result[2 * i] = temp[0];
			result[2 * i + 1] = temp[1];
		}
		return result;
	}
};