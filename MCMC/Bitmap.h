#pragma once
#include "Vector.h"
#include "Utils.h"
#include <vector>
#include "CImg.h"

using Rgb = Vector<3>;

INLINE Rgb makeRgb(const Float r, const Float g, const Float b) {
	Rgb out;
	out[0] = r;
	out[1] = g;
	out[2] = b;
	return out;
}

class Bitmap {
	cimg_library::CImg<unsigned char> mImage;
public:
	INLINE Bitmap(const uint32_t width, const uint32_t height):mImage(width,height,1,3) {
	}

	INLINE const Rgb get(const uint32_t x, const uint32_t y) const {
		return makeRgb(mImage(x, y, 0, 0) / 255.f, mImage(x, y, 0, 1) / 255.f, mImage(x, y, 0, 2) / 255.f);
	}

	INLINE void set(const uint32_t x, const uint32_t y, const Rgb& value) {
		mImage(x, y, 0, 0) = unsigned char(value[0] * 255.f);
		mImage(x, y, 0, 1) = unsigned char(value[1] * 255.f);
		mImage(x, y, 0, 2) = unsigned char(value[2] * 255.f);
	}

	INLINE void save(const char * filePath) const {
		mImage.save(filePath);
	}
private:
};