#pragma once
#include "Vector.h"

template<uint32_t TDimension>
struct SampleAndPdf {
	Vector<TDimension> sample;
	Float pdf;
};