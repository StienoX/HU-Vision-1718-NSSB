#pragma once
#include <math.h>
#include <array>
#include <opencv2/imgproc/imgproc.hpp>


class Deriche {
private:
	double k, a1, a2, a3, a4, a5, a6, a7, a8, b1, b2, c1, c2;
	double sigma;
	double expNeqSigma, expNeq2Sigma;
public:
	Deriche(const double & a) {
		setSigma(a);
	};

	void setSigma(const double & new_sigma) {
		sigma = new_sigma;
		// calculate using the deriche coefficients. see https://en.wikipedia.org/wiki/Deriche_edge_detector and http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.476.5736&rep=rep1&type=pdf for more information.
		expNeqSigma = exp(-sigma);
		expNeqSigma = exp(-2 * sigma);
		k = pow(1 - expNeqSigma, 2.0) / (1 + 2 * sigma*expNeqSigma - expNeq2Sigma);
		b1 = 2 * expNeq2Sigma;
		b2 = -1 * expNeq2Sigma;
	}

	void smooth(cv::Mat & imageMatrix) {

		//asserting if its not CV_8U (unsigned char type matrix)
		CV_Assert(imageMatrix.depth() == CV_8U);

		// setting up smoothing variables
		a1 = k;
		a2 = k * expNeqSigma * (sigma - 1);
		a3 = k * expNeqSigma * (sigma + 1);
		a4 = -k * expNeq2Sigma;

		size_t rows = imageMatrix.rows;
		size_t cols = imageMatrix.cols;
		size_t i, j;
		uchar* pixelPointer;
		std::array<uchar,cols> y1;
		uchar y2[cols];

		for (i = 0; i < rows; ++i) {
			pixelPointer = imageMatrix.ptr<uchar>(i);

			//do stuff for first 2 pixels

			for (j = 2; j < cols; ++j) {

			}
		}
	}
};