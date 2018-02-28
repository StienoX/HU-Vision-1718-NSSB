#pragma once
#include <math.h>
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

	void smooth(const cv::Mat & imageMatrix) {
		// setting up smoothing variables
		a1 = k;
		a2 = k * expNeqSigma * (sigma - 1);
		a3 = k * expNeqSigma * (sigma + 1);
		a4 = -k * expNeq2Sigma;

		for (size_t m = 0; m < imageMatrix.cols; ++m) {
			for (size_t n = 2; n < imageMatrix.rows; ++n) {

			}
			for (size_t n = imageMatrix.rows - 3; n >= 0; --n) {

			}
			for (size_t n = 0; n < imageMatrix.rows; ++n) {

			}
		}

	}
};