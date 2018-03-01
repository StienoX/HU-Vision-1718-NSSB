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
	Deriche(const double & new_sigma = 1) {
		setSigma(new_sigma);
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

	cv::Mat smooth(cv::Mat & imageMatrix) {

		//asserting if its not CV_8U (unsigned char type matrix)
		CV_Assert(imageMatrix.depth() == CV_8U);

		// setting up smoothing variables
		a1 = k;
		a2 = k * expNeqSigma * (sigma - 1);
		a3 = k * expNeqSigma * (sigma + 1);
		a4 = -k * expNeq2Sigma;

		const int rows = imageMatrix.rows;
		const int cols = imageMatrix.cols;

		cv::Mat outputMatrix = cv::Mat(rows,cols,CV_8U);

		int i, j;

		//uchar* pixelPointer;
		uchar* pixelPointer;
		std::vector<double> y1(cols);
		std::vector<double> y2(cols);
		

		for (i = 0; i < rows; ++i) {
			pixelPointer = imageMatrix.ptr<uchar>(i);

			//do stuff for first 2 pixels NOTE: might need improvement.
			y1[0] = a1 * pixelPointer[0];
			y1[1] = a1 * pixelPointer[1] + a2 * pixelPointer[0] + b1 * y1[0];

			for (j = 2; j < cols; ++j) {
				y1[j] = a1 * pixelPointer[j] + a2 * pixelPointer[j - 1] + b1 * y1[j - 1] + b2 * y1[j - 2];
			}

			//do other stuff for first 2 pixels..
			y2[cols-1] = a3 * pixelPointer[cols-1];
			y2[cols-2] = a3 * pixelPointer[cols-2] + a4 * pixelPointer[0] + b1 * y2[cols - 1];

			for (j = cols - 3; j >= 0; --j) {
				y2[j] = a3 * pixelPointer[j+1] + a4 * pixelPointer[j + 2] + b1 * y1[j + 1] + b2 * y1[j + 2];
			}

			//Merge both y1 and y2 and store them.
			for (j = 0; j < cols; ++j) {
				outputMatrix.at<uchar>(i,j) = (uchar)(y1[j] + y2[j]);
			}
		}

		return outputMatrix;
		uchar pixel;

		for (j = 0; j < cols; ++j) {
			
			y1[0] = a1 * pixelPointer[0];
			y1[1] = a1 * pixelPointer[1] + a2 * pixelPointer[0] + b1 * y1[0];
			for (i = 2; i < rows; ++i) {
				y1[i] = a1 * pixelPointer[i] + a2 * pixelPointer[i - 1] + b1 * y1[i - 1] + b2 * y1[i - 2];
			}
		}
	}
};