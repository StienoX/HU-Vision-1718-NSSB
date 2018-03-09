#pragma once
#include <math.h>
#include <array>
#include <opencv2/imgproc/imgproc.hpp>


class Deriche {
private:
	double k, a1, a2, a3, a4, b1, b2, c;
	double sigma;
	double expNeqSigma, expNeq2Sigma;
	uchar* pixelPointer;


	void dericheIIR(cv::Mat & src, cv::Mat & dst, const int & rows, const int & cols, std::vector<double> & y1, std::vector<double> & y2) {
		int i, j;
		
		for (i = 0; i < rows; ++i) {
			pixelPointer = src.ptr<uchar>(i);

			//do stuff for first 2 pixels NOTE: might need improvement.
			y1[0] = (a1 + a2) * pixelPointer[0];
			y1[1] = a1 * pixelPointer[1] + a2 * pixelPointer[1] + b1 * y1[0];

			for (j = 2; j < cols; ++j) {
				y1[j] = a1 * pixelPointer[j] + a2 * pixelPointer[j - 1] + b1 * y1[j - 1] + b2 * y1[j - 2];
			}

			//do other stuff for first 2 pixels..
			y2[cols - 1] = (a3 + a4) * pixelPointer[cols - 2];
			y2[cols - 2] = a3 * pixelPointer[cols - 3] + a4 * pixelPointer[cols - 4] + b1 * y2[cols - 1];

			for (j = cols - 3; j >= 0; --j) {
				y2[j] = a3 * pixelPointer[j + 1] + a4 * pixelPointer[j + 2] + b1 * y2[j + 1] + b2 * y2[j + 2];
			}

			//Merge both y1 and y2 and store them.
			for (j = 0; j < cols; ++j) {

				//we transpose the image here!
				dst.at<uchar>(j, i) = (uchar)(y1[j] + y2[j]);
			}
		}
	}

	void dericheIIR2(cv::Mat & src, cv::Mat & dst, const int & rows, const int & cols, std::vector<double> & y1, std::vector<double> & y2) {
		int i, j;

		for (i = 0; i < rows; ++i) {
			pixelPointer = src.ptr<uchar>(i);

			//do stuff for first 2 pixels NOTE: might need improvement.
			y1[0] = pixelPointer[0];
			y1[1] = pixelPointer[1] + b1 * y1[0];

			for (j = 2; j < cols; ++j) {
				y1[j] = pixelPointer[j - 1] + b1 * y1[j - 1] + b2 * y1[j - 2];
			}

			//do other stuff for first 2 pixels..
			y2[cols - 1] = -pixelPointer[cols - 2];
			y2[cols - 2] = -pixelPointer[cols - 3] + b1 * y2[cols - 1];

			for (j = cols - 3; j >= 0; --j) {
				y2[j] = -pixelPointer[j + 1] + b1 * y2[j + 1] + b2 * y2[j + 2];
			}

			//Merge both y1 and y2 and store them.
			for (j = 0; j < cols; ++j) {

				//we transpose the image here!
				dst.at<uchar>(j, i) = c * (uchar)(y1[j] + y2[j]);
			}
		}
	}

public:
	Deriche(const double & new_sigma = 1) {
		setSigma(new_sigma);
	};

	void setSigma(const double & new_sigma) {
		sigma = new_sigma;
		// calculate using the deriche coefficients. see https://en.wikipedia.org/wiki/Deriche_edge_detector and http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.476.5736&rep=rep1&type=pdf for more information.
		expNeqSigma = exp(-sigma);
		expNeq2Sigma = exp(-2 * sigma);
		k = pow(1 - expNeqSigma, 2.0) / (1 + 2 * sigma*expNeqSigma - expNeq2Sigma);
		b1 = 2 * expNeqSigma;
		b2 = -expNeq2Sigma;
		a1 = k;
		a2 = k * expNeqSigma * (sigma - 1);
		a3 = k * expNeqSigma * (sigma + 1);
		a4 = -k * expNeq2Sigma;
		c = -pow(1 - expNeqSigma, 2.0);

		std::cout << "DEBUG\n" <<
			"Sigma: " << sigma << "\n " <<
			"expNeqSigma: " << expNeqSigma << "\n " <<
			"expNeq2Sigma: " << expNeq2Sigma << "\n " <<
			"k: " << k << "\n " <<
			"b1: " << b1 << "\n " <<
			"b2: " << b2 << "\n ";
		}

	void smooth(cv::Mat & matrix) {

		//asserting if its not CV_8U (unsigned char type matrix)
		CV_Assert(matrix.depth() == CV_8U);

		const int rows = matrix.rows;
		const int cols = matrix.cols;

		cv::Mat tempMatrix = cv::Mat(cols, rows, CV_8U);

		std::vector<double> y1((cols >= rows) ? cols : rows);
		std::vector<double> y2((cols >= rows) ? cols : rows);
		
		dericheIIR(matrix, tempMatrix, rows, cols, y1, y2);
		dericheIIR(tempMatrix, matrix, cols, rows, y1, y2);
	}

	void derivativeX(cv::Mat & matrix) {

		//asserting if its not CV_8U (unsigned char type matrix)
		CV_Assert(matrix.depth() == CV_8U);

		const int rows = matrix.rows;
		const int cols = matrix.cols;

		cv::Mat tempMatrix = cv::Mat(cols, rows, CV_8U);

		std::vector<double> y1((cols >= rows) ? cols : rows);
		std::vector<double> y2((cols >= rows) ? cols : rows);

		dericheIIR2(matrix, tempMatrix, rows, cols, y1, y2);
		dericheIIR(tempMatrix, matrix, cols, rows, y1, y2);
		dericheIIR(matrix, tempMatrix, rows, cols, y1, y2);
		dericheIIR2(tempMatrix, matrix, cols, rows, y1, y2);


		

	}

};