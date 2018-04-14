#pragma once
#include <math.h>
#include <array>
#include <opencv2/imgproc/imgproc.hpp>


class Deriche {
private:
	float k, a1, a2, a3, a4, b1, b2, c;
	float sigma;
	float expNeqSigma, expNeq2Sigma;

	template<typename T>
	void dericheIIR(cv::Mat & src, cv::Mat & dst, const int & rows, const int & cols, std::vector<float> & y1, std::vector<float> & y2) {
		int i, j;
		T * pixelPointer;

		for (i = 0; i < rows; ++i) {
			pixelPointer = src.ptr<T>(i);

			//do stuff for first 2 pixels NOTE: might need improvement.
			y1[0] = a1 * pixelPointer[0];
			y1[1] = a1 * pixelPointer[1] + a2 * pixelPointer[0] + b1 * y1[0];

			for (j = 2; j < cols; ++j) {
				y1[j] = a1 * pixelPointer[j] + a2 * pixelPointer[j - 1] + b1 * y1[j - 1] + b2 * y1[j - 2];
			}

			//do other stuff for first 2 pixels..
			y2[cols - 1] = a3 * pixelPointer[cols - 2];
			y2[cols - 2] = a3 * pixelPointer[cols - 3] + a4 * pixelPointer[cols - 4] + b1 * y2[cols - 1];

			for (j = cols - 3; j >= 0; --j) {
				y2[j] = a3 * pixelPointer[j + 1] + a4 * pixelPointer[j + 2] + b1 * y2[j + 1] + b2 * y2[j + 2];
			}

			//Merge both y1 and y2 and store them.
			for (j = 0; j < cols; ++j) {

				//we transpose the image here!
				dst.at<T>(j, i) = y1[j] + y2[j];
			}
		}
	}

	template<typename T>
	void dericheIIR2(cv::Mat & src, cv::Mat & dst, const int & rows, const int & cols, std::vector<float> & y1, std::vector<float> & y2) {
		int i, j;
		T * pixelPointer;

		for (i = 0; i < rows; ++i) {
			pixelPointer = src.ptr<T>(i);

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
				dst.at<T>(j, i) = c * (y1[j] + y2[j]);
			}
		}
	}

public:
	Deriche(const float & new_sigma = 1.f) {
		setSigma(new_sigma);
	};

	void setSigma(const float & new_sigma) {
		sigma = new_sigma;
		// calculate using the deriche coefficients. see https://en.wikipedia.org/wiki/Deriche_edge_detector and http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.476.5736&rep=rep1&type=pdf for more information.
		expNeqSigma = expf(-sigma);
		expNeq2Sigma = expf(-2.f * sigma);
		k = powf(1.f - expNeqSigma, 2.f) / (1.f + 2.f * sigma * expNeqSigma - expNeq2Sigma);
		b1 = 2.f * expNeqSigma;
		b2 = -expNeq2Sigma;
		a1 = k;
		a2 = k * expNeqSigma * (sigma - 1.f);
		a3 = k * expNeqSigma * (sigma + 1.f);
		a4 = -k * expNeq2Sigma;
		c = -powf(1.f - expNeqSigma, 2.f);

		std::cout << "DEBUG\n" <<
			"Sigma: " << sigma << "\n" <<
			"expNeqSigma: " << expNeqSigma << "\n" <<
			"expNeq2Sigma: " << expNeq2Sigma << "\n" <<
			"k: " << k << "\n" <<
			"b1: " << b1 << "\n" <<
			"b2: " << b2 << "\n";
		}

	void smooth(cv::Mat & matrix) {

		//asserting if its not CV_8U (unsigned char type matrix)
		CV_Assert(matrix.depth() == CV_8UC1 || matrix.depth() == CV_32FC1);

		if (matrix.depth() == CV_8UC1) {
			matrix.convertTo(matrix, CV_32FC1); // Convert matrix to float.
		}

		const int rows = matrix.rows;
		const int cols = matrix.cols;

		cv::Mat tempMatrix = cv::Mat(cols, rows, CV_32FC1);

		std::vector<float> y1((cols >= rows) ? cols : rows);
		std::vector<float> y2((cols >= rows) ? cols : rows);

		dericheIIR<float>(matrix, tempMatrix, rows, cols, y1, y2);
		dericheIIR<float>(tempMatrix, matrix, cols, rows, y1, y2);
	}

	void derivativeX(cv::Mat & matrix) {

		//asserting if its not CV_8U (unsigned char type matrix)
		CV_Assert(matrix.depth() == CV_32FC1);

		const int rows = matrix.rows;
		const int cols = matrix.cols;

		cv::Mat tempMatrix = cv::Mat(cols, rows, CV_32FC1);

		std::vector<float> y1((cols >= rows) ? cols : rows);
		std::vector<float> y2((cols >= rows) ? cols : rows);

		dericheIIR2<float>(matrix, tempMatrix, rows, cols, y1, y2);
		dericheIIR<float>(tempMatrix, matrix, cols, rows, y1, y2);
	}

	void derivativeY(cv::Mat & matrix) {

		//asserting if its not CV_8U (unsigned char type matrix)
		CV_Assert(matrix.depth() == CV_32FC1);

		const int rows = matrix.rows;
		const int cols = matrix.cols;

		cv::Mat tempMatrix = cv::Mat(cols, rows, CV_32FC1);

		std::vector<float> y1((cols >= rows) ? cols : rows);
		std::vector<float> y2((cols >= rows) ? cols : rows);

		dericheIIR<float>(matrix, tempMatrix, rows, cols, y1, y2);
		dericheIIR2<float>(tempMatrix, matrix, cols, rows, y1, y2);
	}

	cv::Mat deriche(cv::Mat & matrix) {
		if (matrix.depth() != CV_32FC1) {
			matrix.convertTo(matrix, CV_32FC1);
		}

		const int rows = matrix.rows;
		const int cols = matrix.cols;

		cv::Mat tempMatrix = cv::Mat(cols, rows, CV_32FC1);

		std::vector<float> y1((cols >= rows) ? cols : rows);
		std::vector<float> y2((cols >= rows) ? cols : rows);

		// Smooth
		dericheIIR<float>(matrix, tempMatrix, rows, cols, y1, y2);
		dericheIIR<float>(tempMatrix, matrix, cols, rows, y1, y2);

		// Prepare derivative matrices.
		cv::Mat imageMatrixX, imageMatrixY;
		matrix.convertTo(imageMatrixX, CV_32FC1);
		matrix.convertTo(imageMatrixY, CV_32FC1);

		// Gradient derivative X.
		dericheIIR2<float>(imageMatrixX, tempMatrix, rows, cols, y1, y2);
		dericheIIR<float>(tempMatrix, imageMatrixX, cols, rows, y1, y2);

		// Gradient derivative Y.
		dericheIIR<float>(imageMatrixY, tempMatrix, rows, cols, y1, y2);
		dericheIIR2<float>(tempMatrix, imageMatrixY, cols, rows, y1, y2);

		// Get gradient edges.
		cv::Mat edge_gradients = cv::Mat(imageMatrixX.rows, imageMatrixX.cols, CV_32FC1);
		for (int x = 0; x < imageMatrixX.rows; ++x) {
			for (int y = 0; y < imageMatrixX.cols; ++y) {
				edge_gradients.at<float>(x, y) = hypotf(imageMatrixX.at<float>(x, y), imageMatrixY.at<float>(x, y));
			}
		}

		// Get gradient directions.
		cv::Mat angles = cv::Mat(imageMatrixX.rows, imageMatrixX.cols, CV_8UC1);
		for (int y = 0; y < imageMatrixY.rows; ++y) {
			for (int x = 0; x < imageMatrixX.cols; ++x) {
				float Ypixel = imageMatrixY.at<float>(x, y);
				float Xpixel = imageMatrixX.at<float>(x, y); // Access violation on some images.

				float degrees = atan2f(Ypixel, Xpixel) * (180.f / 3.14159265358979f); // Convert radians to degrees using the float version of pi.

				if (degrees < 0.f) degrees = -degrees;

				angles.at<uchar>(x, y) = (uchar)(roundf(degrees / 45.f) * 45.f) % 180; // Round to nearest 45 for non-max-suppression.
			}
		}

		return nonMaxSuppression(edge_gradients, angles);
	}

	cv::Mat nonMaxSuppression(cv::Mat & gradient_edges, cv::Mat & gradient_directions) {
		cv::Mat output;
		gradient_edges.convertTo(output, CV_8UC1);

		float * gradPointer;
		uchar * dirPointer, *outputPointer;

		// All other pixels
		for (int y = 1; y < gradient_edges.rows - 1; ++y) {
			gradPointer = gradient_edges.ptr<float>(y);
			dirPointer = gradient_directions.ptr<uchar>(y);
			outputPointer = output.ptr<uchar>(y);

			for (int x = 1; x < gradient_edges.cols - 1; ++x) {
				switch (dirPointer[x]) {
					case 0:
						if (gradient_edges.at<float>(y, x - 1) > gradPointer[x] || 
							gradient_edges.at<float>(y, x + 1) > gradPointer[x]) {
							outputPointer[x] = 0;
						} else {
							outputPointer[x] = (uchar)gradPointer[x];
						}
						break;
					case 45:
						if (gradient_edges.at<float>(y - 1, x + 1) > gradPointer[x] ||
							gradient_edges.at<float>(y + 1, x - 1) > gradPointer[x]) {
							outputPointer[x] = 0;
						} else {
							outputPointer[x] = (uchar)gradPointer[x];
						}
						break;
					case 90:
						if (gradient_edges.at<float>(y - 1, x) > gradPointer[x] ||
							gradient_edges.at<float>(y + 1, x) > gradPointer[x]) {
							outputPointer[x] = 0;
						} else {
							outputPointer[x] = (uchar)gradPointer[x];
						}
						break;
					case 135:
						if (gradient_edges.at<float>(y - 1, x - 1) > gradPointer[x] || 
							gradient_edges.at<float>(y + 1, x + 1) > gradPointer[x]) {
							outputPointer[x] = 0;
						} else {
							outputPointer[x] = (uchar)gradPointer[x];
						}
						break;
				}
			}
		}

		return output;
	}
};