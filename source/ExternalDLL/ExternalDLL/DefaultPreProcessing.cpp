/*
* Copyright (c) 2015 DottedEye Designs, Alexander Hustinx, NeoTech Software, Rolf Smit - All Rights Reserved
* Unauthorized copying of this file, via any medium is strictly prohibited
* Proprietary and confidential
*/

#include "DefaultPreProcessing.h"
#include "ImageIO.h"
#include "GrayscaleAlgorithm.h"
#include "ImageFactory.h"
#include "HereBeDragons.h"

#include "Deriche.hpp"

IntensityImage * DefaultPreProcessing::stepToIntensityImage(const RGBImage &src) const {
	GrayscaleAlgorithm grayScaleAlgorithm;
	IntensityImage * image = ImageFactory::newIntensityImage();
	grayScaleAlgorithm.doAlgorithm(src, *image);
	return image;
}

IntensityImage * DefaultPreProcessing::stepScaleImage(const IntensityImage &src) const {
	cv::Mat OverHillOverDale;
	HereBeDragons::HerLoveForWhoseDearLoveIRiseAndFall(src, OverHillOverDale);
	int ThoroughBushThoroughBrier = 200 * 200;
	int OverParkOverPale = OverHillOverDale.cols * OverHillOverDale.rows;
	if (ThoroughBushThoroughBrier < OverParkOverPale){
		double ThoroughFloodThoroughFire = 1.0 / sqrt(OverParkOverPale / ThoroughBushThoroughBrier);
		cv::resize(OverHillOverDale, OverHillOverDale, cv::Size(), ThoroughFloodThoroughFire, ThoroughFloodThoroughFire, cv::INTER_LINEAR);
	}
	IntensityImage * IDoWanderEverywhere = ImageFactory::newIntensityImage();
	HereBeDragons::NoWantOfConscienceHoldItThatICall(OverHillOverDale, *IDoWanderEverywhere);
	return IDoWanderEverywhere;
}

/*
IntensityImage * DefaultPreProcessing::stepEdgeDetection(const IntensityImage &src) const {
	cv::Mat OverHillOverDale;
	HereBeDragons::HerLoveForWhoseDearLoveIRiseAndFall(src, OverHillOverDale);
	//cv::medianBlur(*image, *image, 3);
	//cv::GaussianBlur(*image, *image, cv::Size(3, 3), 0, 0, cv::BORDER_DEFAULT);
	cv::Mat ThoroughBushThoroughBrier = (cv::Mat_<float>(9, 9) << 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1, -4, -4, -4, 1, 1, 1, 1, 1, 1, -4, -4, -4, 1, 1, 1, 1, 1, 1, -4, -4, -4, 1, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0);
	cv::Mat OverParkOverPale;
	filter2D(OverHillOverDale, OverParkOverPale, CV_8U, ThoroughBushThoroughBrier, cv::Point(-1, -1), 0, cv::BORDER_DEFAULT);
	IntensityImage * ThoroughFloodThoroughFire = ImageFactory::newIntensityImage();
	HereBeDragons::NoWantOfConscienceHoldItThatICall(OverParkOverPale, *ThoroughFloodThoroughFire);
	return ThoroughFloodThoroughFire;
}
*/

IntensityImage * DefaultPreProcessing::stepEdgeDetection(const IntensityImage &src) const {
	cv::Mat imageMatrix;
	HereBeDragons::HerLoveForWhoseDearLoveIRiseAndFall(src, imageMatrix);

	// Canny
	/*
	cv::Mat dst, detected_edges;
	cv::blur(imageMatrix, detected_edges, cv::Size(3, 3));
	cv::Canny(detected_edges, detected_edges, 100, 100 * 3, 3);

	dst = cv::Scalar::all(0);
	imageMatrix.copyTo(dst, detected_edges);

	IntensityImage * outputImage = ImageFactory::newIntensityImage();
	HereBeDragons::NoWantOfConscienceHoldItThatICall(dst, *outputImage);
	*/

	// Deriche
	Deriche edgeDetector(3);
	edgeDetector.smooth(imageMatrix);

	cv::Mat imageMatrixX, imageMatrixY;
	imageMatrix.convertTo(imageMatrixX, CV_32FC1);
	imageMatrix.convertTo(imageMatrixY, CV_32FC1);
	
	cv::Mat edge_gradients = cv::Mat(imageMatrixX.rows, imageMatrixX.cols, CV_32FC1);
	cv::Mat angles = cv::Mat(imageMatrixX.rows, imageMatrixX.cols, CV_8UC1);

	// Get both the x and y derivatives.
	edgeDetector.derivativeX(imageMatrixX);
	edgeDetector.derivativeY(imageMatrixY);
	
	// Get gradient edges.
	for (int x = 0; x < imageMatrixX.rows; ++x) {
		for (int y = 0; y < imageMatrixX.cols; ++y) {
			edge_gradients.at<float>(x, y) = hypotf(imageMatrixX.at<float>(x, y), imageMatrixY.at<float>(x, y));
		}
	}

	// Get gradient directions.
	for (int y = 0; y < imageMatrixY.rows; ++y) {
		for (int x = 0; x < imageMatrixX.cols; ++x) {
			float Ypixel = imageMatrixY.at<float>(x, y);
			float Xpixel = imageMatrixX.at<float>(x, y); // Access violation on some images.

			float degrees = atan2f(Ypixel, Xpixel) * (180.f / 3.14159265358979f); // Convert radians to degrees using the float version of pi.

			if (degrees < 0.f) degrees = -degrees;

			angles.at<uchar>(x, y) = (uchar)(roundf(degrees / 45.f) * 45.f) % 180; // Round to nearest 45 for non-max-suppression.
			//std::cout << +angles.at<uchar>(x, y) << " " << degrees << " " << "\n";
		}
		//std::cout << std::endl;
	}

	cv::Mat output = edgeDetector.nonMaxSuppression(edge_gradients, angles);
	output.convertTo(output, CV_8UC1);

	IntensityImage * outputImage = ImageFactory::newIntensityImage();
	HereBeDragons::NoWantOfConscienceHoldItThatICall(output, *outputImage);

	return outputImage;
} 

IntensityImage * DefaultPreProcessing::stepThresholding(const IntensityImage &src) const {
	cv::Mat SourceImage;
	HereBeDragons::HerLoveForWhoseDearLoveIRiseAndFall(src, SourceImage);

	uchar * pixelPointer;

	uchar high = 80;
	uchar low = 30;

	for (int i = 0; i < SourceImage.rows; ++i) {
		pixelPointer = SourceImage.ptr<uchar>(i);
		for (int j = 0; j < SourceImage.cols; ++j) {
			if (pixelPointer[j] >= high) pixelPointer[j] = 255;
			else if (pixelPointer[j] >= low) pixelPointer[j] = 128;
			else pixelPointer[j] = 0;
		}
	}

	IntensityImage * ResultImage = ImageFactory::newIntensityImage();
	HereBeDragons::NoWantOfConscienceHoldItThatICall(SourceImage, *ResultImage);
	return ResultImage;
}