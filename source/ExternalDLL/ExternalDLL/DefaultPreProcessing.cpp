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
	cv::Mat output = edgeDetector.deriche(imageMatrix);	
	IntensityImage * outputImage = ImageFactory::newIntensityImage();
	HereBeDragons::NoWantOfConscienceHoldItThatICall(output, *outputImage);

	return outputImage;
} 

IntensityImage * DefaultPreProcessing::stepThresholding(const IntensityImage &src) const {
	cv::Mat SourceImage;
	HereBeDragons::HerLoveForWhoseDearLoveIRiseAndFall(src, SourceImage);

	uchar * pixelPointer;

	float highThresholdRatio = 0.15;
	float lowThresholdRatio = 0.7;

	// Determine highest pixel value
	uchar max_pixel = 0;
	for (int i = 0; i < SourceImage.rows; ++i) {
		pixelPointer = SourceImage.ptr<uchar>(i);
		for (int j = 0; j < SourceImage.cols; ++j) {
			if (pixelPointer[j] > max_pixel) {
				max_pixel = pixelPointer[j];
			}
			else if (max_pixel == 255) {
				break;
			}
		}
	}

	uchar highThreshold = max_pixel * highThresholdRatio;
	uchar lowThreshold = highThreshold * lowThresholdRatio;

	// Find high and weak edges (double thresholding)
	for (int i = 0; i < SourceImage.rows; ++i) {
		pixelPointer = SourceImage.ptr<uchar>(i);
		for (int j = 0; j < SourceImage.cols; ++j) {
			if (pixelPointer[j] >= highThreshold) pixelPointer[j] = 255;
			else if (pixelPointer[j] >= lowThreshold) pixelPointer[j] = 128;
			else pixelPointer[j] = 0;
		}
	}
	
	// Follow weak edges (edge tracking by hysteresis)
	for (int i = 1; i < SourceImage.rows - 1; ++i) {
		pixelPointer = SourceImage.ptr<uchar>(i);
		for (int j = 1; j < SourceImage.cols - 1; ++j) {
			bool connected = false;
			if (pixelPointer[j] == 128) {
				for (int a = i - 1; a < i + 1; ++a) {
					if (connected) break;
					for (int b = j - 1; b < i + 1; ++b) {
						if (SourceImage.at<uchar>(a, b) == 255) connected = true;
					}
				}
				if (!connected) pixelPointer[j] = 0;
				else pixelPointer[j] = 255;
			}
		}
	}
	

	IntensityImage * ResultImage = ImageFactory::newIntensityImage();
	HereBeDragons::NoWantOfConscienceHoldItThatICall(SourceImage, *ResultImage);
	return ResultImage;
}