#include "StudentPreProcessing.h"
#include "HereBeDragons.h"
#include "ImageIO.h"
#include "Deriche.hpp"



IntensityImage * StudentPreProcessing::stepToIntensityImage(const RGBImage &image) const {
	return nullptr;
}

IntensityImage * StudentPreProcessing::stepScaleImage(const IntensityImage &image) const {
	return nullptr;
}

IntensityImage * StudentPreProcessing::stepEdgeDetection(const IntensityImage &image) const {
	/*
	cv::Mat inputImageMatrix;
	HereBeDragons::HerLoveForWhoseDearLoveIRiseAndFall(image,inputImageMatrix);

	Deriche edgeDetector(0.5);
	cv::Mat outputImageMatrix = edgeDetector.smooth(inputImageMatrix);
	std::cout << outputImageMatrix;
	IntensityImage * outputImage;
	//return outputImage;
	*/
	return nullptr;
}

IntensityImage * StudentPreProcessing::stepThresholding(const IntensityImage &image) const {
	return nullptr;
}