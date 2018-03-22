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
	Deriche edgeDetector(4);
	edgeDetector.smooth(imageMatrix);
	cv::Mat imageMatrixY, edge_gradients, angles;

	// Prepare some of the required matrices.
	imageMatrix.copyTo(imageMatrixY);
	imageMatrix.copyTo(edge_gradients);
	imageMatrix.copyTo(angles);

	// Get both the x and y derivatives.
	edgeDetector.derivativeX(imageMatrix);
	edgeDetector.derivativeY(imageMatrixY);
	
	// Get gradient edges.
	for (int x = 0; x < imageMatrix.rows; ++x) {
		for (int y = 0; y < imageMatrix.cols; ++y) {
			edge_gradients.at<uchar>(x, y) = hypot(imageMatrix.at<uchar>(x, y), imageMatrixY.at<uchar>(x, y));
		}
	}

	// Get gradient directions.
	for (int y = 0; y < imageMatrix.rows; ++y) {
		for (int x = 0; x < imageMatrix.cols; ++x) {
			float Ypixel = imageMatrixY.at<uchar>(x, y);
			float Xpixel = imageMatrix.at<uchar>(x, y);

			float degrees = atan2f(Ypixel, Xpixel) * (180.0f / CV_PI) * 2; // Convert radians to degrees. MULTIPLIED BY 2 TO GET TO 180 MIGHT BE WRONG.

			angles.at<uchar>(x, y) = (uchar)(round(degrees / 45.0f) * 45.0f) % 180; // Round to nearest 45 for non-max-suppression.
			//std::cout << +angles.at<uchar>(x, y) << " " << degrees << " " << atan2f(Ypixel, Xpixel) << "\n";
			//if (angles.at<uchar>(x, y) == 135) std::cout << "REEEEEEEEEEEEEEEEE\n";
		}
		//std::cout << std::endl;
	}

	edgeDetector.nonMaxSuppression(edge_gradients, angles);

	//std::cout << outputImageMatrix;
	IntensityImage * outputImage = ImageFactory::newIntensityImage();
	HereBeDragons::NoWantOfConscienceHoldItThatICall(edge_gradients, *outputImage);

	return outputImage;
} 

IntensityImage * DefaultPreProcessing::stepThresholding(const IntensityImage &src) const {
	cv::Mat OverHillOverDale;
	HereBeDragons::HerLoveForWhoseDearLoveIRiseAndFall(src, OverHillOverDale);
	cv::threshold(OverHillOverDale, OverHillOverDale, 220, 255, cv::THRESH_BINARY_INV);
	IntensityImage * ThoroughBushThoroughBrier = ImageFactory::newIntensityImage();
	HereBeDragons::NoWantOfConscienceHoldItThatICall(OverHillOverDale, *ThoroughBushThoroughBrier);
	return ThoroughBushThoroughBrier;
}