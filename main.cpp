#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <algorithm>
#include <climits>
#include <complex>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <random>

#include <math.h>
#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr
#define M_PI 3.14159265358979323846 /* pi */

#include "image.h"

using namespace std;

void readImage(char fname[], ImageType &image);
void writeImage(const char fname[], ImageType &image);
void fft(float data[], unsigned long nn, int isign);

void test2dfft();
void fft2d(int N, int M, double **real_fuv, double **image_fuv, int isign);
void imgtodata(ImageType &image, float **data);
void datatoimg(ImageType &image, float **data);

void experiment1(char fname[]);
void experiment2(char fname[]);
void degrade(char fname[], float mu, float sigma, double **real_fuv, double **image_fuv);
void inversefilter(char fname[], float mu, float sigma, float radius);
void wienerfilter(char fname[], float mu, float sigma, float k);
void experiment4(char fname[], float gh, float gl);

void visualizespectrum(double **real, double **image, int M, int N, const char name[]);

float box_muller(float m, float s);
double generateGaussianNoise(const double &variance);

const float sobel_mask[3][3] = {{-1, 0, 1}, {-2, 0, 2}, {-1, 0, 1}};

static const int mask_7x7[7][7] = {

	{1, 1, 2, 2, 2, 1, 1},
	{1, 2, 2, 4, 2, 2, 1},
	{2, 2, 4, 8, 4, 2, 2},
	{2, 4, 8, 16, 8, 4, 2},
	{2, 2, 4, 8, 4, 2, 2},
	{1, 2, 2, 4, 2, 2, 1},
	{1, 1, 2, 2, 2, 1, 1}

};

int main(int argc, char *argv[])
{
	char girl[] = "girl.pgm";
	char boy[] = "boy_noisy.pgm";
	char lenna[] = "lenna.pgm";

	//test2dfft();

	// experiment1(boy);
	experiment2(lenna);

	// inversefilter(lenna, 0, 1, 40);
	// wienerfilter(lenna, 0, 1, 0.0001);
	// experiment4(girl, 1.5, 0.5);


	return 0;
}

void visualizespectrum(double **real, double **image, int M, int N, const char name[]){
	double **spec = new double *[M];
	for (int i = 0; i < M; i++)
	{
		spec[i] = new double[N];
	}
	ImageType spectrum(M, N, 255);
	for (int i = 0; i < M; i++)
	{
		for (int j = 0; j < N; j++)
		{
			spec[i][j] = sqrt(real[i][j] * real[i][j] + image[i][j] * image[i][j]);
			spec[i][j] = log(1 + spec[i][j]);
		}
	}
	// Normalization
	float max = spec[0][0];
	float min = spec[0][0];
	for (int i = 0; i < M; i++)
	{
		for (int j = 0; j < N; j++)
		{
			if (spec[i][j] > max)
			{
				max = spec[i][j];
			}
			if (spec[i][j] < min)
			{
				min = spec[i][j];
			}
		}
	}
	for (int i = 0; i < M; i++)
	{
		for (int j = 0; j < N; j++)
		{
			spec[i][j] = 255 * (spec[i][j] - min) / (max - min);
			spectrum.setPixelVal(i, j, spec[i][j]);
		}
	}
	writeImage(name, spectrum);

	for (int i = 0; i < M; ++i)
	{
		delete[] spec[i];
	}
	delete[] spec;
}

void experiment1(char fname[])
{
	ImageType baseImage(512, 512, 255);
	ImageType paddedImage(1024, 1024, 255);
	readImage(fname, baseImage);

	// Part A-------------------------------------------------------------------------------------------------------------------------
	//  Padding base image
	int temp;
	for (int i = 0; i < 512; i++)
	{
		for (int j = 0; j < 512; j++)
		{
			baseImage.getPixelVal(i, j, temp);
			paddedImage.setPixelVal(i, j, temp);
		}
	}

	// Step 1, FT
	double **real_fuv = new double *[1024];
	for (int i = 0; i < 1024; i++)
	{
		real_fuv[i] = new double[1024];
	}
	double **image_fuv = new double *[1024];
	for (int i = 0; i < 1024; i++)
	{
		image_fuv[i] = new double[1024];
	}

	// center
	for (int i = 0; i < 1024; i++)
	{
		for (int j = 0; j < 1024; j++)
		{
			paddedImage.getPixelVal(i, j, temp);
			real_fuv[i][j] = temp * pow(-1, i + j);
			image_fuv[i][j] = 0;
		}
	}

	fft2d(1024, 1024, real_fuv, image_fuv, -1);

	double **spectrum = new double *[1024];
	for (int i = 0; i < 1024; i++)
	{
		spectrum[i] = new double[1024];
	}
	for (int i = 0; i < 1024; i++)
	{
		for (int j = 0; j < 1024; j++)
		{
			// Calculate spectrum
			spectrum[i][j] = sqrt(real_fuv[i][j] * real_fuv[i][j] + image_fuv[i][j] + image_fuv[i][j]);
			spectrum[i][j] = log(1 + spectrum[i][j]);
		}
		
	}
	// Print spectrum
	char boySpectrum[] = "boySpec.pgm";
	ImageType boySpec(1024, 1024, 255);
	for (int i = 0; i < 1024; i++)
	{
		for (int j = 0; j < 1024; j++)
		{
			temp = spectrum[i][j] * pow(-1, i + j);
			boySpec.setPixelVal(i, j, temp);
		}
	}
	writeImage(boySpectrum, boySpec);
	
	// complex division by H(u,v)
	for (int i = 0; i < 1024; i++)
	{
		for (int j = 0; j < 1024; j++)
		{
			// Band-Reject
			// H(u,v) = 1 / 1 + [ DW / D^2 - D0^2 ]^2n
			int i_adj = i - 1024 / 2, j_adj = j - 1024 / 2;
			double D = sqrt(i_adj * i_adj + j_adj * j_adj);
			double DW = D * 1;
			double DD0 = D * D - 35.5 * 35.5;
			double denom = 1 + pow(DW / DD0, 2 * 4);
			double H = 1 / denom;

			real_fuv[i][j] = real_fuv[i][j] * H;
			image_fuv[i][j] = image_fuv[i][j] * H;

			// cout << H << " " << endl;
		}
	}

	// Step 4 Inverse FT
	fft2d(1024, 1024, real_fuv, image_fuv, 1);

	// Step 5 uncenter
	ImageType finalImage(512, 512, 255);
	for (int i = 0; i < 512; i++)
	{
		for (int j = 0; j < 512; j++)
		{
			temp = real_fuv[i][j] * pow(-1, i + j);
			finalImage.setPixelVal(i, j, temp);
		}
	}

	// Normalization
	finalImage.getPixelVal(0, 0, temp);
	float rmax = temp;
	float rmin = temp;
	for (int i = 0; i < 512; i++)
	{
		for (int j = 0; j < 512; j++)
		{
			finalImage.getPixelVal(i, j, temp);
			if (temp > rmax)
			{
				rmax = temp;
			}
			if (temp < rmin)
			{
				rmin = temp;
			}
		}
	}
	for (int i = 0; i < 512; i++)
	{
		for (int j = 0; j < 512; j++)
		{
			finalImage.getPixelVal(i, j, temp);
			temp = 255 * (temp - rmin) / (rmax - rmin);
			finalImage.setPixelVal(i, j, temp);
		}
	}

	// print and clean up
	char finImage_reject[] = "reject.pgm";
	writeImage(finImage_reject, finalImage);

	for (int i = 0; i < 1024; ++i)
	{
		delete[] real_fuv[i];
	}
	delete[] real_fuv;
	for (int i = 0; i < 1024; ++i)
	{
		delete[] image_fuv[i];
	}
	delete[] image_fuv;
	for (int i = 0; i < 1024; ++i)
	{
		delete[] spectrum[i];
	}
	delete[] spectrum;
}

void experiment2(char fname[])
{
	ImageType baseImage(256, 256, 255);
	readImage(fname, baseImage);

	// Spatial Filter start ///////////////////////////////////////////////////////////////
	int temp;
	float pixels[256 * 256];
	for (int i = 0; i < 256; i++)
	{
		for (int j = 0; j < 256; j++)
		{
			float result = 0;

			for (int k = -1; k < 2; k++)
			{
				for (int l = -1; l < 2; l++)
				{
					if (i + k < 0 || i + k >= 256 || j + l < 0 || j + l >= 256)
					{
						result += 0;
					}

					else
					{
						baseImage.getPixelVal(i + k, j + l, temp);
						result += temp * sobel_mask[3 / 2 - k][3 / 2 - l];
					}
				}
			}
			pixels[i * 256 + j] = result;
		}
	}

	// Normalization
	float max = INT_MIN;
	float min = INT_MAX;

	for (int i = 0; i < 256; i++)
	{
		for (int j = 0; j < 256; j++)
		{
			if (pixels[i * 256 + j] > max)
			{
				max = pixels[i * 256 + j];
			}
			if (pixels[i * 256 + j] < min)
			{
				min = pixels[i * 256 + j];
			}
		}
	}

	ImageType finalImage_spatial(256, 256, 255);
	for (int i = 0; i < 256; i++)
	{
		for (int j = 0; j < 256; j++)
		{
			temp = (int)(255.0 * (((double)pixels[i * 256 + j] - min) / (double)(max - min)));
			finalImage_spatial.setPixelVal(i, j, temp);
		}
	}

	// Getting Spectrum
	double **real_fuv = new double *[256];
	for (int i = 0; i < 256; i++)
	{
		real_fuv[i] = new double[256];
	}
	double **image_fuv = new double *[256];
	for (int i = 0; i < 256; i++)
	{
		image_fuv[i] = new double[256];
	}

	double **real_fuv2 = new double *[256];
	for (int i = 0; i < 256; i++)
	{
		real_fuv2[i] = new double[256];
	}
	double **image_fuv2 = new double *[256];
	for (int i = 0; i < 256; i++)
	{
		image_fuv2[i] = new double[256];
	}

	for (int i = 0; i < 256; i++)
	{
		for (int j = 0; j < 256; j++)
		{
			baseImage.getPixelVal(i, j, temp);
			real_fuv[i][j] = temp;
			image_fuv[i][j] = 0;

			finalImage_spatial.getPixelVal(i, j, temp);
			real_fuv2[i][j] = temp;
			image_fuv2[i][j] = 0;
		}
	}

	// Shift the spectrum
	for (int i = 0; i < 256; i++)
	{
		for (int j = 0; j < 256; j++)
		{
			real_fuv[i][j] = real_fuv[i][j] * pow(-1, i + j);
			image_fuv[i][j] = image_fuv[i][j] * pow(-1, i + j);

			real_fuv2[i][j] = real_fuv2[i][j] * pow(-1, i + j);
			image_fuv2[i][j] = image_fuv2[i][j] * pow(-1, i + j);
		}
	}
	// Call 2DFFT
	fft2d(256, 256, real_fuv, image_fuv, -1);
	fft2d(256, 256, real_fuv2, image_fuv2, -1);

	visualizespectrum(real_fuv, image_fuv, 256, 256, "spatial_spectrum_before.pgm");
	visualizespectrum(real_fuv2, image_fuv2, 256, 256, "spatial_spectrum_after.pgm");

	char finImage_spatial[] = "part2_spatial.pgm";
	writeImage(finImage_spatial, finalImage_spatial);

	for (int i = 0; i < 256; ++i)
	{
		delete[] real_fuv[i];
		delete[] real_fuv2[i];
	}
	delete[] real_fuv;
	delete[] real_fuv2;
	for (int i = 0; i < 256; ++i)
	{
		delete[] image_fuv[i];
		delete[] image_fuv2[i];
	}
	delete[] image_fuv;
	delete[] image_fuv2;

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	// Frequency Filter, image itself
	double **real_fuv3 = new double *[512];
	for (int i = 0; i < 512; i++)
	{
		real_fuv3[i] = new double[512];
	}
	double **image_fuv3 = new double *[512];
	for (int i = 0; i < 512; i++)
	{
		image_fuv3[i] = new double[512];
	}

	for (int i = 0; i < 512; i++)
	{
		for (int j = 0; j < 512; j++)
		{
			if(i < 256 && j < 256){
				baseImage.getPixelVal(i, j, temp);
				real_fuv3[i][j] = temp * pow(-1, i + j);
				image_fuv3[i][j] = 0;
			}
			else {
				real_fuv3[i][j] = 0;
				image_fuv3[i][j] = 0;
			}
		}
	}

	fft2d(512, 512, real_fuv3, image_fuv3, -1);

	visualizespectrum(real_fuv3, image_fuv3, 512, 512, "frequency_spectrum_before.pgm");

	double **real_huv = new double *[512];
	for (int i = 0; i < 512; i++)
	{
		real_huv[i] = new double[512];
	}
	double **image_huv = new double *[512];
	for (int i = 0; i < 512; i++)
	{
		image_huv[i] = new double[512];
	}

	for (int i = 0; i < 512; i++){
		for (int j = 0; j < 512; j++){
			real_huv[i][j] = 0;
			image_huv[i][j] = 0;
		}
	}
		

	// Initialize mask inside center
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			int x = 512 / 2 - 1 + i;
			int y = 512 / 2 - 1 + j;

			real_huv[x][y] = sobel_mask[i][j];
			image_huv[x][y] = 0;
		}
	}
	for (int i = 0; i < 512; i++){
		for (int j = 0; j < 512; j++){
			real_huv[i][j] *= pow(-1, i + j);
		}
	}

	ImageType finImage_frequency(256, 256, 255);
	char finImage_freq[] = "part2_frequency.pgm";

	fft2d(512, 512, real_huv, image_huv, -1);
	for (int i = 0; i < 512; i++){
		for (int j = 0; j < 512; j++){
			real_huv[i][j] = 0;
			image_huv[i][j] *= pow(-1, i + j);
		}
	}

	visualizespectrum(real_huv, image_huv, 512, 512, "sobel_spectrum.pgm");

	// Element-wise complex multiplication
	for (int i = 0; i < 512; i++)
	{
		for (int j = 0; j < 512; j++)
		{
			// set real part to zero, undo cenetering of sobel mask
			// mask_transform[index] = std::complex<float>(0, mask_transform[index].imag()) * (float)pow(-1, i + j);

			// multiply F(x,y) and H(x,y)
			real_fuv3[i][j] = real_huv[i][j] * real_fuv3[i][j] - image_huv[i][j] * image_fuv3[i][j];
			image_fuv3[i][j] = image_huv[i][j] * real_fuv3[i][j] + real_huv[i][j] * image_fuv3[i][j];
		}
	}
	visualizespectrum(real_fuv3, image_fuv3, 512, 512, "frequency_spectrum_after.pgm");

	fft2d(512, 512, real_fuv3, image_fuv3, 1);

	for (int i = 0; i < 512; i++){
		for (int j = 0; j < 512; j++){
			real_fuv3[i][j] *= pow(-1, i + j);
		}
	}

	double rmax2 = real_fuv3[0][0];
	double rmin2 = real_fuv3[0][0];
	for (int i = 0; i < 256; i++)
	{
		for (int j = 0; j < 256; j++)
		{
			if (real_fuv3[i][j] > rmax2)
			{
				rmax2 = real_fuv3[i][j];
			}
			if (real_fuv3[i][j] < rmin2)
			{
				rmin2 = real_fuv3[i][j];
			}
		}
	}
	for (int i = 0; i < 256; i++)
	{
		for (int j = 0; j < 256; j++)
		{
			real_fuv3[i][j] = 255 * (real_fuv3[i][j] - rmin2) / (rmax2 - rmin2);
		}
	}

	for (int i = 0; i < 256; i++)
	{
		for (int j = 0; j < 256; j++)
		{
			temp = real_fuv3[i][j];
			finImage_frequency.setPixelVal(i, j, temp);
		}
	}

	writeImage(finImage_freq, finImage_frequency);

	for (int i = 0; i < 512; ++i)
	{
		delete[] real_fuv3[i];
		delete[] real_huv[i];
	}
	delete[] real_fuv3;
	delete[] real_huv;
	for (int i = 0; i < 512; ++i)
	{
		delete[] image_fuv3[i];
		delete[] image_huv[i];
	}
	delete[] image_fuv3;
	delete[] image_huv;
}

void degrade(char fname[], float mu, float sigma, double ** real_fuv, double ** image_fuv){
	// Step 1, FT
	double **copyr = new double * [512];
	for (int i = 0; i < 512; i++)
	{
		copyr[i] = new double [512];
	}
	double **copyi = new double * [512];
	for (int i = 0; i < 512; i++)
	{
		copyi[i] = new double [512];
	}

	double **noiser = new double * [512];
	for (int i = 0; i < 512; i++)
	{
		noiser[i] = new double [512];
	}
	double **noisei = new double * [512];
	for (int i = 0; i < 512; i++)
	{
		noisei[i] = new double [512];
	}

	//noise
	for (int i = 0; i < 512; i++)
	{
		for (int j = 0; j < 512; j++)
		{
			noiser[i][j] = generateGaussianNoise(sigma * sigma);
			noisei[i][j] = generateGaussianNoise(sigma * sigma);
		}
	}

	fft2d (512, 512, real_fuv, image_fuv, -1);
	// fft2d (512, 512, noiser, noisei, -1);

	visualizespectrum(real_fuv, image_fuv, 512, 512, "lenna_spec_before.pgm");

	// Step 3 Apply H(u,v) and N(u,v)
	double a = .1, b = .1, T = 1;
	for (int i = 0; i < 512; i++)
	{
		for (int j = 0; j < 512; j++)
		{
			int i_adj = i - 512 / 2, j_adj = j - 512 / 2;
			double uavb = (a * i_adj + b * j_adj) * M_PI;
			//H(u,v) = (T / uavb) * sin(uavb) * (cos(uavb) - jsin(uavb)), therefor:
			double real_huv = (T / uavb) * sin(uavb) * cos(uavb);
			double image_huv = (T / uavb) * sin(uavb) * -1 * sin(uavb);
			//none of this works at (256, 256), use limit as uavb approaches 0
			if (i_adj == -j_adj){
				real_huv = 1;
				image_huv = 0;
			}
			complex<double> f(real_fuv[i][j], image_fuv[i][j]);
			complex<double> h(real_huv, image_huv);
			complex<double> fh = f * h;
			real_fuv[i][j] = fh.real();
			image_fuv[i][j] = fh.imag();
			real_fuv[i][j] = real_fuv[i][j] + noiser[i][j];
			image_fuv[i][j] = image_fuv[i][j] + noisei[i][j];
			
		}
	}

	visualizespectrum(real_fuv, image_fuv, 512, 512, "lenna_spec_after.pgm");

	for (int i = 0; i < 512; i++)
	{
		for (int j = 0; j < 512; j++)
		{
			copyr[i][j] = real_fuv[i][j];
			copyi[i][j] = image_fuv[i][j];
		}
	}

	// Step 4 Inverse FT
	fft2d (512, 512, copyr, copyi, 1);

	//Step 5 uncenter
	ImageType finalImage(256, 256, 255);
	for (int i = 0; i < 256; i++)
	{
		for (int j = 0; j < 256; j++)
		{
			copyr[i][j] *= pow(-1, i + j);
		}
	}

	int max = copyr[0][0];
	int min = copyr[0][0];
	for (int i = 0; i < 256; i++)
	{
		for (int j = 0; j < 256; j++)
		{
			if (copyr[i][j] > max)
			{
				max = copyr[i][j];
			}
			if (copyr[i][j] < min)
			{
				min = copyr[i][j];
			}
		}
	}
	for (int i = 0; i < 256; i++)
	{
		for (int j = 0; j < 256; j++)
		{
			copyr[i][j] = 255 * (copyr[i][j] - min) / (max - min);
		}
	}

	for (int i = 0; i < 256; i++)
	{
		for (int j = 0; j < 256; j++)
		{
			finalImage.setPixelVal(i, j, copyr[i][j]);
		}
	}

	//print and clean up
	char finImage[] = "degraded.pgm";
	writeImage(finImage, finalImage);

	for (int i = 0; i < 512; ++i)
	{
		delete[] copyr[i];
	}
	delete[] copyr;
	for (int i = 0; i < 512; ++i)
	{
		delete[] copyi[i];
	}
	delete[] copyi;
	for (int i = 0; i < 512; ++i)
	{
		delete[] noiser[i];
	}
	delete[] noiser;
	for (int i = 0; i < 512; ++i)
	{
		delete[] noisei[i];
	}
	delete[] noisei;
}
void inversefilter(char fname[], float mu, float sigma, float radius){
	char degImage[] = "degraded.pgm";
	ImageType degradedImage(256, 256, 255);
	ImageType paddedImage(512, 512, 255);
	readImage(fname, degradedImage);

	// Padding base image
	int temp;
	for (int i = 0; i < 256; i++)
	{
		for (int j = 0; j < 256; j++)
		{

			degradedImage.getPixelVal(i, j, temp);
			paddedImage.setPixelVal(i, j, temp);
		}
	}
	
	// Step 1, FT
	double **real_fuv = new double * [512];
	for (int i = 0; i < 512; i++)
	{
		real_fuv[i] = new double [512];
	}
	double **image_fuv = new double * [512];
	for (int i = 0; i < 512; i++)
	{
		image_fuv[i] = new double [512];
	}
	double **real_huv = new double * [512];
	for (int i = 0; i < 512; i++)
	{
		real_huv[i] = new double [512];
	}
	double **image_huv = new double * [512];
	for (int i = 0; i < 512; i++)
	{
		image_huv[i] = new double [512];
	}

	//center
	for (int i = 0; i < 512; i++)
	{
		for (int j = 0; j < 512; j++)
		{
			paddedImage.getPixelVal(i, j, temp);
			real_fuv[i][j] = temp * pow(-1, i + j);
			image_fuv[i][j] = 0;
		}
	}
	degrade(fname, mu, sigma, real_fuv, image_fuv);

	// fft2d (512, 512, real_fuv, image_fuv, -1);

	// complex division by H(u,v)
	double a = .1, b = .1, T = 1;
	for (int i = 0; i < 512; i++)
	{
		for (int j = 0; j < 512; j++)
		{
			int i_adj = i - 512 / 2, j_adj = j - 512 / 2;
			double uavb = (a * i_adj + b * j_adj) * M_PI;
			//H(u,v) = (T / uavb) * sin(uavb) * (cos(uavb) - jsin(uavb)), therefor:
			real_huv[i][j] = (T / uavb) * sin(uavb) * cos(uavb);
			image_huv[i][j] = (T / uavb) * sin(uavb) * -1 * sin(uavb);
			//none of this works at (256, 256), use limit as uavb approaches 0
			if (i_adj == -j_adj){
				real_huv[i][j] = 1;
				image_huv[i][j] = 0;
			}
			complex<double> f(real_fuv[i][j], image_fuv[i][j]);
			complex<double> h(real_huv[i][j], image_huv[i][j]);
			complex<double> fh = f / h;
			if(sqrt(i_adj * i_adj + j_adj * j_adj) <= radius){
				// real_fuv[i][j] = (real_fuv[i][j] * real_huv[i][j] + image_fuv[i][j] * image_huv[i][j]) / (real_huv[i][j] * real_huv[i][j] + image_huv[i][j] * image_huv[i][j]);
				// image_fuv[i][j] = (image_fuv[i][j] * real_huv[i][j] - real_fuv[i][j] * image_huv[i][j]) / (real_huv[i][j] * real_huv[i][j] + image_huv[i][j] * image_huv[i][j]);
				real_fuv[i][j] = fh.real();
				image_fuv[i][j] = fh.imag();
			}
		}
	}
	visualizespectrum(real_huv, image_huv, 512, 512, "h_spectrum.pgm");
	visualizespectrum(real_fuv, image_fuv, 512, 512, "lennainversespec.pgm");

	// Step 4 Inverse FT
	fft2d (512, 512, real_fuv, image_fuv, 1);

	// Step 5 uncenter
	ImageType finalImage(256, 256, 255);
	for (int i = 0; i < 256; i++)
	{
		for (int j = 0; j < 256; j++)
		{
			real_fuv[i][j] = real_fuv[i][j] * pow(-1, i + j);
		}
	}

	// Normalization
	float rmax = real_fuv[0][0];
	float rmin = real_fuv[0][0];
	for (int i = 0; i < 256; i++)
	{
		for (int j = 0; j < 256; j++)
		{
			if (real_fuv[i][j] > rmax)
			{
				rmax = real_fuv[i][j];
			}
			if (real_fuv[i][j] < rmin)
			{
				rmin = real_fuv[i][j];
			}
		}
	}
	for (int i = 0; i < 256; i++)
	{
		for (int j = 0; j < 256; j++)
		{
			real_fuv[i][j] = 255 * (real_fuv[i][j] - rmin) / (rmax - rmin);
			finalImage.setPixelVal(i, j, real_fuv[i][j]);
		}
	}


	// for (int i = 0; i < 512; i++)
	// {
	// 	for (int j = 0; j < 512; j++)
	// 	{
	// 		cout << real_fuv[i][j] << " " << image_fuv[i][j] << endl;
	// 	}
	// }

	//print and clean up
	char finImage[] = "inversefiltered.pgm";
	writeImage(finImage, finalImage);

	for (int i = 0; i < 512; ++i)
	{
		delete[] real_fuv[i];
	}
	delete[] real_fuv;
	for (int i = 0; i < 512; ++i)
	{
		delete[] image_fuv[i];
	}
	delete[] image_fuv;
	for (int i = 0; i < 512; ++i)
	{
		delete[] real_huv[i];
	}
	delete[] real_huv;
	for (int i = 0; i < 512; ++i)
	{
		delete[] image_huv[i];
	}
	delete[] image_huv;

}
void wienerfilter(char fname[], float mu, float sigma, float k){
char degImage[] = "degraded.pgm";
	ImageType degradedImage(256, 256, 255);
	ImageType paddedImage(512, 512, 255);
	readImage(fname, degradedImage);

	// Padding base image
	int temp;
	for (int i = 0; i < 256; i++)
	{
		for (int j = 0; j < 256; j++)
		{

			degradedImage.getPixelVal(i, j, temp);
			paddedImage.setPixelVal(i, j, temp);
		}
	}
	
	// Step 1, FT
	double **real_fuv = new double * [512];
	for (int i = 0; i < 512; i++)
	{
		real_fuv[i] = new double [512];
	}
	double **image_fuv = new double * [512];
	for (int i = 0; i < 512; i++)
	{
		image_fuv[i] = new double [512];
	}
	double **real_huv = new double * [512];
	for (int i = 0; i < 512; i++)
	{
		real_huv[i] = new double [512];
	}
	double **image_huv = new double * [512];
	for (int i = 0; i < 512; i++)
	{
		image_huv[i] = new double [512];
	}

	//center
	for (int i = 0; i < 512; i++)
	{
		for (int j = 0; j < 512; j++)
		{
			paddedImage.getPixelVal(i, j, temp);
			real_fuv[i][j] = temp * pow(-1, i + j);
			image_fuv[i][j] = 0;
		}
	}
	degrade(fname, mu, sigma, real_fuv, image_fuv);

	// fft2d (512, 512, real_fuv, image_fuv, -1);

	// complex division by H(u,v)
	double a = .1, b = .1, T = 1;
	for (int i = 0; i < 512; i++)
	{
		for (int j = 0; j < 512; j++)
		{
			int i_adj = i - 512 / 2, j_adj = j - 512 / 2;
			double uavb = (a * i_adj + b * j_adj) * M_PI;
			//H(u,v) = (T / uavb) * sin(uavb) * (cos(uavb) - jsin(uavb)), therefor:
			real_huv[i][j] = (T / uavb) * sin(uavb) * cos(uavb);
			image_huv[i][j] = (T / uavb) * sin(uavb) * -1 * sin(uavb);
			//none of this works at (256, 256), use limit as uavb approaches 0
			if (i_adj == -j_adj){
				real_huv[i][j] = 1;
				image_huv[i][j] = 0;
			}
			double hpower = (real_huv[i][j] * real_huv[i][j] + image_huv[i][j] * image_huv[i][j]);
			hpower = hpower / (hpower + k);
			complex<double> f(real_fuv[i][j], image_fuv[i][j]);
			complex<double> h(real_huv[i][j], image_huv[i][j]);
			complex<double> fh = f / h;
				// real_fuv[i][j] = (real_fuv[i][j] * real_huv[i][j] + image_fuv[i][j] * image_huv[i][j]) / (real_huv[i][j] * real_huv[i][j] + image_huv[i][j] * image_huv[i][j]);
				// image_fuv[i][j] = (image_fuv[i][j] * real_huv[i][j] - real_fuv[i][j] * image_huv[i][j]) / (real_huv[i][j] * real_huv[i][j] + image_huv[i][j] * image_huv[i][j]);
			real_fuv[i][j] = fh.real();
			image_fuv[i][j] = fh.imag();
			real_fuv[i][j] *= hpower;
			image_fuv[i][j] *= hpower;

		}
	}
	visualizespectrum(real_huv, image_huv, 512, 512, "h_spectrum.pgm");
	visualizespectrum(real_fuv, image_fuv, 512, 512, "lennainversespec.pgm");

	// Step 4 Inverse FT
	fft2d (512, 512, real_fuv, image_fuv, 1);

	// Step 5 uncenter
	ImageType finalImage(256, 256, 255);
	for (int i = 0; i < 256; i++)
	{
		for (int j = 0; j < 256; j++)
		{
			real_fuv[i][j] = real_fuv[i][j] * pow(-1, i + j);
		}
	}

	// Normalization
	float rmax = real_fuv[0][0];
	float rmin = real_fuv[0][0];
	for (int i = 0; i < 256; i++)
	{
		for (int j = 0; j < 256; j++)
		{
			if (real_fuv[i][j] > rmax)
			{
				rmax = real_fuv[i][j];
			}
			if (real_fuv[i][j] < rmin)
			{
				rmin = real_fuv[i][j];
			}
		}
	}
	for (int i = 0; i < 256; i++)
	{
		for (int j = 0; j < 256; j++)
		{
			real_fuv[i][j] = 255 * (real_fuv[i][j] - rmin) / (rmax - rmin);
			finalImage.setPixelVal(i, j, real_fuv[i][j]);
		}
	}


	// for (int i = 0; i < 512; i++)
	// {
	// 	for (int j = 0; j < 512; j++)
	// 	{
	// 		cout << real_fuv[i][j] << " " << image_fuv[i][j] << endl;
	// 	}
	// }

	//print and clean up
	char finImage[] = "wienerfiltered.pgm";
	writeImage(finImage, finalImage);

	for (int i = 0; i < 512; ++i)
	{
		delete[] real_fuv[i];
	}
	delete[] real_fuv;
	for (int i = 0; i < 512; ++i)
	{
		delete[] image_fuv[i];
	}
	delete[] image_fuv;
	for (int i = 0; i < 512; ++i)
	{
		delete[] real_huv[i];
	}
	delete[] real_huv;
	for (int i = 0; i < 512; ++i)
	{
		delete[] image_huv[i];
	}
	delete[] image_huv;

}

void experiment4(char fname[], float gh, float gl)
{
	ImageType baseImage(512, 512, 255);
	ImageType paddedImage(1024, 1024, 255);
	readImage(fname, baseImage);

	// Padding base image
	int temp;
	double temp2;
	for (int i = 0; i < 512; i++)
	{
		for (int j = 0; j < 512; j++)
		{

			baseImage.getPixelVal(i, j, temp);
			paddedImage.setPixelVal(i, j, temp);
		}
	}
	for (int i = 512; i < 1024; i++)
	{
		for (int j = 512; j < 1024; j++)
		{
			paddedImage.setPixelVal(i, j, 0);
		}
	}

	// Step 1 and 2, ln and FT
	double **real_fuv = new double *[1024];
	for (int i = 0; i < 1024; i++)
	{
		real_fuv[i] = new double[1024];
	}
	double **image_fuv = new double *[1024];
	for (int i = 0; i < 1024; i++)
	{
		image_fuv[i] = new double[1024];
	}

	for (int i = 0; i < 1024; i++)
	{
		for (int j = 0; j < 1024; j++)
		{
			paddedImage.getPixelVal(i, j, temp);
			temp2 = log(temp + 1);
			real_fuv[i][j] = temp2;
			image_fuv[i][j] = 0;
		}
	}

	// Shift the spectrum
	for (int i = 0; i < 1024; i++)
	{
		for (int j = 0; j < 1024; j++)
		{
			real_fuv[i][j] = real_fuv[i][j] * pow(-1, i + j);
			image_fuv[i][j] = image_fuv[i][j] * pow(-1, i + j);
		}
	}

	fft2d(1024, 1024, real_fuv, image_fuv, -1);

	// Step 3 Apply H(u,v)
	float c = 1, D_0 = 1.8;
	float coef = -c / (D_0 * D_0);
	float gammaDiff = gh - gl;
	for (int i = 0; i < 1024; i++)
	{
		for (int j = 0; j < 1024; j++)
		{
			int i_adj = i - 1024 / 2, j_adj = j - 1024 / 2;
			real_fuv[i][j] *= gammaDiff * (1 - exp(coef * (i_adj * i_adj + j_adj * j_adj))) + gl;
			image_fuv[i][j] *= gammaDiff * (1 - exp(coef * (i_adj * i_adj + j_adj * j_adj))) + gl;
		}
	}

	// Step 4 Inverse FT
	fft2d(1024, 1024, real_fuv, image_fuv, 1);

	// Step 5 Take exp
	ImageType finalImage(512, 512, 255);
	for (int i = 0; i < 512; i++)
	{
		for (int j = 0; j < 512; j++)
		{
			temp2 = real_fuv[i][j] * pow(-1, i + j);
			temp2 = exp(temp2) - 1;
			finalImage.setPixelVal(i, j, temp2);
		}
	}

	// Normalization
	finalImage.getPixelVal(0, 0, temp);
	float rmax = temp;
	float rmin = temp;
	for (int i = 0; i < 512; i++)
	{
		for (int j = 0; j < 512; j++)
		{
			finalImage.getPixelVal(i, j, temp);
			if (temp > rmax)
			{
				rmax = temp;
			}
			if (temp < rmin)
			{
				rmin = temp;
			}
		}
	}
	for (int i = 0; i < 512; i++)
	{
		for (int j = 0; j < 512; j++)
		{
			finalImage.getPixelVal(i, j, temp);
			temp = 255 * (temp - rmin) / (rmax - rmin);
			finalImage.setPixelVal(i, j, temp);
		}
	}

	//print and clean up
	char finImage[] = "part4_blurred.pgm";
	writeImage(finImage, finalImage);

	for (int i = 0; i < 1024; ++i)
	{
		delete[] real_fuv[i];
	}
	delete[] real_fuv;
	for (int i = 0; i < 1024; ++i)
	{
		delete[] image_fuv[i];
	}
	delete[] image_fuv;
}

void test2dfft()
{
	char lenna[] = "lenna.pgm";
	ImageType testimage(256, 256, 255);
	readImage(lenna, testimage);
	// double real_fuv[256][256];
	// double image_fuv[256][256];

	double **real_fuv = new double *[256];
	for (int i = 0; i < 256; i++)
	{
		real_fuv[i] = new double[256];
	}
	double **image_fuv = new double *[256];
	for (int i = 0; i < 256; i++)
	{
		image_fuv[i] = new double[256];
	}

	int temp;
	for (int i = 0; i < 256; i++)
	{
		for (int j = 0; j < 256; j++)
		{
			testimage.getPixelVal(i, j, temp);
			real_fuv[i][j] = temp;
			image_fuv[i][j] = 0;
		}
	}

	// Call 2DFFT
	fft2d(256, 256, real_fuv, image_fuv, -1);
	fft2d(256, 256, real_fuv, image_fuv, 1);

	for (int i = 0; i < 256; ++i)
	{
		delete[] real_fuv[i];
	}
	delete[] real_fuv;
	for (int i = 0; i < 256; ++i)
	{
		delete[] image_fuv[i];
	}
	delete[] image_fuv;
}

void imgtodata(ImageType &image, float **data)
{
	int temp = 0;
	int N, M, L;
	image.getImageInfo(N, M, L);
	for (int i = 0; i < N; i++)
	{
		for (int j = 1; j < M * 2 + 1; j += 2)
		{
			image.getPixelVal(i, (j - 1) / 2, temp);
			data[i][j] = temp;
			data[i][j + 1] = 0;
		}
	}
}

void datatoimg(ImageType &image, float **data)
{
	int temp = 0;
	int N, M, L;
	image.getImageInfo(N, M, L);
	for (int i = 0; i < N; i++)
	{
		for (int j = 1; j < M * 2 + 1; j += 2)
		{
			temp = data[i][j];
			image.setPixelVal(i, (j - 1) / 2, temp);
		}
	}
}

/*

   The real part is stored in the odd locations of the array
   (data[1], data[3], data[5]. etc) and the imaginary part
   in the even locations (data[2], data[4], data[6], etc.

   The elements in the array data must be stored from data[1]
   to data[2*nn] -  data[0] is not used!

   nn is the length of the input which should be power of 2.
   Warning: the program does not check this condition.

   isign: -1 Forward FFT, isign: 1  Inverse FFT (based on our definition)

   Warning: the FFT routine provided does not multiply by the normalization
   factor 1/N that appears in the forward DFT equation; you should do this
   yourself (see page 506 from the fft handout).

*/

void fft(float data[], unsigned long nn, int isign)
{
	unsigned long n, mmax, m, j, istep, i;
	double wtemp, wr, wpr, wpi, wi, theta;
	float tempr, tempi;

	n = nn << 1;
	j = 1;
	for (i = 1; i < n; i += 2)
	{
		if (j > i)
		{
			SWAP(data[j], data[i]);
			SWAP(data[j + 1], data[i + 1]);
		}
		m = n >> 1;
		while (m >= 2 && j > m)
		{
			j -= m;
			m >>= 1;
		}
		j += m;
	}
	mmax = 2;
	while (n > mmax)
	{
		istep = mmax << 1;
		theta = isign * (6.28318530717959 / mmax);
		wtemp = sin(0.5 * theta);
		wpr = -2.0 * wtemp * wtemp;
		wpi = sin(theta);
		wr = 1.0;
		wi = 0.0;
		for (m = 1; m < mmax; m += 2)
		{
			for (i = m; i <= n; i += istep)
			{
				j = i + mmax;
				tempr = wr * data[j] - wi * data[j + 1];
				tempi = wr * data[j + 1] + wi * data[j];
				data[j] = data[i] - tempr;
				data[j + 1] = data[i + 1] - tempi;
				data[i] += tempr;
				data[i + 1] += tempi;
			}
			wr = (wtemp = wr) * wpr - wi * wpi + wr;
			wi = wi * wpr + wtemp * wpi + wi;
		}
		mmax = istep;
	}
}
#undef SWAP

// void fft2d(float **data, int size, unsigned long nn, int isign)
void fft2d(int N, int M, double ** real_fuv, double ** image_fuv, int isign)
{
	float **data = new float *[N];

	for (int i = 0; i < N; i++)
	{
		data[i] = new float[M * 2 + 1];
	}

	// Initialize data with values
	for (int i = 0; i < N; i++)
	{
		for (int j = 1; j < M * 2 + 1; j += 2)
		{
			data[i][j] = real_fuv[i][(j - 1) / 2];
			data[i][j + 1] = image_fuv[i][(j - 1) / 2];
		}
	}

	//create dynamic 2d arrays
	float** colsdata = new float*[M];
	for(int i = 0; i < M; i++){
    	colsdata[i] = new float[N * 2 + 1];
	}
	float** rowsdata = new float*[N];
	for(int i = 0; i < N; i++){
    	rowsdata[i] = new float[M * 2 + 1];
	}

	// copy the image into a 2d array, flipped so the columns become the rows and vice versa
	for (int i = 0; i < M; i++)
	{
		for (int j = 1; j < N * 2 + 1; j = j + 2)
		{
			colsdata[i][j] = data[(j - 1) / 2][2 * i + 1];
			colsdata[i][j + 1] = data[(j - 1) / 2][2 * i + 1 + 1];
		}
	}

	// perform 1d fft on each column
	for (int i = 0; i < M; i++)
	{
		fft(colsdata[i], M, isign);
	}

	// flip back to the way it should be so we can transform the rows
	for (int i = 0; i < N; i++)
	{
		for (int j = 1; j < M * 2 + 1; j = j + 2)
		{
			rowsdata[i][j] = colsdata[(j - 1) / 2][2 * i + 1];
			rowsdata[i][j + 1] = colsdata[(j - 1) / 2][2 * i + 1 + 1];
		}
	}

	// perform id fft on each row
	for (int i = 0; i < N; i++)
	{
		fft(rowsdata[i], N, isign);
	}

	// copy to original array and do 1/N
	for (int i = 0; i < N; i++)
	{
		for (int j = 1; j < M * 2 + 1; j = j + 2)
		{
			data[i][j] = rowsdata[i][j] / N;
			data[i][j + 1] = rowsdata[i][j + 1] / N;
			// cout << data[i][j] << " " << data[i][j + 1] << endl;
		}
	}

	// copy from data to original arrays
	for (int i = 0; i < N; i++)
	{
		for (int j = 1; j < N * 2 + 1; j = j + 2)
		{
			real_fuv[i][(j - 1) / 2] = data[i][j];
			image_fuv[i][(j - 1) / 2] = data[i][j + 1];
			// cout << data[i][j] << " " << data[i][j + 1] << endl;
		}
	}

	// Checking if the inverse worked
	/*if (isign == -1)
	{
		ImageType final(N, M, 255);
		datatoimg(final, data);
		char fft2dtest[] = "test2dfft.pgm";
		writeImage(fft2dtest, final);
	}*/

	//delete the dynamic 2d arrays
	for(int i = 0; i < N; ++i){
    	delete[] colsdata[i];
	}
	delete[] colsdata;
	for(int i = 0; i < M; ++i){
    	delete[] rowsdata[i];
	}
	delete[] rowsdata;
	for(int i = 0; i < N; ++i){
    	delete[] data[i];
	}
	delete[] data;
}

float box_muller(float m, float s) /* normal random variate generator */
{ /* mean m, standard deviation s */
default_random_engine generator;
uniform_real_distribution<double> distribution(0.0,1.0);
float x1, x2, w, y1;
static float y2;
static int use_last = 0;

if (use_last) /* use value from previous call */
{
y1 = y2;
use_last = 0;
}
else
{
do {
x1 = 2.0 * distribution(generator) - 1.0;
x2 = 2.0 * distribution(generator) - 1.0;
w = x1 * x1 + x2 * x2;
} while ( w >= 1.0 );

w = sqrt( (-2.0 * log( w ) ) / w );
y1 = x1 * w;
y2 = x2 * w;
use_last = 1;
}

return( m + y1 * s );
}

double generateGaussianNoise(const double &variance)
{
    static bool haveSpare = false;
    static double rand1, rand2;
  
    if(haveSpare)
    {
        haveSpare = false;
        return sqrt(variance * rand1) * sin(rand2);
    }
  
    haveSpare = true;
  
    rand1 = rand() / ((double) RAND_MAX);
    if(rand1 < 1e-100) rand1 = 1e-100;
    rand1 = -2 * log(rand1);
    rand2 = (rand() / ((double) RAND_MAX)) * 2 * M_PI;
  
    return sqrt(variance * rand1) * cos(rand2);
}