#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <algorithm>
#include <climits>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include <math.h>
#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr
#define M_PI 3.14159265358979323846 /* pi */

#include "image.h"

using namespace std;

void readImage(char fname[], ImageType &image);
void writeImage(char fname[], ImageType &image);
void fft(float data[], unsigned long nn, int isign);

void test2dfft();
void fft2d(int N, int M, double **real_fuv, double **image_fuv, int isign);
void imgtodata(ImageType &image, float **data);
void datatoimg(ImageType &image, float **data);

void experiment2(char fname[]);
void experiment4(char fname[]);

const float sobel_mask[3][3] = {{-1, 0, 1}, {-2, 0, 2}, {-1, 0, 1}};

int main(int argc, char *argv[])
{
	char girl[] = "girl.pgm";
	char boy[] = "boy_noisy.pgm";
	char lenna[] = "lenna.pgm";

	experiment2(lenna);
	// experiment4(girl);

	return 0;
}

void experiment2(char fname[])
{
	ImageType baseImage(512, 512, 255);
	readImage(fname, baseImage);

	// Spatial Filter
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
	char finImage_spatial[] = "part2_spatial.pgm";
	for (int i = 0; i < 256; i++)
	{
		for (int j = 0; j < 256; j++)
		{
			temp = (int)(255.0 * (((double)pixels[i * 256 + j] - min) / (double)(max - min)));
			finalImage_spatial.setPixelVal(i, j, temp);
		}
	}
	writeImage(finImage_spatial, finalImage_spatial);

	// Frequency Filter
	// double **real_fuv = new double *[256];
	// for (int i = 0; i < 256; i++)
	// {
	// 	real_fuv[i] = new double[256 * 2 + 1];
	// }
	// double **image_fuv = new double *[256];
	// for (int i = 0; i < 256; i++)
	// {
	// 	image_fuv[i] = new double[256 * 2 + 1];
	// }

	// for (int i = 0; i < 256; i++)
	// {
	// 	for (int j = 0; j < 256; j++)
	// 	{
	// 		baseImage.getPixelVal(i, j, temp);
	// 		real_fuv[i][j] = temp;
	// 		image_fuv[i][j] = 0;
	// 	}
	// }

	// fft2d(256, 256, real_fuv, image_fuv, 1);

	// double **real_fuv2 = new double *[256];
	// for (int i = 0; i < 256; i++)
	// {
	// 	real_fuv2[i] = new double[256 * 2 + 1];
	// }
	// double **image_fuv2 = new double *[256];
	// for (int i = 0; i < 256; i++)
	// {
	// 	image_fuv2[i] = new double[256 * 2 + 1];
	// }

	// 	//Initialize mask inside
	// for (int i = 0; i < 3; i++)
	// {
	// 	for (int j = 0; j < 3; j++)
	// 	{
	// 		int x = 256 / 2 - 1 + i;
	// 		int y = 256 / 2 - 1 + j;

	// 		real_fuv2[i][j] = sobel_mask[i][j];
	// 		image_fuv2[i][j] = 0;
	// 	}
	// }

	// fft2d(256, 256, real_fuv2, image_fuv2, 1);
}

void experiment4(char fname[])
{
	ImageType baseImage(512, 512, 255);
	ImageType paddedImage(1024, 1024, 255);
	readImage(fname, baseImage);

	// Padding base image
	int temp;
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

	// Step 1 Ln f(x,y)
	for (int i = 0; i < 1024; i++)
	{
		for (int j = 0; j < 1024; j++)
		{
			paddedImage.getPixelVal(i, j, temp);
			temp = log(temp);
			paddedImage.setPixelVal(i, j, temp);
		}
	}

	// Step 2 FT
	double **real_fuv = new double *[1024];
	for (int i = 0; i < 1024; i++)
	{
		real_fuv[i] = new double[1024 * 2 + 1];
	}
	double **image_fuv = new double *[1024];
	for (int i = 0; i < 1024; i++)
	{
		image_fuv[i] = new double[1024 * 2 + 1];
	}

	for (int i = 0; i < 1024; i++)
	{
		for (int j = 0; j < 1024; j++)
		{
			paddedImage.getPixelVal(i, j, temp);
			real_fuv[i][j] = temp;
			image_fuv[i][j] = 0;
		}
	}
	fft2d(1024, 1024, real_fuv, image_fuv, 1);

	// Shift the spectrum
	for (int i = 0; i < 1024; i++)
	{
		for (int j = 0; j < 1024; j++)
		{
			real_fuv[i][j] = real_fuv[i][j] * pow(-1, i + j);
			image_fuv[i][j] = image_fuv[i][j] * pow(-1, i + j);
		}
	}

	// Step 3 Apply H(u,v)
	float gamma_H = 1.5, gamma_L = .5, c = 1, D_0 = 1.8;
	float coef = -c / D_0 / D_0;
	float gammaDiff = gamma_H - gamma_L;
	for (int i = 0; i < 1024; i++)
	{
		for (int j = 0; j < 1024; j++)
		{
			int i_adj = i - 1024 / 2, j_adj = j - 1024 / 2;
			real_fuv[i][j] *= gammaDiff * (1 - exp(coef * (i_adj * i_adj + j_adj * j_adj))) + gamma_L;
			image_fuv[i][j] *= gammaDiff * (1 - exp(coef * (i_adj * i_adj + j_adj * j_adj))) + gamma_L;
		}
	}

	// Step 4 Inverse FT
	fft2d(1024, 1024, real_fuv, image_fuv, -1);

	// Normalization
	float rmax = real_fuv[0][0];
	float rmin = real_fuv[0][0];
	for (int i = 0; i < 1024; i++)
	{
		for (int j = 0; j < 1024; j++)
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
	for (int i = 0; i < 1024; i++)
	{
		for (int j = 0; j < 1024; j++)
		{
			real_fuv[i][j] = 255 * (real_fuv[i][j] - rmin) / (rmax - rmin);
		}
	}

	// Step 5 Take exp
	ImageType finalImage(512, 512, 255);
	int temp2;
	for (int i = 0; i < 512; i++)
	{
		for (int j = 0; j < 512; j++)
		{
			temp2 = real_fuv[i][j];
			temp2 = exp(temp2);
			finalImage.setPixelVal(i, j, temp2);
		}
	}

	char finImage[] = "part4.pgm";
	writeImage(finImage, finalImage);
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
		real_fuv[i] = new double[256 * 2 + 1];
	}
	double **image_fuv = new double *[256];
	for (int i = 0; i < 256; i++)
	{
		image_fuv[i] = new double[256 * 2 + 1];
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
	if (isign == 1)
	{
		ImageType final(N, M, 255);
		datatoimg(final, data);
		char fft2dtest[] = "test2dfft.pgm";
		writeImage(fft2dtest, final);
	}

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