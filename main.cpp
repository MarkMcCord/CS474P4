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
// void fft2d(float **data, int size, unsigned long nn, int isign);
void imgtodata(ImageType &image, float **data);
void datatoimg(ImageType &image, float **data);

int main(int argc, char *argv[])
{
	return 0;
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

/*
void part2(char fname[], int size)
{
	ImageType baseImage(512, 512, 255);
	readImage(fname, baseImage);

	double **real_fuv = new double *[512];
	for (int i = 0; i < 512; i++)
	{
		real_fuv[i] = new double[512 * 2 + 1];
	}
	double **image_fuv = new double *[512];
	for (int i = 0; i < 512; i++)
	{
		image_fuv[i] = new double[512 * 2 + 1];
	}

	int temp;
	for (int i = 0; i < 512; i++)
	{
		for (int j = 0; j < 512; j++)
		{
			baseImage.getPixelVal(i, j, temp);
			real_fuv[i][j] = temp;
			image_fuv[i][j] = 0;
		}
	}

	// FFT2D for no shift
	fft2d(512, 512, real_fuv, image_fuv, -1);
	// Spectrum no shift
	double spectrum[512][512];
	for (int i = 0; i < 512; i++)
	{
		for (int j = 0; j < 512; j++)
		{
			spectrum[i][j] = sqrt(real_fuv[i][j] * real_fuv[i][j] + image_fuv[i][j] + image_fuv[i][j]);
		}
	}
	// Inverse for shift
	fft2d(512, 512, real_fuv, image_fuv, 1);

	// Shift the spectrum
	for (int i = 0; i < 512; i++)
	{
		for (int j = 0; j < 512; j++)
		{
			real_fuv[i][j] = real_fuv[i][j] * pow(-1, i + j);
			image_fuv[i][j] = image_fuv[i][j] * pow(-1, i + j);
		}
	}

	// Call 2DFFT
	fft2d(512, 512, real_fuv, image_fuv, -1);

	double magnitude[512][512];
	// Calculate magnitude
	for (int i = 0; i < 512; i++)
	{
		for (int j = 0; j < 512; j++)
		{
			magnitude[i][j] = sqrt(real_fuv[i][j] * real_fuv[i][j] + image_fuv[i][j] + image_fuv[i][j]);
		}
	}

	// Apply Log for better visualization
	for (int i = 0; i < 512; i++)
	{
		for (int j = 0; j < 512; j++)
		{
			magnitude[i][j] = log(1 + magnitude[i][j]);
			spectrum[i][j] = log(1 + spectrum[i][j]);
		}
	}

	if (size == 32)
	{
		char tempname[] = "fft_32.pgm";
		char tempname2[] = "fft_unshifted_32.pgm";
		ImageType final(512, 512, 255);
		ImageType final2(512, 512, 255);

		// Normalization
		float rmax = magnitude[0][0];
		float rmin = magnitude[0][0];
		float rmax2 = spectrum[0][0];
		float rmin2 = spectrum[0][0];
		for (int i = 0; i < 512; i++)
		{
			for (int j = 0; j < 512; j++)
			{
				if (magnitude[i][j] > rmax)
				{
					rmax = magnitude[i][j];
				}
				if (magnitude[i][j] < rmin)
				{
					rmin = magnitude[i][j];
				}
				if (spectrum[i][j] > rmax2)
				{
					rmax2 = spectrum[i][j];
				}
				if (spectrum[i][j] < rmin2)
				{
					rmin2 = spectrum[i][j];
				}
			}
		}
		for (int i = 0; i < 512; i++)
		{
			for (int j = 0; j < 512; j++)
			{
				magnitude[i][j] = 255 * (magnitude[i][j] - rmin) / (rmax - rmin);
				spectrum[i][j] = 255 * (spectrum[i][j] - rmin2) / (rmax2 - rmin2);
			}
		}

		int temp2;
		int temp3;
		// Write to image
		for (int i = 0; i < 512; i++)
		{
			for (int j = 0; j < 512; j++)
			{
				temp3 = spectrum[i][j];
				temp2 = magnitude[i][j];
				final.setPixelVal(i, j, temp2);
				final2.setPixelVal(i, j, temp3);
			}
		}

		writeImage(tempname, final);
		writeImage(tempname2, final2);
	}
	else if (size == 64)
	{
		char tempname[] = "fft_64.pgm";
		char tempname2[] = "fft_unshifted_64.pgm";
		ImageType final(512, 512, 255);
		ImageType final2(512, 512, 255);

		// Normalization
		float rmax = magnitude[0][0];
		float rmin = magnitude[0][0];
		float rmax2 = spectrum[0][0];
		float rmin2 = spectrum[0][0];
		for (int i = 0; i < 512; i++)
		{
			for (int j = 0; j < 512; j++)
			{
				if (magnitude[i][j] > rmax)
				{
					rmax = magnitude[i][j];
				}
				if (magnitude[i][j] < rmin)
				{
					rmin = magnitude[i][j];
				}
				if (spectrum[i][j] > rmax2)
				{
					rmax2 = spectrum[i][j];
				}
				if (spectrum[i][j] < rmin2)
				{
					rmin2 = spectrum[i][j];
				}
			}
		}
		for (int i = 0; i < 512; i++)
		{
			for (int j = 0; j < 512; j++)
			{
				magnitude[i][j] = 255 * (magnitude[i][j] - rmin) / (rmax - rmin);
				spectrum[i][j] = 255 * (spectrum[i][j] - rmin2) / (rmax2 - rmin2);
			}
		}

		int temp2;
		int temp3;
		// Write to image
		for (int i = 0; i < 512; i++)
		{
			for (int j = 0; j < 512; j++)
			{
				temp3 = spectrum[i][j];
				temp2 = magnitude[i][j];
				final.setPixelVal(i, j, temp2);
				final2.setPixelVal(i, j, temp3);
			}
		}

		writeImage(tempname, final);
		writeImage(tempname2, final2);
	}
	if (size == 128)
	{
		char tempname[] = "fft_128.pgm";
		char tempname2[] = "fft_unshifted_128.pgm";
		ImageType final(512, 512, 255);
		ImageType final2(512, 512, 255);

		// Normalization
		float rmax = magnitude[0][0];
		float rmin = magnitude[0][0];
		float rmax2 = spectrum[0][0];
		float rmin2 = spectrum[0][0];
		for (int i = 0; i < 512; i++)
		{
			for (int j = 0; j < 512; j++)
			{
				if (magnitude[i][j] > rmax)
				{
					rmax = magnitude[i][j];
				}
				if (magnitude[i][j] < rmin)
				{
					rmin = magnitude[i][j];
				}
				if (spectrum[i][j] > rmax2)
				{
					rmax2 = spectrum[i][j];
				}
				if (spectrum[i][j] < rmin2)
				{
					rmin2 = spectrum[i][j];
				}
			}
		}
		for (int i = 0; i < 512; i++)
		{
			for (int j = 0; j < 512; j++)
			{
				magnitude[i][j] = 255 * (magnitude[i][j] - rmin) / (rmax - rmin);
				spectrum[i][j] = 255 * (spectrum[i][j] - rmin2) / (rmax2 - rmin2);
			}
		}

		int temp2;
		int temp3;
		// Write to image
		for (int i = 0; i < 512; i++)
		{
			for (int j = 0; j < 512; j++)
			{
				temp3 = spectrum[i][j];
				temp2 = magnitude[i][j];
				final.setPixelVal(i, j, temp2);
				final2.setPixelVal(i, j, temp3);
			}
		}

		writeImage(tempname, final);
		writeImage(tempname2, final2);
	}

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
}*/

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