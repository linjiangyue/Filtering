#include "pfm.h"
#include "Vector.hpp"

#include <iostream>
#include <math.h>
#include <cstring>
#include <algorithm>

using namespace std;

// indexing 1d image array
int imIdx(int i, int j, int width, int channels = 3) {
	return channels * width * i + channels * j;
}

// clamp function
float clamp(float num, float low, float high){
	if (num > high){
		return high;
	}
	else if (num < low){
		return low;
	}
	else{
		return num;
	}
}

// Find the distance between pixels
float distanceSqr(int i_x, int i_y, int j_x, int j_y){
	return float(pow((i_x - j_x), 2) + pow((i_y - j_y), 2));
}

// Gaussian filter kernel
float G_kernel(float distanceSqr, float sigma_p){
	return exp(-(distanceSqr / (2 * pow(sigma_p, 2))));
}

// Bilateral filter kernel
float B_kernel(float distanceSqr, float colorDistSqr, float sigma_p, float sigma_c){
	return exp(-(distanceSqr / (2 * pow(sigma_p, 2))) - (colorDistSqr / (2 * pow(sigma_c, 2))));
}

// Joint bilateral filter kernel
float J_kernel(float distanceSqr, float colorDistSqr, float normalDistSqr, float positionDistSqr, 
			   float sigma_p, float sigma_c, float sigma_n, float sigma_d){

	return exp(-(distanceSqr / (2 * pow(sigma_p, 2))) - (colorDistSqr / (2 * pow(sigma_c, 2))) 
			   - (normalDistSqr / (2 * pow(sigma_n, 2))) - (positionDistSqr / (2 * pow(sigma_d, 2))));
}

float* GaussianFilter(float* image, unsigned width, unsigned height) {
	// C_output
	float* gauss = new float[width*height*3];
	memcpy(gauss, image, sizeof(float)*width*height*3);

	// Parameters
	float sigma_p = 5;

	// TODO 1
	// Implement Gaussian Filtering

	int i = 0;

	// For each pixel i
	for (int i_x = 0; i_x < width; i_x++){
		for (int i_y = 0; i_y < height; i_y++){

			// Weights and values
			float sow = 0.0;

			float sow_v_R = 0.0;
			float sow_v_G = 0.0;
			float sow_v_B = 0.0;

			// For each pixel j around i
			for (int j_x = i_x - 10; j_x <= i_x + 10; j_x++){
				for (int j_y = i_y - 10; j_y <= i_y + 10; j_y++){

					float w_ij = 0.0;
					
					// Out of bounds
					if (j_x < 0 || j_x >= width || j_y < 0 || j_y >= height){
						sow_v_R += 0.0;
						sow_v_G += 0.0;
						sow_v_B += 0.0;
					}else{
						// Calculate the weight
						w_ij = G_kernel(distanceSqr(i_x, i_y, j_x, j_y), sigma_p);

						// Add to values
						sow_v_R += w_ij * image[imIdx(j_y, j_x, width) + 0];
						sow_v_G += w_ij * image[imIdx(j_y, j_x, width) + 1];
						sow_v_B += w_ij * image[imIdx(j_y, j_x, width) + 2];
					}

					// Add to weights
					sow += w_ij;
				}
			}

			i++;

			// Output values
			gauss[imIdx(i_y, i_x, width) + 0] = sow_v_R / sow;
			gauss[imIdx(i_y, i_x, width) + 1] = sow_v_G / sow;
			gauss[imIdx(i_y, i_x, width) + 2] = sow_v_B / sow;
		}
	}

	cout << i << " TOTAL\n";
	
	cout << "Gaussian Filtering -- DONE\n";
	return gauss;
}

float* BilateralFilter(float* image, unsigned width, unsigned height) {
	// C_output
	float* bilateral = new float[width*height*3];
	memcpy(bilateral, image, sizeof(float)*width*height*3);

	// Parameters
	float sigma_p = 5;
	float sigma_c = 0.1;

	// TODO 2
	// Implement Bilateral Filtering 

	int i = 0;

	for (int i_x = 0; i_x < width; i_x++){
		for (int i_y = 0; i_y < height; i_y++){

			float sow = 0.0;

			float sow_v_R = 0.0;
			float sow_v_G = 0.0;
			float sow_v_B = 0.0;


			for (int j_x = i_x - 10; j_x <= i_x + 10; j_x++){
				for (int j_y = i_y - 10; j_y <= i_y + 10; j_y++){

					float w_ij = 0.0;
					
					if (j_x < 0 || j_x >= width || j_y < 0 || j_y >= height){
						sow_v_R += 0.0;
						sow_v_G += 0.0;
						sow_v_B += 0.0;
					}else{

						// Calculate distance between colors
						float colorDistSqr = float(pow((image[imIdx(i_y, i_x, width) + 0] - image[imIdx(j_y, j_x, width) + 0]), 2)
							+ pow((image[imIdx(i_y, i_x, width) + 1] - image[imIdx(j_y, j_x, width) + 1]), 2)
							+ pow((image[imIdx(i_y, i_x, width) + 2] - image[imIdx(j_y, j_x, width) + 2]), 2));

						w_ij = B_kernel(distanceSqr(i_x, i_y, j_x, j_y), colorDistSqr, sigma_p, sigma_c);
						
						sow_v_R += w_ij * image[imIdx(j_y, j_x, width) + 0];
						sow_v_G += w_ij * image[imIdx(j_y, j_x, width) + 1];
						sow_v_B += w_ij * image[imIdx(j_y, j_x, width) + 2];
					}

					sow += w_ij;
				}
			}

			i++;

			bilateral[imIdx(i_y, i_x, width) + 0] = sow_v_R / sow;
			bilateral[imIdx(i_y, i_x, width) + 1] = sow_v_G / sow;
			bilateral[imIdx(i_y, i_x, width) + 2] = sow_v_B / sow;
		}
	}

	cout << i << " TOTAL\n";
	
	cout << "Bilateral Filtering -- DONE\n";
	return bilateral;
}

float* JointBilateralFilter(float* image, float* normal, float* position, unsigned width, unsigned height) {
	// C_output
	float* joint = new float[width*height*3];
	memcpy(joint, image, sizeof(float)*width*height*3);

	// Parameters
	float sigma_p = 5;
	float sigma_c = 0.3;
	float sigma_n = 0.1;
	float sigma_d = 0.1;

	// TODO 3
	// Implement Joint Bilateral Filtering

	int i = 0;

	for (int i_x = 0; i_x < width; i_x++){
		for (int i_y = 0; i_y < height; i_y++){

			float sow = 0.0;

			float sow_v_R = 0.0;
			float sow_v_G = 0.0;
			float sow_v_B = 0.0;


			for (int j_x = i_x - 10; j_x <= i_x + 10; j_x++){
				for (int j_y = i_y - 10; j_y <= i_y + 10; j_y++){

					float w_ij = 0.0;
				
					if (j_x < 0 || j_x >= width || j_y < 0 || j_y >= height){
						sow_v_R += 0.0;
						sow_v_G += 0.0;
						sow_v_B += 0.0;
					}else{

						float colorDistSqr = float(pow((image[imIdx(i_y, i_x, width) + 0] - image[imIdx(j_y, j_x, width) + 0]), 2)
							+ pow((image[imIdx(i_y, i_x, width) + 1] - image[imIdx(j_y, j_x, width) + 1]), 2)
							+ pow((image[imIdx(i_y, i_x, width) + 2] - image[imIdx(j_y, j_x, width) + 2]), 2));
			
						// Calculate distance between surface normals
						Vector3f N_i = Vector3f(normal[imIdx(i_y, i_x, width) + 0], normal[imIdx(i_y, i_x, width) + 1], normal[imIdx(i_y, i_x, width) + 2]);
						Vector3f N_j = Vector3f(normal[imIdx(j_y, j_x, width) + 0], normal[imIdx(j_y, j_x, width) + 1], normal[imIdx(j_y, j_x, width) + 2]);

						float normalDot = float(clamp(dotProduct(N_i, N_j), 0, 1));
						float normalDistSqr = float(clamp(pow(acos(normalDot), 2), 0, 1));

						// Non-planarity measurement
						Vector3f P_i = Vector3f(position[imIdx(i_y, i_x, width) + 0], position[imIdx(i_y, i_x, width) + 1], position[imIdx(i_y, i_x, width) + 2]);
						Vector3f P_j = Vector3f(position[imIdx(j_y, j_x, width) + 0], position[imIdx(j_y, j_x, width) + 1], position[imIdx(j_y, j_x, width) + 2]);

						Vector3f P_i_j = P_j - P_i;
						P_i_j = normalize(P_i_j);

						float positionDot = float(clamp(dotProduct(N_i, P_i_j), 0, 1));
						float positionDistSqr = float(clamp(pow(dotProduct(N_i, P_i_j), 2), 0, 1));

						w_ij = J_kernel(distanceSqr(i_x, i_y, j_x, j_y), colorDistSqr, normalDistSqr, positionDistSqr, sigma_p, sigma_c, sigma_n, sigma_d);
						
						sow_v_R += w_ij * image[imIdx(j_y, j_x, width) + 0];
						sow_v_G += w_ij * image[imIdx(j_y, j_x, width) + 1];
						sow_v_B += w_ij * image[imIdx(j_y, j_x, width) + 2];
					}

					sow += w_ij;
				}
			}

			i++;

			joint[imIdx(i_y, i_x, width) + 0] = sow_v_R / sow;
			joint[imIdx(i_y, i_x, width) + 1] = sow_v_G / sow;
			joint[imIdx(i_y, i_x, width) + 2] = sow_v_B / sow;
		}
	}

	cout << i << " TOTAL\n";

	cout << "Joint Bilateral Filtering -- DONE\n";
	return joint;
}

int main() {
    unsigned *w = new unsigned;
    unsigned *h = new unsigned;

	// Input buffers
	float* imageBuffer = read_pfm_file3("im/s2.pfm", w, h);      		 // load image buffer - 3 channels
    float* normalBuffer = read_pfm_file3("im/s2_normal.pfm", w, h);      // load normal buffer - 3 channels
	float* positionBuffer = read_pfm_file3("im/s2_position.pfm", w, h);  // load position buffer - 3 channels
	

	float* gausssianFiltered = GaussianFilter(imageBuffer, *w, *h);
	write_pfm_file3("im/s2_gaussian.pfm", gausssianFiltered, *w, *h);
	delete[] gausssianFiltered;
	
	float* bilateralFiltered = BilateralFilter(imageBuffer, *w, *h);
	write_pfm_file3("im/s2_bilateral.pfm", bilateralFiltered, *w, *h);
	delete[] bilateralFiltered;
	
	float* jointBilateralFiltered = JointBilateralFilter(imageBuffer, normalBuffer, positionBuffer, *w, *h);
	write_pfm_file3("im/s2_jointBilateral.pfm", jointBilateralFiltered, *w, *h);
	delete[] jointBilateralFiltered;


	// 
	delete[] imageBuffer;
	delete[] normalBuffer;
	delete[] positionBuffer;
	
	return 0;
}