/*
  LevMar.h - Library for Levenberg-Marquardt Algorithm
  Created by Bryan Chen, May 8, 2017
  Released into the public domain.
*/
#ifndef LevMar_h
#define LevMar_h
#include "Arduino.h"
#include <MatrixMath.h>

class LevMar
{
	public:
		LevMar(float mdist);
		float norm(float* vector, uint8_t vector_size);
		float* Run(float* measurements);
		
/********************************************************************/
/*
 * Global Variables
 */

// 2D Microphone Array
#define N  (2)
float c = 34.3;
//uint8_t mdist = 20;
float m1[N]; // vector of mic 1
float m2[N]; // vector of mic 2
float m3[N]; // vector of mic 3

// LevMar Parameters
uint8_t n_iters = 10; // # of iterations for LM
float lambda= 0.01; // initial damping factor 

//LevMar Variables
float state_est[N]; //state estimate (maybe can be optimized)
float state_lm[N]; //Lev_Marstate estimate
float dp[N]; //step
float y_est [3]; //y estimate, derived from state_est
float y_est_lm [3]; //y estimate, derived from LevMar
float d[3]; //Residual
float d_lm[3]; //Residual after LevMar step
float e; //Error (residual^2)
float e_lm; //Error (residual^2) after LevMar step
float J [3][2]; //Jacobian
float J_transpose [2][3]; //Transpose of Jacobian
float temp[2]; //temp matrix to help J'*d
float H[2][2]; //Hessian
float H_lm[2][2];
float res1[N]; //result vector for Jacobian math
float res2[N]; //result vector for Jacobian math
float res3[N]; //result vector for Jacobian math
float eye[N][N];

/********************************************************************/
};
#endif
