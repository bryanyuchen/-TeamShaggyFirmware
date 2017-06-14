/*
  LevMar.h - Library for Levenberg-Marquardt Algorithm
  Created by Bryan Chen, May 8, 2017
  Released into the public domain.
*/
#ifndef LevMar3D_h
#define LevMar3D_h
#include "Arduino.h"
#include <MatrixMath.h>

class LevMar3D
{
	public:
		LevMar3D();
		float norm(float* vector, uint8_t vector_size);
		float* Run(float* measurements);
		
/********************************************************************/
/*
 * Global Variables
 */

// 2D Microphone Array
#define N  (3) // num state variables
#define M  (4) // num measurements
float c = 34.3;
//uint8_t mdist = 20;
float m1[N]; // vector of mic 1 position
float m2[N]; // vector of mic 2 position
float m3[N]; // vector of mic 3 position
float m4[N]; // vector of mic 4 position
// LevMar Parameters
uint8_t n_iters = 5; // # of iterations for LM - 5 is about halfway to convergence
float lambda = 0.01f; // initial damping factor 

//LevMar Variables
float state_est[N]; //state estimate (maybe can be optimized)
float state_lm[N]; //Lev_Marstate estimate
float dp[N]; //step
float y_est [M]; //y estimate, derived from state_est
float y_est_lm [M]; //y estimate, derived from LevMar
float d[M]; //Residual
float d_lm[M]; //Residual after LevMar step
float e; //Error (residual^2)
float e_lm; //Error (residual^2) after LevMar step
float J [M][N]; //Jacobian
float J_transpose [N][M]; //Transpose of Jacobian
float LevMarTemp[N]; //temp matrix to help J'*d
float H[N][N]; //Hessian
float H_lm[N][N];
float res1[N]; //result vector for Jacobian math [norm(state_est-m)]
float res2[N]; //result vector for Jacobian math [norm(state_est-m)] 
float res3[N]; //result vector for Jacobian math [norm(state_est-m)]
float res4[N]; //result vector for Jacobian math [norm(state_est-m)]
float eye[N][N];

/********************************************************************/
};
#endif
