/*
  Kalmf.h - Library for Basic Kalman Filtering
  Created by Bryan Chen, May 15, 2017
  Released into the public domain.
*/
#ifndef Kalmf_h
#define Kalmf_h
#include "Arduino.h"
#include <MatrixMath.h>

class Kalmf
{
  public:
    Kalmf();
    float* Run(float* measurements);

/********************************************************************/
/*
 * Global Variables
 */
#define N  (2) //how many state variables?
#define M  (2) //number measurements

//define noise
float process_noise = 0;
float meas_noise1 = 1.1f;
float meas_noise2 = 1.1f;

//State-space Model
float A[N][N];
float A_transpose[N][N];
//float B
float v1[2][2]; //process noise matrix
float C[M][N];
float C_transpose[N][M];
//float D
float v2[2][2]; //measurement noise matrix

//Kalmf variables
float Kalmf_state[N];
float next_state[N];
float new_state[N];
float new_state_temp[N];
float cov[N][N];
float cov_temp1[N][N];
float cov_temp2[N][N];
float Y[N];
float Y_temp[N];
float S[N][N];
float S_temp1[N][N];
float S_temp2[N][N];
float K[N][N];
float K_temp[N][N];
float identity[N][N];

/********************************************************************/
};
#endif
