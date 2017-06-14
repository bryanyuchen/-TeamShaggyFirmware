#include "Arduino.h"
#include "Kalmf.h"

/********************************************************************/
/*
 * Kalmf Constructor - initializes state space model
 */

Kalmf::Kalmf()
{
//init state
Kalmf_state[0] = 0;
Kalmf_state[1] = 0;

//init matrices
A[0][0] = 1;
A[0][1] = 0;
A[1][0] = 0;
A[1][1] = 1;

v1[0][0] = process_noise*process_noise;
v1[0][1] = 0;
v1[1][0] = 0;
v1[1][1] = process_noise*process_noise;

C[0][0] = 1;
C[0][1] = 0;
C[1][0] = 0;
C[1][1] = 1;

v2[0][0] = (meas_noise1*meas_noise1);
v2[0][1] = 0;
v2[1][0] = 0;
v2[1][1] = (meas_noise2*meas_noise2);

identity[0][0] = 1;
identity[0][1] = 0;
identity[1][0] = 0;
identity[1][1] = 1;

Matrix.Transpose((float*)A,N,N,(float*)A_transpose);
Matrix.Transpose((float*)C,N,N,(float*)C_transpose);

//init cov as identity
cov[0][0] = 1;
cov[0][1] = 0;
cov[1][0] = 0;
cov[1][1] = 1;
}
/********************************************************************/
/*
 * Kalman Filter Algorithm
 */
float* Kalmf::Run(float* curr_measurements) {
  //predict next state  
  Matrix.Multiply((float*)A,Kalmf_state,N,N,1,(float*)next_state);
  
  Matrix.Multiply((float*)cov,(float*)A_transpose,N,N,N,(float*)cov_temp1);
  Matrix.Multiply((float*)A,(float*)cov_temp1,N,N,N,(float*)cov_temp2);
  Matrix.Add((float*)cov_temp2,(float*)v1,N,N,(float*)cov);

  //Innovation
  Matrix.Multiply((float*)C,(float*)next_state,N,N,1,(float*)Y_temp);
  Matrix.Subtract(curr_measurements,(float*)Y_temp,N,1,(float*)Y);
  
  Matrix.Multiply((float*)cov,(float*)C_transpose,N,N,M,(float*)S_temp1);
  Matrix.Multiply((float*)C,(float*)S_temp1,M,N,N,(float*)S_temp2);
  Matrix.Add((float*)S_temp2,(float*)v2,M,M,(float*)S);
  
  //Kalman Gain
  Matrix.Invert((float*)S,N);
  Matrix.Multiply((float*)C_transpose,(float*)S,N,M,N,(float*)K_temp);
  Matrix.Multiply((float*)cov,(float*)K_temp,N,N,N,(float*)K);

  //Update
  Matrix.Multiply((float*)K,(float*)Y,N,N,1,(float*)new_state_temp);
  Matrix.Add((float*)next_state,(float*)new_state_temp,N,1,(float*)new_state);

  Matrix.Multiply((float*)K,(float*)C,N,N,N,(float*)S_temp1); //reusing old matrices to save space
  Matrix.Subtract((float*)identity,(float*)S_temp1,N,N,(float*)cov_temp1);
  Matrix.Multiply((float*)cov_temp1,(float*)cov,N,N,N,(float*)cov_temp2);
  Matrix.Copy((float*)cov_temp2,N,N,(float*)cov);
  
  //Return update
  Matrix.Copy((float*)new_state,N,1,(float*)Kalmf_state);
  return (float*)Kalmf_state;
}
/********************************************************************/
