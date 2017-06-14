#include "Arduino.h"
#include "LevMar3D.h"

/********************************************************************/
/*
 * LevMar3D Constructor - all that is needed at this point is microphone position
 * and identity matrix
 */
LevMar3D::LevMar3D()
{
// 3D Microphone Array ..........^(y axis)..........
m1[0] = 0; //....................|..................
m1[1] = 11.62f; //..............,m1.................
m1[2] = 0; //.................../|\.................
m2[0] = -10.0f; //............./.|.\..(z axis perp).
m2[1] = -5.7f; //............../.|.\................
m2[2] = 0; //................./..|..\...............
m3[0] = 10.0f; //............./..m4-\-->(x axis)....
m3[1] = -5.7f; //............/../.\..\..............
m3[2] = 0; //..............././.....\.\.............
m4[0] = 0; //.............././.......\.\............
m4[1] = 0; //.............//...........\\...........
m4[2] = 16.28f; //.......m2-------------m3..........

// lambda diagonal matrix init
for (int i = 0; i < N; i++) {
	for (int j = 0; j < N; j++){
		if (i == j){
			eye[i][j] = lambda;
		}
		else {
			eye[i][j] = 0;
		}
	}
}
}
/********************************************************************/
/*
 * Norm Function
 */
float LevMar3D::norm(float* vector, uint8_t vector_size) {
  float total = 0;
  for (int i = 0; i < vector_size; i++) {
    total += (vector[i]*vector[i]);
  }
  return sqrt(total);
}

/********************************************************************/
/*
 * Levenberg-Marquardt Algorithm
 */
float* LevMar3D::Run(float* measurements) {
uint8_t updateJ=1; //update?
state_est[0] = 5;
state_est[1] = 5;
state_est[2] = 5;
for (uint8_t it = 1; it <= n_iters; it++) {
  // if updated, calculate new parameters
  if (updateJ == 1) {
    // Evaluate the Jacobian matrix (4x3) at the current parameters 
    Matrix.Subtract((float*)state_est,(float*)m1,N,1,(float*)res1);
    Matrix.Subtract((float*)state_est,(float*)m2,N,1,(float*)res2);
    Matrix.Subtract((float*)state_est,(float*)m3,N,1,(float*)res3);
	Matrix.Subtract((float*)state_est,(float*)m4,N,1,(float*)res4);

    J[0][0] = (1/c) * ((state_est[0]-m1[0])/norm(res1,N) - (state_est[0]-m2[0])/norm(res2,N));
    J[0][1] = (1/c) * ((state_est[1]-m1[1])/norm(res1,N) - (state_est[1]-m2[1])/norm(res2,N));
	J[0][2] = (1/c) * ((state_est[2]-m1[2])/norm(res1,N) - (state_est[2]-m2[2])/norm(res2,N));
    J[1][0] = (1/c) * ((state_est[0]-m2[0])/norm(res2,N) - (state_est[0]-m3[0])/norm(res3,N));
    J[1][1] = (1/c) * ((state_est[1]-m2[1])/norm(res2,N) - (state_est[1]-m3[1])/norm(res3,N));
	J[1][2] = (1/c) * ((state_est[2]-m2[2])/norm(res2,N) - (state_est[2]-m3[2])/norm(res3,N));
    J[2][0] = (1/c) * ((state_est[0]-m3[0])/norm(res3,N) - (state_est[0]-m4[0])/norm(res4,N));
    J[2][1] = (1/c) * ((state_est[1]-m3[1])/norm(res3,N) - (state_est[1]-m4[1])/norm(res4,N));
	J[2][2] = (1/c) * ((state_est[2]-m3[2])/norm(res3,N) - (state_est[2]-m4[2])/norm(res4,N));
	J[3][0] = (1/c) * ((state_est[0]-m1[0])/norm(res1,N) - (state_est[0]-m4[0])/norm(res4,N));
	J[3][1] = (1/c) * ((state_est[1]-m1[1])/norm(res1,N) - (state_est[1]-m4[1])/norm(res4,N));
	J[3][2] = (1/c) * ((state_est[2]-m1[2])/norm(res1,N) - (state_est[2]-m4[2])/norm(res4,N));

    // Evaluate the distance error at the current parameters
    y_est[0] = (1/c) * (norm(res1,N) - norm(res2,N));
    y_est[1] = (1/c) * (norm(res2,N) - norm(res3,N));
    y_est[2] = (1/c) * (norm(res3,N) - norm(res4,N));
	y_est[3] = (1/c) * (norm(res1,N) - norm(res4,N));

    // Evaluate residual matrix
    Matrix.Subtract(measurements,(float*)y_est,M,1,(float*)d);

    // compute the approximated Hessian matrix, J’ is the transpose of J
    Matrix.Transpose((float*)J,M,N,(float*)J_transpose);
    Matrix.Multiply((float*)J_transpose,(float*)J,N,M,N,(float*)H);

    //if the first iteration, compute the total error    
    if (it==1) { 
      //slight modification to error calculation using abs() to save memory and improve precision
      //e = (d[0]*d[0] + d[1]*d[1] + d[2]*d[2])..;
	  e = (abs(d[0]) + abs(d[1]) + abs(d[2]) + abs(d[3]));
      Serial.print("e = ");
      Serial.println(e,7);
    }
  }

  Matrix.Add((float*)H,(float*)eye,N,N,(float*)H_lm);

  // Compute the updated parameters
  Matrix.Invert((float*)H_lm,N);
  Matrix.Multiply((float*)J_transpose,(float*)d,N,M,1,(float*)LevMarTemp);
  Matrix.Multiply((float*)H_lm,(float*)LevMarTemp,N,N,1,(float*)dp);

  Matrix.Add((float*)state_est,(float*)dp,N,1,(float*)state_lm);
  
  //going towards direction of wrong minima
  //if (state_lm[1] < 1.0f) {
  //    state_lm[1] = state_est[1]-(10*dp[1]);
  //}
  
  //Evaluate the total distance error at the updated parameters
  Matrix.Subtract((float*)state_lm,(float*)m1,N,1,(float*)res1);
  Matrix.Subtract((float*)state_lm,(float*)m2,N,1,(float*)res2);
  Matrix.Subtract((float*)state_lm,(float*)m3,N,1,(float*)res3);
  Matrix.Subtract((float*)state_lm,(float*)m4,N,1,(float*)res4);
  
  //Evaluate new y estimate with levmar step
  y_est_lm[0] = (1/c) * (norm(res1,N) - norm(res2,N));
  y_est_lm[1] = (1/c) * (norm(res2,N) - norm(res3,N));
  y_est_lm[2] = (1/c) * (norm(res3,N) - norm(res4,N));
  y_est_lm[3] = (1/c) * (norm(res1,N) - norm(res4,N));
  
  Matrix.Subtract((float*)measurements,(float*)y_est_lm,M,1,d_lm);
  
  //Evaluate new total distance error with levmar step
  e_lm = abs(d_lm[0]) + abs(d_lm[1]) + abs(d_lm[2]) + abs(d_lm[3]);
  //Serial.print("e_lm = ");
  //Serial.println(e_lm,7);
  Matrix.Print((float*)state_est,N,1,"state_est");
  
  //compare total distance errors
  
  //if less error, continue to step in that direction
  if (e_lm < e) { 
   Serial.println("update step");
   Matrix.Scale((float*)eye,N,N,(1/10));
   Matrix.Copy((float*)state_lm,1,N,(float*)state_est);
   e = e_lm;
   Serial.print("error: ");
   Serial.println(e,7);
   updateJ = 1;
  }
  
  //if more error, increase damping factor
  else {
    Serial.println("increase damping");
    updateJ = 0;
	Matrix.Scale((float*)eye,N,N,10);
  }
}
return state_est;
}
/********************************************************************/

