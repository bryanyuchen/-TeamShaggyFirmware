#include <Kalmf.h>
#include <LevMar3D.h>
#include <MatrixMath.h>
LevMar3D levmar;
Kalmf kalmf;
float y_estimate[4];

void setup() {
 Serial.begin(9600); 
}

// norm function for time delay calculation
float norm(float* vector, uint8_t vector_size) {
  float total = 0;
  for (int i = 0; i < vector_size; i++) {
    total += (vector[i]*vector[i]);
  }
  return sqrt(total);
}

// state -> time delay given the following microphone matrix where origin is the centroid of the bottom surface
// 3D Microphone Array ....^(y axis)..........
//.........................|..................
//........................,m1.................
//......................../|\.................
//......................./.|.\..(z axis perp).
//......................./.|.\................
//....................../..|..\...............
//....................../..m4-\-->(x axis)....
//...................../../.\..\..............
//...................././..o..\.\.............
//..................././.......\.\............
//..................//...........\\...........
//.................m2-------------m3..........
// m1 = position of microphone 1
// m2 = position of microphone 2
// m3 = position of microphone 3
// m4 = position of microphone 4
// o = origin
// x axis is horizontal
// y axis is vertical
// z axis is perpendicular to the computer screen
void Y_EST(float* state) {
  //microphone variables
  float c = 34.3;
  float m1[3] = {0,11.62f,0};
  float m2[3] = {-10.0f,-5.7f,0};
  float m3[3] = {10.0f,-5.7f,0};
  float m4[3] = {0,0,16.28f};

  //y estimate variables
  float y_est_temp1[3];
  float y_est_temp2[3];
  float y_est_temp3[3];
  float y_est_temp4[3];

  Matrix.Subtract((float*)state,(float*)m1,3,1,(float*)y_est_temp1);
  Matrix.Subtract((float*)state,(float*)m2,3,1,(float*)y_est_temp2);
  Matrix.Subtract((float*)state,(float*)m3,3,1,(float*)y_est_temp3);
  Matrix.Subtract((float*)state,(float*)m4,3,1,(float*)y_est_temp4);

  y_estimate[0] = (1/c) * (norm(y_est_temp1,3) - norm(y_est_temp2,3));
  y_estimate[1] = (1/c) * (norm(y_est_temp2,3) - norm(y_est_temp3,3));
  y_estimate[2] = (1/c) * (norm(y_est_temp3,3) - norm(y_est_temp4,3));
  y_estimate[3] = (1/c) * (norm(y_est_temp1,3) - norm(y_est_temp4,3));
  Matrix.Print((float*)y_estimate,4,1,"y_est");
  //OPTIONAL add noise here
  float var = 0.001;
}

void loop() {
float state_set[3] = {25.0f, 25.0f, 25.0f};

//convert to time delay {-0.44f,0.31f,0.17f,0.04f}
Y_EST((float*)state_set);
float measurements[4];
Matrix.Copy((float*)y_estimate,4,1,(float*)measurements);

//Levenberg Marquardt State Estimation
float* state_estimate = levmar.Run(measurements);
Matrix.Print((float*)state_estimate,3,1,"LevMar State Estimation");
delay(100000);
}

