#include <Kalmf.h>
#include <LevMar3D.h>
#include <MatrixMath.h>
LevMar3D levmar;
Kalmf kalmf;

void setup() {
 Serial.begin(9600); 
}
float norm(float* vector, uint8_t vector_size) {
  float total = 0;
  for (int i = 0; i < vector_size; i++) {
    total += (vector[i]*vector[i]);
  }
  return sqrt(total);
}

float* Y_EST(float* state) {
  //microphone variables
  float c = 34.3;
  float m1[3] = {0,11.62f,0};
  float m2[3] = {-10.0f,-5.7f,0};
  float m3[3] = {10.0f,-5.7f,0};
  float m4[3] = {0,0,16.28f};

  //y estimate procedure
  float y_estimate[4];
  float temp1[N];
  float temp2[N];
  float temp3[N];
  float temp4[N];
  Matrix.Subtract((float*)state,(float*)m1,3,1,(float*)temp1);
  Matrix.Subtract((float*)state,(float*)m2,3,1,(float*)temp2);
  Matrix.Subtract((float*)state,(float*)m3,3,1,(float*)temp3);
  Matrix.Subtract((float*)state,(float*)m4,3,1,(float*)temp4);

  y_estimate[0] = (1/c) * (norm(temp1,3) - norm(temp2,3));
  y_estimate[1] = (1/c) * (norm(temp2,3) - norm(temp3,3));
  y_estimate[2] = (1/c) * (norm(temp3,3) - norm(temp4,3));
  y_estimate[3] = (1/c) * (norm(temp1,3) - norm(temp4,3));
  Matrix.Print((float*)y_estimate,4,1,"y_est");
  //OPTIONAL add noise here
  float var = 0.001;

  return y_estimate;
}

void loop() {
//float measurements[4] = {-0.4592,0.2977,0.1331,-0.0283}; // ideal measurement
float state_set[3] = {25.0f, 25.0f, 25.0f};
float* measurements = Y_EST(state_set);
Matrix.Print((float*)state_set,3,1,"state_set");
float* state_estimate = levmar.Run(measurements);
Matrix.Print((float*)state_estimate,3,1,"LevMar State Estimation");
delay(100000);
}

