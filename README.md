# Team Shaggy TDOA Firmware by Bryan Chen
 June 13, 2017 Los Angeles, CA.

## This repository contains all of the firmware code for acoustic source localization
### Arduino code contains all code necessary to run on the Arduino. 
The algorithms utilize the Arduino matrix math library: https://playground.arduino.cc/Code/MatrixMath
This library is also included in this repo for convenience but all copyrights go to Arduino

-2D code contains both Levenberg-Marquardt Algorithm and Kalman Filter - Explanation is in the paper.

-3D code only utilizes Levenberg-Marquardt Algorithm because it assumes higher sampling precision

### Matlab code contains all code necessary to run simulations on Matlab and recreate the figures in the report.
The plots are the .jpg files and the animations are the .avi files. Animations will need to be downloaded in order to view.

-2D algorithm accuracy comparison compares the accuracy between Levenberg-Marquardt and IEKF given real data

-2D algorithm efficiency comparison compares the processing speed between all three algorithms

-3D Simulation contains all of the 3D simulation and animation code
