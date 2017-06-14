function [state, cov] = Kalmf(state, cov, measurements)
%Input: 
%-state:        state variables x,y in centimeters( 2x1 matrix )
%-prev_state:   previous state (2x1), since function holds no global memory
%-cov:          feedback error covariance matrix (initialized to eye(2)), unless otherwise
%-measurements: Time delays between the 3 microphones in an equilateral triangle in ms (3x1 matrix) 
%Output: 
%-state:        updated state prediction in centimeters
%-prev_state:   pass back current state as prev_state
%- cov:         updated cov matrix

% look at http://robotsforroboticists.com/kalman-filtering/
% http://blog.tkjelectronics.dk/2012/09/a-practical-approach-to-kalman-filter-and-how-to-implement-it/
%% Global variables
c = 34.3; %speed of sound in dry air cm/ms
dt = 1; %time step - not yet needed until model for matrix A is produced
u = 0.005;
%% Microphone position
 m1 = [0;11.62;0];
 m2 = [-10;-5.7;0];
 m3 = [10;-5.7;0];
 m4 = [0;0;16.28];

%% Covariances
process_noise = 0.1; % process noise (variance of stationary speaker model in cm)
measurement_noise_TD1 = 0.1;  % LevMar measurement noise in x direction 
measurement_noise_TD2 = 0.1;  % LevMar measurement noise in y direction 
measurement_noise_TD3 = 0.1;  % LevMar measurement noise in z direction 

%% State Space Model assuming stationary speaker

% next state is the same as current state (no B matrix)
    A = eye(3)
     
    B = [((dt.^2) / 2); ...
         ((dt.^2) / 2); ...
         ((dt.^2) / 2); ...
         dt; ...
         dt; ...
         dt];

% process noise
    v1 = eye(3).*process_noise^2;
    
% function of measured time delay = some function of position estimate

C = eye(3)
             
%measurement noise    
    v2 = [measurement_noise_TD1^2 0 0; ...
            0 measurement_noise_TD2^2 0; ...
            0 0 measurement_noise_TD3^2];

%% KF process
    %predict next state
    next_state = A * state;% + B.*u; % predict next state (basically does not change for now)
    cov = A * cov * A' + v1; % predict next cov (basically does not change for now)

    %Innovation
    Y = measurements - C*next_state; %Innovation in milliseconds
    S = C * cov * C' + v2; %Innovation covariance
    K = cov * C' / (S); %Kalman Gain (other words compare innovation covariance with covariance
    
    %update
    new_state = next_state + K * Y; %new predicted state in cm
    cov = (eye(3) - K * C) * cov;    

    %return updated values
    state = new_state;
end