function [state, cov] = TDOA(state, cov, measurements)
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
%dt = 1; %time step - not yet needed until model for matrix A is produced

%% Microphone position
% diagram of microphones and positioning
%....................[m1](0,0)...................
%................................................
%................................................
%................................................
%................................................
%................................................
%................................................
%................................................
%................................................
%.......[m2]........................[m3]...........

mdist = 20; %20 cm separation between microphones in an equilateral triangle
m1 = [0;0];
m2 = [mdist*-0.5; mdist*-0.8660254];
m3 = [mdist*0.5; mdist*-0.8660254];

%plot out the positions of the microphones for visualization
%mplot = horzcat(m1, m2, m3);
%scatter(mplot(1,:),mplot(2,:));
%% Covariances
process_noise = 3; % process noise (variance of stationary speaker model in cm)
measurement_noise_TD1 = 1;  % measurement noise between microphones 1 and 2 (refers to time delay in milliseconds)
measurement_noise_TD2 = 1;  % measurement noise between microphones 2 and 3 (refers to time delay in milliseconds)
measurement_noise_TD3 = 1;  % measurement noise between microphones 1 and 3 (refers to time delay in milliseconds)
%% State Space Model assuming stationary speaker

% next state is the same as current state (no B matrix)
    A = [1 0 ; ... 
         0 1 ]; 
     
%    B = [cos(theta_actual) * ((t.^2) / 2); ...
%         sin(theta_actual) * ((t.^2) / 2); ...
%         cos(theta_actual) * t; ...
%         sin(theta_actual) * t;
%         0; ...
%         0];

% process noise
    v1 = eye(2).*process_noise^2;
    
% function of measured time delay = some function of position estimate
% NOTE: (C(t,x(t))is NONLINEAR

C = (1/c) * [norm(state-m1) - norm(state-m2); ...
                norm(state-m2) - norm(state-m3); ...
                norm(state-m1) - norm(state-m3)];

            x = state(1,1);
            y = state(2,1);
            
C_jacobian = (1/c) * [(x-m1(1,1))/norm(state-m1) - (x-m2(1,1))/norm(state-m2), ...
                    (y-m1(2,1))/norm(state-m1) - (y-m2(2,1))/norm(state-m2); ...
                    (x-m2(1,1))/norm(state-m2) - (x-m3(1,1))/norm(state-m3), ...
                    (y-m2(2,1))/norm(state-m2) - (y-m3(2,1))/norm(state-m3); ...
                    (x-m1(1,1))/norm(state-m1) - (x-m3(1,1))/norm(state-m3), ...
                    (y-m1(2,1))/norm(state-m1) - (y-m3(2,1))/norm(state-m3)];
                
%     C = transpose((1/c)*[(state-m1)/norm(state-m1) - (state-m2)/norm(state-m2), ...
%                  (state-m2)/norm(state-m2) - (state-m3)/norm(state-m3), ...
%                  (state-m1)/norm(state-m1) - (state-m3)/norm(state-m3)]);
             
%measurement noise    
    v2 = [measurement_noise_TD1^2 0 0; ...
            0 measurement_noise_TD2^2 0; ...
            0 0 measurement_noise_TD3^2];

%% IEKF process 
    
    %predict next state
    next_state = A * state; % predict next state (basically does not change for now)
    cov = A * cov * A' + v1; % predict next cov (basically does not change for now)

    %Innovation
    Y = measurements - C; %Innovation in milliseconds
    S = C_jacobian * cov * C_jacobian' + v2; %Innovation covariance
    K = cov * C_jacobian' / (S); %Kalman Gain (other words compare innovation covariance with covariance
    
    %update
    new_state = next_state + K * Y; %new predicted state in cm
    cov = (eye(2) - K * C_jacobian) * cov;
    
% % Measurement adjustment
%     T_prev = (1/c) * [(norm(state-m1) - norm(state-m2)); ...
%                (norm(state-m2) - norm(state-m3)); ...
%                (norm(state-m1) - norm(state-m3))]; % intermediate calculation of previous estimated delay
% 
%     measurements_macron = measurements - (T_prev - C * state) ; % adjusted measurement in ms
    
    %iterated
%      n_next = [0;0];
%      n_curr = next_state;
%      n_prev = next_state;
%      for i = 1:10
%      C_curr = (1/c) * [norm(n_curr-m1) - norm(n_curr-m2); ...
%                          norm(n_curr-m2) - norm(n_curr-m3); ...
%                          norm(n_curr-m1) - norm(n_curr-m3)];
%      x = n_curr(1,1);
%      y = n_curr(2,1);
%      C_jacobian_curr = (1/c) * [(x-m1(1,1))/norm(n_curr-m1) - (x-m2(1,1))/norm(n_curr-m2), ...
%                     (y-m1(2,1))/norm(n_curr-m1) - (y-m2(2,1))/norm(n_curr-m2); ...
%                     (x-m2(1,1))/norm(n_curr-m2) - (x-m3(1,1))/norm(n_curr-m3), ...
%                     (y-m2(2,1))/norm(n_curr-m2) - (y-m3(2,1))/norm(n_curr-m3); ...
%                     (x-m1(1,1))/norm(n_curr-m1) - (x-m3(1,1))/norm(n_curr-m3), ...
%                     (y-m1(2,1))/norm(n_curr-m1) - (y-m3(2,1))/norm(n_curr-m3)];
%                 
%      G_curr = cov * C_jacobian_curr' / (C_jacobian_curr * cov * C_jacobian_curr' + v2);
%  
%      alpha = measurements - C_curr;
%      z = alpha - C_jacobian_curr*(n_prev-n_curr);
%      n_next = n_prev + G_curr * z;
%      n_prev = n_curr;
%      n_curr = n_next;
%      end

    %return updated values
    state = new_state;
    %state = n_curr;
end

