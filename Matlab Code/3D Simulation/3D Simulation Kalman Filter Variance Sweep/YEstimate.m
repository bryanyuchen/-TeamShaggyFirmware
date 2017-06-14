function [ y_est ] = YEstimate( state_est )

%% Microphone variables
 c = 34.3;
 
 m1 = [0;11.62;0];
 m2 = [-10;-5.7;0];
 m3 = [10;-5.7;0];
 m4 = [0;0;16.28];
 
 %% y estimate
        y_est = (1/c) * [norm(state_est-m1) - norm(state_est-m2); ...
                         norm(state_est-m2) - norm(state_est-m3); ...
                         norm(state_est-m3) - norm(state_est-m4); ...
                         norm(state_est-m1) - norm(state_est-m4)];
%% add noise
std = 0.001; %15 us standard deviation
noise = std * randi([-3 3],4,1);
y_est = y_est+noise;
              
end

