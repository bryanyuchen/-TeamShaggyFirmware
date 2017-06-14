clear all;
close all;

for i = 1:50
    path(:,i) = [2*log(i)+20;4*log(i)+11;5*sin(i)+20];
    B(:,i) = YEstimate(path(:,i));
end

lengthB = length(B);
random = 0.001*randn(4,lengthB);
B = B + random;
 
%LPF
% for i = 1:lengthB-8
%       temp = B(:,i:i+8); 
%       B(:,i+4) = mean(temp,2);
% end 

length = 1000;

for a = 1:length
mnoise = (10.1-a*(10/length))
cov = eye(3);
state = [1;1;1];
state_est = LevMarFunc(B(:,i),[1;1;1]);
x(1, 1) = state_est(1, 1);
y(1, 1) = state_est(2, 1); 
z(1, 1) = state_est(3, 1);
state(1:3,1) = state_est;
[state, cov] = Kalmf(state,cov,[x(1);y(1);z(1)],mnoise);
x2(1, 1) = state(1, 1);
y2(1, 1) = state(2, 1); 
z2(1, 1) = state(3, 1);
    for i = 2:lengthB 
      state_est = LevMarFunc(B(:,i),state(1:3));
      x(1, i) = state_est(1, 1);
      y(1, i) = state_est(2, 1); 
      z(1, i) = state_est(3, 1);
      
      [state, cov] = Kalmf(state,cov,[x(i);y(i);z(i)],mnoise);
      x2(1, i) = state(1, 1);
      y2(1, i) = state(2, 1); 
      z2(1, i) = state(3, 1);
    end
path2 = vertcat(x,y,z);
path3 = vertcat(x2,y2,z2);
path2 = path2-path;
path3 = path3-path;
path2 = path2(:,2:50);
path3 = path3(:,2:50);
varLM(:,a) = var(path2,0,2);
varKF(:,a) = var(path3,0,2);
t(1,a) = mnoise;
end
varKF2 = varKF(1,:).*varKF(2,:).*varKF(3,:);
varLM2 = varLM(1,:).*varLM(2,:).*varLM(3,:);
plot(t,varLM2(1,:),'b');
hold on
plot(t,varKF2(1,:),'r');
title('Aggregate KF Variance over Measurement Noise');
xlabel('Measurement Noise (cm)');
ylabel('Variance (cm)');
legend('LM Variance', 'KF Variance');
% 
% for i=1:length(B) 
%     %LevMar Estimation
%     scatter3(x(i), y(i), z(i),  'b');
%     hold on
%     %Kalman Filtered Estimation
%     scatter3(x2(i),y2(i),z2(i),'r');
%     hold on
%     %Actual Position
%     scatter3(path(1,i), path(2,i), path(3,i),'m');
%     hold on
%     %plot mic matrix
%     scatter3(0,11.62,0,'g');
%     hold on
%     scatter3(-10,-5.7,0,'g');
%     hold on
%     scatter3(10,-5.7,0,'g');
%     hold on
%     scatter3(0,0,16.28,'g');
% 
%     set(gcf,'units','normalized','outerposition',[0 0 1 1]);
%     
%     view(-18,24)
%     axis([-11 39 -11 39 -11 39])
%     title('3d Source Localization Simulation, acceleration input state-space KF (3 mnoise), 1us noise (1Mhz precision), 50 samples');
%     xlabel('x axis (cm)');
%     ylabel('y axis (cm)');
%     zlabel('z axis (cm)');
%     legend('Levenberg Marquardt Estimation','LM after Kalman Filtering', 'Actual Position', 'Microphone Matrix');
%    
%     pause(.1)
%     hold off
% end

% scatter3(x, y, z,  'b');
%     hold on
%     %Kalman Filtered Estimation
%     scatter3(x2,y2,z2,'r');
%     hold on
%     %plot mic matrix
%     scatter3(0,11.62,0,'g');
%     hold on
%     scatter3(-10,-5.7,0,'g');
%     hold on
%     scatter3(10,-5.7,0,'g');
%     hold on
%     scatter3(0,0,16.28,'g');
% 
%     view(-18,24)
%     axis([-11 39 -11 39 -11 39])
%     title('3d Source Localization Simulation (acceleration input state-space), 15us noise, 50 samples');
%     xlabel('x axis (cm)');
%     ylabel('y axis (cm)');
%     zlabel('z axis (cm)');
%     legend('Levenberg Marquardt Estimation','LM after Kalman Filtering', 'Actual Position', 'Microphone Matrix');

