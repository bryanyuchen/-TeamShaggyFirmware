function [ state_est ] = LevMarFunc( measurements, kalmf_est )


%% Microphone variables
 c = 34.3;
 %mdist = 20;
 %m1 = [0;mdist*0.43301270189];
 %m2 = [mdist*-0.5; mdist*-0.43301270189];
 %m3 = [mdist*0.5; mdist*-0.43301270189];  
  m1 = [0;11.62;0];
 m2 = [-10;-5.7;0];
 m3 = [10;-5.7;0];
 m4 = [0;0;16.28];


%measurements = [-0.4592;0.2977;0.1331;-0.0283];
%% LevMar variables
%Nparams=2; % defines size of lambda matrix
%%%%%%%%%%%%?? L - M ??%%%%%%%%%%%%%%%
Nparams=3;
%%%%%%%%%%%%?? L - M ??%%%%%%%%%%%%%%%
n_iters = 9; % set # of iterations for the LM
lamda= 0.01; % set an initial value of the damping factor for the LM
updateJ=1;
a_est= kalmf_est(1);
b_est= kalmf_est(2); 
c_est= kalmf_est(3);
state_est = kalmf_est; % intial state estimate (maybe can be optimized)
%%%%%%%%%%%%?? L - M ??%%%%%%%%%%%%%%%

%% LevMar Algorithm
for it=1:n_iters
    if updateJ==1
     
 % Evaluate the Jacobian matrix at the current parameters (a_est, b_est)
%%%%%%%%%%%%?? L - M ??%%%%%%%%%%%%%%%
        J = (1/c) * [(a_est-m1(1,1))/norm(state_est-m1) - (a_est-m2(1,1))/norm(state_est-m2), ...
                    (b_est-m1(2,1))/norm(state_est-m1) - (b_est-m2(2,1))/norm(state_est-m2), ...
                    (c_est-m1(3,1))/norm(state_est-m1) - (c_est-m2(3,1))/norm(state_est-m2); ...
                    (a_est-m2(1,1))/norm(state_est-m2) - (a_est-m3(1,1))/norm(state_est-m3), ...
                    (b_est-m2(2,1))/norm(state_est-m2) - (b_est-m3(2,1))/norm(state_est-m3), ...
                    (c_est-m2(3,1))/norm(state_est-m2) - (c_est-m3(3,1))/norm(state_est-m3); ...
                    (a_est-m3(1,1))/norm(state_est-m3) - (a_est-m4(1,1))/norm(state_est-m4), ...
                    (b_est-m3(2,1))/norm(state_est-m3) - (b_est-m4(2,1))/norm(state_est-m4), ...
                    (c_est-m3(3,1))/norm(state_est-m3) - (c_est-m4(3,1))/norm(state_est-m4); ...
                    (a_est-m1(1,1))/norm(state_est-m1) - (a_est-m4(1,1))/norm(state_est-m4), ...
                    (b_est-m1(2,1))/norm(state_est-m1) - (b_est-m4(2,1))/norm(state_est-m4), ...
                    (c_est-m1(3,1))/norm(state_est-m1) - (c_est-m4(3,1))/norm(state_est-m4)];
 
 % Evaluate the distance error at the current parameters
        y_est = (1/c) * [norm(state_est-m1) - norm(state_est-m2); ...
                         norm(state_est-m2) - norm(state_est-m3); ...
                         norm(state_est-m3) - norm(state_est-m4); ...
                         norm(state_est-m1) - norm(state_est-m4)];
%%%%%%%%%%%%?? L - M ??%%%%%%%%%%%%%%%

        d = measurements-y_est;
 
 % compute the approximated Hessian matrix, J’ is the transpose of J
        H=J'*J;
 

        if it==1 % the first iteration : compute the total error
            e=abs(d(1))+abs(d(2))+abs(d(3))+abs(d(4));
            disp(e);
        end
    end
 
 % Apply the damping factor to the Hessian matrix
    H_lm=H+(lamda*eye(Nparams,Nparams));

 % Compute the updated parameters
    dp = (H_lm)\(J'*d);
    a_lm=a_est+dp(1,1);
    b_lm=b_est+dp(2,1);
    c_lm=c_est+dp(3,1);
    
%if b_lm < 1 %going towards direction of wrong minima
%    b_lm=b_est-10*dp(2,1);
%end
    state_lm = [a_lm;b_lm;c_lm];

 % Evaluate the total distance error at the updated parameters
    y_est_lm = (1/c) * [norm(state_lm-m1) - norm(state_lm-m2); ...
                         norm(state_lm-m2) - norm(state_lm-m3); ...
                         norm(state_lm-m3) - norm(state_lm-m4); ...
                         norm(state_lm-m1) - norm(state_lm-m4)];
    d_lm= measurements-y_est_lm;
    e_lm=dot(d_lm,d_lm);
    
 % check total distance error
 if e_lm<e
 disp(state_est);
 disp('update step');
 lamda=lamda/10;
 a_est=a_lm;
 b_est=b_lm;
 c_est=c_lm;
 state_est = [a_est; b_est; c_est];
 e=e_lm;
 disp(e);
 updateJ=1;
 else % otherwise increases the value of the damping factor
 disp('inc damping ');
 updateJ=0;
 lamda=lamda*10;
 end
end 

end

