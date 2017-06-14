function [ state_est, x ] = LevMarFunc( measurements, x )
%% Microphone variables
 c = 34.3;
 %c = 0.034;
 mdist = 20;
 m1 = [0;0];
 m2 = [mdist*-0.5; mdist*-0.8660254];
 m3 = [mdist*0.5; mdist*-0.8660254];   


%% LevMar variables
Nparams=2; % defines size of lambda matrix
n_iters = 100; % set # of iterations for the LM
lamda= 0.01; % set an initial value of the damping factor for the LM
updateJ=1;
a_est= 1;
b_est= 1; 
state_est = [a_est;b_est]; % intial state estimate (maybe can be optimized)

%% LevMar Algorithm
for it=1:n_iters
    if updateJ==1
     
 % Evaluate the Jacobian matrix at the current parameters (a_est, b_est)

        J = (1/c) * [(a_est-m1(1,1))/norm(state_est-m1) - (a_est-m2(1,1))/norm(state_est-m2), ...
                    (b_est-m1(2,1))/norm(state_est-m1) - (b_est-m2(2,1))/norm(state_est-m2); ...
                    (a_est-m2(1,1))/norm(state_est-m2) - (a_est-m3(1,1))/norm(state_est-m3), ...
                    (b_est-m2(2,1))/norm(state_est-m2) - (b_est-m3(2,1))/norm(state_est-m3); ...
                    (a_est-m1(1,1))/norm(state_est-m1) - (a_est-m3(1,1))/norm(state_est-m3), ...
                    (b_est-m1(2,1))/norm(state_est-m1) - (b_est-m3(2,1))/norm(state_est-m3)];
 
 % Evaluate the distance error at the current parameters
        y_est = (1/c) * [norm(state_est-m1) - norm(state_est-m2); ...
                         norm(state_est-m2) - norm(state_est-m3); ...
                         norm(state_est-m1) - norm(state_est-m3)];

        d = measurements - y_est;
 
 % compute the approximated Hessian matrix, J’ is the transpose of J
        H=J'*J;
 
        if it==1 % the first iteration : compute the total error
            e=dot(d,d);
        end
    end
 
 % Apply the damping factor to the Hessian matrix
    H_lm=H+(lamda*eye(Nparams,Nparams));

 % Compute the updated parameters
    dp = (H_lm)\(J'*d);
    a_lm=a_est+dp(1,1);
    b_lm=b_est+dp(2,1);
    
if b_lm < 1 %going towards direction of wrong minima
    b_lm=b_est-10*dp(2,1);
end
    state_lm = [a_lm;b_lm];

 % Evaluate the total distance error at the updated parameters
    y_est_lm = (1/c) * [norm(state_lm-m1) - norm(state_lm-m2); ...
                         norm(state_lm-m2) - norm(state_lm-m3); ...
                         norm(state_lm-m1) - norm(state_lm-m3)];
    d_lm= measurements-y_est_lm;
    e_lm=dot(d_lm,d_lm);
    
 % check total distance error
 x(:,it)=state_est;
 if e_lm<e
 disp('update step');
 lamda=lamda/10;
 a_est=a_lm;
 b_est=b_lm;
 state_est = [a_est; b_est];
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

