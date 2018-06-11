function [x,y, F, F_real] = simulate_linear(p,T, alpha, ro, b, phi)
%Simulate data for linear forecasting (fan et al 2017 chapter 4.1)

%Variables declared outside of simulation
K = 5;

%Empty array declaration
F = zeros(T,K);
u = zeros(T,p);
x = zeros(T,p);
y = zeros(T,1);

F_real = zeros(T,K);

%Possible change in Sigma_y
%sigma_y = sqrt(sum(1./(1 - alpha.^2)));
sigma_y = 1;

for t=1:T-1
    %Compute the AR part of the factors and idiosyncratic part
    if(t == 1)
        F(t,:) = normrnd(0,1,[1,K]);
        u(t,:) = normrnd(0,1,[1,p]);
        
        F_real(t,:) = F(t,:);
    else
        F(t,:) = alpha' .* F(t-1,:) +  normrnd(0,1,[1,K]);
        u(t,:) = ro' .* u(t-1,:) + normrnd(0,1,[1,p]);
        
        F_real(t,:) = alpha' .* F(t-1,:);
    end
    
    %Compute the value of the predictors
    for j=1:p
        x(t,j) = b(j,:) * F(t,:)' + normrnd(0,1);
    end
    
    y(t+1) = phi' * F(t,:)' + sigma_y * normrnd(0,1);
end

x = x - mean(x);

% %test plot
% hold off
% indices = (phi' * f')';
% scatter(indices(1:end-1),y(2:end));
% legend('Real predictive indices');

end

