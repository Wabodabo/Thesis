function [x,y, f_real] = simulate_linear(p,T)
%Simulate data for linear forecasting (fan et al 2017 chapter 4.1)

%Variables declared outside of simulation
K = 5;

%Empty array declaration
f = zeros(T,K);
u = zeros(T,p);
x = zeros(T,p);
y = zeros(T,1);

f_real = zeros(T,K);

%Parameter values 
alpha = 0.2 + 0.6 * rand(K,1);
ro = 0.2 + 0.6 * rand(p,1);
phi = [0.8 0.5 0.3 0 0]';
b = normrnd(0,1, [p,K]);

%Possible change in Sigma_y
%sigma_y = sqrt(sum(1./(1 - alpha.^2)));
sigma_y = 1;

for t=1:T-1
    %Compute the AR part of the factors and idiosyncratic part
    if(t == 1)
        f(t,:) = normrnd(0,1,[1,K]);
        u(t,:) = normrnd(0,1,[1,p]);
        
        f_real(t,:) = f(t,:);
    else
        f(t,:) = alpha' .* f(t-1,:) +  normrnd(0,1,[1,K]);
        u(t,:) = ro' .* u(t-1,:) + normrnd(0,1,[1,p]);
        
        f_real(t,:) = alpha' .* f(t-1,:);
    end
    
    %Compute the value of the predictors
    for j=1:p
        x(t,j) = b(j,:) * f(t,:)' + normrnd(0,1);
    end
    
    y(t+1) = phi' * f(t,:)' + sigma_y * normrnd(0,1);
end

x = x - mean(x);

% %test plot
% hold off
% indices = (phi' * f')';
% scatter(indices(1:end-1),y(2:end));
% legend('Real predictive indices');

end

