function [x,y] = simulate_linear(p,T)
%Simulate data for linear forecasting (fan et al 2017 chapter 4.1)

%Variables declared outside of simulation
K = 5;

%Empty array declaration
f = zeros(T,K);
u = zeros(T,p);
x = zeros(T,p);
y = zeros(T,1);

%Parameter values 
alpha = 0.2 + 0.6 * rand(K,1);
ro = 0.2 + 0.6 * rand(p,1);
phi = [0.8 0.5 0.3 0 0]';
b = normrnd(0,1, [p,K]);
sigma_y = 1;            %sigma checken!!!!!

for t=1:T
    %Compute the AR part of the factors and idiosyncratic part
    if(t == 1)
        f(t,:) = normrnd(0,1,[1,K]);
        u(t,:) = normrnd(0,1,[1,p]);
    else
        f(t,:) = alpha' .* f(t-1,:) +  normrnd(0,1,[1,K]);
        u(t,:) = ro' .* u(t-1,:) + normrnd(0,1,[1,p]);
    end
    
    %Compute the value of the predictors
    for j=1:p
        x(t,j) = b(j,:) * f(t,:)' + normrnd(0,1);
    end
    
    y(t) = phi' * f(t,:)' + sigma_y * normrnd(0,1);
end

end

