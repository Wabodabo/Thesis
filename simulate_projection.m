function [x, y, F] = simulate_projection(T, alpha, ro, phi, g)
%This function simulates a time series as outlined in fan et al. 2017
%section 4.2

%alpha and ro should be size p and K, indicating the value of the AR(1)
%parameter of the factors and error terms respectively
%b should be of size [p,K] and indicates the value of the factor loadings

%Variables declared outside of simulation
K = size(alpha,1);
p = size(ro,1);

%Empty array declaration
F = zeros(T,K);
u = zeros(T,p);
x = zeros(T,p);
y = zeros(T,1);

for t=1:T-1
    %Compute the AR part of the factors and idiosyncratic part
    if(t == 1)
        F(t,:) = normrnd(0,1,[1,K]);
        u(t,:) = normrnd(0,1,[1,p]);      
    else
        F(t,:) = alpha' .* F(t-1,:) +  normrnd(0,1,[1,K]);
        u(t,:) = ro' .* u(t-1,:) + normrnd(0,1,[1,p]);             
    end
    
    %Compute the value of the predictors
    
    for j=1:p
        z = normrnd(0,1);
        b = g(z);
        x(t,j) = b * F(t,:)' + normrnd(0,1);
    end
    
    y(t+1) = F(t,1)* (F(t,2) + F(t,3) + 1) +  normrnd(0,1);
end

x = x - mean(x);


%test plot
hold off
t = 1:T;
plot(t,y);


end