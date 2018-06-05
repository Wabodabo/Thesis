function [X, y, f_real] = generate_interaction(p,T)
%This function simulates a time series as outlined in fan et al. 2017
%section 4.2

%Variables declared outside of simulation
K = 7;

%Empty array declaration
f = zeros(T,K);
u = zeros(T,p);
x = zeros(T,p);
y = zeros(T,1);

f_real = zeros(T,K);

%Parameter values 
alpha = 0.2 + 0.6 * rand(K,1);
ro = 0.2 + 0.6 * rand(p,1);
phi = [1 0 0 0 0 0 0; 0 1/sqrt(2) 1/sqrt(2) 0 0 0 0]';
b = normrnd(0,1, [p,K]);


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
    
    y(t+1) = f(t,1)* (f(t,2) + f(t,3) + 1) +  normrnd(0,1);
end

x = x - mean(x);


%test plot
hold off
t = 1:200;
plot(t,y);

end