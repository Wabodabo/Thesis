function [x, y, F, DGP] = simulate_interaction(T, alpha, ro, b, phi)
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

error = zeros(T,1);

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
        x(t,j) = b(j,:) * F(t,:)' + normrnd(0,1);
    end
    
    error(t) = normrnd(0,1);
    y(t+1) = F(t,1)* (F(t,2) + F(t,3) + 1) +  error(t);
end

x = x - mean(x);


%test plot
% hold off
% t = 1:T;
% plot(t,y);

%Scatter plot of predictive indices hold off
hold off
pred_ind = zeros(T,2);
pred_ind(:,1) = (phi(:,1)' * F')';
pred_ind(:,2) = (phi(:,2)' * F')';
DGP = pred_ind(:,1) .* pred_ind(:,2) + pred_ind(:,1);

%scatter(y(2:end),pred(1:end-1,:));

end