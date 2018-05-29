function [R_squared] = PCR_factors()

p = 50;
T = 200;
H = 10;
K = 5;
L = 1;

%Simulate Data
[X,y, ~] = simulate_linear(p,T); 

%Apply principal component analysis
[coeff, score,latent] = pca(X);

factors = score(1:end-1,1:5);
y = y(2:end);

beta_hat = factors \ y;

y_pred = (beta_hat' * factors')';

e = y - y_pred;
y_bar = mean(y);

SS_tot = sum((y - y_bar).^2);
SS_res = sum(e.^2);

R_squared = 1 - SS_res/SS_tot;


t = 1:T-1;
hold off
plot(t,y);
hold on
plot(t,y_pred);
legend('Real data', 'in-sample prediction');

end