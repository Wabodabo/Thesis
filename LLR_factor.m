function [forecast, R_squared] = LLR_factor(pred_ind, y)

T = length(y);
L = size(pred_ind,2);

% Epanechnikov kernel function
%kerf=@(z) 3/4 * (1-z).^2;

%Gaussian kernel funcetion
kerf=@(z)exp(-z.*z/2)/sqrt(2*pi);

% bandwidth selection, bowman and azzalini 1997 p.31
hx=median(abs(pred_ind-median(pred_ind)))/0.6745*(4/3/T)^0.2;
hy=median(abs(y-median(y)))/0.6745*(4/3/T)^0.2;
h=sqrt(hy*hx);


%initializing arrays
regr = ones(T - 1, L + 1);
forecast = nan(T-1,1);
y_pred = y(2:end);

for i = 1:T-1
    W = kerf((pred_ind(1:end-1,:) - pred_ind(i))/h);
    W = diag(W);
    
    regr(:,2:end) = pred_ind(1:end-1,:) - pred_ind(i); 
    
    beta_hat = inv(regr' * W * regr) * regr' * W * y_pred;
    
    forecast(i) = beta_hat' * [1; 0]; 
    %f(i+1) = [1,0] * inv(regr' * W * regr) * regr' * W * y(i+1);
end

e = y_pred - forecast;
y_bar = mean(y_pred);

SS_tot = sum((y_pred - y_bar).^2);
SS_res = sum(e.^2);

R_squared = 1 - SS_res/SS_tot;

%test plots
hold off
scatter(pred_ind(1:end-1),y_pred);
hold on
tmp = [pred_ind(1:end-1) forecast];
tmp = sort(tmp);
plot(tmp(:,1),tmp(:,2));
scatter(pred_ind(1:end-1), forecast);
end