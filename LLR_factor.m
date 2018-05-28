function [f, R_squared] = LLR_factor(pred_ind, y)

T = length(y);
L = size(pred_ind,2);

% Gaussian kernel function
kerf=@(z)exp(-z.*z/2)/sqrt(2*pi);

%bowman and azzalini 1997 p.31
hx=median(abs(pred_ind-median(pred_ind)))/0.6745*(4/3/T)^0.2;
hy=median(abs(y-median(y)))/0.6745*(4/3/T)^0.2;
h=sqrt(hy*hx);


regr = ones(T, L + 1);
f = nan(T+1,1);

for i = 1:T
    W = kerf((pred_ind - pred_ind(i))/h);
    W = diag(W);
    
    regr(:,2:end) = pred_ind - pred_ind(i);
    
    f(i+1) = [1,0] * inv(regr' * W * regr) * regr' * W * y;
end

e = y(2:end) - f(2:end-1);
y_bar = mean(y(2:end));

SS_tot = sum((y(2:end) - y_bar).^2);
SS_res = sum(e.^2);

R_squared = 1 - SS_res/SS_tot;

end