function [R_OOS] = R_sq_oos(y_hat, y)
%Computes and returns the R squared out of sample

R_OOS = 1 - sum((y - y_hat).^2)/sum((y - mean(y)).^2);
end