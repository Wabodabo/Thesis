function [forecast_in_sample, forecast_out_sample, R_squared] = LLR_factor(pred_ind, y, in_sample)
%forecast_in_sample:    Returns a vector of in sample forecasts
%forecast_out_sample:   Returns a single value of the one step ahead forecast
%R_squared:             Returns the R_squared of the in_sample forecasts

forecast_in_sample = 0;
R_squared = 0;

T = length(pred_ind);
L = size(pred_ind,2);

% Epanechnikov kernel function
%kerf=@(z) (abs(z) <= 1) .* 3/4 .* (1-z).^2 + (abs(z) > 1) .* 0;

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


%Does the local linear regression for every point if in_sample, otherwise
%just the last
if(in_sample)
    for i = 1:T
        W = kerf((pred_ind(1:end-1,:) - pred_ind(i))/h);
        W = diag(W);

        regr(:,2:end) = pred_ind(1:end-1,:) - pred_ind(i); 

        beta_hat = (regr' * W * regr) \ regr' * W * y_pred;

        forecast(i) = beta_hat(1);
    end

    forecast_in_sample = forecast(1:end-1);
    forecast_out_sample = forecast(end);

    R_squared = R_sq(forecast_in_sample, y(2:end));
else
    W = kerf((pred_ind(1:end-1,:)  - pred_ind(end))/h);
    W = diag(W);

    regr(:,2:end) = pred_ind(1:end-1,:) - pred_ind(end); 

    beta_hat = (regr' * W * regr) \ regr' * W * y_pred;

    forecast_out_sample = beta_hat(1); 
    
end

end