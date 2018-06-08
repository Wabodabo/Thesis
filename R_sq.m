function [R_2] = R_sq(forecast, y)
%Computes the in sample R squared between the forecasts and y

if(size(forecast,1) ~= size(y,1))
   error('Size of the forecast and y are not the same');
end

e = y - forecast;
y_bar = mean(y);

SS_tot = sum((y - y_bar).^2);
SS_res = sum(e.^2);

R_2 = 1 - SS_res/SS_tot;

end