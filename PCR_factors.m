function [R_squared,R_squared_PC1, R_OOS_PCR, R_OOS_PC1] = PCR_factors(p,T)

H = 10;
K = 5;
L = 1;
sample = T/2;

%Simulate Data
[X,y, ~] = simulate_linear(p,T); 

fc_PCR = zeros(sample,1);
fc_PC1 = zeros(sample,1);

for t = 1:T/2
    %Apply principal component analysis
    [~, factors,~] = pca(X(1:sample,:), 'NumComponents', K);
    
    %Calculating the forecasts using PCR
    %Select relevant nr of factors (assumed to be known
    y_sample = y(2:sample);

    beta_hat = factors(1:sample-1,:) \ y_sample;

    y_pred_PCR = (beta_hat' * factors')';
    fc_PCR(t) = y_pred_PCR(end);

    %Estimate Beta when only using the first principal component
    factors_PC1 = factors(:,1);

    beta_PC1 = factors_PC1(1:end-1,1) \ y_sample;

    y_pred_PC1 = (beta_PC1' * factors_PC1')';
    fc_PC1(t) = y_pred_PC1(end);
    
    if(t == 1)
        %Calculating the R_squared of the PCR part
        e = y(2:T/2) - y_pred_PCR(1:end-1);
        y_bar = mean(y(1:T/2));

        SS_tot = sum((y(1:T/2) - y_bar).^2);
        SS_res = sum(e.^2);

        R_squared = 1 - SS_res/SS_tot;

        %Calculating the R_squared of the single PC
        e = y(2:T/2) - y_pred_PC1(1:end-1);
        SS_res = sum(e.^2);

        R_squared_PC1 = 1 - SS_res/SS_tot;
    end    
    
    sample = sample + 1;
end

%Calculate the R squared OOS for the PCR model
R_OOS_PCR = 1 - (sum((y(T/2 + 1:end) - fc_PCR).^2))/sum((y(T/2 + 1:end)- mean(y(T/2 + 1:end))).^2);

%Calculate the R squared IS for the PC1 model
R_OOS_PC1 = 1 - (sum((y(T/2+1:end) - fc_PC1).^2))/sum((y(T/2 + 1:end)- mean(y(T/2 + 1:end))).^2);

% t = 1:T-1;
% hold off
% plot(t,y);
% hold on
% plot(t,y_pred_PCR);
% legend('Real data', 'in-sample prediction');

end