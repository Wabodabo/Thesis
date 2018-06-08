function [R_IS, R_OOS] = suff_forecast_nonlin(X, y, K, F)
%Performs the sufficient forecast method as in fan et al 2017 
%X and y should be the same size
%Returns the out of sample and in sample R squared

% X should be of size p * T
if(size(X,2) ~= size(y,1))
    disp('ERROR X and y not the same size');
    return;
end

%Variable definition
T_full = size(y,1);
T = T_full/2;
H = 10;
L = 2;

OOS_forecast = nan(T,1);

%Repeats the procedure T/2 times to compute all 1 step ahead forecasts
for i =1:T
    %Step 1
    %calculate the eigenvectors
    [eigenvectors, ~] = eigs(X(:,1:T + i - 1)' * X(:,1:T + i -1), K);

    %Calculate F_hat and B_hat
    F_hat = eigenvectors * sqrt(T + i -1);
    B_hat = 1/(T + i -1) * X(:,1:T + i -1) * F_hat;

    %Step 2
    c = ceil((T + i - 2)/H);

    %Order the factors
    order_f = zeros(T + i - 2, K + 1);
    order_f(:,1) = y(2:T + i - 1);
    order_f(:,2:end) = F_hat(1:T + i - 2, :);

    order_f = sortrows(order_f);

    sigma_hat_1 = zeros(K);

    %Calculate Sigma_1 based on factors
    cont = true;
    h = 1;
    while(cont)
        if(h*c < size(order_f,1))
            temp = order_f((h-1) * c + 1: h*c,2:end);
            temp = mean(temp);
            sigma_hat_1 = sigma_hat_1 + temp'  * temp;
            h = h +1;
        else
            temp = order_f(h:end,2:end);
            temp = mean(temp);
            sigma_hat_1 = sigma_hat_1 + temp' * temp;
            cont = false;
            H = h;
        end
    end

    % %Calculate Sigma_2 based on original data
    % order_x = zeros(T-1, p + 1);
    % order_x(:,1) = y(2:end);
    % order_x(:,2:end) = X(:,1:T-1)';
    % 
    % order_x = sortrows(order_x);
    % 
    % Lambda_hat = B_hat' * B_hat \ B_hat';
    % 
    % temp_2 = zeros(p);
    % 
    % for h = 1:H
    %     if(h ~= H)
    %         temp = order_x((h-1) * 10 + 1: h*c,2:end);
    %         temp = mean(temp);
    %         temp_2 = temp' * temp;
    %     else
    %         temp = order_x(h:end,2:end);
    %         temp = mean(temp);
    %         temp_2 = temp' * temp;
    %     end
    % end

    sigma_hat_1 = sigma_hat_1/H;
    H = 10;
    %sigma_hat_2 = Lambda_hat * temp_2/H * Lambda_hat';

    %Step 3
    [psi,~] = eigs(sigma_hat_1, L);

    %Step 4
    pred_ind = (psi' * F_hat')';

    %Step 5
%     if(i == 1)
%         [forecast,OOS_forecast(i), R_IS] = LLR_factor(pred_ind,y(1:T + i -1), true);
%     else
%         [~,OOS_forecast(i), ~] = LLR_factor(pred_ind, y(1:T+i-1), false);
%     end
    
end

R_OOS = 0;

end