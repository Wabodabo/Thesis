function [R_squared, R_OOS] = suff_forecast(p,T_full)
%excecutes the sufficient forecasting procedure

T = T_full/2;    %This size of the out of sample and in-sample periods
H = 10;
K = 5;
L = 1;

[X,y, real_factors] = simulate_linear(p,T_full); 
X = X';

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
    if(i == 1)
        [forecast,OOS_forecast(i), R_squared] = LLR_factor(pred_ind,y(1:T + i -1), real_factors, true);
    else
        [~,OOS_forecast(i), ~] = LLR_factor(pred_ind, y(1:T+i-1), real_factors, false);
    end
    
end

%Compute the out of sample R_squared
R_OOS = 1 - (sum((y(T+1:end) - OOS_forecast).^2))/sum((y(T + 1:end)- mean(y(T + 1:end))).^2);

%plot for test
t_full = 2:T_full;
t_in = 2:T;
t_out = T+1:T_full;

hold off
plot(t_full, y(2:end));
hold on

plot(t_in,forecast);
plot(t_out,OOS_forecast);
legend('DGP', 'in\_sample forecast', 'out\_sample forecast');


end