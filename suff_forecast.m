function [R_squared] = suff_forecast()
%excecutes the sufficient forecasting procedure

p = 50;
T = 500;
H = 10;
K = 5;
L = 1;

[X,y, ~] = simulate_linear(p,T);

%%%%%%%%
%STEP 1%
%%%%%%%%

%Calculate the factors using the description of fan 
X = X';

%calculate the eigenvectors
[eigenvectors, ~] = eigs(X' * X, K);

%Calculate F_hat and B_hat
F_hat = eigenvectors * sqrt(T);
B_hat = 1/T * X * F_hat;

%Check if the constraints still hold, throw error otherwise
if(~isdiag(round(B_hat' * B_hat)) || ~isdiag(round(F_hat' * F_hat)))
    disp('Constraints might not hold');
    return;
end

%%%%%%%%
%STEP 2%
%%%%%%%%

c = ceil((T-1)/H);

%Order the factors
order_f = zeros(T-1, K + 1);
order_f(:,1) = y(2:T);
order_f(:,2:end) = F_hat(1:T-1, :);

order_f = sortrows(order_f);

sigma_hat_1 = zeros(K);

%Calculate Sigma_1 based on factors
for h = 1:H
    if(h ~= H)
        temp = order_f((h-1) * 10 + 1: h*c,2:end);
        temp = mean(temp);
        sigma_hat_1 = sigma_hat_1 + temp'  * temp;
    else
        temp = order_f(h:end,2:end);
        temp = mean(temp);
        sigma_hat_1 = sigma_hat_1 + temp' * temp;
    end
end

%Calculate Sigma_2 based on original data
order_x = zeros(T-1, p + 1);
order_x(:,1) = y(2:end);
order_x(:,2:end) = X(:,1:T-1)';

order_x = sortrows(order_x);

Lambda_hat = B_hat' * B_hat \ B_hat';

temp_2 = zeros(p);

for h = 1:H
    if(h ~= H)
        temp = order_x((h-1) * 10 + 1: h*c,2:end);
        temp = mean(temp);
        temp_2 = temp' * temp;
    else
        temp = order_x(h:end,2:end);
        temp = mean(temp);
        temp_2 = temp' * temp;
    end
end

sigma_hat_1 = sigma_hat_1/H;
sigma_hat_2 = Lambda_hat * temp_2/H * Lambda_hat';


%%%%%%%%
%STEP 3%
%%%%%%%%
[psi,~] = eigs(sigma_hat_1, L);

%%%%%%%%
%STEP 4%
%%%%%%%%
pred_ind = (psi' * F_hat')';

%%%%%%%%
%STEP 5%
%%%%%%%%

% Gaussian kernel function
kerf=@(z)exp(-z.*z/2)/sqrt(2*pi);

h = 10;

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

%plot for test
t = 1:T;
hold off
plot(t,y);
hold on

scatter(t,pred_ind);
plot(t,f(1:T))
legend('DGP', 'predictive indices', 'fitted');


end