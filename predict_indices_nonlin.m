function [F_hat, psi] = predict_indices_nonlin(X, y, K, B_hat)
%This function calculates the predictive indices and returns them

% X should be of size p * T
if(size(X,2) ~= size(y,1))
    disp('ERROR X and y not the same size');
    return;
end

%Variable definition
T = size(y,1);
H = 10;
L = 2;
p = size(X,1);


%Step 1
%calculate the eigenvectors
[eigenvectors, ~] = eigs(X' * X, K);

%Calculate F_hat and B_hat
F_hat = eigenvectors * sqrt(T);
B_hat = (1/T) *  X * F_hat;

%Step 2
c = ceil((T-1)/H);

%Order the factors
order_f = zeros(T-1, K + 1);
order_f(:,1) = y(2:T);
order_f(:,2:end) = F_hat(1:T-1, :);

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
        temp = order_f((h-1)*c+1:end,2:end);
        temp = mean(temp);
        sigma_hat_1 = sigma_hat_1 + temp' * temp;
        cont = false;
        H1 = h;
    end
end

% %Calculate Sigma_2 based on original data
order_x = zeros(T-1, p + 1);
order_x(:,1) = y(2:end);
order_x(:,2:end) = X(:,1:T-1)';

order_x = sortrows(order_x);

Lambda_hat = B_hat' * B_hat \ B_hat';

temp_2 = zeros(p);
cont = true;
h = 1;
while(cont)
    if(h*c < size(order_x,1))
        temp = order_x((h-1) * 10 + 1: h*c,2:end);
        temp = mean(temp);
        temp_2 = temp_2 + temp' * temp;
        h = h+1;
    else
        temp = order_x((h-1) * c + 1:end,2:end);
        temp = mean(temp);
        temp_2 = temp_2 + temp' * temp;
        cont = false;
        H2 = h;
    end
end

sigma_hat_1 = sigma_hat_1/H1;
sigma_hat_2 = Lambda_hat * (temp_2/H2) * Lambda_hat';

%Step 3
[psi,~] = eigs(sigma_hat_1, L);

end
