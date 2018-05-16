function [] = suff_forecast()
%excecutes the sufficient forecasting procedure

p = 50;
T = 100;
H = 10;
K = 5;
L = 3;

[X,y] = simulate_linear(p,T);

%%%%%%%%
%STEP 1%
%%%%%%%%

%Calculate the factors using the description of fan 
X = X';

%calculate the eigenvectors
[V, D] = eig(X' * X);

[~,ind] = sort(diag(D));
ind = ind(end - K + 1:end);
eigenvectors = V(:,ind);

%Calculate F_hat and B_hat
F_hat = eigenvectors * sqrt(T);
B_hat = 1/T * X * F_hat;

%Check if the constraints still hold, throw error otherwise
if(~isdiag(round(B_hat' * B_hat)) || ~isdiag(round(F_hat' * F_hat)))
    disp("Constraints might not hold");
    return;
end

%%%%%%%%
%STEP 2%
%%%%%%%%

order_stats = zeros(T-1, K + 1);
order_stats(:,1) = y(2:T);
order_stats(:,2:end) = F_hat(1:T-1, :);

order_stats = sortrows(order_stats);

c = ceil((T-1)/H);
sigma_hat = zeros(K);

for h = 1:H
    if(h ~= H)
        temp = order_stats((h-1) * 10 + 1: h*c,2:end);
        temp = mean(temp);
        sigma_hat = sigma_hat + temp'  * temp;
    else
        temp = order_stats(h:end,2:end);
        temp = mean(temp);
        sigma_hat = sigma_hat + temp' * temp;
    end
end

sigma_hat = sigma_hat/H;


%%%%%%%%
%STEP 3%
%%%%%%%%
[V,D] = eig(sigma_hat);
[~,ind] = sort(diag(D));

ind = ind(end - L + 1: end);
psi = V(ind);

%%%%%%%%
%STEP 4%
%%%%%%%%
pred_ind = psi' * F_hat(T);

end