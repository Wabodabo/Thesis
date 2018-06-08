function [b, F_hat] = PCR(X, y, K)
%This function applies PCA on X and calculates the PCR directions, it takes
%into account K factors. 

%The model assumes regresses y_t+1 on f_t

X = X';
[~, F_hat, ~] = pca(X, 'NumComponents', K);
b = F_hat(1:end-1,:) \ y(2:end);

end