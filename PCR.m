function [b, y_hat] = PCR(F_hat, y, K)
%This function applies PCA on X and calculates the PCR directions, it takes
%into account K factors. 

%The model assumes regresses y_t+1 on F_t

b = F_hat(1:end-1,:) \ y(2:end);
y_hat = F_hat(1:end-1,:) * b;

end