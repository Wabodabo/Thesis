function [R] = corr_coeff(phi, phi_hat, H)
%This function calculates the squared multiple correlation coefficient to
%evaluate the effectiveness of phi_hat
%R is a constant value

num_phi = size(phi,2);
R = zeros(num_phi,1);

%Computes the rotated subspace
subspace = H \ phi;
subspace = subspace / norm(subspace);

for i = 1:num_phi
    R(i) = (phi_hat' * subspace(:,i))^2;
end

R = max(R);

end