function [H] = compute_H(F_hat, F, B, X)
%Computes an invertible matrix H such that the identififiability conditions
%are met

T = size(F_hat,1);

V = diag(eigs(X' * X / T, size(B,2)));
H = inv(V) * F_hat' * F * B' * B / T;
end