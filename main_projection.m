%This function runs the sufficient forecasting as described in section 4.3

p = 100;
T = 100;
K = 3;
num_sim = 100;

%Parameters used for the DGP simulation
alpha = 0.2 + 0.6 * rand(K,1);
ro = 0.2 + 0.6 * rand(p,1);
phi = [1 0 0 0 0 0 0; 0 1/sqrt(2) 1/sqrt(2) 0 0 0 0]';

g = @(z) [z, z^2-1, z^3-2*z];


for i = 1:num_sim
    
    [X,y, F] = simulate_projection(T, alpha, ro, phi, g);
    
end

