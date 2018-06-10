%This function performs a monte carlo simulation of the sufficient
%forecasting method on simulated non-linear data

p = 100;
T = 100;
K = 7;
num_sim = 1000;

%Parameters used for the DGP simulation
alpha = 0.2 + 0.6 * rand(K,1);
ro = 0.2 + 0.6 * rand(p,1);
B = normrnd(0,1, [p,K]);
phi = [1 0 0 0 0 0 0; 0 1/sqrt(2) 1/sqrt(2) 0 0 0 0]';

%Matrices to save data each simulation
R_corrcoeff = zeros(num_sim,1);
R_corrcoeff_PCR = zeros(num_sim,1);
R_IS = zeros(num_sim,1);
R_IS_PCR = zeros(num_sim,1);

for i = 1:num_sim
    [X,y, F, DGP] = simulate_interaction(T,alpha,ro,B, phi);
    X = X';

    [F_hat, psi] = predict_indices_nonlin(X,y,K);
    
    [phi_PCR, y_hat_PCR] = PCR(F_hat, y,K);    
    R_IS_PCR(i) = R_sq(y_hat_PCR, y(2:end));

    %Computing the correlation coef ficient
    H = compute_H(F_hat, F, B, X);
    R_corrcoeff_PCR(i) =  corr_coeff(phi_PCR, b, H);    
    R_corrcoeff(i) = corr_coeff(phi, psi(:,1), H);
    
    %Regress y on predictive indices
    pred_ind = (psi' * F_hat')';
    
    pred_ind = pred_ind(1:end-1,:);
    y = y(2:end);
    regr = [ones(T-1,1) pred_ind(:,1) pred_ind(:,2) pred_ind(:,1) .* pred_ind(:,2)];
    
    b =  regr \ y;
    y_hat = regr * b;
    
    hold off
    scatter(y,DGP(1:end-1));
    hold on
    pred = pred_ind(:,1) + pred_ind(:,1) .* pred_ind(:,2);
    pred2 = pred_ind(:,2) + pred_ind(:,2) .* pred_ind(:,1);
    scatter(y,pred2);
    scatter(y,pred);
   
    R_IS(i) = R_sq(y_hat, y);

    %test plot
    hold off
    t = 1:T-1;
    plot(t,y);
    hold on
    plot(t,y_hat_PCR);
    
      
    disp(i);
end

