function [R_IS, R_OS] = main_nonlin(p,T,num_sim)

%This function performs a monte carlo simulation of the sufficient
%forecasting method on simulated non-linear data

K = 7;
test_sample = T/2;

%Parameters used for the DGP simulation
alpha = 0.2 + 0.6 * rand(K,1);
ro = 0.2 + 0.6 * rand(p,1);
phi = [1 0 0 0 0 0 0; 0 1/sqrt(2) 1/sqrt(2) 0 0 0 0]';

%Matrices to save data each simulation
R_corrcoeff = zeros(num_sim,2);
R_corrcoeff_PCR = zeros(num_sim,1);
R_IS_SFi = zeros(num_sim,1);
R_IS_PCR = zeros(num_sim,1);
R_IS_PCRi = zeros(num_sim,1);

R_OS_SFi = zeros(num_sim,1);
R_OS_PCR = zeros(num_sim,1);
R_OS_PCRi = zeros(num_sim,1);

for i = 1:num_sim
    %empty vectors to save out of sample predictions
    y_hat_os_SFi = zeros(test_sample,1);
    y_hat_os_PCR = zeros(test_sample,1);
    y_hat_os_PCRi = zeros(test_sample,1);
    
    B = normrnd(0,1, [p,K]);
    [X,y, F, ~] = simulate_interaction(T,alpha,ro,B, phi);
    X = X';
    
    for t = 1:test_sample
        X_sample = X(:,1:test_sample + t -1);
        y_sample = y(1:test_sample + t -1);
        
        %Compute the predictive indices
        [F_hat, psi] = predict_indices_nonlin(X_sample,y_sample,K);
        
        %Regress y on predictive indices
        pred_ind = (psi' * F_hat')';
        regr = [pred_ind(:,1) pred_ind(:,2) pred_ind(:,1) .* pred_ind(:,2)];
        
        %make in sample predictions for SFi, PCR and PCi
        b =  regr(1:end-1,:) \ y_sample(2:end);
        y_hat_SFi = regr(1:end-1,:) * b;
        [phi_PCR, y_hat_PCR] = PCR(F_hat, y_sample,K); 
        regr_PCi = [F_hat, F_hat(:,1) .* F_hat(:,2)];
        [phi_PCi, y_hat_PCRi] = PCR(regr_PCi, y_sample, K);
        if(t == 1)        
            R_IS_SFi(i) = R_sq(y_hat_SFi, y_sample(2:end)); 
            R_IS_PCR(i) = R_sq(y_hat_PCR, y_sample(2:end));
            R_IS_PCRi(i) = R_sq(y_hat_PCRi, y_sample(2:end));
        end
        
        %Make out of sample predictions for SFi, PCR and PCi
        y_hat_os_SFi(t) = regr(end,:) * b;    
        y_hat_os_PCR(t) = F_hat(end,:) * phi_PCR;
        y_hat_os_PCRi(t) = regr_PCi(end,:) * phi_PCi;
        
        %Computing the correlation coefficient
%         H = compute_H(F_hat, F, B, X);
%         R_corrcoeff_PCR(i) =  corr_coeff(phi_PCR, phi_PCR, H);    
%         R_corrcoeff(i,1) = corr_coeff(phi, psi(:,1), H);
%         R_corrcoeff(i,2) = corr_coeff(phi, psi(:,2), H);

     

        %Plot the predictive indices
%         hold off
%         scatter(y,DGP(1:end-1));
%         hold on
%         pred = pred_ind(:,1) + pred_ind(:,1) .* pred_ind(:,2);
%         pred2 = pred_ind(:,2) + pred_ind(:,2) .* pred_ind(:,1);
%         scatter(y,pred2);
%         scatter(y,pred);

    end
    
    R_OS_SFi(i) = R_sq_oos(y_hat_os_SFi, y(test_sample + 1:end));
    R_OS_PCR(i) = R_sq_oos(y_hat_os_PCR, y(test_sample + 1:end));
    R_OS_PCRi(i) = R_sq_oos(y_hat_os_PCRi, y(test_sample + 1:end));
    
%     %test plot
    hold off
    t_in = 2:test_sample;
    t_out = test_sample + 1: T;
    plot([t_in t_out], y(2:end));
    hold on
    plot(t_in, y_hat_SFi(1:test_sample - 1));
    plot(t_out, y_hat_os_SFi);
          
    disp(i);
end

R_IS = [mean(R_IS_SFi), mean(R_IS_PCR), mean(R_IS_PCRi)];
R_OS = [mean(R_OS_SFi), mean(R_OS_PCR), mean(R_OS_PCRi)];

end