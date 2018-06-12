function [R_IS, R_OS] = main_linear(p,T, num_sim)
%This function returns the mean squared error both in sample and out of
%sample for SF, PCR and PC1

K = 5;
L = 1;
test_sample = T/2;

%Parameters used for the DGP simulation
alpha = 0.2 + 0.6 * rand(K,1);
ro = 0.2 + 0.6 * rand(p,1);
phi = [0.8; 0.5; 0.3; 0; 0];

%Matrices to keep track
R_SF_IS = zeros(num_sim,1);
R_PCR_IS = zeros(num_sim,1);
R_PC1_IS = zeros(num_sim,1);
R_SF_OS = zeros(num_sim,1);
R_PCR_OS = zeros(num_sim,1);
R_PC1_OS = zeros(num_sim,1);

for i=1:num_sim
     %empty matrices to save out of sample prediction
     y_hat_SF_OS = zeros(test_sample,1);    
    y_hat_PCR_OS = zeros(test_sample,1);
    y_hat_PC1_OS = zeros(test_sample,1);
    
    
     B = normrnd(0,1, [p,K]);
     [X, y, ~, ~] = simulate_linear(p, T, alpha, ro, B, phi);
     X = X';
     for t = 1:test_sample
        X_sample = X(:, 1:test_sample + t -1);
        y_sample = y(1:test_sample + t -1);
        %Obtain factors
        [eigenvectors, ~] = eigs(X_sample' * X_sample / (test_sample + t -1), K);
        F_hat = eigenvectors * sqrt(test_sample + t -1);
        B_hat = (1/(test_sample + t -1)) * X_sample * F_hat;

        %Obtain the sliced covariance matrix
        Sigma = sliced_covariance(F_hat, X_sample, y_sample, B_hat);

        %Obtain psi
        [psi, ~] = eigs(Sigma, L);

        %Compute predictive indices
        pred_ind = (psi' * F_hat')';

        %Perform local linear regerssion
        if(t ==1)
            [~, y_hat_SF_OS(t), R_SF_IS(i)] = LLR_factor(pred_ind, y_sample, true);
        else
            [~, y_hat_SF_OS(t), ~] = LLR_factor(pred_ind, y_sample, false);
        end    

        %Perform PCR as benchmark
        [b_pcr, y_hat_pcr] = PCR(F_hat, y_sample, K);
        if(t ==1)
            R_PCR_IS(i) = R_sq(y_hat_pcr, y_sample(2:end));
        end
        y_hat_PCR_OS(t) = F_hat(end,:) * b_pcr;
        
        %Compute PC1
        b_pc1 = F_hat(1:end-1,1) \ y_sample(2:end);
        if(t ==1)            
            R_PC1_IS(i) = R_sq(F_hat(1:end-1,1) * b_pc1, y_sample(2:end));
        end
        y_hat_PC1_OS(t) = F_hat(end,1) * b_pc1;
     end
    
     R_SF_OS(i) = R_sq_oos(y_hat_SF_OS, y(test_sample + 1:end));
     R_PCR_OS(i) = R_sq_oos(y_hat_PCR_OS, y(test_sample + 1:end));
     R_PC1_OS(i) = R_sq_oos(y_hat_PC1_OS, y(test_sample + 1: end));
    
    %test plots
%     hold off
%     scatter(pred_ind(1:end-1),y(2:end));
%     hold on
%     tmp = [pred_ind(1:end-1) y_hat_IS];
%     tmp = sortrows(tmp);
%     plot(tmp(:,1),tmp(:,2));
%     real_ind = phi' * F_real';
%     scatter(pred_ind(1:end-1), real_ind(1:end-1));
%     legend('Estimated predictive indices', 'Estimated link function', 'Real predictive indices');
%     xlabel('Value of predictive index');
%     ylabel('Value of target variable');

    disp(i);
end

R_IS = [mean(R_SF_IS), mean(R_PCR_IS), mean(R_PC1_IS)];
R_OS = [mean(R_SF_OS), mean(R_PCR_OS), mean(R_PC1_OS)];

end


    
