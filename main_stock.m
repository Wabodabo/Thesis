function [R_IS, R_OS] = main_stock(X,y)
%This function returns the mean squared error both in sample and out of
%sample for SF, PCR and PC1

K = 7;
L = 2;
T = size(y,1);
test_sample = ceil(T/2);

X = X';

%empty matrices to save out of sample prediction
y_hat_SF_OS = zeros(test_sample,1);    
y_hat_SF2_OS = zeros(test_sample,1);    
y_hat_PCR_OS = zeros(test_sample,1);
y_hat_PC1_OS = zeros(test_sample,1);

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
    pred_ind_2 = (psi' * F_hat')';
    pred_ind_1 = pred_ind_2(:,1);

    %Perform local linear regerssion
    if(t ==1)
        [y_hat_SF_IS, y_hat_SF_OS(t), R_SF_IS] = LLR_factor(pred_ind_1, y_sample, true);
    else
        [~, y_hat_SF_OS(t), ~] = LLR_factor(pred_ind_1, y_sample, false);
    end    
    

    %Perform LOWESS regression
    f = fit(pred_ind_2, y_sample, 'Lowess', 'Normalize', 'on');
    if(t == 1)
       y_hat_SF2_IS = f(pred_ind_2(1:end-1, :));
       y_hat_SF2_OS(t) =  f(pred_ind_2(end,:));
       R_SF2_IS = R_sq(y_hat_SF2_IS, y_sample(2:end));
    else
       y_hat_SF2_OS(t) = f(pred_ind_2(end,:));         
    end
    
    %Perform PCR as benchmark
    [b_pcr, y_hat_pcr] = PCR(F_hat, y_sample, K);
    if(t ==1)
        R_PCR_IS = R_sq(y_hat_pcr, y_sample(2:end));
    end
    y_hat_PCR_OS(t) = F_hat(end,:) * b_pcr;

    %Compute PC1
    b_pc1 = F_hat(1:end-1,1) \ y_sample(2:end);
    if(t ==1)            
        R_PC1_IS = R_sq(F_hat(1:end-1,1) * b_pc1, y_sample(2:end));
    end
    y_hat_PC1_OS(t) = F_hat(end,1) * b_pc1;
 end

 R_SF_OS = R_sq_oos(y_hat_SF_OS, y((T - test_sample) + 1:end));
 R_PCR_OS = R_sq_oos(y_hat_PCR_OS, y((T -test_sample) + 1:end));
 R_PC1_OS = R_sq_oos(y_hat_PC1_OS, y((T - test_sample) + 1: end));
 R_SF2_OS = R_sq_oos(y_hat_SF2_OS, y((T - test_sample) + 1:end));

hold off
t = 1:size(y,1);
plot(t,y);
hold on
t_in = 2:test_sample;
t_out = ((size(y,1) - test_sample)+1):size(y,1);
plot(t_in, y_hat_SF_IS);
plot(t_out, y_hat_SF_OS);

R_IS = [R_SF_IS, R_SF2_IS, R_PCR_IS, R_PC1_IS];
R_OS = [R_SF_OS, R_SF2_OS, R_PCR_OS, R_PC1_OS];

end


    
