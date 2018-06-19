%This function runs main_stock for different target variables. The index of
%the target variable will be given in a vector, all variables are part of
%stockwatson.mat

%load data
load('stockwatson.mat');

%save data to file
file_name = 'Simulation_output.xlsx';
sheet_name = 'stock';

%parameters
targets = [5, 30, 40, 52, 104, 79, 54, 96, 101, 102, 74];
num_targets = length(targets);

%empty arrays to save
R_IS = zeros(num_targets, 4);
R_OOS = zeros(num_targets, 4);

for i = 1:num_targets
    y = stockwatson(:,targets(i));
    X = stockwatson(:, [1:targets(i)-1, targets(i) + 1:end]);
   [R_IS(i,:), R_OOS(i,:)] = main_stock(X, y);
   
   disp(num2str(i));
end

xlswrite(file_name, [R_IS, R_OOS], sheet_name, 'C3');
