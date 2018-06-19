%This function runs the main_goyal method with different split dates of the
%sample, it returns in sample R_squared and out of sample R_squared for
%every split date

%load data
load('goyalwelch.mat');

%save data to file
file_name = 'Simulation_output.xlsx';
sheet_name = 'goyal';

%parameters
start_date = 168;
end_date = 504;
num_obs = end_date - start_date;

%empty arrays to save
R_IS = zeros(num_obs, 4);
R_OOS = zeros(num_obs, 4);

for i = start_date:171
   [R_IS(i - 167,:), R_OOS(i - 167,:)] = main_goyal(X, y, i, 3);
   
   disp(num2str(i));
end

xlswrite(file_name, [R_IS, R_OOS], 'goyal', strcat('K', num2str(start_date-165)));
  