%This script simulates the sufficient forecasting a certain amount of times
%for every combination of p and T and saves the R_squared to an excel file

p_T = [50,100;50,200;100,100;100,500;500,100;500,500];
num_sim = 1000;

%Create the file where the data will be stored
file_name = 'Simulation_output.xlsx';
sheet_name = 'Linear Simulation';

for j = 1:6
    p = p_T(j,1);
    T = p_T(j,2);
    
    R_IS = zeros(num_sim,1);
    R_OOS = zeros(num_sim,1);

    for i =1:num_sim
       [R_IS(i), R_OOS(i)] = suff_forecast(p,T);
       disp(i);
    end

    output = [p,T;mean(R_IS),mean(R_OOS)];
    xlswrite(file_name, output, sheet_name, strcat('B', int2str(2*(j) - 1)));
        
    p
    T
    disp('The mean R\_squared is:');
    mean(R_IS)
    disp('The mean R\_squared out of sample is: ');
    mean(R_OOS)
    
end