%This script simulates the PCR and PC1 regression 
%for every combination of p and T and saves the R_squared to an excel file

p_T = [50,100;50,200;100,100;100,500;500,100;500,500];
num_sim = 1000;

%Create the file where the data will be stored
file_name = 'Simulation_output.xlsx';
sheet_name = 'Linear Simulation';


for j=1:3
   
    p = p_T(j,1);
    T = p_T(j,2);
    
    R = zeros(num_sim,4);

    for i =1:num_sim
       [R(i,1), R(i,2), R(i,3), R(i,4)] = PCR_factors(p,T);
       disp([j, i]);
    end
    
    R = mean(R);
    
    xlswrite(file_name, R(:,1:2), sheet_name, strcat('D', num2str(j+1)));
    xlswrite(file_name, R(:,3:4), sheet_name, strcat('G', num2str(j+1)));
    
end

