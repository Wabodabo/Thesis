
p = [100,100,100,500,500,500];
T = [100,200,500,100,200,500];
num_sim = 5;

file_name = 'simulation_output.xlsx';
sheet_name = 'Non-linear Simulation';


for i = 1:6
    disp(strcat('STARTING Non-Linear RUN: ', num2str(i)));
    [R_IS, R_OS] = main_nonlin(p(i),T(i),num_sim);

    xlswrite(file_name, R_IS, sheet_name, strcat('c', num2str(i+2)));
    xlswrite(file_name, R_OS, sheet_name, strcat('f', num2str(i+2)));
end
