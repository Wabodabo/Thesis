
p = [50,50,100,100,500,500];
T = [100,200,100,500,100,500];
num_sim = 5;

file_name = 'simulation_output.xlsx';
sheet_name = 'Linear Simulation';


for i = 5:6
    disp(strcat('STARTING Linear RUN: ', num2str(i)));
    [R_IS, R_OS] = main_linear(p(i),T(i),num_sim);

    xlswrite(file_name, R_IS, sheet_name, strcat('c', num2str(i+2)));
    xlswrite(file_name, R_OS, sheet_name, strcat('f', num2str(i+2)));
end
