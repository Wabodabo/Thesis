num_sim = 500;

R = zeros(num_sim,1);

for i =1:num_sim
   R(i) = PCR_factors();
   disp(i);
end

mean(R)