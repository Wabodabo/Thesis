num_sim = 500;

R = zeros(num_sim,1);

for i =1:num_sim
   R(i) = suff_forecast();
   disp(i);
end

mean(R)