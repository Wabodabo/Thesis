%function handle
dgp = @(z) (2 * z)./(1 + z.^2);

T = 200;

x = zeros(T,1);
x(1) = 2 * rand() -1;
t = 1:T;
    

for i = 2:T
   x(i) = dgp(x(i-1)) + (2 * rand -1);
end

forecast = LLR_factor(x,x);

hold off
ax1 = subplot(2,1,1);
plot(t,x);


ax2 = subplot(2,1,2);
y = x(2:end);
x = x(1:end-1);
scatter(x,y);
t_dgp = linspace(-2,2,T);
hold on
plot(t_dgp,dgp(t_dgp));
scatter(x,forecast);

legend('Data points', 'True DGP', 'nonparametric estimation', 'location', 'northwest');
