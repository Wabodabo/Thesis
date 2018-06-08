function [] = RESET(x,y)
%This test performs the Ramsey reset test as described in Econometric
%methods by Christian Heij chapter 5.2.2

if(size(x,1) ~= size(y,1))
    disp('x and y do not have the same dimensions');
    return;
end

t = size(y,1);

b = x\y;
polynomials = zeros(t,2);
polynomials(:,1) = (x * b).^2;
polynomials(:,2) = (x * b).^3;

gamma = [x, polynomials] \ y;
gamma = gamma(end-1:end);

end
