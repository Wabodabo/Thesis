

for i = 1:108
   if(~adftest(stockwatson(:,i)))
       if(adftest(log(stockwatson(:,i))))
          stockwatson(:,i) = log(stockwatson(:,i)); 
          disp(strcat(num2str(i), ' log '));
       else
           stockwatson(2:end, i) = stockwatson(2:end,i) - stockwatson(1:end-1,i);
           stockwatson(1,i) = NaN;
           
           disp(strcat(num2str(i), ' firstdiff '));      
       end
   end
end