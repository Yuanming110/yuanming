function [x1 x2]=boundary(x0,B1,B2,range)
i=fix(x0);
            if i-range<=B1 
               x1=B1;x2=i+range;
            elseif i+1>=B2
               x1=i-range;x2=B2;
            else
               x1=i-range;x2=i+range;
            end
            
end