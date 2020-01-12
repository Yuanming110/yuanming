function F=fon(di,r1,r2)

for i=1:length(di)
    if di(i)>=0 & di(i) <= r2
        F(i)=1;
    end
    if di(i)>r2 & di(i) <= r1+r2
        F(i)=0.1^((di(i)-r2)/r1);
    end
    if di(i) > r1+r2
        F(i)=0;
    end
end
end