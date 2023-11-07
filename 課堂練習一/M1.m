function y=M1(x)
y1=1/4000*sum(x.^2);
y2=1;
for i=1:size(x,2)
    y2=y2*cos(x(i)/sqrt(i));
end
y=y1-y2+1;

    
    