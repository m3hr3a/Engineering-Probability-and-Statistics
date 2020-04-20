function [C]=RandomCode54(n);
p=0.5;
for i = 1 : 54
C(i,:)=Bernoulli(n,p);
end
end