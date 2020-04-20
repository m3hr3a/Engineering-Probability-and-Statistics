function [C]=RandomCode(n);
p=0.5;
for i = 1 : 26
C(i,:)=Bernoulli(n,p);
end
end