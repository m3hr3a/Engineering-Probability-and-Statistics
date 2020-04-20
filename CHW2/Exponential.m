function [Output]=Exponential(landa,n);
U=rand(1,n); % unif 
Output=-landa*log(1-U); %exp dist , proof in report
end