function [Output]=Rayleigh(sigma,n);
expdist=Exponential(1,n);  % exp dist
Output=sigma*(2*expdist).^(1/2); %rayleigh dist , proof in report
end