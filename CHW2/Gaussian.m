function [Output]=Gaussian(sigma,n);
raydist=Rayleigh(sigma,n); %rayleigh dist
unidist=2*pi*rand(1,n);    %uniform dist
Output=raydist.*cos(unidist); % gaussian dist , proof in reprt
end