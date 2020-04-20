%%  
%% Probability and Statistics - Dr.Mirmohseni Spring 2018
%% Course Project /Computer Assignmet phase 1 
%% Name                  Student ID
%% Mehrsa Pourya         95101247
%% 1.1 Bernoulli Disturbution
% functions needed : Bernoulli
clear all
close all
clc 
% test Bernoulli function 
% you can change n and p and check output
p = 0.3    % value between 0 and 1 
n= 10    % positive integer
ourbernoulliseries=Bernoulli(n,p)
%% 2 Exponential Disturbution
% functions needed : Exponential
clear all
close all
clc 
landa=1;    % you can choose landa here 1 : default
n1=10^3;    % first given n in hw
n2=10^5;    % second given n 
n3=10^7;    % third given n 
x1=Exponential(landa,n1);      % disturbution for n1
x2=Exponential(landa,n2);      % disturbution for n1
x3=Exponential(landa,n3);      % disturbution for n1
m = 100;    % number of numbins for histogram
t = linspace(min([x1 x2 x3]),max([x1 x2 x3]),m);   % proper t vector 
                                                   % between min and max
                                                   % outputs
dt=t(2)-t(1);  
f_x1=hist(x1,t)/(n1*dt);  % pdf exp for n1
f_x2=hist(x2,t)/(n2*dt);  % pdf exp for n2
f_x3=hist(x3,t)/(n3*dt);  % pdf exp for n3
pdfteory=(1/landa)*exp(-t/landa);   %theoritical pdf with dist formula
figure % plot wantededs
plot(t,f_x1,'b','linewidth',1)
hold on 
plot(t,f_x2,'r','linewidth',1)
hold on 
plot(t,f_x3,'g','linewidth',1)
hold on 
plot(t,pdfteory,'black','linewidth',1)
grid on
legend(['n=',num2str(n1)],['n=',num2str(n2)],['n=',num2str(n3)], ...
    'Theoretical pdf')
title('Theoretical and Generated Exponential RV \lambda=1')
xlabel('x')
ylabel('f_{X}(x)~ (1/\lambda)e^{-x/\lambda}')
%% 3 Rayleigh Disturbution
% functions needed : Exponential , Rayleigh
clear all
close all
clc 
sigma=1;    % you can choose sigma here 1 : default
% define given n values
n1=10^3;
n2=10^5;
n3=10^7;
x1=Rayleigh(sigma,n1);   % Rayleigh dist generator for n1
x2=Rayleigh(sigma,n2);   % Rayleigh dist generator for n1
x3=Rayleigh(sigma,n3);   % Rayleigh dist generator for n1
m = 100;
t = linspace(min([x1 x2 x3]),max([x1 x2 x3]),m);  % proper t vector 
                                                   % between min and max
                                                   % outputs
dt=t(2)-t(1);
f_x1=hist(x1,t)/(n1*dt);   % pdf rayleigh for n1
f_x2=hist(x2,t)/(n2*dt);   % pdf rayleigh for n2
f_x3=hist(x3,t)/(n3*dt);   % pdf rayleigh for n3
pdfteory=(t./(sigma^2)).*exp(-(t.^2)./(2*sigma^2));  %theoritical pdf with
                                                     %dist formula
figure % plot wanteds
plot(t,f_x1,'b','linewidth',1)
hold on 
plot(t,f_x2,'r','linewidth',1)
hold on 
plot(t,f_x3,'g','linewidth',1)
hold on 
plot(t,pdfteory,'black','linewidth',1)
grid on
legend(['n=',num2str(n1)],['n=',num2str(n2)],['n=',num2str(n3)], ...
    'Theoretical pdf')
title('Theoretical and Generated Rayleigh RV \sigma=1')
xlabel('y')
ylabel('f_{Y}(y)~ (y/\sigma^{2})e^{-y/2\sigma^{2}}')
%% 4 Rayleigh Disturbution
% functions needed : Exponential , Rayleigh , Gaussian
clear all
close all
clc 
sigma=1;    % you can choose sigma here 1 : default
n1=10^3;
n2=10^5;
n3=10^7;
x1=Gaussian(sigma,n1);  % gaussian generator n1
x2=Gaussian(sigma,n2);  % gaussian generator n2
x3=Gaussian(sigma,n3);  % gaussian generator n3
m = 100;
t = linspace(min([x1 x2 x3]),max([x1 x2 x3]),m);
dt=t(2)-t(1);
f_x1=hist(x1,t)/(n1*dt); % pdf gaussain for n1
f_x2=hist(x2,t)/(n2*dt); % pdf gaussian for n2
f_x3=hist(x3,t)/(n3*dt); % pdf gaussian for n3
pdfteory=(1/sqrt(2*pi*sigma^2)).*exp(-(t.^2)./(2*sigma^2));%theoritical pdf 
figure %plot wanteds
plot(t,f_x1,'b','linewidth',1)
hold on 
plot(t,f_x2,'r','linewidth',1)
hold on 
plot(t,f_x3,'g','linewidth',1)
hold on 
plot(t,pdfteory,'black','linewidth',1)
grid on
legend(['n=',num2str(n1)],['n=',num2str(n2)],['n=',num2str(n3)], ...
    'Theoretical pdf')
title('Theoretical and Generated Gaussian RV \sigma=1')
xlabel('y')
ylabel('f_{Y}(y)~ (1/sqrt(2\pi\sigma^{2})e^{-y/2\sigma^{2}}')
