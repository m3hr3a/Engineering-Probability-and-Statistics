%% Probability and Statistic 
%% Dr Mirmohseni
%% Spring 2018
%% Mehrsa Pourya 95101247
%% ------------------------------------------------------
%% please add to path all included functions before running the code
%% BSC 
%% part alef 
% no code needed
%% part b
% BSC funtion , you need Bernoulli function too for this ;
% Both functions are attached
%% part p simple coding 
%initialize
clc 
clear
close
for run = 1 : 10   % repeat 10 times
message='Information is the resolution of uncertainty. C.Shannon'; % given message
% convert from string to binary series
messagevector=reshape(uint8(dec2bin(message)'-'0'),55*7,1); 
e=0.01; % channel e , given
outputm=BSC(messagevector,e); % send message using channel
% convert output from binary series to string
out=cellstr(char(bin2dec(char(reshape(outputm,7,55)'+'0')))'); 
disp(['Run ',num2str(run)])
oo=out{1,1};   % output message , string format
sqerror(run)=mean((outputm-messagevector).^2); % calculate square error
disp(['square error : ',num2str(sqerror(run))])
disp(oo)
end
disp(['mean square error : ',num2str(mean(sqerror))])
%% BSC t repetition coding
clc
clear
close
th=0:0.1:1;
% given message
message='Information is the resolution of uncertainty. C.Shannon';
% convert from string to binary series
messagevector=reshape(uint8(dec2bin(message)'-'0'),55*7,1); 
for i = 1 : 10
repeatedmessage(i,:)=messagevector;         % repeat each bit 10 times
end
reshapedmessage=reshape(repeatedmessage,1,[]);   % reshape repeated matrix to vector
e=0.01; % channel e , given
for i = 1 : 11 % tresholds
for send = 1 : 100  %  10 times we repeat
output=BSC(reshapedmessage,e); % send message using channel
outputm=mean(reshape(output,10,[]),1);
meanout=mean(outputm,1);
index1=find(meanout>th(i));
index0=find(meanout<=th(i));
meanmessage(index1)=1;
meanmessage(index0)=0;
error(send)=mean((meanmessage'-double(messagevector)).^2);
end
sqerrorr(i)=mean(error); % calculate square error
end
plot(th(2:10),sqerrorr(2:10),'b','lineWidth',1.5)
hold on 
scatter(th(2:10),sqerrorr(2:10),'fill')
title('Square Error based on chosen thresohld')
xlabel('threshold')
ylabel('Square error')
grid on
%% BSC t disturbution of channel output for each bit in repetition coding
clear
clc
close  
th = 0 : 0.1 : 1 ; 
e= 0.01;
% 0 is send 
for j = 0 : 10
Pzi0(j+1)=factorial(10)/(factorial(10-j)*factorial(j))*e.^(j)*(1-e).^(10-j);
Pzi1(j+1)=factorial(10)/(factorial(10-j)*factorial(j))*e^(10-j)*(1-e)^(j);
end
plot(th,Pzi0,'b','lineWidth',1.5)
hold on 
plot(th,Pzi1,'r','lineWidth',1.5)
grid on 
title('Probability of each value of Z')
xlabel('Z')
ylabel('Probability')
legend('0 is sent' , '1 is sent')
%% BSC s Random coding matrix generator
clc
n=10;
Rmat=RandomCode(n);
%% BSC j distrubution of hamming of two raws
clear
clc
t=0;
n=100; % chose n
for run = 1 : 100 % repeat 100 time
    Rmat=RandomCode(n); % generate random coding mat
    for i = 1 : 26
        for j = i+1 : 26
        t= t+1;
        hdij(t)=sum(abs(Rmat(i,:)-Rmat(j,:))); % calulating all possible hammings
        end                                    
    end
end
bin=0 : 1 : n+1;
h=histogram(hdij,bin);
values=h.Values;
close
stem(bin(1:n+1),values/t)
grid on 
for k =  0 : n  % calculate Pn(k)
teovalues(k+1)=factorial(n)/(factorial(bin(k+1))*factorial(n-bin(k+1)))*(0.5)^n;
end
hold on 
plot(bin(1:n+1),teovalues)
ylabel('Probability')
legend('Simulation','teoretical')
xlabel('hamming distanse')
title('Hamming distance disturbution')
%% BSC ch min hamming distance disturbution
clear
n=100;
t= 0;
for run = 1 : 10000 % repeat 10000 times
Rmat=RandomCode(n);
for i = 1 : 26
    for j = 1 + i :26
        t= t + 1;
       hdij(run,t)=sum(abs(Rmat(i,:)-Rmat(j,:)));
    end
end
t=0;
end
hdmin=min(hdij');  % chose min hammings
bin=0 : 1 : n+1;
h2=histogram(hdmin,bin);
values2=h2.Values;
close
stem(bin(1:n+1),values2/10000)
grid on 
hold on
% teoretical disturbution
bin=0 : 1 : n+1;
for k =  0 : n % caculate Pn(k)
teovalues(k+1)=factorial(n)/(factorial(bin(k+1))*factorial(n-bin(k+1)))*(0.5)^n;
end
for i = 1 : n+1  % bionomial cdf
    mycdf(i)=sum(teovalues(1:i));
end 
mm=[0,mycdf];
for k = 0 : n
    p(k+1)=(1-mm(k+1)).^325-(1-mm(k+2)).^325;
end
plot(bin(1:n+1),p)
legend('Simulation','theoritical')
xlabel('min hamming')
ylabel('Probability')
title(['Min Hamming Disturbution for n = ',num2str(n)])
%% BSC Error Calculation
clear 
clc
close 
ev=[0.1 0.01];
for g = 1 : 2
  e=ev(g)  ;
for n = 1 : 100
bin=0 : 1 : n+1;
for k =  0 : n % caculate Pn(k)
teovalues(k+1)=factorial(n)/(factorial(bin(k+1))*factorial(n-bin(k+1)))*(0.5)^n;
end
for k =  0 : n % caculate Pn(k)
pn(k+1)=factorial(n)/(factorial(bin(k+1))*factorial(n-bin(k+1)));
end
for i = 1 : n+1  % bionomial cdf
    mycdf(i)=sum(teovalues(1:i));
end 
mm=[0,mycdf];
for k = 0 : n
    p(k+1)=(1-mm(k+1)).^325-(1-mm(k+2)).^325;
end
for m = 0 : n
    t=ceil(m/2)-1;
    if t==-1
    t=0;
    end
    x=0;
    for i = t+1 : n;
        x=x+1;
    p2(x)=e^i*(1-e)^(n-i)*pn(i);
    end
    p3(m+1)=sum(p2);
    %clear p2 t
end
for m = 0 : n
    p4(m+1)=p3(m+1)*p(m+1);
end
perror(n)=sum(p4);
clear p4 p3 p mycdf teovalues bin
end
plot(1:n,perror)
hold on
end
legend('e=0.1','e=0.01')
xlabel('n')
ylabel('Error Probability')
title('teoretical error of hamming coding for one character')
%% BSC h Random coding and decoding 
clc 
clear
n=100;
e=0.01;
disp('for e = 0.01 and n = 100')
alphabet='abcdefghijklmnopqrstuvwxyz .ABCDEFGHIJKLMNOPQRSRUVWXYZ';
alchar=double(dec2bin(alphabet))-'0';
mycoding=RandomCode54(n);
message='Information is the resolution of uncertainty. C.Shannon';
messagemat=double(dec2bin(message))-'0';
for i = 1 : 55
    for  j = 1 : 54
        if isequal(messagemat(i,:),alchar(j,:))
            index(i)=j;
            break 
        end
    end
end
for i = 1 : 55
codedmessage(i,:)=mycoding(index(i),:);
end
codedmessagevector=reshape(codedmessage',1,[]);
outputmessage=BSC(codedmessagevector,e);
squareerrorbeforenearesthamming=sum(sum((outputmessage-codedmessagevector).^2))
outputmat=reshape(outputmessage,[],55)';
for i = 1 :55
    for j = 1 : 54
    hamming(j)=sum(abs(outputmat(i,:)-mycoding(j,:)));
    end
    chosenindex=min(find(hamming==min(hamming)));
    decodedmat(i,:)=alchar(chosenindex,:);
end
t=0;
for i = 1 : 54
    for j = i+1 : 54
        t=t+1;
        hamming(t)=sum(abs(mycoding(i,:)-mycoding(j,:)));
    end
end
b=cellstr(char(bin2dec(char(decodedmat+'0')))');
detectedmessage=b{1,1}
squareerrorafternearesthamming=sum(sum((decodedmat-messagemat).^2))
minhamming=min(hamming)
%% BSC kh Random coding and decoding 
clc 
clear
n=200;
e=0.1;
disp('for e = 0.1 and n = 200')
alphabet='abcdefghijklmnopqrstuvwxyz .ABCDEFGHIJKLMNOPQRSRUVWXYZ';
alchar=double(dec2bin(alphabet))-'0';
mycoding=RandomCode54(n);
message='Information is the resolution of uncertainty. C.Shannon';
messagemat=double(dec2bin(message))-'0';
for i = 1 : 55
    for  j = 1 : 54
        if isequal(messagemat(i,:),alchar(j,:))
            index(i)=j;
            break 
        end
    end
end
for i = 1 : 55
codedmessage(i,:)=mycoding(index(i),:);
end
codedmessagevector=reshape(codedmessage',1,[]);
outputmessage=BSC(codedmessagevector,e);
squareerrorbeforenearesthamming=sum(sum((outputmessage-codedmessagevector).^2))
outputmat=reshape(outputmessage,[],55)';
for i = 1 :55
    for j = 1 : 54
    hamming(j)=sum(abs(outputmat(i,:)-mycoding(j,:)));
    end
    chosenindex=min(find(hamming==min(hamming)));
    decodedmat(i,:)=alchar(chosenindex,:);
end
t=0;
for i = 1 : 54
    for j = i+1 : 54
        t=t+1;
        hamming(t)=sum(abs(mycoding(i,:)-mycoding(j,:)));
    end
end
b=cellstr(char(bin2dec(char(decodedmat+'0')))');
detectedmessage=b{1,1}
squareerrorafternearesthamming=sum(sum((decodedmat-messagemat).^2))
minhamming=min(hamming)
%% 2 PAM Part 
close 
clc
clear 
message='Information is the resolution of uncertainty. C.Shannon'; % given message
% convert from string to binary series
messagevector=reshape(uint8(dec2bin(message)'-'0'),55*7,1);
clear 
for N = 1 : 20

Eg=100;
n=10^5;
mynoise=Gaussian(sqrt(N),n);
mymessage=Bernoulli(n,0.5);
message1=find(mymessage==1);
sendmessage(message1)=sqrt(Eg);
message0=find(mymessage==0);
sendmessage(message0)=-sqrt(Eg);
sendmessage=sendmessage+mynoise;
r1=find(sendmessage>=0);
r0=find(sendmessage<0);
detmessage(r1)=1;
detmessage(r0)=0;
error1=sum(abs(detmessage-mymessage));
error(N)=error1/10^5;
teo(N)=0.5-0.5*normcdf(sqrt(Eg)/sqrt(N),0,1)+0.5*normcdf(-sqrt(Eg)/sqrt(N),0,1);
end
plot(1:N,error,'b','LineWidth',1.5)
hold on
plot(1:N,teo,'r','LineWidth',1.5)
grid on 
title('error based on noise variance E_{g}=100 p=0.5')
legend('simulation','teoritical')
xlabel('N')
ylabel('error')
%% 2 PAM Part s | repetation coding
%initialize
clear 
clc
close
n=10^5;       % given n = 10^5
p=0.5;        % each bit has a Bernouli(p = 0.5) disturbution 
Eg=4;         % given Eg = 4
N=1;          % noise variance
message=Bernoulli(n,p);               % generate Bernoulli message
for i = 1 : 10
repeatedmessage(i,:)=message;         % repeat each bit 10 times
end
reshapedmessage=reshape(repeatedmessage,1,[]);   % reshape repeated matrix to vector   
message1=find(reshapedmessage==1);      % message 1's 
sendmessage(message1)=sqrt(Eg);         % replace one with Eg
message0=find(reshapedmessage==0);      % message 0's
sendmessage(message0)=-sqrt(Eg);        % replace 0 with -Eg
mynoise=Gaussian(sqrt(N),10*n);         % generate white gaussain noise 
sendmessage=sendmessage+mynoise;        % add noise
% repeat above procedure but this time we do not reapet each bit 10 times
message1nor=find(message==1);
sendmessagenor(message1nor)=sqrt(Eg);
message0nor=find(message==0);
sendmessagenor(message0nor)=-sqrt(Eg);
mynoise2=Gaussian(sqrt(N),n);
sendmessagenor=sendmessagenor+mynoise2;
% restruct output messages due to our decision rulre
reshapedoutput=mean(reshape(sendmessage,10,[]),1);  % average and reshape
r1=find(reshapedoutput>0); 
r0=find(reshapedoutput<=0);
detmessage(r1)=1;
detmessage(r0)=0;
% output for nonrepeated 
r1nor=find(sendmessagenor>=0);
r0nor=find(sendmessagenor<0);
detmessagenor(r1nor)=1;
detmessagenor(r0nor)=0; 
% calculate errors
error=sum((detmessage-message).^2);      % square error for reapeted coding
errornor=sum((detmessagenor-message).^2); % square error for simple coding
disp('square error for reapeted coding')
disp(error)
disp('square error for simple coding')
disp(errornor)








