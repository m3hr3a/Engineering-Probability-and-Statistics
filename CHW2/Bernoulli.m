function [Output]=Bernoulli(n,p);
if (n<1)  % check n , it should be positive
    disp ('error : n should be a positive integer.')
    Output='Not assigend because of error';
else if (~(0<=p<=1)) % p should be in [0,1] interval
    disp ('error : p should be in [0,1] interval.')
    Output='Not assigend because of error';
    else
    Output=rand(1,n) < p ;   % generate bernoulli disturbution by 
                             %comparing each vector element with p
                             % Output(i)=1 i th elemet of rand generated
                             %vector < p and 0 if >p
    Output=double(Output);   %change Output from logical to double
    end
end