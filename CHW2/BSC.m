function Output = BSC(Input , e);
n=length(Input);  % length of input vector
probv= Bernoulli(n,e); % probability vector
indexchanged=find(probv==1); % elements change with prob e
Output=Input;   % first copy input into output and then make changes
Output(indexchanged)=1-Input(indexchanged); % 1-0 = 1
                                            % 1-1 = 0
end