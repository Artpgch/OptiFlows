function MiddleQueue = TargetFunction(LambdaFlows)
global Mu nFlows
fS=zeros(1, length(Mu));
MiddleQueue=0;
for i = 1 : length(Mu)
    for j = 0 : nFlows - 1
        fS(i) = fS(i) + LambdaFlows(length(Mu) * j + i); 
    end
end
for k = 1 : length(Mu)
    MiddleQueue = MiddleQueue + (fS(k) / Mu(k)) / (1-fS(k) / Mu(k));
end
end