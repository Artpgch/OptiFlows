nM = readmatrix('NetMatrix.xlsx');
fM = readmatrix('FlowMatrix.xlsx');
qn = length(nM);
global Mu nFlows
nEdges = nnz(nM);
target = zeros( 1, nEdges );
source = target;
Mu = source;
u = 1;
for p = 1 : qn
    for q = 1 : qn
        if nM(p, q) > 0
                source(u) = p;
                target(u) = q;
                Mu(u) = nM(p, q);
                u = u + 1;
        end
     end
end
iFlows = zeros(qn, nEdges);
for l = 1 : nEdges
    iFlows(source(l), l) = -1;
    iFlows(target(l), l) =  1;
end
nFlows = nnz(fM);
ftarget = zeros(1, nFlows );
fsource = ftarget;
flows = fsource;
u = 1;
for p = 1 : qn
        for q = 1 : qn
            if fM(p, q) > 0
                    fsource(u) = p;
                    ftarget(u) = q;
                    flows(u) = fM(p, q);
                    u = u + 1;
            end
        end
end
fST = zeros(qn, nFlows);
for l=1 : nFlows
    fST(fsource(l), l) = -flows(l);
    fST(ftarget(l), l) =  flows(l);
end
flowSourceTarget = zeros(qn + 2, nFlows);
flowSourceTarget(1 : qn, :) = fST;
stLimit = zeros(2, nEdges, nFlows);
for v = 1 : nFlows
    for w = 1 : nEdges
        if source(w) == ftarget(v)
            stLimit(1, w, v) = 1;
        elseif target(w) == fsource(v)
            stLimit(2, w, v) = 1;
        end
    end
end
QFlows = zeros(qn + 2, nEdges, nFlows);
for o = 1 : nFlows
QFlows(1 : qn, :, o)=iFlows;
QFlows(qn + 1 : qn + 2, :, o) = stLimit(:, :, o);
end
for k = 1 : nFlows
    for i = 1 : qn
        switch sign(flowSourceTarget(i, k))
            case -1
                for j = 1 : nEdges
                    if QFlows(i, j, k) == 1
                        QFlows(i, j, k) = 0;
                    end
                end
            case 1                                  
                for j = 1 : nEdges
                    if QFlows(i, j, k) == -1
                        QFlows(i, j, k) = 0;
                    end
                end
        end
    end
end
Aeq = zeros((qn + 2) * nFlows, nEdges * nFlows);
Beq = zeros(( qn + 2) * nFlows, 1);
for j = 0 : nFlows - 1
    Aeq((qn + 2) * j + 1 : (qn + 2) * (j + 1), nEdges * j + 1 : nEdges * (j + 1)) = QFlows(:, :, j + 1);
    Beq((qn + 2) * j + 1 : (qn + 2) * (j + 1), 1) = flowSourceTarget(:, j + 1);
end
Lambda0 = zeros(1, nFlows * nEdges);
A = zeros(nEdges, nEdges*nFlows );
for z=1:nEdges
    for y = 0 : nFlows - 1
        A(z, z + y * nEdges) = 1; % Формирование левой части ограничительного неравенства для не превышения пропускной способности.
    end
end
B(:,1)=Mu; % Формирование столбца-массива пропускных способностей
options = optimset('Algorithm', 'sqp', 'TolCon', 2^-32, 'MaxFunEvals', 2^18); %interior-point
[Lambda, N] = fmincon(@TargetFunction, Lambda0, A, B, Aeq, Beq, Lambda0, [], [], options);
oFlows=zeros(nFlows,nEdges);
for n = 0 : nFlows - 1
    oFlows(n + 1, :) = Lambda(n * nEdges + 1 : nEdges * (n + 1));
end
for i = 1 : nFlows
    for j = 1 : nEdges
        if oFlows(i, j) < 10^-10
                        oFlows(i, j) = 0;
        end
    end
end
OFQM = zeros(qn, qn, nFlows);
for k = 1 : nFlows
    for x = 1 : nEdges
        OFQM(source(x), target(x), k) = oFlows(k, x); % Формирование трёхмерного массива оптимальных потоков
    end
end
qL = sum(oFlows)./Mu;  % Коэффициенты загрузки каналов связи.
qLM = zeros(qn);
for v = 1 : nEdges
     qLM(source(v), target(v)) = qL(v); % Формирование матрицы загрузки каналов связи.
end
if N > 0
    disp('Optimization completed.');
    for i=1:nFlows
    writematrix(OFQM(:, :, i), 'OptiFlows.xlsx', 'Sheet', i);
    end
else disp("Optimization impossible. Need to increase the bandwidth of the channels.");
end 
writematrix(qLM, 'LoadCoeffs.xlsx');
