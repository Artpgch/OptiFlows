nM = readmatrix('NetMatrix.xlsx');
fM = readmatrix('FlowMatrix.xlsx');
[N, nFlows, OFQM, qLM] = OptiFlows(nM, fM);
nl = flowDir(fM);
clc
if N > 0
    fprintf('Optimization completed.\n');
    for n = 1 : nFlows
        writematrix(OFQM(:, :, n), 'OptiFlows.xlsx', 'Sheet', strcat(num2str(nl(1, n)), '⟶', num2str(nl(2, n))));
    end
else 
    fprintf('Optimization impossible. Need to increase the bandwidth of the channels. \n');
end
    writematrix(qLM, 'LoadFactors.xlsx');

function nl = flowDir(nM)

qN = length(nM);
nfl = nnz(nM);
nl = zeros(1,2*nfl);
c = 1;
for a = 1 : qN
    for b = 1 : qN
        if nM(a, b) > 0
            nl(2*c - 1) = a;
            nl(2*c) = b;
            c = c + 1;
        end
    end
end

nl = reshape(nl, [2, nfl]);

end
