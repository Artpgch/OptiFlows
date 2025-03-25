nM = readmatrix('NetMatrix.xlsx');
fM = readmatrix('FlowMatrix.xlsx');
lM = readmatrix('LoadFactors.xlsx');
[nfl, nl, fG] = flowGraph(fM);
nG = netGraph(nM);
sG = simpGraph(nM);
lG = netGraph(lM);
fnums = reshape(nl, [2, nfl]);
oM = zeros(length(nM), length(nM), nfl);
for n = 1 : nfl
    oM(:, :, n) = cell2mat(readcell('OptiFlows.xlsx','Sheet', n));
end

figure('Name','Simplified network graph','NumberTitle','off','MenuBar','none')
plot(sG,'NodeColor', [245/256 64/256 33/256], 'EdgeFontSize', 10, 'EdgeColor', [127/256 255/256 212/256], 'LineWidth',1.5,'EdgeAlpha', 1, 'NodeFontSize',12,'NodeLabelColor',[1 1 1])
title('Simplified network graph')
set(gca, 'color', [62/256 95/256 138/256])
exportgraphics(gcf, 'OptiFlows.pdf','Append', true)

figure('Name','Flow directions','NumberTitle','off','MenuBar','none')
plot(fG, 'NodeLabel', nl, 'EdgeLabel', fG.Edges.Weight,'NodeColor', [245/256 64/256 33/256],'ArrowSize', 13, 'EdgeFontSize', 12, 'EdgeColor', [127/256 255/256 212/256], 'LineWidth',1.5,'EdgeAlpha', 1, 'NodeFontSize',12,'NodeLabelColor',[1 1 1], 'EdgeFontSize', 10, 'EdgeLabelColor',[0 1 0],'EdgeAlpha', 1, 'EdgeFontSize', 11, 'EdgeLabelColor', [1 1 0]);
title('Flow directions')
set(gca, 'color', [62/256 95/256 138/256])
exportgraphics(gcf, 'OptiFlows.pdf','Append', true)

figure('Name','Network graph','NumberTitle','off','MenuBar','none')
plot(nG, 'Layout', 'force', 'EdgeLabel', nG.Edges.Weight, 'NodeColor', [245/256 64/256 33/256], 'EdgeFontSize', 10, 'EdgeColor', [127/256 255/256 212/256], 'ArrowSize', 10,'NodeFontSize',12,'NodeLabelColor',[1 1 1],'LineWidth',1.5,'MarkerSize', 4,'EdgeLabelColor',[0 1 0],'EdgeAlpha',1);
title('Network graph')
set(gca, 'color', [62/256 95/256 138/256])
exportgraphics(gcf, 'OptiFlows.pdf','Append', true)

figure('Name','Load factors of communication lines','NumberTitle','off','MenuBar','none')
plot(lG, 'Layout', 'force', 'EdgeLabel', lG.Edges.Weight, 'NodeColor', [245/256 64/256 33/256], 'EdgeFontSize', 10, 'EdgeColor', [127/256 255/256 212/256], 'ArrowSize', 10,'NodeFontSize',12,'NodeLabelColor',[1 1 1],'LineWidth',1.5,'MarkerSize', 4,'EdgeLabelColor',[0 1 0],'EdgeAlpha',1);
title('Load factors of communication lines')
set(gca, 'color', [62/256 95/256 138/256])
exportgraphics(gcf, 'OptiFlows.pdf','Append', true)

for m = nfl:-1:1
    optiname = strcat('Optimal flows from node', 32, num2str(fnums(1, m)), 32, 'to node', 32, num2str(fnums(2, m)));
    oG = netGraph(oM(:, :, m));
    figure('Name', optiname,'NumberTitle','off','MenuBar','none')
    plot(oG, 'Layout', 'force', 'EdgeLabel', oG.Edges.Weight, 'NodeColor', [245/256 64/256 33/256], 'EdgeFontSize', 10, 'EdgeColor', [127/256 255/256 212/256], 'ArrowSize', 10,'NodeFontSize',12,'NodeLabelColor',[1 1 1],'LineWidth',1.5,'MarkerSize', 4,'EdgeLabelColor',[1 1 0],'EdgeAlpha',1);
    title(optiname)
    set(gca, 'color', [62/256 95/256 138/256])
    exportgraphics(gcf, 'OptiFlows.pdf','Append', true)
end

function nG = netGraph(nM)
qN = length(nM);
t = zeros(1, nnz(nM));
s = t;
bW = s;
u = 1;
for p = 1 : qN
    for q = 1 : qN
        if nM(p, q) > 0
                s(u) = p;
                t(u) = q;       
                bW(u) = nM(p, q);   
                u = u + 1;
        end
     end
end

nG = digraph(s, t, bW);

end

function sG = simpGraph(nM)
qN = length(nM);
rM = nM ./ nM;

for a = 1 : qN
    for b = 1 : qN
        if isnan(rM(a, b))
            rM(a, b) = 0;
        end
    end
end

ruM = triu(rM);
trM = rM- ruM;
rrM = zeros(qN);
for p= 1:qN
    for q = 1:qN
        rrM(q, p) = trM(p,q);
    end
end
grM = rrM + ruM;

t = zeros(1, nnz(grM));
s = t;
bW = s;
u = 1;
for p = 1 : qN
    for q = 1 : qN
        if grM(p, q) > 0
                s(u) = p;
                t(u) = q;        
                bW(u) = grM(p, q);  
                u = u + 1;
        end
     end
end
sG = graph(s,t);

end

function [nfl, nl, fG] = flowGraph(nM)

qN = length(nM);
nfl = nnz(nM);
nl = zeros(1,2*nfl);
flG(1,:) = 1:2:2*nfl-1;
flG(2,:) = 2:2:2*nfl;
flI = zeros(1,nfl);
c = 1;
for a = 1 : qN
    for b = 1 : qN
        if nM(a, b) > 0
            nl(2*c - 1) = a;
            nl(2*c) = b;
            flI(c) = nM(a, b);
            c = c + 1;
        end
    end
end
fG = digraph(flG(1,:), flG(2,:), flI);

end
