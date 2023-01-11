% Данная функция генерирует случайную матрицу структуры сети, 
% пропускной способности и величины информационных потоков от узла к узлу. 
% Цель - моделирование случайных структур сетей для дальнейших
% математических исследований.
% Исходные данные: количество узлов и максимальная пропускная способность.
% Результат: массив случайной матрицы сети с пропускными способностями, массив случайных потоков потоков. 

function [nMatrix, bandWidthMatrix, flowMatrix] = NetMatrix(qNodes, maxBandWidth)

% Исходные переменные:
% qNodes - количество узлов
% maxBandWidth - максимальная пропускная способность
 
% Выходные переменные:
% nMatrix - матрица сети
% bandWidthMatrix - матрица сети с пропускными способностями
% flowMatrix - матрица потоков

diagZeros = eye(qNodes);

for i = 1 : qNodes
    for j = 1 : qNodes
        if diagZeros(i, j) == 1
                diagZeros(i, j) = 0;  % Создание матрицы с нулевой диагональю
        else                    % для исключения зацикливания потоков на одном узле.
                diagZeros(i, j) = 1;
        end
    end
end

nMatrix = diagZeros .* round( rand(qNodes) ); % Создание случайной матрицы сети

for n = 1 : qNodes
    for m = n : qNodes
        nMatrix(m, n) = nMatrix(n, m); % Отражение матрицы по диагонали для обеспечения двусторонней связи узлов.
    end
end

bandWidthMatrix = maxBandWidth .* rand(qNodes) .* nMatrix; % Создание матрицы случайных пропускных способностей.

flowMatrix = maxBandWidth .* diagZeros .* rand(qNodes); % Создание матрицы случайных потоков такой, чтобы сумма потоков
flowMatrix = flowMatrix ./ ( sum(sum(flowMatrix)) / sum(sum(bandWidthMatrix)) ); % была не выше суммы пропускных способностей.

end