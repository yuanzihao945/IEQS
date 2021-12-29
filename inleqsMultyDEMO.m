clearvars; close all; clc

N = 200;
nc = [10, 20, 50, 100];
Nall = zeros(N, 4, 4);
errorX = zeros(N, 4, 4);
options.Display = 0;
options.PlotFcns = 'off';
options.Smethod = 'Two';
options.leqsM = 0;
options.s = 2;
figure(1);
for i = 1 : 4
    n = nc(i);
    options.X0 = ones(n, 1);

    for iter = 1 : N
        [nlef, nlefd, Xtrue] = genNLE(n);
        [Xstar1, X_iter1, errorFun1] = inleqs(n, nlef, nlefd, 'Newton', options);
        [Xstar2, X_iter2, errorFun2] = inleqs(n, nlef, nlefd, 'MNewton', options);
        [Xstar3, X_iter3, errorFun3] = inleqs(n, nlef, nlefd, 'Newton', options);
        [Xstar4, X_iter4, errorFun4] = inleqs(n, nlef, nlefd, 'QNewton', options);
        Nall(iter, :, i) = [length(errorFun1), length(errorFun2), ...
            length(errorFun3), length(errorFun4)];
        errorX(iter, :, i) = [norm(Xstar1 - Xtrue, 2), norm(Xstar2 - Xtrue, 2), ...
            norm(Xstar3 - Xtrue, 2), norm(Xstar4 - Xtrue, 2)];
    end

    figure(1)
    subplot(2,2,i)
    fig1 = boxplot(Nall(:, :, i), {'Newton', 'MNewton', 'Secant', 'QNewton'}, 'Widths',0.3);
    set(fig1, 'Linewidth', 2);
    subtitle(['n = ' num2str(n)])
    ylabel('Stop Number of Iteration')
end
