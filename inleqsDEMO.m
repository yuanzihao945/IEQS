clearvars; close all; clc

nc = [10, 20, 50, 100];
options.Display = 1;
options.PlotFcns = 'on';
options.Smethod = 'Two';
options.leqsM = 0;
options.s = 2;

for i = 1 : 4
    n = nc(i);
    options.X0 = ones(n, 1);
    [nlef, nlefd, Xtrue] = genNLE(n);

    subplot(2,2,i)
    hold on
    [Xstar1, X_iter1, errorFun1] = inleqs(n, nlef, nlefd, 'Newton', options);
    [Xstar2, X_iter2, errorFun2] = inleqs(n, nlef, nlefd, 'MNewton', options);
    [Xstar3, X_iter3, errorFun3] = inleqs(n, nlef, nlefd, 'Secant', options);
    [Xstar4, X_iter4, errorFun4] = inleqs(n, nlef, nlefd, 'QNewton', options);
    legend({'Newton', 'MNewton', 'Secant', 'QNewton'})
end


