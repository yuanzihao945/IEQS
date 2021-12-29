clearvars; close all; clc

n = 10;
options.X0 = ones(n, 1);
options.Display = 1;
options.PlotFcns = 'on';
options.Smethod = 'Two';
options.leqsM = 0;
options.s = 2;

[nlef, nlefd, Xtrue] = genNLE(n);

subplot(2,2,1)
[Xstar1, X_iter1, errorFun1] = inleqs(n, nlef, nlefd, 'Newton', options);
subtitle('by Newton')
subplot(2,2,2)
[Xstar2, X_iter2, errorFun2] = inleqs(n, nlef, nlefd, 'MNewton', options);
subtitle('by MNewton')
subplot(2,2,3)
[Xstar3, X_iter3, errorFun3] = inleqs(n, nlef, nlefd, 'Secant', options);
subtitle('by Secant')
subplot(2,2,4)
[Xstar4, X_iter4, errorFun4] = inleqs(n, nlef, nlefd, 'QNewton', options);
subtitle('by QNewton')