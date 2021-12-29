clearvars; close all; clc;
n = 10;
A = ones(n) + (n - 1) * diag(ones(n, 1));
x0 = ones(n, 1);
f1 = @(x) sum(A(1, :) * reshape(x, numel(x), 1))^2;
f2 = @(x) A(2, :) * reshape(x, numel(x), 1).^2 ;
fx = @(x) [f1(x) - f1(x0) ; f2(x) - f2(x0) ;...
    A(3 : end, :) * (reshape(x, numel(x), 1) - x0)];
fxd = @(x) [2 * sum(A(1, :) * reshape(x, numel(x), 1)) * A(1, :); ...
    2 * A(2, :) * reshape(x, numel(x), 1) * A(2, :); A(3 : end, :)];

options.X0 = rand(n, 1);
options.PlotFcns = 'on';
options.Smethod = 'Two';
options.leqsM = 0;
options.s = 2;

subplot(2,2,1)
[Xstar1, X_iter1, errorFun1] = inleqs(n, fx, fxd, 'Newton', options);
subtitle('by Newton')
subplot(2,2,2)
[Xstar2, X_iter2, errorFun2] = inleqs(n, fx, fxd, 'MNewton', options);
subtitle('by MNewton')
subplot(2,2,3)
[Xstar3, X_iter3, errorFun3] = inleqs(n, fx, fxd, 'Secant', options);
subtitle('by Secant')
subplot(2,2,4)
[Xstar4, X_iter4, errorFun4] = inleqs(n, fx, fxd, 'QNewton', options);
subtitle('by QNewton')


