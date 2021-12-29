function [Xstar, X_iter, errorFun] = inleqs(n, nlef, nlefd, method, options)
% ILEQS - Iterative optimization solution of Nonlinear equations by Newton
% methods or some variant method
%
% Input:
%   n        - Dimension of NLEQS
%   nlef     - n equations of nonlinear funcions
%   nlefd    - diff of nonlinear funcions
%   method   - Method of iteration: 'Newton', 'QNewton', 'MNewton', 'Secant'
%   options  - Options of algorithm(Struct data), include: 
%            - options.p: use p-norm
%            - options.X0: Intial value of algorithm
%            - options.TolFun: tolerance of equation
%            - options.TolX: tolerance of solution
%            - options.Maxiter: Maximum number of iterations
%            - options.Display: If display the final iterations number
%            - options.PlotFcns: Draw the error value of each iteration 
%                                during the execution of the algorithm
%            - options.Smethod: method for Secant
%
% Output:
%   Xstar    - Iterative solutions of equations
%   X_iter   - Iterative solutions in every iteration
%   errorFun - Error of linear equations in each iteration
%
% Usage:
%   [] = ILEQS(nlef, nlefd) uses the default settings soluting nonlinear equations
%   [] = ILEQS(nlef, nlefd, method) uses the input method to iterate nonlinear equatioins
%   [] = ILEQS(nlef, nlefd, [], options) uses the options' settings soluting nonlinear equatioins
%   Xstar = ILEQS( ... ) returns the iterative solution of nonlinear equations
%   [~, X_iter] = ILEQS( ... ) returns the solutions in each iterations
%   [~, ~, error] = ILEQS( ... ) returns the error of nonlinear equations
%       in each iterations

% Author  : ZH.Yuan
% Update  : 2021/12/26 (First Version: 2021/12/26)
% Email   : zihaoyuan@whut.edu.cn (If any suggestions or questions)

% Set default value of value
if ~exist('method', 'var') || isempty(method)
    method = 'Newton';
end

% Set default value of options
if ~exist('options', 'var') || isempty(options)
    options.p = 2;
    options.X0 = zeros(n, 1);
    options.TolFun = 1e-04;
    options.TolX = 1e-04;
    options.Maxiter = max([1000 n]);
    options.Display = 0;
    options.omega = [];
    options.PlotFcns = 'off';
    options.leqsM = 0;
    if strcmp(method, 'Secant')
        options.Smethod = 'N';
    elseif strcmp(method, 'MNewton')
        options.s = 3;
    end
end

% Set default value of options.p
if ~isfield(options, 'p')
    options.p = 2;
elseif isempty(options.p)
    options.p = 2;
end

% Set default value of options.X0
if ~isfield(options, 'X0')
    options.X0 = zeros(n, 1);
elseif isempty(options.X0)
    options.X0 = zeros(n, 1);
end

% Set default value of options.TolFun
if ~isfield(options, 'TolFun')
    options.TolFun = 1e-04;
elseif isempty(options.TolFun)
    options.TolFun = 1e-04;
end

% Set default value of options.TolX
if ~isfield(options, 'TolX')
    options.TolX = 1e-04;
elseif isempty(options.TolX)
    options.TolX = 1e-04;
end

% Set default value of options.Maxiter
if ~isfield(options, 'Maxiter')
    options.Maxiter = max([1000 n]);
elseif isempty(options.Maxiter)
    options.Maxiter = max([1000 n]);
end

% Set default value of options.Display
if ~isfield(options, 'Display')
    options.Display = 0;
elseif isempty(options.Display)
    options.Display = 0;
end

% Set default value of options.PlotFcns
if ~isfield(options, 'PlotFcns')
    options.PlotFcns = 'off';
elseif isempty(options.PlotFcns)
    options.PlotFcns = 'off';
end

% Set default value of options.leqsM
if ~isfield(options, 'leqsM')
    options.leqsM = 0;
elseif isempty(options.leqsM)
    options.leqsM = 0;
end

% Set default value of options.Smethod for 'Secant' method
if strcmp(method, 'Secant')
    if ~isfield(options, 'Smethod')
        options.Smethod = 'N';
    elseif isempty(options.Smethod)
        options.Smethod = 'N';
    end
end

% Set default value of options.s for 'MNewton' method
if strcmp(method, 'MNewton')
    if ~isfield(options, 's')
        options.s = 3;
    elseif isempty(options.s)
        options.s = 3;
    end
end

iter = 0;
iter_m = 0;
X_old = reshape(options.X0, n, 1);
nlefd_m = nlefd(X_old);

while iter < options.Maxiter
    iter = iter + 1;

    switch method
        case 'Newton'
            if options.leqsM == 1
                DX = reshape(ileqs([nlefd(X_old), -nlef(X_old)]), n, 1);
            else
                DX = - (nlefd(X_old))^(-1) * nlef(X_old);
            end

        case 'MNewton'
            iter_m_new = floor(iter / options.s);
            if iter_m_new > iter_m
                nlefd_m = nlefd(X_old);
                iter_m = iter_m_new;
            end
            if options.leqsM == 1
                DX = reshape(ileqs([nlefd_m, -nlef(X_old)]), n, 1);
            else
                DX = - nlefd_m^(-1) * nlef(X_old);
            end

        case 'Secant'
            switch options.Smethod
                case 'Two'
                    if iter == 1
                        DX = 0.3 * randn(1, n);
                    end
                    DXNM = diag(DX);
                    XN = DXNM + X_old;
                    Ak = zeros(n, n);
                    for iAk = 1 : n
                        Ak(:, iAk) = (nlef(XN(:, iAk)) - nlef(X_old))/ DX(iAk);
                    end
                    if options.leqsM == 1
                        DX = reshape(ileqs([Ak, -nlef(X_old)]), n, 1);
                    else
                        DX = - Ak^(-1) * nlef(X_old);
                    end
                case 'N'
                    if iter == 1
                        P = eye(n) - diag(ones(1, n-1), 1);
                        XN = 0.2 * randn(n, n) + X_old;
                        Gammak = zeros(n, n);
                        for iGamma = 1 : n
                            Gammak(:, iGamma) = nlef(XN(:, iGamma)) - nlef(X_old);
                        end
                        Hk = XN - X_old;
                        Hkbar = Hk * P;
                        Gammakbar = Gammak * P;
                    end
                    if options.leqsM == 1
                        DX = Hkbar * reshape(ileqs([Gammakbar, -nlef(X_old)]), n, 1);
                    else
                        DX = Hkbar * Gammakbar^(-1) * -nlef(X_old);
                    end
                    Hkbar = [-DX, Hkbar(:, 1 : end - 1)];
                    Gammakbar = [nlef(X_old) - nlef(X_old + DX), Gammakbar(:, 1 : end - 1)];
            end

        case 'QNewton'
            if iter == 1
                Bk = nlefd(X_old)^(-1);
            end
            DX = - Bk * nlef(X_old);
            Bk = Bk - (Bk * nlef(X_old + DX) * DX' * Bk) / ...
                (DX' * Bk * (nlef(X_old + DX) - nlef(X_old)));
            
    end

    X_new = X_old + DX;
    errorX = (norm(DX, options.p))^(1/options.p);
    X_iter(:, iter) = X_new;
    errorFun(iter) = (norm(nlef(X_old), options.p))^(1/options.p);

    if errorFun(iter) <= options.TolFun || errorX <= options.TolX
        break
    end
    X_old = X_new;
end

if options.Display
    fprintf(['Algorithm-' mfilename ' stop at the %d-th iteration by ' ...
        method ' method.\n'], iter);
end

Xstar = X_new;

if strcmp(options.PlotFcns, 'on')
    plot(errorFun, 'LineWidth', 2)
    title('Error of equations in each iteration')
    xlabel('Iteration number')
    ylabel('Error of equations')
end

