function structure = image_struct_admm_1pass(f, lambda, disc_fun, ax_fun, px_fun, method)
%
% Compute the level-set curvature of the input image and the curvature
% filtered by a low pass filter, then use the difference for computing
% a(x) and p(x).
% Minimize the weighted p-Dirichlet energy
% with L2 or L1 constraint fitting term.
%

if ~exist('method', 'var') || isempty(method)
    method = 'l2';
end

H = fspecial('gaussian',[7,7],1);
%g = imfilter(f,H);
g = imfilter(f,H,'circular');

D = disc_fun(f);
%D = disc_fun(g);

a = ax_fun(D);
p = px_fun(D);

% Setup parameters
opts.rho_r   = 4;
opts.rho_o   = 100;
opts.alpha   = 0.7;
opts.max_itr  =  20;

switch method
    case 'l2'
        out = admm_solver_l2(g, H, lambda, p, a, opts);
    case 'l1'
        out = admm_solver_l1(g, H, lambda, p, a, opts);
    otherwise
        % use 'l2' by default
        out = admm_solver_l2(g, H, lambda, p, a, opts);
end

structure = out.f;

end
