function out = admm_solver_l2(g, H, mu, p, a, opts)

[rows, cols, channels] = size(g); 

max_itr = getoptions(opts, 'max_itr', 20);
tol = getoptions(opts, 'tol', 1e-3);
alpha = getoptions(opts, 'alpha', 0.7);
gamma = getoptions(opts, 'gamma', 2);
rho = getoptions(opts, 'rho_r', 2);

f = g; 
y1 = zeros(rows, cols, channels); 
y2 = zeros(rows, cols, channels);
u1 = zeros(rows, cols, channels); 
u2 = zeros(rows, cols, channels);

HtH = abs(fftn(H, [rows cols channels])).^2;
DtD = abs(fftn([1 -1], [rows cols channels])).^2 + abs(fftn([1 -1]', [rows cols channels])).^2;
Htg = imfilter(g, H, 'circular');

[Df1, Df2] = Grad_(f);

r_norm = sqrt(norm(Df1(:))^2 + norm(Df2(:))^2);

for itr=1:max_itr    
    f_old = f;
    rhs = fftn((mu/rho).*Htg + Dive_(u1-(1/rho)*y1, u2-(1/rho)*y2));
    A = (mu/rho).*HtH + DtD; 
    f = real(ifftn(rhs./A));
    
    [Df1, Df2] = Grad_(f);
    v1 = Df1 + (1/rho)*y1;
    v2 = Df2 + (1/rho)*y2;
    v = sqrt(v1.^2 + v2.^2);
    v(v==0) = 1;
    Df = sqrt(Df1.^2 + Df2.^2);
    v = max(v - (1/rho).*Df.^(p-1).*a, 0)./v;
    u1 = v1.*v;
    u2 = v2.*v;
    
    y1 = y1 - rho*(u1 - Df1);
    y2 = y2 - rho*(u2 - Df2);
    
    r_norm_old = r_norm;
    r_norm = sqrt(norm(Df1(:)-u1(:), 'fro')^2 + norm(Df2(:)-u2(:), 'fro')^2);
    
    if r_norm > (alpha*r_norm_old)
        rho = rho * gamma;
    end
    
    rel_chg = norm(f(:)-f_old(:))/norm(f_old(:));
    if rel_chg < tol
        break
    end
    
end

out.f = f;
out.itr = itr;

end

function [Dux,Duy] = Grad_(U)
Dux = [diff(U,1,2), U(:,1,:) - U(:,end,:)];
Duy = [diff(U,1,1); U(1,:,:) - U(end,:,:)];
end

function DtXY = Dive_(X,Y)
DtXY = [X(:,end,:) - X(:, 1,:), -diff(X,1,2)];
DtXY = DtXY + [Y(end,:,:) - Y(1, :,:); -diff(Y,1,1)];
end
