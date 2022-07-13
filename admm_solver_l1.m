function out = admm_solver_l1(g, H, mu, p, a, opts)

[rows, cols, channels] = size(g);

max_itr = getoptions(opts, 'max_itr', 20);
tol = getoptions(opts, 'tol', 1e-3);
alpha = getoptions(opts, 'alpha', 0.7);
gamma = getoptions(opts, 'gamma', 2);
rho_r = getoptions(opts, 'rho_r', 2);
rho_o = getoptions(opts, 'rho_o', 50);

f = g;
y1 = zeros(rows, cols, channels);
y2 = zeros(rows, cols, channels);
z = zeros(rows, cols, channels);

HtH = abs(fftn(H, [rows cols channels])).^2;
DtD = abs(fftn([1 -1], [rows cols channels])).^2 + abs(fftn([1 -1]', [rows cols channels])).^2;
Htg = imfilter(g, H, 'circular');

[Df1, Df2] = Grad_(f);
w = imfilter(f, H, 'circular') - g;

r_norm = sqrt(norm(Df1(:))^2 + norm(Df2(:))^2);

for itr=1:max_itr
    v1 = Df1 + (1/rho_r)*y1;
    v2 = Df2 + (1/rho_r)*y2;
    v = sqrt(v1.^2 + v2.^2);
    v(v==0) = 1e-6; 
    Df = sqrt(Df1.^2 + Df2.^2);
    v = max(v - (1/rho_r).*Df.^(p-1).*a, 0)./v;
    u1 = v1.*v;
    u2 = v2.*v;
    
    r = max(abs(w + 1/rho_o*z)-mu/rho_o, 0).*sign(w+1/rho_o*z);
    
    f_old = f;
    rhs = rho_o*Htg + imfilter(rho_o*r-z, H, 'circular') + Dive_(rho_r*u1-y1, rho_r*u2-y2);
    A = rho_o*HtH + rho_r*DtD;
    f = real(ifftn(fftn(rhs)./A));
    
    [Df1, Df2] = Grad_(f);
    w = imfilter(f, H, 'circular') - g;
    
    y1 = y1 - rho_r*(u1 - Df1);
    y2 = y2 - rho_r*(u2 - Df2);
    z = z - rho_o*(r - w);
    
    r_norm_old  = r_norm;
    r_norm = sqrt(norm(Df1(:)-u1(:), 'fro')^2 + norm(Df2(:)-u2(:), 'fro')^2);
    
    if r_norm > (alpha*r_norm_old)
        rho_r  = rho_r * gamma;
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
