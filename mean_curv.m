function [mcurv] = mean_curv(img)

% The image is interpreted as a height surface: z = f(x, y) 
% 2 H = div((\grad f - \grad z)/|\sqrt(1 + |\grad f|^2)|)

dx = [ -1  0  1]./2;   
dy = [ -1 ; 0 ; 1 ]./2;

fx = imfilter(img, dx, 'replicate');
fy = imfilter(img, dy, 'replicate');

gradf_norm2 = fx.*fx + fy.*fy;
denom = sqrt(1.0 + gradf_norm2);

twoH = imfilter(fx./denom,dx,'replicate') + imfilter(fy./denom,dy,'replicate');
mcurv = 0.5 .* twoH;

end
