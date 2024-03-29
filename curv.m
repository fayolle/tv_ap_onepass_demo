function Icurv = curv(I)
dx = [ -1  0  1]./2;
dy = [ -1 ; 0 ; 1 ]./2;

a = 0.001;

Ix = imfilter(I,dx,'replicate');
Iy = imfilter(I,dy,'replicate');

Igrad2 = Ix.*Ix + Iy.*Iy;
Igrad1 = sqrt(Igrad2);

Icurv = imfilter(Ix./(Igrad1+a),dx,'replicate') + imfilter(Iy./(Igrad1+a),dy,'replicate');
end

