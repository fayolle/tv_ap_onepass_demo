function [ D2 ] = soft_step1( D, T )

if nargin == 1 || isempty(T)
    T = min(0.8*mean(D(:)),0.4);
end

D2 = (atan(10*(D-T))/pi+0.5).*D;

end
