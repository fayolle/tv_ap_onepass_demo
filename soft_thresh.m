function [ D ] = soft_thresh( D, T )

if nargin == 1 || isempty(T)
    T = median(D(:));
end

beta = 100.0;
D = softplus(D, T, beta);

end
