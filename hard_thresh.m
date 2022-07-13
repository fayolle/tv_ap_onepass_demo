function [ D ] = hard_thresh( D, T )

if nargin == 1 || isempty(T)
    T = median(D(:));
end

%D(D>T)=1;
D(D<=T)=0;

end
