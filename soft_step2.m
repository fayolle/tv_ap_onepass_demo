function [ D2 ] = soft_step2( D, T )

if nargin == 1 || isempty(T)
    T = min(0.8*mean(D(:)),0.4);
end

st = 5.0; %10.0;
D2 = 1./(1+exp(-st.*(D-T)));

end
