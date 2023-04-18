function [ D ] = discrepancy_curv( f, b )

K0 = curv(f); % clean curvature

butter = lbutter(f,b,10); % low-pass Butterworth filter
for c=1:size(f,3)
    K(:,:,c) = ifft2(ifftshift(fftshift(fft2(K0(:,:,c))).*butter));
end % smoothed curvature

D = abs(K-K0);
D1 = D; 
for c=1:size(f,3)
    D(:,:,c) = sum(D1,3)/size(D1,3);
end

H3 = fspecial('gaussian',[7,7],3);
D1 = imfilter(D,H3);

D = D1./max(D1(:));

end
