function out=lbutter(im,d,n)
% LBUTTER(IM,D,N) creates a low-pass Butterworth filter 
% of the same size as image IM, with cutoff D, and order N
%
% Use:
%   x=imread('cameraman.tif');
%   l=lbutter(x,25,2);
%
height=size(im,1);
width=size(im,2);
[x,y]=meshgrid(-floor(width/2):floor((width-1)/2),-floor(height/2): ...
	       floor((height-1)/2));
% z = sqrt(x.^2+y.^2);
% out = (z<d);
out=1./(1+(sqrt(2)-1)*((x.^2+y.^2)/d^2).^n);
