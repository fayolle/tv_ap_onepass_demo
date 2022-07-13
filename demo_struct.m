clear
close all
clc

f = im2double(imread('barbara_color.png')); 
b = 80; 
lambda = 15; 
T=0.4; 
alpha = 1.0; 

figure,imshow(f),title('original');

disc_curv = @(img) discrepancy_curv(img, b);
disc_mcurv = @(img) discrepancy_mean_curv(img, b);
ax = @(D) alpha*D; % a(x,y)
px = @(D) (1-0.5*D); % p(x,y)
ht = @(img) hard_thresh(disc_curv(img));
ss = @(img) soft_step1(disc_curv(img)); 
st = @(img) soft_thresh(disc_curv(img));
ht_c = @(img) hard_thresh(disc_curv(img), T);
ss_c = @(img) soft_step1(disc_curv(img), T);
st_c = @(img) soft_thresh(disc_curv(img), T);

%u = image_struct_admm_1pass(f, lambda, disc_curv, ax, px);
%u = image_struct_admm_1pass(f, lambda, disc_mcurv, ax, px);
%u = image_struct_admm_1pass(f, lambda, ht, ax, px);
u = image_struct_admm_1pass(f, lambda, ss, ax, px);
%u = image_struct_admm_1pass(f, lambda, st, ax, px);
%u = image_struct_admm_1pass(f, lambda, ht_c, ax, px);
%u = image_struct_admm_1pass(f, lambda, ss_c, ax, px);
%u = image_struct_admm_1pass(f, lambda, st_c, ax, px);


figure,imshow(u),title('structure');
figure,imshow(f-u),title('texture');


s = size(f);
hline = 330;
figure,
subplot(3,1,1); plot(f(hline,:)); axis([0 s(2) 0 1]); title('original signal');
subplot(3,1,2); plot(u(hline,:)); axis([0 s(2) 0 1]); title(['smoothed signal, \lambda=',num2str(lambda),', b=',num2str(b)]);
subplot(3,1,3); plot(f(hline,:)-u(hline,:)); axis([0 s(2) -0.5 0.5]); title('difference');

