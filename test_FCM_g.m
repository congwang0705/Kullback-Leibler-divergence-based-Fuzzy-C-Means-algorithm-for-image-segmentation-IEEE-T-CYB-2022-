%%
getd = @(p)path(p,path);
%%
clc; clear all; close all;
Io = imread('sy2.bmp');
rng('default');
[m, n, p] = size(Io);
if p~=1
    Io=rgb2gray(Io);
    p=1;
end

figure(1)
imshow(Io);

%%
I=imnoise(uint8(Io),'gaussian',0,30^2/255^2);
I=imnoise(uint8(I),'salt & pepper',0.15);
I=uint8(I);
figure(2)
imshow(I);
%%
[counts,x]=imhist(I);%Grayscale statistics
figure(3)
imhist(I);
%%
se=3;
If=w_recons_CO(double(I),strel('square',se));
If=uint8(If);
figure(4)
imshow(If);
%%
X = reshape(double(I), m*n, p);
Xbar = reshape(double(If), m*n, p);
k = 4; b = 1;alpha=15000;beta=3;
tic;[KL, eta, P,C, dist, J,P2] = KLFCM(X,Xbar,k,b,m,n,alpha,beta);toc;

%%
[~, label] = min(dist, [], 2);
label=reshape(label, m, n, p);
Is=reshape(C(label, :), m, n, p);
figure(5)
imshow(uint8(Is),'border','tight')
figure(6)
plot(1:length(J), J, 'r-*'), xlabel('#iterations'), ylabel('objective function')


Ce=reshape(double(C), size(C,1)*size(C,2), p);
t = P.^b;
R = (t*Ce)./(sum(t')'*ones(1, p));
RR=reshape(R, m, n, p);

figure(7)
imshow(uint8(RR),'border','tight')

RR=Optivalue(RR,k);
figure(8)
imshow(uint8(RR),'border','tight')

