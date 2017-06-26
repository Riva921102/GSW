
% ----------------------GSW algorithm-------------------------- 
% Generate 2017-06-26 
% Last modified 2017-06-26
% Copyright (c) Chengliang Chang & Yijun Qi
clear;
k=0;
m=3;
mm=512;
nn=512;
w=[1,1,1]; % weight 
z=[0.0,0.06,0.02];
r=532e-9; % wave length
f=0.3;
dx=8e-6;
dfx=r*f/dx/nn;
dy=8e-6;
dfy=r*f/dy/mm;
KK=20; % iteration
x=linspace(-dx*nn/2,dx*nn/2,nn);
y=linspace(-dy*mm/2,dy*mm/2,mm);
[x,y]=meshgrid(x,y);
x1m=0; y1m=0;
x2m=200*dfx; y2m=200*dfy;
x3m=-200*dfx; y3m=-200*dfy;

%------random phase theta---------+
theta1=2*pi*rand(1,1);           %|
theta2=2*pi*rand(1,1);           %|
theta3=2*pi*rand(1,1);           %|
%------random phase theta---------+

%------delta-----------------------------------------------------------------------------------+
                                                                                              %|

delta1=z(1)*pi/r/f/f.*(x.^2+y.^2)+2*pi/r/f.*(x1m.*x+y1m.*y);
delta2=z(2)*pi/r/f/f.*(x.^2+y.^2)+2*pi/r/f.*(x2m.*x+y2m.*y);
delta3=z(3)*pi/r/f/f.*(x.^2+y.^2)+2*pi/r/f.*(x3m.*x+y3m.*y);

                                                                                              %|
%------delta-----------------------------------------------------------------------------------+


%/////////iteration/////////////////////////////////////////////////////////////////////////////////////
for k=1:KK                                                                                            %/
	phij=angle(w(1).*exp(1i.*(delta1+theta1))+w(2).*exp(1i.*(delta2+theta2))+w(3).*exp(1i.*(delta3+theta3)));
	V1=mean(mean(exp(1i.*(phij-delta1))));
    V2=mean(mean(exp(1i.*(phij-delta2))));
    V3=mean(mean(exp(1i.*(phij-delta3))));
	w(1)=w(1)*(abs(V1)+abs(V2)+abs(V3))/3/abs(V1);
    w(2)=w(2)*(abs(V1)+abs(V2)+abs(V3))/3/abs(V2);
    w(3)=w(3)*(abs(V1)+abs(V2)+abs(V3))/3/abs(V3);
	theta1=angle(V1);
    theta2=angle(V2);
    theta3=angle(V3);
    
    k
end                                                                                                  %/
%///////iteration//////////////////////////////////////////////////////////////////////////////////////

%///////verification///////////////////////////////////////////////////////
img=double(rgb2gray(imread('pic/B512(3).jpg')))/255;                     %/
fft_img=fftshift(fft2(fftshift(img.*exp(1i*rand(512,512)))));  
fft_out=fft_img.*exp(1i*phij);
output=fftshift(ifft2(fftshift(fft_out)));                            
output=output/max(max(abs(output)));
figure;
imshow(abs(output));
figure;
imshow((angle(output)+pi)/1/pi);
                                                                          %/
%///////verification///////////////////////////////////////////////////////

