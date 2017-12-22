clear; close all; clc;
addpath(genpath('./'));


%% ********************* PART 1 EDGE DETECTION ****************************

% ********* chose either J10 or J20 for different noise levels ************ 
I=imread('edgetest_11.png');        %read the image
I=double(I)/255;                    %normalize values from 0 to 1
Imax=max(I(:));                     %compute the maximum value
Imin=min(I(:));                     %compute the minimum value
PSNR=10;                            %set PSNR=10dB
sigma=(Imax-Imin)/(10^(PSNR/20));   %compute optimal sigma
var=sigma^2;                        %compute variance
J10=imnoise(I,'gaussian',0,var);    %add noise
PSNR=20;                            %set PSNR=20dB
sigma=(Imax-Imin)/(10^(PSNR/20));   %compute sigma
var=sigma^2;                        %compute variance
J20=imnoise(I,'gaussian',0,var);    %add noise

figure;
set(gcf, 'Position', get(0, 'Screensize'));
subplot(1,3,1); imshow(I);   title('Default image without noise');
subplot(1,3,2); imshow(J10); title('Image with gaussian noise (PSNR=10dB)');
subplot(1,3,3); imshow(J20); title('image with gaussian noise (PSNR=20dB)');

%******************** parameters for EdgeDetect **************************
sigma=1.4;                      %standard deviation for EdgeDetecet
th_ed=0.16;                     %threshold levels
Lap=1;                          %linear laplacian for L=0 (else non-linear)
J=J20;                          %J10 for image with PSNR-10dB else J20
%*************************************************************************%

%************** edge detection with morphological operator ***************%
I=imread('edgetest_11.png');    %image to read
SE = strel('diamond', 1);       %morphological operator
Idil = imdilate(I,SE);          %dilation
Iero = imerode(I,SE);           %erosion
M = Idil-Iero;                  %difference between dilation and erosion
T=(M > th_ed);                  %true edge detection
%*************************************************************************%

%*************** edge detection with laplace *****************************%
D = EdgeDetect(J,sigma,th_ed,Lap);  %EdgeDetection with laplacian
Z=(D&T);                            %image of true positives
PrDT=sum(Z(:))/sum(T(:));           %Pr(D|T)
PrTD=sum(Z(:))/sum(D(:));           %Pr(T|D)
C=(PrDT+PrTD)/2;                    %quality of edge detection
%*************************************************************************%

%************************* show results ***********************************
figure;
set(gcf, 'Position', get(0, 'Screensize'));
subplot(1,2,1); imshow(T);
title('True Edges');		
subplot(1,2,2); imshow(D);
title('Edge detection with Laplacian');
%**************************************************************************





%% ******************** PART 2 - CORNERS DETECTION ************************
pause(1.0); clear; clc;


%************************** initial image ********************************%
I = imread('blox.png');                   %read image
I = double(I);                            %convert to double
%*************************************************************************%

s=1.1;                                             %parameter sigma
p=1.8;                                             %parameter p 
Gs = fspecial('gaussian', ceil(3*s)*2+1 ,s);       %Kernel Gs
Gp = fspecial('gaussian', ceil(3*p)*2+1 ,p);       %Kernel Gp
Is = imfilter(I,Gs,'symmetric');                   %filtered image Is
[Dx,Dy] = gradient(Is);                            %calculate gradient
J1 = imfilter(Dx.*Dx,Gp);                          %calculate J1
J2 = imfilter(Dx.*Dy,Gp);                          %calculate J2
J3 = imfilter(Dy.*Dy,Gp);                          %calculate J3

lamda1=(1/2)*(J1+J3+((J1-J3).^2+4*(J2.^2)).^(1/2));   %lambda+ eigenvalue
lamda2=(1/2)*(J1+J3-((J1-J3).^2+4*(J2.^2)).^(1/2));   %lambda- eigenvalue

k = 0.03;                                       %parameter k
R = (lamda1.*lamda2)-k*((lamda2+lamda1).^2);    %"angularity" computation

th_corn = 0.009;                    %parameter th_corn
B_sq = strel('square',3);           %create a 3x3 window
Cond1 = ( R==imdilate(R,B_sq));     %first condition
Rmax = max(R(:));                   %maximum value
Cond2 = ( R > th_corn*Rmax);        %second condition
Corners = (Cond1&Cond2);            %final image of corners

Ic = markCorn(I/255,Corners,[0,0,1],2); 

% ************************** visualizations ******************************
figure;
set(gcf, 'Position', get(0, 'Screensize'));
subplot(2,3,1); imshow(I,[]); title('initial image')
subplot(2,3,2); imshow(lamda1,[]); title('\lambda +');
subplot(2,3,3); imshow(lamda2,[]); title('\lambda -');
subplot(2,3,4); imshow(R,[]); title('angularity matrix')
subplot(2,3,5); imshow(Corners); title('Corners')
subplot(2,3,6); imshow(Ic,[]); title('Corners Marked on Image')
%*************************************************************************











