close all; clear; clc;

addpath(genpath('./'))

%% *************** PART 1 - MORPHOLOGICAL SEGMENTATION *******************%
f = imread('stemmcells.tif');
f=double(f);
% Gs = fspecial('gaussian', 3, 0.1);
% f = imfilter(f,Gs);

B = strel('disk',4,8);
op = imopen(f,B);
g = f - op;
sat = 0.05;
gmax = max(g(:));
G = (g > gmax*sat);

h = 170;
m = f - h;
frec = imreconstruct(m,f);
g2 = f - frec;
G2 = (g2 > 20);
[L1,num1] = bwlabel(G,8);
[L2,num2] = bwlabel(G2,8);



figure; imshow(zeros(1,500)); title('A. Minkowski Opening Method');
xlabel('press any key')
pause;
figure; imshow(f,[]); title('Initial Image');
xlabel('press any key')
pause;
figure; imshow(op,[]); title('Minkowski Opening of Image');
xlabel('press any key')
pause;
figure; imshow(g,[]); title('top-hat');
xlabel('press any key')
pause;
figure; imshow(G,[]); title('Thresholded top-hat');
xlabel('press any key')
pause;
figure; imshow(L1,[]); title(['Connected Components (=',num2str(num1),')']);
pause;
figure; imshow(zeros(1,500)); title('B. Method with Reconstruction Opening');
xlabel('press any key')
pause;
figure; imshow(f,[]); title('Initial Image');
xlabel('press any key')
pause;
figure; imshow(frec,[]); title('Reconstruction Opening of Image');
xlabel('press any key')
pause;
figure; imshow(g2,[]); title('top-hat');
xlabel('press any key')
pause;
figure; imshow(G2,[]); title('Thresholded top-hat');
xlabel('press any key')
pause;
figure; imshow(L2,[]); title(['Connected Components (=',num2str(num2),')']);
xlabel('press any key')
pause;




%% ******************* PART 2 - WATERSHED ********************************%

clear; close all; clc

g = imread('prostate.tif');
g = double(g);

figure; imshow(zeros(1,500)); title('C. Watershed Method')
xlabel('press any key')
pause;

G = g;
for r = 1:3
    B = strel('disk',r,8);
    m1 = imerode(G,B);
    G = imreconstruct(m1,G);

    m2 = imdilate(G,B);
    G = imreconstructclosing(m2,G);
end
figure; imshow(g,[]); title('Initial Image');
xlabel('press any key')
pause;
figure; imshow(G,[]); title('Image after ASF');
xlabel('press any key');
pause;

f = G;
h = 20;
a = 0.3;
Regmin = imreconstructclosing(f+h,f) - f;
Min = (Regmin > a*h);
figure; imshow(Regmin,[]); title('Regional minima');
xlabel('press any key')
pause;
figure; imshow(Min,[]); title('Inside Markers');
xlabel('press any key')
pause;



I = imimposemin(G,Min);
figure; imshow(I,[]); title('Change of topology on filtered image');
xlabel('press any key')
pause;

M = watershed(I);
figure; imshow(M,[]); title('Watershed of the previous image');
xlabel('press any key')
pause;

Mout = (M == 0);
figure; imshow(Mout); title('Outside Markers');
xlabel('press any key')
pause;

markers = Min|Mout;
imshow(markers); title('Total Markers');
xlabel('press any key')
pause;
[Gx ,Gy] = gradient(G);
DG = (Gx.^2+Gy.^2).^(1/2);
% Gs = fspecial('gaussian', 18, 2); %see report
% DG = imfilter(DG,Gs);
GF = imimposemin(DG,markers);
figure; imshow(GF,[]); title('Gradient magnitude after topological change');
xlabel('press any key')
pause;
W = watershed(GF);
figure; imshow(W,[]); title('Watershed on the previous image');
xlabel('press any key')
pause;
F = (W==0);
figure; imshow(F); title('Boundaries of Basins');
xlabel('press any key')
pause;
g2 = g/max(g(:));
final = imoverlay(g2,F,[0 0 0]);
figure; imshow(final,[]); title('Final Segmentation');













