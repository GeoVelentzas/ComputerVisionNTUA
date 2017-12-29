clear all; close all; clc;
addpath(genpath('./'));
%% Part1

%Load skin samples
load('skinSamplesRGB.mat');

%convert RGB values of samples to YCbCr values
skinSamplesYCbCr = rgb2ycbcr(skinSamplesRGB);
Cb = im2double(skinSamplesYCbCr(:,:,2));
Cr = im2double(skinSamplesYCbCr(:,:,3));

%calculate mean mu = [mu_cb mu_cr] and covariance matrix S
mu = [mean(Cb(:)) mean(Cr(:))];
S = cov([Cb(:) Cr(:)]);

%Plot skin texture propability density distribution
%Create Cb/Cr 2D space
X=linspace(0,1,255);
Y=linspace(0,1,255);
[mX, mY]=meshgrid(X,Y);
%Calculate propability density for all space
Z=mvnpdf([mX(:) mY(:)],mu,S);
Z=reshape(Z,length(Y),length(X));
%Plot it
figure
subplot(1,2,1)
mesh(X,Y,Z)
xlabel 'Cb'
ylabel 'Cr'
zlabel 'Propability Density'
subplot(1,2,2)
X=linspace(mu(1)-0.1,mu(1)+0.1,100);
Y=linspace(mu(2)-0.1,mu(2)+0.1,100);
[mX, mY]=meshgrid(X,Y);
Z=mvnpdf([mX(:) mY(:)],mu,S);
Z=reshape(Z,length(Y),length(X));
subplot(1,2,2)
mesh(X,Y,Z)
xlabel 'Cb'
ylabel 'Cr'
zlabel 'Propability Density'

%Read first frame
Irgb = im2double(imread('1.png'));

%calculate bounding box from first frame
Irgb = im2double(imread('1.png'));
boxInit = fd(Irgb,mu,S);

%put images in a 3d matrix for gray images in time
I_all = zeros(size(Irgb,1),size(Irgb,2),72);
for i = 1:72;
    Irgb = im2double(imread([num2str(i),'.png']));
    I_all(:,:,i) = rgb2gray(Irgb);

end

%% Part 2

%initialization for optical flow vector at zero
dx0 = zeros(boxInit(4),boxInit(3));
dy0 = zeros(boxInit(4),boxInit(3));

%Setting undersampling factor for optical flow field
res = 0.3;

%Initializing bounding box bank
box=zeros(4,72);
box(1,1)=boxInit(1);
box(2,1)=boxInit(2);
box(3,1)=boxInit(3);
box(4,1)=boxInit(4);

%Initializing displacement bank
DISPL = zeros(2,72-1);

%Performing optical flow tracking for all video frames
for i = 1:72-1;
    %Calculating bounding box edges
    row_start = round(box(2,i));
    row_end = row_start+box(4,i)-1;
    col_start = round(box(1,i));
    col_end = col_start+box(3,i)-1;
    
    %Cropping frame
    I1 = I_all(row_start:row_end,col_start:col_end,i);
    I2 = I_all(row_start:row_end,col_start:col_end,i+1);
    
    %Performing Lucas-Kanade optical flow extraction
    [dx,dy] = lk(I1,I2,1,0.06,dx0,dy0);
    %Thersholding based on flow energy
    [DISPL(1,i) DISPL(2,i)] = displ(dx,dy);
    %Constructing display optical field
    dx_r = imresize(dx,res);
    dy_r = imresize(dy,res);
    DX(:,:,i) = dx_r(:,:);
    DY(:,:,i) = dy_r(:,:);
    %Constructing energy field
    ENERGY(:,:,i) = flipud((dx(:,:).^2+dy(:,:).^2));
    
    %Updating bouding box position
    box(1,i+1)=box(1,i)-DISPL(1,i);
    box(2,i+1)=box(2,i)-DISPL(2,i);
    box(3,i+1)=box(3,i);
    box(4,i+1)=box(4,i);
end

%make a vector for representation;
max_displ_x = max(abs(DISPL(1,:)));
max_displ_y = max(abs(DISPL(2,:)));

%Displaying video with overlay and side-panel information
scrsz = get(0,'ScreenSize');
figure('Position',[0 0 scrsz(3) scrsz(4)]);
subplot(2,4,[1 2 5 6]); imshow(I_all(:,:,1)); title('Image');
pause(0.1);
for i = 1:72-1
    subplot(2,4,[1 2 5 6]); imshow(I_all(:,:,i+1)); title('Image');
    rectangle('EdgeColor',[1 1 1],'Position',round(box(:,i+1)));
    subplot(2,4,3); quiver(-DX(:,:,i),DY(:,:,i));  title('Optical Flow'); axis equal; 
    subplot(2,4,4); imshow(ENERGY(:,:,i),[]); title('Optical Flow Energy');
    subplot(2,4,[7 8]); quiver(-DISPL(1,i),DISPL(2,i),'LineWidth',3);
    title(['Estimated Move of the Face. Pixels Down : ', num2str(-DISPL(2,i)), ', Pixels Left : ',num2str(DISPL(1,i))]);
    xlim([-max_displ_x+1 max_displ_x+1]); ylim([-max_displ_y+1 max_displ_y+1]); 
    pause(0.1);
end


