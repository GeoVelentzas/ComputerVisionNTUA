close all; clear; clc;

addpath(genpath('./'));

plithos = 43;  %choose number of images (max 43)


for i = 1:plithos
    if i<10
        text = ['Frame_no.00000', num2str(i), '.jpg'];
    else
        text = ['Frame_no.0000', num2str(i), '.jpg'];
    end
    I{i} = double(rgb2gray(imread(text)))/255;
    F(:,:,i) = I{i};
end


% number of lines and columns
[lin,col] = size(I{1});

% initial condition
dx0 = zeros(lin,col);
dy0 = zeros(lin,col);
res =0.3;

h = waitbar(0,'Please wait...');
for i = 1:plithos-1
    [dx,dy] = LucasKanade(F(:,:,i),F(:,:,i+1),1.5,0.1,dx0,dy0);
    dx_r = imresize(dx,res);
    dy_r = imresize(dy,res);
    DX(:,:,i) = dx_r(:,:);
    DY(:,:,i) = dy_r(:,:);
    MAGN(:,:,i) = (dx(:,:).^2+dy(:,:).^2).^(1/2);
    waitbar(i/(plithos-1))
end
close(h)

figure;
set(gcf, 'Position', get(0, 'Screensize'));
for i = 1:plithos-1
	subplot(1,3,1);
    imshow(F(:,:,i));
    title('original image');
    subplot(1,3,2);
    quiver(-DX(:,:,i),-DY(:,:,i)); axis square;
    title('optical flow vectors');
    subplot(1,3,3);
  	imshow(flipud(MAGN(:,:,i)),[]);
    title('optical flow magnitude');
    pause(0.1);
end


