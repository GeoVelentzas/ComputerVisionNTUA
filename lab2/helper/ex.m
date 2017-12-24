close all; clear all; clc;

%% 1.1
f = im2double(imread('D17D55.tiff'));
f_ref = imread('D17D55_ref.gif');
Fs = max(size(f)); %sampling frequency

figure; 
subplot(1,2,1); imshow(f);      title('initial image');
subplot(1,2,2); imshow(f_ref);  title('segmentation');

%% 1.2 gabor filters (an example)
U0 = 10; V0 = 10; a = 25; b = 25; %parameters
h = gabor2D(U0,V0,a,b,Fs);        %create the filter

N = floor(ceil((3*Fs/a)+1)/2);    %we need these to show the result
[x,y] = meshgrid(-N:1:N,-N:1:N);  %although they are also computed in function

H = freqz2(h,Fs,Fs);              %compute the frequency response
[fx,fy]= meshgrid(-Fs/2+1:Fs/2,-Fs/2+1:Fs/2); %grid for both positive frequencies 

figure;
mesh(x,y,real(h)); title('U_0 = 10 , V_0 = 10, a = b = 25, space domain'); 
axis square;
figure; 
mesh(fx,fy,fliplr(flipud(abs(H)))); title('U_0 = 10 , V_0 = 10, a = b = 25, frequency domain');
axis square;

%% 1.3
[M,N] = size(f);         %image size to M,N
u = (0:N-1)-floor(N/2);  %create vector for x-frequency
v = (0:M-1)-floor(N/2);  %create vector for y-frequency
F = fft2(f,M,N);         %calculate fft

figure; 
mesh(u*pi/(N/2),v*pi/(M/2),log(1+abs(fftshift(F)))); %show FFT in [-pi,pi]
title('log(abs(F)+1)'); %to be positive
colormap(jet);

F_sh = fftshift(F);
F_14 = F_sh(:,(N/2+1 : end));
[Fmax, Imax, Fmin, Imin] = extrema2(abs(F_14));  %find maximums
[I,J] = ind2sub(size(F_14),Imax);                %in row,col
max_coor = [J-1,I-N/2-1];    %indices of maximum frequencies in x,y 


n = 1;
i = 1;
K = 10;
K_max_points = [];            %K maximum frequencies in u,v couples
while n<K+1
    if max_coor(i,1)==0
        if max_coor(i,2)>0
            K_max_points = [K_max_points ; max_coor(i,:)];
            n = n+1;
        end
    else K_max_points = [K_max_points ; max_coor(i,:)];
         n = n+1;
    end
    i = i+1;
end

K_max_freq = K_max_points*2*pi/N;        %K maximum frequencies in rad


%% 1.4
a=20;
b=20;
g = zeros(M,N,K); %initializations
Omega = zeros(1,K);
theta = zeros(1,K);

for i = 1:K
    h = gabor2D(K_max_points(i,1),K_max_points(i,2),a,b,Fs);
    g(:,:,i) = conv2(f,h,'same');  %compute convolution for every filter
    Omega(i) = norm(K_max_points(i,:),2);
    theta(i) = atan2(K_max_points(i,2),K_max_points(i,1));
end


for i= 1:K
    figure;
    subplot(1,2,1); imshow(abs(g(:,:,i)),[]);
    title(['E(x,y) for Ω = ', num2str(Omega(i)), ' and θ = ', num2str(theta(i))])
    subplot(1,2,2); imshow(angle(g(:,:,i)),[]);
    title(['phase(g(x,y)) for Ω = ', num2str(Omega(i)), ' and θ = ', num2str(theta(i))])
end


%% 1.5
P1 = im2double(imread('pattern1.png'));
P2 = im2double(imread('pattern2.png'));
figure;
subplot(1,2,1); imshow(P1); title('πρότυπο 1')
subplot(1,2,2); imshow(P2); title('πρότυπο 2')


g1 = zeros(M,N,K); %initializations
Omega1 = zeros(1,K);
theta1 = zeros(1,K);
g2 = g1; 
mean1 = zeros(1,K);
mean2 = zeros(1,K);
std1 = zeros(1,K);
std2 =zeros(1,K);

for i = 1:K
    h = gabor2D(K_max_points(i,1),K_max_points(i,2),a,b,Fs);
    g1(:,:,i) = conv2(P1,h,'same');  %compute convolution for every filter
    g2(:,:,i) = conv2(P2,h,'same');  %compute convolution for every filter
    Omega(i) = norm(K_max_points(i,:),2);
    theta(i) = atan2(K_max_points(i,2),K_max_points(i,1));
    mean1(i) = mean(reshape(abs(g1(:,:,i)),1,[]));
    mean2(i) = mean(reshape(abs(g2(:,:,i)),1,[]));
    std1(i) = std(reshape(abs(g1(:,:,i)),1,[]));
    std2(i) = std(reshape(abs(g2(:,:,i)),1,[]));
end

figure;
hold on;
for i = 1:K
    x = mean1(i)-4*abs(std1(i)): 0.01: mean1(i)+4*abs(std1(i));
    y = gaussmf(x,[std1(i) , mean1(i)]);
    plot(x,y); 
end
title('Γκαουσιανή προσέγγιση τιμών ενέργειας για εικόνα πρότυπο 1')


figure;
hold on;
for i = 1:K
    x = mean2(i)-4*abs(std2(i)): 0.01: mean2(i)+4*abs(std2(i));
    y = gaussmf(x,[std2(i) , mean2(i)]);
    plot(x,y); 
end
title('Γκαουσιανή προσέγγιση τιμών ενέργειας για εικόνα πρότυπο 2')



for i = 1:K
    figure; hold on; box on
    x = mean1(i)-4*abs(std1(i)): 0.01: mean1(i)+4*abs(std1(i));
    y = gaussmf(x,[std1(i) , mean1(i)]);
    plot(x,y,'r','linewidth',2);
    x = mean2(i)-4*abs(std2(i)): 0.01: mean2(i)+4*abs(std2(i));
    y = gaussmf(x,[std2(i) , mean2(i)]);
    plot(x,y,'m','linewidth',2) 
    legend('Π1','Π2')
    title(['Φιλτρο ', num2str(i)])
end


%% 1.6

E1 = abs(g1);
E2 = abs(g2);
E = abs(g);

M1 = max(E1,[],3);
M2 = max(E2,[],3);
M = max(E,[],3);

figure; 
subplot(1,2,1); imshow(M1,[]); title('M1(x,y)');
subplot(1,2,2); imshow(M2,[]); title('M2(x,y)');
colormap(jet);

figure;
imshow(M,[]); title('M(x,y)')
colormap(jet);


thresh = 0.5:0.1:0.7;
B = zeros(length(M),length(M),length(thresh));
for i = 1:length(thresh)
    B(:,:,i) = (M > thresh(i)*max(M(:)));
    figure;
    imshow(B(:,:,i)); title(['threshold th_B = ', num2str(thresh(i))])
end    

thb = graythresh(mat2gray(M));


BS1 = 0*E;
BS2 = 0*E;
for i = 1:K
    BS1(:,:,i) = normpdf(E(:,:,i),mean1(i),std1(i));
    BS2(:,:,i) = normpdf(E(:,:,i),mean2(i),std2(i));
    figure;
    subplot(1,2,1); imshow(BS1(:,:,i)); title(['P(E',num2str(i),'|t1)']);
    subplot(1,2,2); imshow(BS2(:,:,i)); title(['P(E',num2str(i),'|t2)']);
    pause;
end


logP = sum(log10((BS1)./(BS2)),3);
logP(logP==-Inf)= min(logP(logP>-Inf));
logP(logP==Inf) = max(logP(logP<Inf));
logP = mat2gray(logP);

thlog = 0.80;
%figure; imshow(logP); title('logP');
logP = logP>thlog;
figure; imshow(logP); title(['thresholded logP with th_{log} = ', num2str(thlog)]);

SE1 = strel('disk',20,0);
SE2 = strel('disk',40,0);
logP = imclose(logP,SE1);
logP = imopen(logP,SE1);
logP = imfill(logP,'holes');
logP = imclose(logP,SE2);
figure; imshow(logP); title('after morphological filtering')














































