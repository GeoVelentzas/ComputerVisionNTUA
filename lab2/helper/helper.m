close all; clear all; clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f = im2double(imread('D17D55.tiff'));
f_ref = imread('D17D55_ref.gif');
Fs = max(size(f));
[M,N] = size(f);         %image size to M,N
u = (0:N-1)-floor(N/2);  %create vector for x-frequency
v = (0:M-1)-floor(N/2);  %create vector for y-frequency
F = fft2(f,M,N);         %calculate fft
F_sh = fftshift(F);
F_14 = F_sh(:,(N/2+1 : end));
[Fmax, Imax, Fmin, Imin] = extrema2(abs(F_14));  %find maximums
[I,J] = ind2sub(size(F_14),Imax);                %in row,col
max_coor = [J-1,I-N/2-1];    %indices of maximum frequencies in x,y 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
K = 10;
a = 1;
energy_smoother = 1.2;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
b = a;
n = 1;
i = 1;
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


g = zeros(N,N,K); %initializations
Omega = zeros(1,K);
theta = zeros(1,K);

for i = 1:K
    h = gabor2D(K_max_points(i,1),K_max_points(i,2),a,b,Fs);
    g(:,:,i) = conv2(f,h,'same');  %compute convolution for every filter
    Omega(i) = norm(K_max_points(i,:),2);
    theta(i) = atan2(K_max_points(i,2),K_max_points(i,1));
end

E = abs(g);

    
% sigma = energy_smoother*a;
% hsize = ceil(3*sigma)*2+1;
% smoother = fspecial('gaussian',hsize,sigma);
% 
% figure;
% for i = 1:K
%     E(:,:,i) = imfilter(E(:,:,i),smoother,'replicate','conv','same');
%     %imshow(E(:,:,i),[]); pause;
% end


P1 = im2double(imread('pattern1.png'));
P2 = im2double(imread('pattern2.png'));
%figure;
%subplot(1,2,1); imshow(P1); title('πρότυπο 1')
%subplot(1,2,2); imshow(P2); title('πρότυπο 2')


g1 = zeros(N,N,K); %initializations
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


M = max(E,[],3);

M = mat2gray(M);

for thB = 0.2:0.1:0.9
B = (M >thB);
figure;
imshow(B);title(['K =',num2str(K), '  a =', num2str(a), '  θ_B = ',num2str(thB)])
pause;
end


BS1 = 0*E;
BS2 = 0*E;
for i = 1:K
    BS1(:,:,i) = normpdf(E(:,:,i),mean1(i),std1(i));
    BS2(:,:,i) = normpdf(E(:,:,i),mean2(i),std2(i));
end


logP = sum(log10((BS1)./(BS2)),3);
logP(logP==-Inf)= min(logP(logP>-Inf));
logP(logP==Inf) = max(logP(logP<Inf));
logP = mat2gray(logP);

for thlog = 0.2:0.1:0.9
figure; imshow(logP>thlog); title(['K =',num2str(K), '  a =', num2str(a), ', θ_{log} = ', num2str(thlog)]);
pause;
end























