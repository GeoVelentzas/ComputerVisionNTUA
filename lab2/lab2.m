close all; clear; clc;

addpath(genpath('./'));
%% 1.1 Reading and plotting input image and ideal segmentation
f = im2double(imread('D17D55.tiff'));
f_ref = imread('D17D55_ref.gif');
Fs = max(size(f)); %sampling frequency

figure;
set(gcf, 'Position', get(0, 'Screensize'));
subplot(1,2,1); imshow(f);      title('initial image');
subplot(1,2,2); imshow(f_ref);  title('segmentation');
suptitle('press any key to proceed');
pause;

%% 1.2 gabor filters (an example)
close all
U0 = 10; V0 = 10; a = 25; b = 25; %parameters
h = gabor2D(U0,V0,a,b,Fs);        %create the filter

N = floor(ceil((3*Fs/a)+1)/2);    %we need these to show the result
[x,y] = meshgrid(-N:1:N,-N:1:N);  %although they are also computed in function

H = freqz2(h,Fs,Fs);              %compute the frequency response
[fx,fy]= meshgrid(-Fs/2+1:Fs/2,-Fs/2+1:Fs/2); %grid for both positive frequencies 

%Plotting filters both in space and frequency
figure;
set(gcf, 'Position', get(0, 'Screensize'));
subplot(1,2,1);
mesh(x,y,real(h)); title('U_0 = 10 , V_0 = 10, a = b = 25, space domain'); 
axis square;
subplot(1,2,2);
mesh(fx,fy,fliplr(flipud(abs(H)))); title('U_0 = 10 , V_0 = 10, a = b = 25, frequency domain');
axis square;
suptitle('press any key to proceed');
pause;
%% 1.3 Image spectrum and spatial frequency selection
close all
[M,N] = size(f);         %image size to M,N
u = (0:N-1)-floor(N/2);  %create vector for x-frequency
v = (0:M-1)-floor(N/2);  %create vector for y-frequency
F = fft2(f,M,N);         %calculate fft

figure;
set(gcf, 'Position', get(0, 'Screensize'));
mesh(u*pi/(N/2),v*pi/(M/2),log(1+abs(fftshift(F)))); %show FFT in [-pi,pi]
title('log(abs(F)+1)'); %to be positive
colormap(jet);

F_sh = fftshift(F);
F_14 = F_sh(:,(N/2+1 : end));
[Fmax, Imax, Fmin, Imin] = extrema2(abs(F_14));  %find maxima
[I,J] = ind2sub(size(F_14),Imax);                %in row,col
max_coor = [J-1,I-N/2-1];    %indices of maximum frequencies in x,y 


n = 1;
i = 1;
K = 10;
K_max_points = [];            %K maximum frequencies number in u,v couples
while n<K+1
    if max_coor(i,1)==0 %If maximum along y axis
        if max_coor(i,2)>0 %Make sure it is on the positive semiaxis
            K_max_points = [K_max_points ; max_coor(i,:)]; 
            n = n+1;
        end
    else K_max_points = [K_max_points ; max_coor(i,:)];
         n = n+1;
    end
    i = i+1;
end

K_max_freq = K_max_points*2*pi/N;        %K maximum frequencies in rad

suptitle('press any key to proceed');
pause;
%% 1.4 Image filtering and energy response
close all
%%%%%%%%%%%Gabor filter coefficients%%%%%%%%%%%%%%%
a=10; %Gabor filter horizontal span coefficient; smaller=longer execution time
b=10; %Gabor filter vertical span coefficient; smaller=longer execution time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['Isomorphic Gabor filter span coefficient a= ',num2str(a)]);
g = zeros(M,N,K); %initializations
Omega = zeros(1,K);
theta = zeros(1,K);

for i = 1:K %For each spectrum peak
    h = gabor2D(K_max_points(i,1),K_max_points(i,2),a,b,Fs); %Create Gabor filter
    g(:,:,i) = conv2(f,h,'same');  %compute convolution for every filter
    Omega(i) = norm(K_max_points(i,:),2); %Compute frequency
    theta(i) = atan2(K_max_points(i,2),K_max_points(i,1)); %Compute angle
end

%Plot filter responses
for i= 1:K
    figure;
    set(gcf, 'Position', get(0, 'Screensize'));
    subplot(1,2,1); imshow(abs(g(:,:,i)),[]);
    title(['E(x,y) for \Omega = ', num2str(Omega(i)), ' and \theta = ', num2str(theta(i))])
    subplot(1,2,2); imshow(angle(g(:,:,i)),[]);
    title(['phase(g(x,y)) for \Omega = ', num2str(Omega(i)), ' and \theta = ', num2str(theta(i))])
end

suptitle('press any key to proceed');
pause;
%% 1.5 Texture pattern training
close all
P1 = im2double(imread('pattern1.png')); %Load pattern 1
P2 = im2double(imread('pattern2.png')); %Load pattern 2
figure;
set(gcf, 'Position', get(0, 'Screensize'));
subplot(1,2,1); imshow(P1); title('pattern 1')
subplot(1,2,2); imshow(P2); title('pattern 2')


g1 = zeros(M,N,K); %initializations
Omega1 = zeros(1,K);  
theta1 = zeros(1,K);
g2 = g1; 
mean1 = zeros(1,K);
mean2 = zeros(1,K);
std1 = zeros(1,K);
std2 =zeros(1,K);

%Extract response characteristics for every filter on each pattern
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
set(gcf, 'Position', get(0, 'Screensize'));
hold on;
for i = 1:K
    x = mean1(i)-4*abs(std1(i)): 0.01: mean1(i)+4*abs(std1(i));
    y = gaussmf(x,[std1(i) , mean1(i)]);
    plot(x,y); 
end
title('gaussian estimation of energy for pattern 1')


figure;
set(gcf, 'Position', get(0, 'Screensize'));
hold on;
for i = 1:K
    x = mean2(i)-4*abs(std2(i)): 0.01: mean2(i)+4*abs(std2(i));
    y = gaussmf(x,[std2(i) , mean2(i)]);
    plot(x,y); 
end
title('gaussian estimation of energy for pattern 2')

for i = 1:K
    figure;
    set(gcf, 'Position', get(0, 'Screensize'));
    hold on; box on;
    x = mean1(i)-4*abs(std1(i)): 0.01: mean1(i)+4*abs(std1(i));
    y = gaussmf(x,[std1(i) , mean1(i)]);
    plot(x,y,'r','linewidth',2);
    x = mean2(i)-4*abs(std2(i)): 0.01: mean2(i)+4*abs(std2(i));
    y = gaussmf(x,[std2(i) , mean2(i)]);
    plot(x,y,'m','linewidth',2) 
    legend('P1','P2')
    title(['filter', num2str(i)])
end

suptitle('press any key to proceed');
pause;
%% 1.6 Segmentation and Sorting
close all
%Create response magnitude
E1 = abs(g1);
E2 = abs(g2);
E = abs(g);
%Find maximum response
M1 = max(E1,[],3);
M2 = max(E2,[],3);
M = max(E,[],3);

figure;
set(gcf, 'Position', get(0, 'Screensize'));
subplot(1,2,1); imshow(M1,[]); title('M1(x,y)');
subplot(1,2,2); imshow(M2,[]); title('M2(x,y)');
colormap(jet);

figure;
set(gcf, 'Position', get(0, 'Screensize'));drawnow;
imshow(M,[]); title('M(x,y)')
colormap(jet);

M = mat2gray(M); %Convert to [0,255] space

%%%%%%%%%%%%%Binary images threshold%%%%%%%%%%%%%%%%%
thB=0.85;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
B=(M>thB);
figure;
set(gcf, 'Position', get(0, 'Screensize'));drawnow;
imshow(B);
title(['Unsupervised method: K =',num2str(K), '  a =', num2str(a), '  \theta_B = ',num2str(thB)])

% %Enable the following code to examine multiple thresholds at once
% for thB = 0.2:0.1:0.9
% B = (M >thB);
% figure;
% imshow(B);title(['K =',num2str(K), '  a =', num2str(a), '  \theta_B = ',num2str(thB)])
% end

%Calculate pattern propability for each filter
BS1 = 0*E;
BS2 = 0*E;
for i = 1:K
    BS1(:,:,i) = normpdf(E(:,:,i),mean1(i),std1(i));
    BS2(:,:,i) = normpdf(E(:,:,i),mean2(i),std2(i));
    figure;
    set(gcf, 'Position', get(0, 'Screensize'));
    subplot(1,2,1); imshow(BS1(:,:,i)); title(['P(E',num2str(i),'|t1)']);
    subplot(1,2,2); imshow(BS2(:,:,i)); title(['P(E',num2str(i),'|t2)']);
end

%Calculate total logarithmic propability
logP = sum(log10((BS1)./(BS2)),3);
logP(logP==-Inf)= min(logP(logP>-Inf));
logP(logP==Inf) = max(logP(logP<Inf));
logP = mat2gray(logP);

%%%%%%%%%%%Supervised method threshold%%%%%%%%%%%%%%
thlog = thB;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
set(gcf, 'Position', get(0, 'Screensize')); drawnow;
imshow(logP); title('logP');
logP = logP>thlog;

figure;
set(gcf, 'Position', get(0, 'Screensize')); drawnow;
imshow(logP); title(['thresholded logP with th_{log} = ', num2str(thlog)]);


%Perform morphological processing
SE = strel('disk',20,0); %Create structuring element (disk)
L2=bwareaopen(logP>thlog,60); %Remove islands
L2=imopen(L2,SE); %Remove peninsulas
L2=imclose(L2,SE); %Rounden shape
L2=imfill(L2,'holes'); %Fill holes

figure;
set(gcf, 'Position', get(0, 'Screensize')); drawnow;
imshow(L2)
title('supervised method after morphological filtering');

suptitle('press any key to proceed');
pause;
%% 1.7 Quantitative result extraction
close all
fScoreTable=zeros(3,1);
sortErrorTable=zeros(3,1);

%Display method comparison visually
figure;
set(gcf, 'Position', get(0, 'Screensize'));
subplot(1,3,1)
imshow(f_ref);
hold on;
perim=bwboundaries(B>thB);
for k = 1:numel(perim)
    plot(perim{k}(:,2), perim{k}(:,1), 'r', 'Linewidth', 2)
end
title(['Unsupervised: a=',num2str(a),', \theta=',num2str(thB)]);

subplot(1,3,2)
imshow(f_ref);
hold on;
perim=bwboundaries(logP>thlog);
for k = 1:numel(perim)
    plot(perim{k}(:,2), perim{k}(:,1), 'r', 'Linewidth', 2)
end
title(['Supervised: a=',num2str(a),', \theta=',num2str(thlog)]);

subplot(1,3,3)
imshow(f_ref);
hold on;
perim=bwboundaries(L2>thlog);
for k = 1:numel(perim)
    plot(perim{k}(:,2), perim{k}(:,1), 'r', 'Linewidth', 2)
end
title(['Supervised+Morphological: a=',num2str(a),', \theta=',num2str(thlog)]);

%Store method performance
fScoreTable(1,1)=fScore(B,f_ref);
fScoreTable(2,1)=fScore(logP,f_ref);
fScoreTable(3,1)=fScore(L2,f_ref);
sortErrorTable(1,1)=mse(B,f_ref);
sortErrorTable(2,1)=mse(logP,f_ref);
sortErrorTable(3,1)=mse(L2,f_ref);

%Display performance metrics
figure;
set(gcf, 'Position', get(0, 'Screensize'));
bar(fScoreTable);
set(gca,'xTickLabel',{'non-supervised','supervised','supervised+morphological'});
xlabel 'method'
ylabel 'F1-Score'
title('F-Score for various segmentation methods')
grid

figure;
set(gcf, 'Position', get(0, 'Screensize'));
bar(sortErrorTable);
set(gca,'xTickLabel',{'non-supervised','supervised','supervised+morphological'});
xlabel 'method'
ylabel 'mean sorting error'
title('mean sorting error for various segmentation methods')
grid

suptitle('press any key to exit');
pause;

close all