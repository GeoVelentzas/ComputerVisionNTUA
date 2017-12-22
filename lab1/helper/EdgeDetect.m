function  D  = EdgeDetect( Im,sig,theta_edge,LaplacType )
% functin which finds the edges of the images Im after gaussian filtering
% with standard deviation sig, and threshold parameter theta_edge.
% Using LaplacType=0 results in linear apprach while a non-zero value
% results in using a non-linear approach

hsize = ceil(3*sig)*2+1;                  %computation of kernel length
Gs = fspecial('gaussian', hsize, sig);    %impulse response of gaussian
h = fspecial('log', hsize, sig);          %Laplacian of Gaussian
Is = imfilter(Im,Gs);                     %filtered with Gaussian kernel
SE = strel('diamond', 1);                 %morhological operator

if LaplacType==0
    L1 = imfilter(Im,h);                  %linear filtering
else
    Idil = imdilate(Is,SE);               %non-linear filtering
    Iero = imerode(Is,SE);
    L1 = Idil + Iero - 2*Is;
end

X=(L1>=0);                                %binary image
Xdil = imdilate(X,SE);                    %dilation
Xero = imerode(X,SE);                     %erosion
Y = Xdil-Xero;                            %find the perigram

[Gx,Gy]=gradient(Is);                     %gradient computation
G=Gx+Gy*1i;                               %for computing norm...
max_grad=max(abs(G(:)));                  %maximum gradient
Z = (abs(G)) > (theta_edge*max_grad);     %condition
D=(Z&Y);                                  %computation of edge matrix

end



