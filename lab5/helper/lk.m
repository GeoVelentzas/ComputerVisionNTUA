function [ d_x , d_y ] = lk(I1,I2,rho,epsilon,d_x0,d_y0)
%Lucas-Kanade algorithm for estimation of the optical flow

%Gaussian for window
n = 2*ceil(3*rho)+1;
Gp = fspecial('gaussian', n, rho);

%partial derivatives
[I1x ,I1y] = gradient(I1);

%initial values of field
dx = d_x0; dy = d_y0;


for i=1:7
    
    %meshgrid for the points of calculation
    [x_0,y_0] = meshgrid(1:size(I1,2),1:size(I1,1));
    %calculate the values of image at those points with interpolation
    I1_i = interp2(I1,x_0+dx,y_0+dy,'linear',0);

    %the same for the x-derivative
    [x_0,y_0] = meshgrid(1:size(I1x,2),1:size(I1x,1));
    A1 = interp2(I1x,x_0+dx,y_0+dy,'linear',0);

    %the same for the y-derivative
    [x_0,y_0] = meshgrid(1:size(I1y,2),1:size(I1y,1));
    A2 = interp2(I1y,x_0+dx,y_0+dy,'linear',0);

    %now we got the partial derivatives we calculate the blocks of the
    %first matrix we need for u
    K11 = imfilter(A1.^2,Gp)+epsilon;
    K12 = imfilter(A1.*A2,Gp);
    K21 = imfilter(A1.*A2,Gp);
    K22 = imfilter(A2.^2,Gp)+epsilon;
    det = 1./(K11.*K22-K12.*K21);
    %and the second matrix
    E = I2-I1_i;
    P1 = imfilter(A1.*E,Gp);
    P2 = imfilter(A2.*E,Gp);
    
    %by looking at the linear algebra relations we now calculate u
    ux = det.*(K22.*P1-K12.*P2);
    uy = det.*(-K21.*P1+K11.*P2);

    %we update the values of d
    dx = dx + ux;
    dy = dy + uy;
    
end

%just to convert to matrix mode from grid mode
d_x = flipud(dx); 
d_y = flipud(dy);


end

