function [ dxf , dyf ] = LucasKanade(I1,I2,p,m,dx0,dy0)
%Lucas-Kanade algorithm for estimation of the optical flow

Gp = fspecial('gaussian', 3, p);
[I1x ,I1y] = gradient(I1);
dx = dx0;
dy = dy0;


for i=1:7
    
    [x_0,y_0] = meshgrid(1:size(I1,2),1:size(I1,1));
    I1_i = interp2(I1,x_0+dx,y_0+dy,'linear',0);

    [x_0,y_0] = meshgrid(1:size(I1x,2),1:size(I1x,1));
    A1 = interp2(I1x,x_0+dx,y_0+dy,'linear',0);

    [x_0,y_0] = meshgrid(1:size(I1y,2),1:size(I1y,1));
    A2 = interp2(I1y,x_0+dx,y_0+dy,'linear',0);

    K11 = imfilter(A1.^2,Gp)+m;
    K12 = imfilter(A1.*A2,Gp);
    K21 = imfilter(A1.*A2,Gp);
    K22 = imfilter(A2.^2,Gp)+m;
    det = 1./(K11.*K22-K12.*K21);

    E = I2-I1_i;
    P1 = imfilter(A1.*E,Gp);
    P2 = imfilter(A2.*E,Gp);

    ux = det.*(K22.*P1-K12.*P2);
    uy = det.*(-K21.*P1+K11.*P2);

    dx = dx + ux;
    dy = dy + uy;
    
end

dxf = flipud(dx); 
dyf = flipud(dy);


end

