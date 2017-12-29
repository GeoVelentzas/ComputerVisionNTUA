function [ displ_x displ_y ] = displ( d_x , d_y )
%Optical field total displacement extraction
%based on thresholded vector values

thresh = 0.9;

%Calculating energy
E = d_x.^2+d_y.^2;
E = mat2gray(E);

%Calculating mean displacements
displ_x = mean(mean(d_x((E>thresh))));
displ_y = mean(mean(d_y((E>thresh))));

end

