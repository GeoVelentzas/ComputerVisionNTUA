function  imgNcorns = markCorn(img,corners,mcolor,mradius)
% Put markers at the image 'img', at the corner points, defined from the
% binary image 'corners'.
% imgNcorns: the image 'img' with the markers
% mcolor: the color used for the markers (a 1-by-3 array describing the
% color in the RGB format) 
% mradius: the radius of the corner markers (which are disks)
%
% Author: Anastasios Roussos, 28 Mar. 2008.

corners = imdilate(corners,strel('disk',mradius));

if size(img,3)==1 % grayscale image
    %create a "color" image, but with gray colors only: 
    img_col = repmat(img,[1,1,3]);
else
    img_col = img;
end

imgNcorns = zeros(size(img_col));

for i=1:3
    %i=1: RED, i=2: GREEN, i=3: BLUE
    imgNcorns_chan = img_col(:,:,i); imgNcorns_chan(corners) = mcolor(i);
    imgNcorns(:,:,i) = imgNcorns_chan;
end
