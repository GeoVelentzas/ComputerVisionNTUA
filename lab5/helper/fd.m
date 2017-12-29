function boundingBox = fd(I,mu,cov)
%Extract bounding box information for speaker face
%For rgb image I and given mu,S of skin colour

%Skin propability density threshold
thresh = 0.15;

%Structuring elements
ropen = 3;
rclose = 15;

%Convert image format
Iycbcr = rgb2ycbcr(I);
Icb = Iycbcr(:,:,2);
Icr = Iycbcr(:,:,3);

Icb_vector = reshape(Icb,[],1);
Icr_vector = reshape(Icr,[],1);

%Calculate propability density for all pixels
P = mvnpdf([Icb_vector(:) Icr_vector(:)],mu,cov);
P = reshape(P,size(Icb,1),size(Icb,2));
P = mat2gray(P);
P = P > thresh;

%Perform morphological filtering
se1 = strel('disk',ropen,0);
se2 = strel('disk',rclose,0);

P = imopen(P,se1);
P = imclose(P,se2);

%Locate skin contiguous areas
[L, num]  = bwlabel(P,8);
STATS = regionprops(L,'Area','BoundingBox');

A = zeros(1,num);
for i = 1:num
    A(i) = STATS(i).Area;
end

%Finding larger area
ind = (A == max(A));

boundingBox = STATS(ind).BoundingBox;

end

