% Function to calculate the shape features of a 2D lesion annotation.
% Input: file with the manual annotation
% Output: 8 shape features, based on Kenneth Gilhuijs' and Keelin Murphy's
% papers.
%
% The features vector will be of size [1,8] with these features in this
% order:
% - circularity
% - irregularity
% - compactness1
% - compactness2
% - ratio max_dim:min_dim
% - Sphericity (sort of DSC/area sphere)
% - ratio sphericity:r
% - volume
%
% r = radius = (dimX+dimY)/4

function features = Shape(manualVolume)
Mask = nii_read_volume(manualVolume);

volume = nnz(Mask);

[i,s] = max(max(max(Mask)));
Maskslice = Mask(:,:,s);

effectiveDiameter = round(2*sqrt(volume/(pi)));

circleVolume = pi * (effectiveDiameter/2)^2;

circularity = circleVolume/volume;

SE1 = [1 1 1];
SE2 = [1;1;1];
surfaceVolume1 = Maskslice - imerode(Maskslice,SE1);
surfaceVolume2 = Maskslice - imerode(Maskslice,SE2);

irregularity = 1 - ((pi * (effectiveDiameter))/(nnz(surfaceVolume1)+nnz(surfaceVolume2)));

stats = regionprops(Maskslice,'BoundingBox');
BB = ceil(stats.BoundingBox);
dimX = BB(4);
dimY = BB(3);

compactness1 = volume/(dimX*dimY);
compactness2 = volume/((max(dimX,dimY))^2);
ratioDim = max(dimX,dimY)/min(dimX,dimY);

radius = (dimX + dimY)/4;

% Create a logical image of a circle with specified
% diameter, center, and image size.
% First create the image.
imageSizeX = size(Maskslice,1);
imageSizeY = size(Maskslice,2);
[columnsInImage rowsInImage] = meshgrid(1:imageSizeX, 1:imageSizeY);
% Next create the circle in the image.
centerY = BB(2)+round(0.5*BB(4));
centerX = BB(1)+round(0.5*BB(3));

circlePixels = (rowsInImage - centerY).^2 ...
    + (columnsInImage - centerX).^2 <= radius.^2;

clear centerY centerX imageSizeX imageSizeY columnsInImage rowsInImage
m= int16(circlePixels).*Maskslice;

sphericity = nnz(m)/(pi*radius^2);
ratioSphericity = sphericity/radius;

clear m stats

features = [circularity,irregularity,compactness1,compactness2, ratioDim, ...
    sphericity, ratioSphericity,volume];
end