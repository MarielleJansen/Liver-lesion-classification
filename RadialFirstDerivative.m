% RadialFirstDerivative is a radial gradient histogram function.
% It calculates the gradient between the center pixel of the mask and all
% the other pixels. It find the maximum gradient on a radial line between
% the center pixel and the outline pixel. The distance at this maximum is
% used for further calculations.
% The values are normalized between 0 and 1, with 1 the maximum distance.
% Features derived from this calculation; are the mean, standard deviation,
% kurtosis, skewness and precentile 10 and 90 radial gradients.
%

function Features = RadialFirstDerivative(image, mask)
% dilate mask with a 3x3 cross strel.
SE = [0,1,0;1,1,1;0,1,0];
mask_d = imdilate(mask,SE);
outlineMask = mask_d-mask;

idx = find(outlineMask==1);
[x,y,z] = ind2sub(size(mask),idx);
clear idx
slice = round(median(z));
clear z

xmin = min(x);
xmax = max(x);

ymin = min(y);
ymax = max(y);

ROI = image(xmin-2:xmax+2, ymin-2:ymax+2, slice);
ROIoutline = outlineMask(xmin-2:xmax+2, ymin-2:ymax+2, slice);
centerx = round((xmin+xmax)/2+.5)-xmin+2;
centery = round((ymin+ymax)/2+.5)-ymin+2;

clear x y xmin xmax ymin ymax

% Compute first derivative wrt center pixel.
centerI = ROI(centerx,centery);
Grad = abs(ROI-centerI);

idx = 1:size(ROI,1)*size(ROI,2);
[x,y] = ind2sub(size(ROI),idx);
distx = abs(x-centerx);
disty = abs(y-centery);
Dist = reshape(sqrt(distx.^2 +disty.^2), size(ROI));
clear x y distx disty idx

idx = find(ROIoutline==1);
[x,y] = ind2sub(size(ROIoutline),idx);
clear idx

% Find the distance from center pixel to peak in first derivative.
for i=1:length(x)
    lineGrad = bresenham(Grad, [y(i),x(i); centery, centerx],0);
    lineDist = bresenham(Dist, [y(i),x(i); centery, centerx],0);
    [M,p] = max(lineGrad);
    DistanceMax(i) = lineDist(p);
    clear M p lineGrad lineDist
end


% Find the maximum distance.
T = double(ROIoutline).*Dist;
T(T==0)=[];
maxDist = max(T);
clear T

% Normalize with the maxDistance.
DMN = DistanceMax./maxDist;
clear maxDist

Average = mean(DMN);
SD = std(DMN);
Skewness = skewness(DMN);
Kurtosis = kurtosis(DMN);

L_sort = sort(DMN);
P90 = L_sort(round((length(DMN)*0.90)));
P10 = L_sort(round((length(DMN)*0.10)+1));


Features = [Average,SD,Skewness,Kurtosis,P10,P90];

end