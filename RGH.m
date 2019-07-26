% RGH is the radial gradient histogram function.
% It calculates the gradient between the center pixel of the mask and all
% the other pixels and multiplies it with the distance from the center
% pixel. The values are then normalized between 0 and 1.
% Features derived from this calculation; are the mean, standard deviation,
% kurtosis, skewness and precentile 10 and 90 radial gradients.
%

function Features = RGH(image, mask)


idx = find(mask==1);
[x,y,z] = ind2sub(size(mask),idx);
clear idx
slice = round(median(z));
clear z

xmin = min(x);
xmax = max(x);

ymin = min(y);
ymax = max(y);

ROI = image(xmin-2:xmax+2, ymin-2:ymax+2, slice);
centerx = round((xmin+xmax)/2+.5)-xmin+2;
centery = round((ymin+ymax)/2+.5)-ymin+2;

clear x y xmin xmax ymin ymax

centerI = ROI(centerx,centery);
Grad = ROI-centerI;

idx = 1:size(ROI,1)*size(ROI,2);
[x,y] = ind2sub(size(ROI),idx);
distx = abs(x-centerx);
disty = abs(y-centery);
Dist = reshape(sqrt(distx.^2 +disty.^2), size(ROI));

RGH = abs(Dist.*double(Grad));
clear Grad Dist
RGHnorm = (RGH - min(min(RGH)))/(max(max(RGH)) - min(min(RGH)));

listRGHnorm = reshape(RGHnorm,size(distx));

Average = mean(listRGHnorm);
SD = std(listRGHnorm);
Skewness = skewness(listRGHnorm);
Kurtosis = kurtosis(listRGHnorm);

L_sort = sort(listRGHnorm);
P90 = L_sort(round((length(listRGHnorm)*0.90)));
P10 = L_sort(round((length(listRGHnorm)*0.10)+1));


Features = [Average,SD,Skewness,Kurtosis,P10,P90];

end