% Function to calculate the basic gray level statistics of the lesion on
% the non-contrast image.
% Input: non-contrast 3D image and the mask of the lesion
% Output: mean, SD, skewness, kurtosis and percentile 1, 10, 90 and 99
% (in this order) of the gray levels in the 
% given image after normalization.
%
% Normalization is done by mapping the intensities between 0 and 1. Where
% the 0th percentile is mapped to 0 and the 99th percentile is mapped to 1.
% The original shape of the histogram is kept the same.

function features = GrayLevel(Image, Mask)
Y=Image(:,:,:,1);

Lesion = double(Y).*double(Mask);
Lesion(Mask==0)=[];

% % set values for linear mapping
% s1=0;   % minimum value mapping
% s2=1;   % maximum value mapping
% pc2=0.998;   % maximum landmark
% X=double(reshape(Y,1,(size(Y,1)*size(Y,2)*size(Y,3)))); % reshape matrix to vector
% X_sort=sort(X); % sort X 
% 
% p1=0;       % corresponds to s1
% p2=X_sort(round((length(X)*pc2)+1));    % find pc2 percentile
% 
% clear X X_sort Y 
% 
% % mapping
% for i=1:length(Lesion);
% 
%     L(i) = s1 + Lesion(i)*(s2-s1)/(p2-p1);
%     
% end
L = Lesion;
Average = mean(L);
%Median = median(L);
SD = std(L);
Skewness = skewness(L);
Kurtosis = kurtosis(L);

L_sort = sort(L);
P90 = L_sort(round((length(L)*0.90)));
% P99 = L_sort(round((length(L)*0.99)));
P10 = L_sort(round((length(L)*0.10)+1));
% P1 = L_sort(round((length(L)*0.1)+1));

features = [Average,SD,Skewness,Kurtosis,P10,P90];
end

