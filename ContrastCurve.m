% Function to calculate features based on the intensity curves in the
% lesion. Features are based on Gilhuijs' and Chen's papers.
% Input: 4D contrast image and lesion mask
% Output: Kinetic curve features in a [1,13] vector
%
% Features:
% - Maximum enhancement
% - Time to peak
% - Uptake rate
% - Washout rate
% - Area under curve
% - Relative maximum enhancement
% - Relative Time to peak
% - Relative Uptake rate
% - Relative Washout rate
% - Variance Relative Maximum enhancement
% - Variance Relative Time to peak
% - Variance Relative Uptake rate
% - Variance Relative Washout rate

function features = ContrastCurve(Image,Mask)
Kinetic = [];
Parenchyma = [];
Image = double(Image);
Mask = double(Mask);

% Extract voxels in ROI defined by the mask
for image = 1:size(Image,4)
    
    K = Image(:,:,:,image).*Mask(:,:,:,1);
    
    K(Mask(:,:,:,1)==0)=[];
    
    Kinetic = [Kinetic;K];
    
    clear K
    
end

% create bounding box around lesion mask for parenchym curve
stats = regionprops(Mask,'BoundingBox');
if nnz(Mask)>0
    BB = ceil(stats.BoundingBox);
    
    BBmask = zeros(size(Mask));
    BBmask(BB(2):(BB(2)+BB(5)-1),BB(1):(BB(1)+BB(4)-1),BB(3):(BB(3)+BB(6)-1))=1;
    BBmask=double(BBmask)-Mask; % inverse of mask
    
    for image = 1:size(Image,4)
        
        P = Image(:,:,:,image).*BBmask;
        
        P(BBmask==0)=[];
        
        Parenchyma = [Parenchyma;P];
        
    end
    
end

KineticCurve = double(Kinetic');
ParenchymaCurve = double(Parenchyma');

volume = size(KineticCurve,1);

if volume <5
    Kinetic=[];
    clear KineticCurve manualVolume Mask Image C F1 F2 F3 F4 F5 CE
    clear CEmean CEvar Image t Mask
    features =[];
    disp( 'Volume too small, cannot calculate curve features.')
    return
end

clear BB BBmask P
%% features 2004 on mean kinetic curve
KineticCurve(KineticCurve==0)=1;

% normalize image by subtracting non-contrast and divide by non-contrast
for t = 2:size(Image,4)
    CE(:,t) = (KineticCurve(:,t)-KineticCurve(:,1))./KineticCurve(:,1);
end
CE(:,1)=0;


if size(CE,1)>1
    CEmean = mean(CE);
else
    CEmean=CE;
end

[F11,F22] = max(CEmean); % first mean, than feature calculation (2004)
% F1 = maximum enhancement, F2 = time to peak
F33 = F11/F22;      % uptake rate

if F22 == size(Image,4)
    F44 = 0;             % washout rate
elseif F22 == size(Image,4)-1
    F44 = 0;
else
    F44 = (F11 - CEmean(size(Image,4)))/(size(Image,4)-F22);
end

F5 = sum(CEmean); % Area under Curve Agliozzo (2012)

% Features on relative curve (Normalization by dividing by mean intensity
% surrounding tissue for each time point)
relativeCurve=mean(KineticCurve)./mean(ParenchymaCurve);
[Fr1,Fr2] = max(relativeCurve);

Fr3 = Fr1/Fr2;

if Fr2 == size(Image,4)
    Fr4 =0;
else
    Fr4 = (Fr1 - relativeCurve(size(Image,4)))/relativeCurve(size(Image,4)-Fr2);
end

% variance
relativeCE = bsxfun(@rdivide, KineticCurve,mean(ParenchymaCurve));
if size(relativeCE,1)>1
    CEvar = var(double(relativeCE));
else
    CEvar(1:16)=0;
end

[Fv1,Fv2] = max(CEvar);
Fv3 = Fv1/Fv2;

if Fv2 == size(Image,4)
    Fv4 =0;
else
    Fv4 = (Fv1 - CEvar(size(Image,4)))/CEvar(size(Image,4)-Fv2);
end

features(1,:) = [F11,F22,F33,F44,F5,Fr1,Fr2,Fr3,Fr4,Fv1,Fv2,Fv3,Fv4];


end