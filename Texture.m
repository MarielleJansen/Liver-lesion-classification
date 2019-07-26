% Function to calculate texture features of dynamic contrast enhanced MR
% images in the liver.
% Input: 4D image (Image) and 3D lesion mask (Mask) and TimeToPeak
% Output: feature vector of size [1,26]
%
% Features are:
%     Angular Second Moment or Energy
%     Contrast
%     Correlation
%     Sum of Squares Variance
%     IDM/Homogeneity
%     Sum Average
%     Sum Variance
%     Sum Entropy
%     Entropy
%     Difference Variance
%     Difference Entropy
%     Information Measure Of Correlation 1
%     Information Measure Of Correlation 2
% calculated on the time to peak image.
%
% Also the variance of these features over the time points are given as
% features.

function features = Texture(Image,Mask,TTP)
Image = double(Image);
Mask = double(Mask);
volume = nnz(Mask);

% Volume of mask should be large enough to calculate texture features on
if volume <5
    Mask = imdilate(Mask,[0,1,0;1,1,1;0,1,0]);
end
clear volume

% % Calculate the time to peak
% Kinetic =[];
% for image = 1:size(Image,4)
%     
%     K = Image(:,:,:,image).*Mask(:,:,:,1);
%     
%     K(Mask(:,:,:,1)==0)=[];
%     
%     Kinetic = [Kinetic;K];
%     
%     clear K
%     
% end
% 
% KineticCurve = double(Kinetic');
% KineticCurve(KineticCurve==0)=1;
% 
% 
% for t = 2:size(Image,4)
%     CE(:,t) = (KineticCurve(:,t)-KineticCurve(:,1))./KineticCurve(:,1);
% end
% 
% [c,TTP] = max(mean(CE));
% 
% clear Kinetic KineticCurve c image

% create bounding box around lesion mask for texture analysis
stats = regionprops(Mask,'BoundingBox');
for j =1:size(Image,4)
    if nnz(Mask)>0
        BB = ceil(stats.BoundingBox);
        
        
        P = Image(BB(2):(BB(2)+BB(5)),BB(1):(BB(1)+BB(4)),BB(3):(BB(3)+BB(6)-1),j);
        
        % due to deformations the boundingbox can
        if size(P,3)>1
            Z = Mask(BB(2):(BB(2)+BB(5)-1),BB(1):(BB(1)+BB(4)-1),BB(3):(BB(3)+BB(6)-1));
            for l = 1:size(Z,3)
                v(1,l)=nnz(Z(:,:,l));
            end
            [z, w] = max(v);
            P = P(:,:,w);
        end
        clear z w v l Z
        
        offsets = [0 1;-1 1;-1 0;-1 -1];
        
        [GLCMS,SI] = graycomatrix(P,'Offset',offsets,'NumLevels',32,'Symmetric',true,'GrayLimits',[]);
        
    end
    
    
    GLCM = sum(GLCMS,3);
    
    clear BB P CE
    
    Features = GLCMFeatures(GLCM); % invariant features for each time point
    
    AMSi(j,1) = Features.energy;% Angular Second Moment or Energy
    Contrasti(j,1) = Features.contrast;
    Correlationi(j,1) = Features.correlation;
    SSVariancei(j,1) = Features.sumOfSquaresVariance; % Sum of Squares Variance
    IDMi(j,1) = Features.homogeneity; % Inverse Difference moment or homogeneity
    SumAveragei(j,1) = Features.sumAverage;
    SumVariancei(j,1) = Features.sumVariance;
    SumEntropyi(j,1) = Features.sumEntropy;
    Entropyi(j,1) = Features.entropy;
    DiffVariancei(j,1) = Features.differenceVariance;
    DiffEntropyi(j,1) = Features.differenceEntropy;
    IMC1i(j,1) = Features.informationMeasureOfCorrelation1;
    IMC2i(j,1) = Features.informationMeasureOfCorrelation2;
    % Maximum Correlation Coefficient not found.
    
    clear Features GLCM GLCMS
end

% Features on Time To Peak
AMS = AMSi(TTP,1);
Contrast = Contrasti(TTP,1);
Correlation = Correlationi(TTP,1);
SSVariance = SSVariancei(TTP,1);
IDM = IDMi(TTP,1);
SumAverage = SumAveragei(TTP,1);
SumVariance = SumVariancei(TTP,1);
SumEntropy = SumEntropyi(TTP,1);
Entropy = Entropyi(TTP,1);
DiffVariance = DiffVariancei(TTP,1);
DiffEntropy = DiffEntropyi(TTP,1);
IMC1 = IMC1i(TTP,1);
IMC2 = IMC2i(TTP,1);

% Variance of features in images over time
vAMS = var(AMSi);
vContrast = var(Contrasti);
vCorrelation = var(Correlationi);
vSSVariance = var(SSVariancei);
vIDM = var(IDMi);
vSumAverage = var(SumAveragei);
vSumVariance = var(SumVariancei);
vSumEntropy = var(SumEntropyi);
vEntropy = var(Entropyi);
vDiffVariance = var(DiffVariancei);
vDiffEntropy = var(DiffEntropyi);
vIMC1 = var(IMC1i);
vIMC2 = var(IMC2i);

% features = [AMS,Contrast,Correlation,SSVariance,IDM,SumAverage,SumVariance ...
%      ,SumEntropy,Entropy,DiffVariance,DiffEntropy,IMC1,IMC2, ...
%      vAMS,vContrast,vCorrelation,vSSVariance,vIDM,vSumAverage,vSumVariance ...
%      ,vSumEntropy,vEntropy,vDiffVariance,vDiffEntropy,vIMC1,vIMC2];
features = [AMS,Contrast,Correlation,SSVariance,IDM,SumAverage,SumVariance ...
     ,SumEntropy,Entropy,DiffVariance,DiffEntropy,IMC1,IMC2];

end