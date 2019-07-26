xlsFile = 'C:\Users\user\Documents\Characterization\LesionList.xlsx';
sheet = 'LesionList';

ATFile  = 'C:\Users\user\Documents\Characterization\AcquisitionTimes.xlsx';
sheetAT = 'Sheet 1';
rangeAT = 'B2:P134';
AT = xlsread(ATFile,sheetAT,rangeAT);
sAT = zeros([size(AT,1),size(AT,2)+1]);
for t = 1:size(AT,1)
    sAT(t,:) = [-10,AT(t,:)];
end
Subject = xlsread(xlsFile, sheet, 'B2:B226');

Ftexture=[];
j=1;
i=1;

% Calculate features
for sub = 1:133
    if exist(['C:\Users\user\Documents\Characterization\DataCharacterization\' num2str(sub) '\e-THRIVE_reg.nii'],'file') == 0
        continue
    end
    Y = nii_read_volume(['C:\Users\user\Documents\Characterization\DataCharacterization\' num2str(sub) '\e-THRIVE_reg.nii']); %_reg
    
    % If not 16 time points, continue
    tp = size(Y,4);
    if tp <16
        continue
    end
    
    if exist(['C:\Users\user\Documents\Characterization\DataCharacterization\' num2str(sub) '\parenchyma.nii'],'file') == 0
        continue
    end
    
    % set values for linear mapping
    s1=0;   % minimum value mapping
    s2=1;   % maximum value mapping
    pc2=0.998;   % maximum landmark
    X=double(reshape(Y,1,(size(Y,1)*size(Y,2)*size(Y,3)*size(Y,4)))); % reshape matrix to vector
    X_sort=sort(X); % sort X
    
    p1=0;       % corresponds to s1
    p2=X_sort(round((length(X)*pc2)+1));    % find pc2 percentile
    
    
    % mapping
    X_mapped = zeros(size(X));
    for b=1:length(X)
        
        X_mapped(b) = s1 + X(b)*(s2-s1)/(p2-p1);
        
    end
    DCEMRI = reshape(X_mapped,size(Y));
    clear X_sort Y X X_mapped
    
    % Get smoothened parenchyma contrast curves
    ROI = nii_read_volume(['C:\Users\user\Documents\Characterization\DataCharacterization\' num2str(sub) '\parenchyma.nii']);
    
    tipsROI = uint8(imdilate(ROI, ones([5,5,5])));
    tipsDCEMRI = permute(DCEMRI, [4,1,2,3]);
     
    sDCEMRI = permute(tipsFilter(single(tipsDCEMRI), [5,5,3] , 20, uint8(tipsROI)), [2,3,4,1]);
     
    clear tipsROI tipsDCEMRI
     
    % calculate mean parenchyma
    for tp = 1:size(sDCEMRI,4)
         parenchyma = double(sDCEMRI(:,:,:,tp)).*double(ROI);
         parenchyma(ROI==0)=[];
         TIC_p(tp) = mean(parenchyma);

         clear parenchyma
         end
     S_parenchyma = TIC_p;

    
    
    %% read lesions per subject
    
    study = xlsread(xlsFile,sheet,['B' num2str(i+1)]);
    if sub ~= study
        l = min(find(Subject==sub));
        if isempty(l) == 1
            continue
        else
            i = l;
        end
        clear l
        study = xlsread(xlsFile,sheet,['B' num2str(i+1)]);
    end

        
    while sub == study
        
        [num,L] = xlsread(xlsFile,sheet,['A' num2str(i+1)]);
        lesion = char(L);
        study = xlsread(xlsFile,sheet,['B' num2str(i+1)]);
        
        lesiontype(j,1) = xlsread(xlsFile,sheet,['C' num2str(i+1)]);
        
        manualVolume = ['C:\Users\user\Documents\Characterization\DataCharacterization\LesionAnnotations\' lesion '.nii']; %without unreg
        
        Mask = nii_read_volume(manualVolume);
        
        if nnz(Mask)<8
            i = i+1;
            study = xlsread(xlsFile,sheet,['B' num2str(i+1)]);
            continue
        end
        
        idx = find(Mask==1);
        [x,y,z] = ind2sub(size(Mask),idx);
        clear x y idx
        slice = round(median(z));
        
      
        DCEMRI(DCEMRI<0.001)=0.001;
        
        % Get contrast curves
        SI_ratio = zeros(size(DCEMRI));
        Sratio_mask=[];
        for s = 1:size(DCEMRI,4)
            SI_ratio(:,:,:,s) = DCEMRI(:,:,:,s)./S_parenchyma(1,s);

            % Get curve only within mask
            K = SI_ratio(:,:,:,s).*double(Mask);
            K(Mask==0)=[];
            Sratio_mask = [Sratio_mask;K];

            clear K 
        end

        
        subjectLesion(j,1) = sub;
        
        Sratio(j,:) = mean(Sratio_mask,2);
        
        clear Kinetic
        
        MaxEnhancement = max(SI_ratio,[],4);
        
        
        % Gray level features on all 6 phases 
        F = GrayLevel(DCEMRI(:,:,:,1),Mask);
        t0GrayLevelHistogram(j,:) = F;
        F = GrayLevel(DCEMRI(:,:,:,2),Mask);
        t1GrayLevelHistogram(j,:) = F;
        F = GrayLevel(DCEMRI(:,:,:,7),Mask);
        t2GrayLevelHistogram(j,:) = F;
        F = GrayLevel(DCEMRI(:,:,:,11),Mask);
        t3GrayLevelHistogram(j,:) = F;
        F = GrayLevel(DCEMRI(:,:,:,13),Mask);
        t4GrayLevelHistogram(j,:) = F;
        F = GrayLevel(DCEMRI(:,:,:,16),Mask);
        t5GrayLevelHistogram(j,:) = F;

        Im = mean(DCEMRI(:,:,:,7:10),4); % Mean late arterial enhancement
        
        F = GrayLevel(Im,Mask);
        LEGrayLevelHistogram(j,:) = F;
        
        clear F
        
        % Texture feature on Late Arterial Enhancement
        Ftexture_LE(j,:) = Texture(Im,Mask,1);
        % Texture feature on T0
        Ftexture_t0(j,:) = Texture(DCEMRI,Mask,1); % Texture features
        
        % Texture feature on T1
        Ftexture_t1(j,:) = Texture(DCEMRI,Mask,2); % Texture features
        % Texture feature on T2
        Ftexture_t2(j,:) = Texture(DCEMRI,Mask,7); % Texture features
        % Texture feature on T3
        Ftexture_t3(j,:) = Texture(DCEMRI,Mask,11); % Texture features
        % Texture feature on T4
        Ftexture_t4(j,:) = Texture(DCEMRI,Mask,13); % Texture features
        % Texture feature on T5
        Ftexture_t5(j,:) = Texture(DCEMRI,Mask,16); % Texture features
        % Texture feature on Max enhancement
        Ftexture_MaxEnh(j,:) = Texture_3D(MaxEnhancement,Mask); % Texture features

        % Radial gradient histogram 
        
        FeaturesRGH(j,:) = RGH(Im, Mask);
        
        clear F1 F2 F3 F4 CErelativeMask CEmask diff mask
        clear L z slice num Mask ROI
        clear TIC_p ttp
        
        clear L lesion manualVolume
        i = i+1;
        study = xlsread(xlsFile,sheet,['B' num2str(i+1)]);
        j = j+1;
        
        clear CE_ratio CE mask Mask dcemri
        
    end
    
    clear DCEMRI SD_p TIC_p p1 p2 pc2 ROI tp z b l L num s1 s2
    clear Kinetic SD_p
    
end



