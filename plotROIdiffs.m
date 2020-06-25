addpath(genpath('~/Documents/MATLAB'));
clear 
close all

root = '/Volumes/Research/fieree01lab/labspace/Projects/RMT-fMRI/testdata/14408200/nii';
data = '20171211_162444BRAINMAPPINGfMRISMSVERBGENs023a001';
%data = '20171211_162444BRAINMAPPINGfMRISMSALTFINGERTAPPINGs013a001';
%roinames = {'brocas','left preCG','left SMA','right BA homolog','right preCG','right WA homolog','wernickes'};
roinames = {'BA-L','BA-R'};
rawf = fullfile(root,[data, '.feat/stats/zstat1.nii']);
dnf = fullfile(root,[data,'_DN.feat/stats/zstat1.nii']);
mask = fullfile(root,[data,'_DN.feat/mask.nii']);
gunzip([mask,'.gz']);
gunzip([rawf,'.gz']);
gunzip([dnf,'.gz']);

nii = load_untouch_nii('/Volumes/labspace/Projects/RMT-fMRI/testdata/14408200/rois/compositeROIs-revised.nii'); roi = nii.img;
nii = load_untouch_nii(mask); mask = logical(nii.img);
nii = load_untouch_nii(rawf); rawz = nii.img;
nii = load_untouch_nii(dnf); dnz = nii.img;

se = strel('sphere',9);
Zthr = -6;

%figure('color','w','Position',[0 0 2000 200])
for i = 1:max(roi(:))
    roi_ = roi;
    roi_(roi_~=i) = 0;
    roi_ = logical(roi_);

    %roi_ = imdilate(roi_,se);
    roi_(rawz==0) = 0;
    dnroi_ = roi_;
    rwroi_ = roi_;
    dnroi_(dnz>Zthr) = 0;
    rwroi_(rawz>Zthr) = 0;
    roi_ = logical(rwroi_+dnroi_);
    [x,y,z] = ind2sub(size(roi_),find(roi_));
    
    box_ = regionprops(roi_, 'BoundingBox');
    c1 = floor(box_(1).BoundingBox);
    mask = zeros(size(roi_));
    mask(c1(2):c1(2)+c1(5),c1(1):c1(1)+c1(4),c1(3):c1(3)+c1(6)) = 1;
    z = floor(mean([c1(3),c1(3)+c1(6)]));
    
    rawz_ = rawz(roi_==1);
    dnz_ = dnz(roi_==1);
    
    [Nr, Edgesr] = histcounts(rawz_,25);
    [Nd, Edgesd] = histcounts(dnz_,25);
    Nres = -(Nd - Nr);
        
%     subplot(1,max(roi(:)),i,'align')
%     bar(Edgesr(2:end),Nr,'FaceAlpha',0.5,'FaceColor','b'); 
%     hold on
%     bar(Edgesr(2:end),Nd,'FaceAlpha',0.5,'FaceColor','r');
%     %bar(Edgesr(2:end),Nres,'FaceAlpha',0.5,'FaceColor','y');
%     hold off
%     legend('raw','denoised','difference')
%     title(roinames{i})
%     xlabel('Z statistic')
%     %ylim([0 1])
%     % xlim([-11 3]);   
%    if i == 1, ylabel('VerbGen voxels'); end
    
    x = round(mean(x)); y = round(mean(y)); z = round(mean(z));
    figure; imagesc(squeeze(rawz(:,:,z)+10.*mask(:,:,z)))
    
    clear dnz_ rawz_ 
end


    