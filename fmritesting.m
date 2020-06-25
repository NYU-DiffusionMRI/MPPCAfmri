close all

cd /mnt/labspace/Projects/RMT-fMRI/PROCESSING_T1_v1/NYU-1/PIANO3_BRAIN_MAPPING_ALT_INDEX_TAP_fMRI_SMS_s9.nii.noproc.feat
addpath(genpath('/mnt/labspace/Projects/RMT-fMRI/scripts'));

mask = niftiread('mask.nii');
raw = niftiread('filtered_func_data.nii');
info = niftiinfo('filtered_func_data.nii');
design = importdata('tsplot/tsplot_zstat1.txt'); 
tstat_fsl = niftiread('stats/tstat1.nii');
zstat_fsl = niftiread('stats/zstat1.nii');

maxpts = 500;
alpha = 0.01;
[Nx, Ny, Nz, Nt] = size(raw);
imtgv = zeros(size(raw));
for i = 1:Nz
    for j = 1:Nt
        input = squeeze(raw(:,:,i,j));
        [utgv2 etgv2] = tgv2_l2_2D_pd(input, input, 2.0*alpha, alpha, maxits);
        imtgv(:,:,i,j) = utgv2;
    end
end
niftiwrite(imtgv, '/mnt/labspace/Projects/RMT-fMRI/PROCESSING_T1_v1/NYU-1/PIANO3_BRAIN_MAPPING_ALT_INDEX_TAP_fMRI_SMS_s9_tgv.mii', info);



%% compute t
square = zeros(1,160);
square([21:40,61:80,101:120,141:160]) = 1;

X1 = raw(:,:,:,square == 1);
X2 = raw(:,:,:,square == 0);

n1 = size(X1,4);
n2 = size(X2,4);

sp2 = (sum((X1-mean(X1,4)).^2, 4) + sum((X2-mean(X2,4)).^2, 4)) ./ (n1 + n2 - 2);
sx = sqrt(sp2/n1 + sp2/n2);

tstat = (mean(X1,4) - mean(X2,4)) ./ sx;
tstat(isnan(tstat)) = 0;

figure; imagesc(tstat(:,:,20))
figure; imagesc(tstat_fsl(:,:,20))
figure; scatter(tstat(:),tstat_fsl(:));

%% lets try to reproduce z
design = design(:,1);
figure; plot(design);

design_ = permute(repmat(design, [1, size(raw,1), size(raw,2), size(raw,3)]), [2 3 4 1]);
r = sum((raw-mean(raw,4)).*(design_-mean(design_,4)), 4)./ ...
    sqrt(sum((raw-mean(raw,4)).^2, 4).*sum((design_-mean(design_,4)).^2, 4));

zstat = atanh(r)*sqrt(size(raw,4) - 3); % z
zstat(isnan(zstat)) = 0;

% %% z using eqs from "software implimentation" I think its wrong ...
% b = sum((raw-mean(raw,4)).*(design_-mean(design_,4)), 4)./ ...
%     sum((design_-mean(design_,4)).^2, 4);
% zstat3 = atanh(b); zstat3 = real(zstat3);
% zstat3(isnan(zstat3)) = 0;

close all
figure; imagesc(zstat(:,:,20));
figure; imagesc(zstat_fsl(:,:,20));
figure; scatter(zstat(:),zstat_fsl(:));
