clear
close all
clc
clear
addpath(genpath('~/Documents/MATLAB'))

root = '/Volumes/Research/fieree01lab/labspace/Projects/RMT-fMRI/PROCESSING_T1_v1/NYU-24';
ser = 'BRAIN_MAPPING_RIGHT_HAND_PIANO3_s53';
ser_ = 'BRAIN_MAPPING_RIGHT_HAND_PIANO_s53';

raw = niftiread(fullfile(root,[ser,'.nii.noproc.feat'],'filtered_func_data.nii')); 
sm = niftiread(fullfile(root,[ser_,'.feat'],'filtered_func_data.nii'))/15.270622000364042; 
ica = niftiread(fullfile(root,[ser_,'_ica.nii.noprocica2.feat'],'filtered_func_data.nii'))/15.2905198777;
mp = niftiread(fullfile(root,[ser_,'_ica.nii.noprocdn2.feat'],'filtered_func_data.nii'))/15.321626035; 
cnn = niftiread(fullfile(root,[ser_,'_ica.nii.noprocdncnn3.feat'],'filtered_func_data.nii'))/15.3374233129; 

mask = niftiread(fullfile(root,[ser_,'_ica.nii.noprocica2.feat'],'mask.nii'));
sigma = niftiread(fullfile(root,[ser_,'_noisemap.nii'])); 

for i = 1:160
r = raw(:,:,:,i);
r(mask==0) = 0;
raw(:,:,:,i) = r;
s = squeeze(sm(:,:,:,i));
s(mask==0) = 0;
sm(:,:,:,i) = s;
%sm(mask==0) = 0;
end
vox = [27, 32, 25];
%%
close all
figure('color','k')
%subtightplot(1,5,1)
%imagesc(imrotate(raw(10:60,:,vox(3)),-90)); caxis([0 1500]); colormap(gray); axis off; axis image; title('Raw','color','w','interpreter','latex','fontsize',20); 
subtightplot(1,4,1)
r1 = squeeze(mean(imrotate(sm(10:60,:,vox(3),:),-90)-imrotate(raw(10:60,:,vox(3),:),-90),4));
imagesc(r1); caxis([-100 100]); colormap(gray); colorbar; axis off; axis image; title('Smoothed','color','w','interpreter','latex','fontsize',20);
subtightplot(1,4,2)
r2 = squeeze(mean(imrotate(ica(10:60,:,vox(3),:),-90)-imrotate(raw(10:60,:,vox(3),:),-90),4));
imagesc(r2); caxis([-1.5 1.5]); colormap(gray); colorbar; axis off; axis image; title('ICA','color','w','interpreter','latex','fontsize',20);
subtightplot(1,4,3)
r3 = squeeze(mean(imrotate(cnn(10:60,:,vox(3),:),-90)-imrotate(raw(10:60,:,vox(3),:),-90),4));
imagesc(r3); caxis([-20 20]); colormap(gray); colorbar; axis off; axis image; title('DNCNN','color','w','interpreter','latex','fontsize',20);
subtightplot(1,4,4)
r4 = squeeze(mean(imrotate(mp(10:60,:,vox(3),:),-90)-imrotate(raw(10:60,:,vox(3),:),-90),4));
imagesc(r4); caxis([-.01 .01]); colormap(gray); colorbar; axis off; axis image; title('MPPCA','color','w','interpreter','latex','fontsize',20);
set(gca,'color','w')
%%
smres = (sm - raw)./sigma; smres = smres(mask==1);
icares = (ica - raw)./sigma; icares = icares(mask==1);
mpres = (mp - raw)./sigma; mpres = mpres(mask==1);
cnnres = (cnn - raw)./sigma; cnnres = cnnres(mask==1);

edges = linspace(0,15,30);
centers = edges(1:end-1)/2 + edges(2:end)/2;
Nsm = histcounts(abs(smres),edges,'Normalization','probability'); %Nsm = Nsm/max(Nsm);
Nica = histcounts(abs(icares),edges,'Normalization','probability'); %Nica = Nica/max(Nica);
Nmp = histcounts(abs(mpres),edges,'Normalization','probability'); %Nmp = Nmp/max(Nmp);
Ncnn = histcounts(abs(cnnres),edges,'Normalization','probability'); %Ncnn = Ncnn/max(Ncnn);
close all
figure('color','w')
plot(centers.^2, Nmp, centers.^2, Nsm, centers.^2, Nica, centers.^2, Ncnn,'LineWidth',1); 
hold on;
plot(centers.^2,1/sqrt(2*pi)*exp(-centers.^2/2),'--k','LineWidth',1);
xlim([0 15]); ylim([1e-3 1])
grid on
set(gca,'LineWidth',1)
set(gca,'yscale','log');
axis square
set(gca,'fontsize',20,'fontstyle','times')
xlabel('$r^2$','interpreter','latex','fontsize',20)
ylabel('$\log p(r)$','interpreter','latex','fontsize',20)
legend('MPPCA','Smoothed','ICA','DnCNN','$$1/\sqrt{2\pi}e^{-r^2/2}$$','interpreter','latex','fontsize',20,'Location','SouthWest'); 
title('Normalized Residuals','interpreter','latex','fontsize',20)
legend box off

%%
close all
raw = niftiread(fullfile(root,[ser,'.nii.noproc.feat'],'stats','zstat1.nii')); 
sm = niftiread(fullfile(root,[ser,'.feat'],'stats','zstat1.nii')); 
ica = niftiread(fullfile(root,[ser,'_ica.nii.noprocica2.feat'],'stats','zstat1.nii.gz'));
mp = niftiread(fullfile(root,[ser,'_ica.nii.noprocdn2.feat'],'stats','zstat1.nii.gz')); 
cnn = niftiread(fullfile(root,[ser,'_ica.nii.noprocdncnn3.feat'],'stats','zstat1.nii')); 
raw(mask==0) = 0; sm(mask==0) = 0; ica(mask==0) = 0; mp(mask==0) = 0; cnn(mask==0) = 0;
raw = abs(raw);
sm = abs(sm);
ica = abs(ica);
mp = abs(mp);
cnn = abs(cnn);
vox = [27, 32, 25];

t1 = niftiread('/Volumes/Research/fieree01lab/labspace/Projects/RMT-fMRI/PROCESSING_T1_v1/NYU-24/BRAIN_MAPPING_RIGHT_HAND_PIANO_s53.feat/reg/highres2func.nii');

C = log([5 20]);
a = imgaussfilt(imrotate(t1(15:55,10:60,vox(3)),-90),.2);
f1 = log(imrotate(raw(15:55,10:60,vox(3)),-90));
f2 = log(imrotate(sm(15:55,10:60,vox(3)),-90));
f3 = log(imrotate(ica(15:55,10:60,vox(3)),-90));
f4 = log(imrotate(cnn(15:55,10:60,vox(3)),-90));
f5 = log(imrotate(mp(15:55,10:60,vox(3)),-90));
acat = cat(2,a,a,a,a,a);
fcat = cat(2,f1,f2,f3,f4,f5);

figure('color','k')
a1 = axes;
imagesc(a1, acat); colormap(a1, gray); axis image
a2 = axes;
c = imagesc(a2, fcat); axis image; colormap(a2,'hot'); 
caxis(C);
falpha = fcat; falpha(falpha<log(5)) = 0; falpha(falpha >= log(5)) = 1;
alpha(c, falpha); a1.Visible = 'off'; a2.Visible = 'off';
title('Raw   Smoothed   ICA   DNCNN   MPPCA','color','w','interpreter','latex','fontsize',20); 


imagesc(a1, a); colormap(a1, gray);  axis image
a2 = axes; f = log(imrotate(raw(15:55,10:60,vox(3)),-90));
c2 = imagesc(a2, f); axis image; colormap(a2,'hot'); caxis(C);
falpha = f; falpha(falpha<log(5)) = 0; falpha(falpha >= log(5)) = 1;
alpha(c2, falpha); a1.Visible = 'off'; a2.Visible = 'off';

imagesc(log(imrotate(raw(15:55,10:60,vox(3)),-90))); caxis(C); colormap(jet); axis off; axis image; title('Raw','color','w','interpreter','latex','fontsize',20); 
subtightplot(1,5,2)
%imagesc(log(imrotate(sm(15:55,10:60,vox(3)),-90))); caxis(C); colormap(hot); axis off; axis image; title('Smoothed','color','w','interpreter','latex','fontsize',20);
subtightplot(1,5,3)
a1 = axes; a = imrotate(t1(15:55,10:60,vox(3)),-90);
c1 = imagesc(a1, a); colormap(a1, gray);  axis image
a2 = axes; f = log(imrotate(ica(15:55,10:60,vox(3)),-90));
c2 = imagesc(a2, f); axis image; colormap(a2,'hot'); caxis(C);
falpha = f; falpha(falpha<log(5)) = 0; falpha(falpha >= log(5)) = 1;
alpha(c2, falpha); a1.Visible = 'off'; a2.Visible = 'off';
%imagesc(log(imrotate(ica(15:55,10:60,vox(3)),-90))); caxis(C); colormap(hot); axis off; axis image; title('ICA','color','w','interpreter','latex','fontsize',20);
subtightplot(1,5,4)
a1 = axes; a = imrotate(t1(15:55,10:60,vox(3)),-90);
c1 = imagesc(a1, a); colormap(a1, gray);  axis image
a2 = axes; f = log(imrotate(cnn(15:55,10:60,vox(3)),-90));
c2 = imagesc(a2, f); axis image; colormap(a2,'hot'); caxis(C);
falpha = f; falpha(falpha<log(5)) = 0; falpha(falpha >= log(5)) = 1;
alpha(c2, falpha); a1.Visible = 'off'; a2.Visible = 'off';
%imagesc(log(imrotate(cnn(15:55,10:60,vox(3)),-90))); caxis(C); colormap(hot); axis off; axis image; title('DNCNN','color','w','interpreter','latex','fontsize',20);
subtightplot(1,5,5)
a1 = axes; a = imrotate(t1(15:55,10:60,vox(3)),-90);
c1 = imagesc(a1, a); colormap(a1, gray);  axis image
a2 = axes; f = log(imrotate(mp(15:55,10:60,vox(3)),-90));
c2 = imagesc(a2, f); axis image; colormap(a2,'hot'); caxis(C);
falpha = f; falpha(falpha<log(5)) = 0; falpha(falpha >= log(5)) = 1;
alpha(c2, falpha); a1.Visible = 'off'; a2.Visible = 'off';
%imagesc(log(imrotate(mp(15:55,10:60,vox(3)),-90))); caxis(C); colormap(hot); axis off; axis image; title('MPPCA','color','w','interpreter','latex','fontsize',20);

% subtightplot(3,5,6)
% imagesc(imrotate(squeeze(raw(10:60,vox(2),:)),90)); caxis(C); colormap(hot); axis off; axis image;
% subtightplot(3,5,7)
% imagesc(imrotate(squeeze(sm(10:60,vox(2),:)),90)); caxis(C); colormap(hot); axis off; axis image; 
% subtightplot(3,5,8)
% imagesc(imrotate(squeeze(ica(10:60,vox(2),:)),90)); caxis(C); colormap(hot); axis off; axis image;
% subtightplot(3,5,9)
% imagesc(imrotate(squeeze(cnn(10:60,vox(2),:)),90)); caxis(C); colormap(hot); axis off; axis image; 
% subtightplot(3,5,10)
% imagesc(imrotate(squeeze(mp(10:60,vox(2),:)),90)); caxis(C); colormap(hot); axis off; axis image; 
% 
% subtightplot(3,5,11)
% imagesc(imrotate(squeeze(raw(vox(1),:,:)),90)); caxis(C); colormap(hot); axis off; axis image;
% subtightplot(3,5,12)
% imagesc(imrotate(squeeze(sm(vox(1),:,:)),90)); caxis(C); colormap(hot); axis off; axis image; 
% subtightplot(3,5,13)
% imagesc(imrotate(squeeze(ica(vox(1),:,:)),90)); caxis(C); colormap(hot); axis off; axis image;
% subtightplot(3,5,14)
% imagesc(imrotate(squeeze(cnn(vox(1),:,:)),90)); caxis(C); colormap(hot); axis off; axis image; 
% subtightplot(3,5,15)
% imagesc(imrotate(squeeze(mp(vox(1),:,:)),90)); caxis(C); colormap(hot); axis off; axis image; 

