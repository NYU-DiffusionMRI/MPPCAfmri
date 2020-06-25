raw = niftiread(fullfile(root,[ser,'.nii.noproc.feat'],'stats','tstat1_acpc.nii.gz')); 
sm = niftiread(fullfile(root,[ser,'.feat'],'stats','tstat1_acpc.nii.gz')); 
ica = niftiread(fullfile(root,[ser,'_ica.nii.noprocica2.feat'],'stats','tstat1_acpc.nii.gz'));
mp = niftiread(fullfile(root,[ser,'_ica.nii.noprocdn2.feat'],'stats','tstat1_acpc.nii.gz')); 
cnn = niftiread(fullfile(root,[ser,'_ica.nii.noprocdncnn3.feat'],'stats','tstat1_acpc.nii.gz')); 
%raw(mask==0) = 0; sm(mask==0) = 0; ica(mask==0) = 0; mp(mask==0) = 0; cnn(mask==0) = 0;
% raw = abs(raw);
% sm = abs(sm);
% ica = abs(ica);
% mp = abs(mp);
% cnn = abs(cnn);

t1f = niftiread('/Volumes/Research/fieree01lab/labspace/Projects/RMT-fMRI/PROCESSING_T1_v1/NYU-24/BRAIN_MAPPING_RIGHT_HAND_PIANO_s53.feat/reg/highres_acpc.nii');
t1f(t1f<0) = 0;
C = log([5 20]);
vox = [120, 90, 171];

a = imrotate(squeeze(t1f(60:206,vox(2),30:190)),-90);
f1 = log(imrotate(squeeze(raw(60:206,vox(2),30:190)),-90));
f2 = log(imrotate(squeeze(sm(60:206,vox(2),30:190)),-90));
f3 = log(imrotate(squeeze(ica(60:206,vox(2),30:190)),-90));
f4 = log(imrotate(squeeze(cnn(60:206,vox(2),30:190)),-90));
f5 = log(imrotate(squeeze(mp(60:206,vox(2),30:190)),-90));
acat1 = cat(2,a,a,a,a,a);
fcat1 = cat(2,f1,f2,f3,f4,f5);

a = squeeze(t1f(vox(3),60:210,40:186));
f1 = log(squeeze(raw(vox(3),60:210,40:186)));
f2 = log(squeeze(sm(vox(3),60:210,40:186)));
f3 = log(squeeze(ica(vox(3),60:210,40:186)));
f4 = log(squeeze(cnn(vox(3),60:210,40:186)));
f5 = log(squeeze(mp(vox(3),60:210,40:186)));
acat2 = cat(2,a,a,a,a,a);
fcat2 = cat(2,f1,f2,f3,f4,f5);

a = imrotate(squeeze(t1f(60:206,60:210,vox(1))),-90);
f1 = log(imrotate(squeeze(raw(60:206,60:210,vox(1))),-90));
f2 = log(imrotate(squeeze(sm(60:206,60:210,vox(1))),-90));
f3 = log(imrotate(squeeze(ica(60:206,60:210,vox(1))),-90));
f4 = log(imrotate(squeeze(cnn(60:206,60:210,vox(1))),-90));
f5 = log(imrotate(squeeze(mp(60:206,60:210,vox(1))),-90));
acat3 = cat(2,a,a,a,a,a);
fcat3 = cat(2,f1,f2,f3,f4,f5);

acat = cat(1,acat1,acat2,acat3);
fcat = real(cat(1,fcat1,fcat2,fcat3));

figure('color','k')
a1 = axes;
imagesc(a1, acat); colormap(a1, gray); axis image
a2 = axes;
c = imagesc(a2, real(fcat)); axis image; colormap(a2,'hot'); 
caxis(C);
falpha = fcat; falpha(falpha<log(5)) = 0; falpha(falpha >= log(5)) = 1;
alpha(c, falpha); a1.Visible = 'off'; a2.Visible = 'off';

title(a1, 'Raw   Smoothed   ICA   DNCNN   MPPCA','color','w','interpreter','latex','fontsize',20); 

