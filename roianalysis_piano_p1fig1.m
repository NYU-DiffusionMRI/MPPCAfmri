addpath(genpath('~/Documents/MATLAB'));
clear 
close all
set(0,'defaultAxesFontSize',20)
set(0,'defaultAxesFontName','Arial')
format short

subjects = [1 2 5 7 8 9 10 11 12 13 14 16 17 18 19 20 21 22 23 24 25 27 28 29 30];
sid = strcat('NYU-',strsplit(num2str(subjects)));

rois = {'CB','CH','CL','CM','CW','IB','IH','IL','IM','IW'};
%rois = {'CH','CM','IB'};
sphere_r = {'.nii','_sphere10_bin.nii','_sphere15_bin.nii','_sphere20_bin.nii'};
spheret = {'single voxel','sphere r=10v','sphere r=15v','sphere r=20v'};

root = '/Volumes/Research/fieree01lab/labspace/Projects/RMT-fMRI';
fdir = fullfile(root,'PROCESSING_T1_v1');
rdir = fullfile(root,'ROI_piano');

series = 'PIANO';
sertit ='Piano Tapping';
comproiraw = zeros(100000,10);
comproidn = zeros(100000,10);

rawdir = dir(fullfile(root,'RAW',sid{20},'BRAIN_MAPPING_RIGHT_HAND_PIANO_s53'));
rawdir(1:2) = [];
dcm0 = dicominfo(fullfile(rawdir(1).folder,rawdir(1).name));
time0 = str2double(dcm0.AcquisitionTime);
% time = [];
% for i = 1:numel(rawdir)
%     dcm = dicominfo(fullfile(rawdir(i).folder,rawdir(i).name));
%     acqt = str2double(dcm.AcquisitionTime);
%     time(i) = acqt - time0;
% end
time = [0, 41:42+157];

%%

dn = niftiread('~/Desktop/fmri_dn.nii');
sm = niftiread('~/Desktop/fmri_sm.nii.gz');
raw = niftiread('~/Desktop/fmri_raw.nii');
ica = niftiread('~/Desktop/fmri_ica.nii.gz');
dncnn = niftiread('~/Desktop/fmri_dncnn.nii');
nii = load_untouch_nii('~/Desktop/testroi.nii'); roi = logical(nii.img);
nii = load_untouch_nii(fullfile(fdir,sid{1},'BRAIN_MAPPING_ALT_INDEX_TAP_fMRI_SMS_s9_noisemap.nii')); sigma = nii.img;

sm = zeros(size(raw));
for v = 1:size(raw,4)
    sm(:,:,:,v) = imgaussfilt3(raw(:,:,:,v), 1.66);
end

%close all
data = cell(1,5);
data{1} = vec(raw,roi);
data{2} = vec(dn,roi);
data{3}= vec(sm,roi);
data{4} = vec(ica,roi);
data{5} = single(vec(dncnn,roi));
labels = {'MP-PCA','Smoothing','ICA','DnCNN'};

col = get(0, 'DefaultAxesColorOrder');

%%
%close all
S=vec(sigma,roi);
sigma_ = mean(vec(sigma,roi),2);
T=75; dt=1; t=-T/2:dt:T/2-dt; 
dw=2*pi/T; w=-pi/dt:dw:pi/dt-dw; 

figure('color','w');
for i = 3:-1:2
    ep = (data{i}(:,:) - data{1}(:,:));
    ep = ep - mean(ep,1);
    
    epf = fft(ep)*dt; 
    Gamma = mean(conj(epf).*epf, 2) / T; 
    Gamma_t = ifft(Gamma)/dt;
    
    subplot(1,2,1); 
    plot(w/2/pi, fftshift(Gamma)/(sigma_.^2),'linewidth',1.5,'color', col(i-1,:)); hold on
    set(gca,'LineWidth',1.5)
    grid on
    xlabel('\omega')
    ylabel('F(\omega)/\sigma^2')

    subplot(1,2,2); 
    plot(t, ifftshift(Gamma_t)/(sigma_.^2), 'linewidth',1.5,'color', col(i-1,:)); hold on
    grid on
    set(gca,'LineWidth',1.5)
    xlabel('\Deltat')
    ylabel('F(\Deltat)/\sigma^2')
    axis([-T/2 T/2 -.2 1.2])
end
legend(flip(labels(1:2))); legend boxoff

%%
%denoised_ts = importdata(fullfile(fdir,sid{1},'PIANO3_BRAIN_MAPPING_ALT_INDEX_TAP_fMRI_SMS_s9_dn.nii.noproc.feat','tsplot','tsplot_zstat1.txt'));
%raw_ts = importdata(fullfile(fdir,sid{1},'PIANO3_BRAIN_MAPPING_ALT_INDEX_TAP_fMRI_SMS_s9.nii.noproc.feat','tsplot','tsplot_zstat1.txt'));
%sm_ts = importdata(fullfile(fdir,sid{1},'tapping_smoothing.feat','tsplot','tsplot_zstat1.txt'));
denoised_ts = importdata(fullfile(fdir,sid{20},'BRAIN_MAPPING_RIGHT_HAND_PIANO_s53_dn.feat','tsplot','tsplot_zstat1.txt'));
sm_ts = importdata(fullfile(fdir,sid{20},'BRAIN_MAPPING_RIGHT_HAND_PIANO_s53.feat','tsplot','tsplot_zstat1.txt'));
raw_ts = importdata(fullfile(fdir,sid{20},'BRAIN_MAPPING_RIGHT_HAND_PIANO3_s53.nii.noproc.feat','tsplot','tsplot_zstat1.txt'));


figure;
plot(1:160,denoised_ts(:,1),1:160,denoised_ts(:,2)); hold on
plot(1:160,raw_ts(:,1),1:160,raw_ts(:,2));

[rho_dn,p_dn] = corr(denoised_ts(:,1),denoised_ts(:,2));
[rho_rw,p_rw] = corr(raw_ts(:,1),raw_ts(:,2));

% mean_dn = 1;
% mean_rw = 1;
% dn_ts = (denoised_ts(:,1) - min(denoised_ts(:,1)))/(max(denoised_ts(:,1)) - min(denoised_ts(:,1)));
% dn_ts = dn_ts + mean_dn;
% dn_mod = (denoised_ts(:,2) - min(denoised_ts(:,1)))/(max(denoised_ts(:,2)) - min(denoised_ts(:,1)));
% dn_mod = dn_mod + mean_dn;
% rw_ts = (raw_ts(:,1) - min(raw_ts(:,1)))/(max(raw_ts(:,1)) - min(raw_ts(:,1)));
% rw_ts = rw_ts + mean_rw;
% rw_mod = (raw_ts(:,2) - min(raw_ts(:,1)))/(max(raw_ts(:,2)) - min(raw_ts(:,1)));
% rw_mod = rw_mod + mean_rw;

dn_ts1 = denoised_ts(:,1) - mean(denoised_ts(:,1));
rw_ts1 = raw_ts(:,1) - mean(raw_ts(:,1));
dn_mod1 = denoised_ts(:,2) - mean(denoised_ts(:,2));
rw_mod1 = raw_ts(:,2) - mean(raw_ts(:,2));

% time = 1:160;
% %dn_ts(62) = 1.49;
% figure('color','w');
% a = plot(time,dn_ts1,'-','LineWidth',1); hold on 
% b = plot(time,rw_ts1,'-','LineWidth',1); 
% plot(time,dn_mod1,'--','color',get(a,'color'),'LineWidth',1); 
% plot(time,rw_mod1,'--','color',get(b,'color'),'LineWidth',1); 
% set(gca,'ytick',[],'linewidth',1)
% legend(['MP-PCA \rho = ',num2str(rho_dn,2)],['original \rho = ',num2str(rho_rw,2)],'model design','location','best'); legend boxoff
% xlabel('Time (s)')
% ylabel('BOLD signal')
% 
dn_ts = normalize(denoised_ts(:,1));
rw_ts = normalize(raw_ts(:,1));
dn_mod = normalize(denoised_ts(:,2));
rw_mod = normalize(raw_ts(:,2));
% 
time = 1:160;
%dn_ts(62) = 1.49;
figure('color','w');
plot(time,dn_ts,'-','LineWidth',2); hold on
plot(time,rw_ts,'-','LineWidth',2); 
plot(time,dn_mod,'k--','LineWidth',2); 
plot(time,rw_mod,'k--','LineWidth',2); 
set(gca,'ytick',[],'linewidth',2)
legend(['MP-PCA \rho = ',num2str(rho_dn,2)],['original \rho = ',num2str(rho_rw,2)],'model design','location','best'); legend boxoff
xlabel('Time (s)')
ylabel('BOLD signal')

timeo = 1:160;
figure('color','w');
subtightplot(2,2,1)
plot(timeo,rw_ts1,'-','LineWidth',1); hold on
plot(timeo,rw_mod1,'k--','LineWidth',1); hold on
set(gca,'ytick',[],'linewidth',1)
xticks([])
xlabel('Time (s)')
ylabel('BOLD signal')
legend(['original \rho = ',num2str(rho_rw,2)],'model design','location','best'); legend boxoff
subtightplot(2,2,3)
plot(timeo,dn_ts1,'-','LineWidth',1); hold on
plot(timeo,dn_mod1,'k--','LineWidth',1); 
set(gca,'ytick',[],'linewidth',1)
legend(['MP-PCA \rho = ',num2str(rho_dn,2)],'model design','location','best'); legend boxoff
xlabel('Time (s)')
ylabel('BOLD signal')



%%
countM = 1;
% for each subject, load their relevant zstats and ROIs
for i = 1:numel(sid)
    % zstats
    zdir = dir(fullfile(fdir,sid{i},['*',series,'*.feat']));
    if isempty(zdir), continue; end
    if ~exist(fullfile(zdir(1).folder,zdir(1).name,'stats','zstat1.nii'),'file'), continue; end
    if ~exist(fullfile(zdir(2).folder,zdir(2).name,'stats','zstat1.nii'),'file'), continue; end
    nii = load_untouch_nii(fullfile(zdir(1).folder,zdir(1).name,'stats','zstat1.nii')); rw = double(nii.img);
    nii = load_untouch_nii(fullfile(zdir(2).folder,zdir(2).name,'stats','zstat1.nii')); dn = double(nii.img);
    nii = load_untouch_nii(fullfile(zdir(2).folder,zdir(2).name,'mask.nii')); mask = logical(nii.img);
    
    compMraw(countM:countM+sum(mask(:))-1) = rw(mask);
    compMdn(countM:countM+sum(mask(:))-1) = dn(mask);
    countM = sum(mask(:)) + 1;
end
 
% compMraw(compMraw==0) = NaN;
% compMdn(compMdn==0) = NaN;

figure('color','w');
histogram(compMraw,50,'Normalization','pdf'); hold on
histogram(compMdn,50,'Normalization','pdf')
%title(rois{i})
ylabel('pdf')
xlabel('Z score')
legend('raw','mppca','location','northwest'); legend boxoff

