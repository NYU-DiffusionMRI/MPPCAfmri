clear
close all

root = '/Volumes/Research/fieree01lab/labspace/Projects/RMT-fMRI/simulate';
addpath(root);
addpath('/Volumes/Research/fieree01lab/labspace/Projects/RMT-fMRI/scripts');
%root = '/mnt/labspace/Projects/RMT-fMRI/simulate';
 Nroi = niftiread('~/Documents/RMT-fMRI/generate_synthetic_brain_data/generate_synthetic_brain_data/N.nii');
 Yroi = niftiread('~/Documents/RMT-fMRI/generate_synthetic_brain_data/generate_synthetic_brain_data/Y.nii');
 Uroi = niftiread('~/Documents/RMT-fMRI/generate_synthetic_brain_data/generate_synthetic_brain_data/U.nii');
data = niftiread(fullfile(root,'simhrf_3.nii'));
cov = niftiread(fullfile(root,'cov_2.nii'));
p = phantom3d(66).*100;
mask = p>0;

covp = p/100+cov*100;

%snr = [1:1:9, 10:10:100, 200:100:1000];
snr = 0.1:1:15;
reali = 1;
mu = mean(p(p>0));
%vox = [27, 32, 19];
vox = [34, 45, 24];

X = load(fullfile(root,'design.mat')); 
design_ = X.X;
%design_ = [-0.462707400000000;-0.462707400000000;-0.462707400000000;-0.462707400000000;-0.462707400000000;-0.462707400000000;-0.462707400000000;-0.462707400000000;-0.462707400000000;-0.462707400000000;-0.462707400000000;-0.462707400000000;-0.462707400000000;-0.462707400000000;-0.462707400000000;-0.462707400000000;-0.462707400000000;-0.462707400000000;-0.462707400000000;-0.462707400000000;-0.462235600000000;-0.442676400000000;-0.371956900000000;-0.251714900000000;-0.106033800000000;0.0398986100000000;0.169366800000000;0.274955600000000;0.355912200000000;0.415119300000000;0.456835200000000;0.485353300000000;0.504369100000000;0.516785900000000;0.524750500000000;0.529781300000000;0.532916800000000;0.534848100000000;0.536025400000000;0.536736400000000;0.536690400000000;0.517261600000000;0.446542100000000;0.326300100000000;0.180619000000000;0.0346865800000000;-0.0947815900000000;-0.200370400000000;-0.281327000000000;-0.340534100000000;-0.382250000000000;-0.410768100000000;-0.429783900000000;-0.442200700000000;-0.450165300000000;-0.455196100000000;-0.458331600000000;-0.460262900000000;-0.461440200000000;-0.462151200000000;-0.462105300000000;-0.442676400000000;-0.371956900000000;-0.251714900000000;-0.106033800000000;0.0398986100000000;0.169366800000000;0.274955600000000;0.355912200000000;0.415119300000000;0.456835200000000;0.485353300000000;0.504369100000000;0.516785900000000;0.524750500000000;0.529781300000000;0.532916800000000;0.534848100000000;0.536025400000000;0.536736400000000;0.536690400000000;0.517261600000000;0.446542100000000;0.326300100000000;0.180619000000000;0.0346865800000000;-0.0947815900000000;-0.200370400000000;-0.281327000000000;-0.340534100000000;-0.382250000000000;-0.410768100000000;-0.429783900000000;-0.442200700000000;-0.450165300000000;-0.455196100000000;-0.458331600000000;-0.460262900000000;-0.461440200000000;-0.462151200000000;-0.462105300000000;-0.442676400000000;-0.371956900000000;-0.251714900000000;-0.106033800000000;0.0398986100000000;0.169366800000000;0.274955600000000;0.355912200000000;0.415119300000000;0.456835200000000;0.485353300000000;0.504369100000000;0.516785900000000;0.524750500000000;0.529781300000000;0.532916800000000;0.534848100000000;0.536025400000000;0.536736400000000;0.536690400000000;0.517261600000000;0.446542100000000;0.326300100000000;0.180619000000000;0.0346865800000000;-0.0947815900000000;-0.200370400000000;-0.281327000000000;-0.340534100000000;-0.382250000000000;-0.410768100000000;-0.429783900000000;-0.442200700000000;-0.450165300000000;-0.455196100000000;-0.458331600000000;-0.460262900000000;-0.461440200000000;-0.462151200000000;-0.462105300000000;-0.442676400000000;-0.371956900000000;-0.251714900000000;-0.106033800000000;0.0398986100000000;0.169366800000000;0.274955600000000;0.355912200000000;0.415119300000000;0.456835200000000;0.485353300000000;0.504369100000000;0.516785900000000;0.524750500000000;0.529781300000000;0.532916800000000;0.534848100000000;0.536025400000000;0.536736400000000];
design = permute(repmat(design_, [1, size(data,1), size(data,2), size(data,3)]), [2 3 4 1]);

zthr = 15;

zgt_ = zstat(data, design);
zgt_ = real(zgt_);

b = [0.00,0.45,0.74];
r = [0.85,0.33,0.10];
y = [0.93,0.69,0.13];

zraw = zeros(size(data,1), size(data,2), size(data,3), reali);
zsm = zraw;
zmp = zsm;

Nmask = zeros(size(mask));
Nmask(21:35,15:25,24) = 1; Nmask = logical(Nmask);
Ymask = zeros(size(mask));
Ymask(21:35,25:35,24) = 1; Ymask = logical(Ymask);
Umask = zeros(size(mask));
Umask(21:35,35:45,24) = 1; Umask = logical(Umask);

Sgt(1) = sum(zgt_(Nmask)>zthr) / sum(Nmask(:)) * 100;
Sgt(2) = sum(zgt_(Ymask)>zthr) / sum(Ymask(:)) * 100;
Sgt(3) = sum(zgt_(Umask)>zthr) / sum(Umask(:)) * 100;

% %% lets choose how many parameters to keep
% count = 1;
% %sigma = 20/3;
% sigma = 20/3;
% design = permute(repmat(design_, [1, 15,31,10]), [2 3 4 1]);
% np = 15;
% S = zeros(3,np);
% zthr = 3;
% td = data(:,:,24,:); 
% td = repmat(td,[1,1,10,1]);
% d = td(21:35,15:45,:,:) + sigma*randn(size(td(21:35,15:45,:,:)));
% figure;
% for p = 0:np
%     %[mp, sig, par] = MP_Loop3_2(data+sigma*randn(size(data)),[5,5,5],p);
%     [mp, sig, par] = MP_Loop3_thr(d,[3,3,3],p);
%     zmp_nf = zstat(mp, design);
%     S(1,count) = sum(sum(zmp_nf(:,1:11,5)>zthr)) / sum(Nmask(:)) * 100;
%     S(2,count) = sum(sum(zmp_nf(:,11:21,5)>zthr)) / sum(Ymask(:)) * 100;
%     S(3,count) = sum(sum(zmp_nf(:,21:31,5)>zthr)) / sum(Umask(:)) * 100;
%     subplot(4,4,count)
%     imagesc(zmp_nf(:,:,5),[0 10]); axis off; title(['snr = ',num2str(snr(i))]);
%     count = count + 1;
% end
% [mp,sig,mppar] = MP_Loop3_3(d,[3,3,3]);
% zmp_n = zstat(mp,design);
% SmpN = sum(sum(zmp_n(:,1:11,5)>zthr)) / sum(Nmask(:)) * 100;
% SmpY = sum(sum(zmp_n(:,11:21,5)>zthr)) / sum(Ymask(:)) * 100;
% SmpU = sum(sum(zmp_n(:,21:31,5)>zthr)) / sum(Umask(:)) * 100;
% 
% figure;
% plot(0:np,S(1,:),'color',b); hold on
% plot(0:np,S(2,:),'color',r);
% plot(0:np,S(3,:),'color',y);
% plot(0:np,sum(sum(Nroi(:,:,24)))/sum(Nmask(:))*100*ones(length(0:np)),'--','color',b);
% plot(0:np,sum(sum(Yroi(:,:,24)))/sum(Ymask(:))*100*ones(length(0:np)),'--','color',r);
% plot(0:np,sum(sum(Uroi(:,:,24)))/sum(Umask(:))*100*ones(length(0:np)),'--','color',y);
% plot(0:np,SmpN*ones(length(0:np)),'.','color',b);
% plot(0:np,SmpY*ones(length(0:np)),'.','color',r);
% plot(0:np,SmpU*ones(length(0:np)),'.','color',y);
%%
S = zeros(3,length(snr));
zthr = 3;
figure;
for i = 1:length(snr)
    sigma = 20/snr(i);
    %d = td(21:35,15:45,:,:) + sigma*randn(size(td(21:35,15:45,:,:)));
    d = sigma*randn(size(td(21:35,15:45,:,:)));
    [mp, sig, par] = MP_Loop3_thr(d,[3,3,3],2);
    zmp_n = zstat(mp, design);
    S(1,i) = sum(sum(zmp_n(:,1:11,5)>zthr)) / sum(Nmask(:)) * 100;
    S(2,i) = sum(sum(zmp_n(:,11:21,5)>zthr)) / sum(Ymask(:)) * 100;
    S(3,i) = sum(sum(zmp_n(:,21:31,5)>zthr)) / sum(Umask(:)) * 100;
    subplot(4,4,i)
    imagesc(zmp_n(:,:,5),[0 3]); axis off; title(['snr = ',num2str(snr(i))]);
end
figure;
plot(snr,S(1,:),'color',b); hold on
plot(snr,S(2,:),'color',r)
plot(snr,S(3,:),'color',y)
plot(snr,Sgt(1)*ones(length(snr)),'-','color',b,'LineWidth',1); hold on
plot(snr,Sgt(2)*ones(length(snr)),'-','color',r,'LineWidth',1); hold on
plot(snr,Sgt(3)*ones(length(snr)),'-','color',y,'LineWidth',1); hold on
xlabel('snr')
ylabel('%z>10')


%%
droiN = d(floor(end/2)-2 : floor(end/2)+2, 3:7, 3:7,:);
droiY = d(floor(end/2)-2 : floor(end/2)+2, 13:17, 3:7,:);
droiU = d(floor(end/2)-2 : floor(end/2)+2, 23:27, 3:7,:);
droiN_ = reshape(droiN, [125,160]);
droiY_ = reshape(droiY, [125,160]);
droiU_ = reshape(droiY, [125,160]);
figure;
[u, vals, v] = svd(droiN_, 'econ');
vals = diag(vals).^2 / 160;  
histogram(vals,100); hold on
[u, vals, v] = svd(droiY_, 'econ');
vals = diag(vals).^2 / 160;  
histogram(vals,100); hold on
[u, vals, v] = svd(droiU_, 'econ');
vals = diag(vals).^2 / 160;  
histogram(vals,100); hold on


%%

td = data(:,:,24,:); 
td = repmat(td,[1,1,10,1]);

figure;
zthr = 20;
for i = 1:length(snr)
    disp(['snr = ',num2str(snr(i))])
    for j = 1:reali
        sigma = mu/snr(i);
        noisy = td(21:35,15:45,:,:) + sigma*randn(size(td(21:35,15:45,:,:)));
        %noisy = data(21:35,15:45,20:30,:) + sigma*randn(size(data(21:35,15:45,20:30,:)));
        %noisy = data + sigma*randn(size(data));   
        [mp, sm, sig, par] = rundn(noisy);
        %P(:,:,:,i) = par;

        zraw_ = zstat(noisy, design);
        zsm_ = zstat(sm, design);
        zmp_ = real(zstat(mp, design));
                      
        zraw(:,:,:,j,i) = zraw_;
        zsm(:,:,:,j,i) = zsm_;
        zmp(:,:,:,j,i) = zmp_;
        
        Srw(i,1) = sum(sum(zraw_(:,1:11,5)>zthr)) / sum(Nmask(:)) * 100;
        Srw(i,2) = sum(sum(zraw_(:,11:21,5)>zthr)) / sum(Ymask(:)) * 100;
        Srw(i,3) = sum(sum(zraw_(:,21:31,5)>zthr)) / sum(Umask(:)) * 100;
        Ssm(i,1) = sum(sum(zsm_(:,1:11,5)>zthr)) / sum(Nmask(:)) * 100;
        Ssm(i,2) = sum(sum(zsm_(:,11:21,5)>zthr)) / sum(Ymask(:)) * 100;
        Ssm(i,3) = sum(sum(zsm_(:,21:31,5)>zthr)) / sum(Umask(:)) * 100;
        Smp(i,1) = sum(sum(zmp_(:,1:11,5)>zthr)) / sum(Nmask(:)) * 100;
        Smp(i,2) = sum(sum(zmp_(:,11:21,5)>zthr)) / sum(Ymask(:)) * 100;
        Smp(i,3) = sum(sum(zmp_(:,21:31,5)>zthr)) / sum(Umask(:)) * 100;        
%         Ssm(i,j) = sum(zsm_(mask)>zthr) / sum(mask(:)) * 100;
%         Smp(i,j) = sum(zmp_(mask)>zthr) / sum(mask(:)) * 100;
        subplot(4,4,i)
        imagesc(zmp_(:,:,5),[0 zthr]); axis off; title(['snr = ',num2str(snr(i))]);
    end
end
% save(fullfile(root,'zraw_MPfull.mat'),'zraw')
% save(fullfile(root,'zsm_MPfull.mat'),'zsm')
% save(fullfile(root,'zmp_MPfull.mat'),'zmp')
%%
zraw_mean = mean(zraw(:,:,:,1,4),4); zraw_mean(mask == 0) = 0;
zmp_mean = real(mean(zmp(:,:,:,1,4),4)); zmp_mean(mask == 0) = 0;
zsm_mean = mean(zsm(:,:,:,1,4),4); zsm_mean(mask == 0) = 0;

p(cov>=.005) = 0;
gtf = real(zgt_ + p/1000);
figure('color','k');
subtightplot(1,4,1)
imagesc(gtf(:,:,vox(3))); colormap(gray); caxis([0 1]); axis off; axis image
title('Ground Truth','color','w','interpreter','latex','fontsize',20); 
vox(3)=5;
subtightplot(1,4,2)
imagesc(zraw_mean(:,:,vox(3))); colormap(hot); caxis([0 zthr]); axis off; axis image
title('No Processing','color','w','interpreter','latex','fontsize',20); 
subtightplot(1,4,3)
imagesc(zsm_mean(:,:,vox(3))); colormap(hot); caxis([0 zthr]); axis off; axis image
title('Smoothing','color','w','interpreter','latex','fontsize',20); 
subtightplot(1,4,4)
imagesc(zmp_mean(:,:,vox(3))); colormap(hot); caxis([0 zthr]); axis off; axis image
title('MP-PCA','color','w','interpreter','latex','fontsize',20); 

figure;
imagesc(P(:,:,24,4)); colorbar;

%%
% figure('color','w');
% plot((squeeze(zraw(:,45,26,1,5))'+squeeze(zsm(34,:,26,1,5)))/2,'color',b,'LineWidth',1); hold on
% plot((squeeze(zsm(:,45,26,1,5))'+squeeze(zraw(34,:,26,1,5)))/2,'color',y','LineWidth',1); hold on
% plot((squeeze(zmp(:,45,26,1,5))'+squeeze(zmp(34,:,26,1,5)))/2,'color',r,'LineWidth',1); hold on
% plot((squeeze(zgt_(:,45,26))'+squeeze(zgt_(34,:,26)))/2,'k-','LineWidth',1); hold on

%%
% close all
% figure('color','w');
% errorbar(snr,mean(Srw,2),std(Srw,[],2)/sqrt(size(Srw,2)),'color',b,'LineWidth',1); hold on
% errorbar(snr,mean(Ssm,2),std(Ssm,[],2)/sqrt(size(Ssm,2)),'color',y,'LineWidth',1); hold on
% errorbar(snr,mean(Smp,2),std(Smp,[],2)/sqrt(size(Smp,2)),'color',r,'LineWidth',1); hold on
% plot(snr,mean(Srw,2),'-',snr,mean(Smp,2),'-',snr,mean(Ssm,2),'-','LineWidth',1)
% plot(snr,Sgt*ones(1,length(snr)),'k--','LineWidth',1);
% L = get(gca,'Children'); 
% legend([L(4),L(3),L(2),L(1)],{'No processing','MP-PCA','Smoothing','Ground Truth'},'Location','NorthWest','Interpreter','latex'); legend boxoff
% xlabel('SNR','Interpreter','latex')
% ylabel('$\%\ z > 3$','interpreter','latex')
% set(gca,'LineWidth',1,'FontSize',20,'FontName','Times New Roman')
% xlim([0 15]);



figure('color','w');
plot(snr,Srw(:,1),'-','color',b,'LineWidth',1); hold on
text(snr,Srw(:,1),'N','color',b);
plot(snr,Smp(:,1),'-','color',r,'LineWidth',1); hold on
text(snr,Smp(:,1),'N','color',r);
plot(snr,Ssm(:,1),'-','color',y,'LineWidth',1); hold on
text(snr,Ssm(:,1),'N','color',y);
plot(snr,Sgt(1)*ones(length(snr)),'-','color','k','LineWidth',1); hold on
plot(snr,Srw(:,2),'-','color',b,'LineWidth',1); hold on
text(snr,Srw(:,2),'Y','color',b);
plot(snr,Smp(:,2),'-','color',r,'LineWidth',1); hold on
text(snr,Smp(:,2),'Y','color',r);
plot(snr,Ssm(:,2),'-','color',y,'LineWidth',1); hold on
text(snr,Ssm(:,2),'Y','color',y);
plot(snr,Sgt(2)*ones(length(snr)),'-','color','k','LineWidth',1); hold on
plot(snr,Srw(:,3),'-','color',b,'LineWidth',1); hold on
text(snr,Srw(:,3),'U','color',b);
plot(snr,Smp(:,3),'-','color',r,'LineWidth',1); hold on
text(snr,Smp(:,3),'U','color',r);
plot(snr,Ssm(:,3),'-','color',y,'LineWidth',1); hold on
text(snr,Ssm(:,3),'U','color',y);
plot(snr,Sgt(3)*ones(length(snr)),'-','color','k','LineWidth',1); hold on

L = get(gca,'Children'); 
legend([L(end),L(end-16),L(end-32),L(5)],{'No processing','MP-PCA','Smoothing','Ground Truth'},'Location','NorthWest','Interpreter','latex'); legend boxoff
xlabel('SNR','Interpreter','latex')
ylabel('$\%\ z > 3$','interpreter','latex')
set(gca,'LineWidth',1,'FontSize',20,'FontName','Times New Roman')
xlim([0 15]);
%%

function z = zstat(input, design)
    r = sum((input-mean(input,4)).*(design-mean(design,4)), 4)./ ...
    sqrt(sum((input-mean(input,4)).^2, 4).*sum((design-mean(design,4)).^2, 4));
    
    z = atanh(r)*sqrt(size(input,4) - 3);
    z(isnan(z)) = 0; z = real(z);
end

function [mp, sm, sig, par] = rundn(noisy)
    [mp, sig, par] = MP_Loop3_3(noisy,[5,5,5]);
    %noisy_ = reshape(noisy,[size(noisy,1)*size(noisy,2)*size(noisy,3), size(noisy,4)]);
    %mp_ = MP_GREG(noisy_);
    %mp = reshape(mp_, [size(noisy,1), size(noisy,2), size(noisy,3), size(noisy,4)]);  
    sm = zeros(size(noisy));
    for v = 1:size(noisy,4)
        sm(:,:,:,v) = imgaussfilt3(noisy(:,:,:,v), 1.66);
    end
    mp(isnan(mp)) = 0;
    sm(isnan(sm)) = 0;
end
