clear 
close all
clc

addpath(genpath('~/Documents/MATLAB'))
set(0,'defaultAxesFontSize',20)
set(0,'defaultAxesFontName','Arial')

subjs={'7','8','9','10','11','12','13','14','16','17','18','19','20','21','22','23','24','25','29','30'};
%subjs={'20','21','24','29','30'};
series = {'verbgen','sentcomp','listening','tapping'};
seriest = {'Verb Generation', 'Sentence Completion','Listening Comprehension','Finger Tapping'};
% rw_v = zeros(numel(subjs),121);
% dn_v = zeros(numel(subjs),121);

for s = 1:numel(subjs)
    root=['/Volumes/Research/fieree01lab/labspace/Projects/RMT-fMRI/PROCESSING_T1_v1/NYU-',subjs{s},'/vary_scan_length'];

    c = 0;
    for i = 40:160
        c = c+1;
        for j = 1:numel(series)
            rw = dir(fullfile(root,num2str(i),[series{j},'.feat'],'stats','zstat1.nii*'));
            dn = dir(fullfile(root,num2str(i),[series{j},'_dn.feat'],'stats','zstat1.nii*'));
            if isempty(rw)
                rw_v(s,c,j) = NaN;
                dn_v(s,c,j) = NaN;
            else
                rw_z = double(niftiread(fullfile(rw.folder,rw.name)));
                dn_z = double(niftiread(fullfile(dn.folder,dn.name)));
                rw_v(s,c,j) = sum(abs(rw_z(:))>3);
                dn_v(s,c,j) = sum(abs(dn_z(:))>3);
            end
        end 
    end
end
%%
close all
mrw_v = squeeze(nanmean(rw_v,1));
mdn_v = squeeze(nanmean(dn_v,1));
x = (40:160)/160*100;

red = [1, 0.37, 0.41];
blue = [0, 0.74, 0.84];

close all
figure('color','w')
count = 1;
for i = [4, 1, 2, 3]
    subtightplot(1,4,count,[.05 .05],.2,.05)
    stdshade(rw_v(:,:,i)/50000*100,0.25,blue,x,.5); hold on
    stdshade(dn_v(:,:,i)/50000*100,0.25,red,x,.5)
    xlim([25 100])
    xlabel('Included Scan %')
    ylabel('% t>3')
    title(seriest{i})
    grid on
    ylim([0 35])
    set(gca,'LineWidth',1.5)
    count = count + 1;
end
g = gcf;
legend([g.Children(1).Children(1),g.Children(1).Children(3)],{'MP-PCA','Original'}); legend boxoff
hgsave(gcf,'~/Desktop/vtime.fig')

%%
%g(1,4).coord_flip;
%g(1,4).set_names('x','','y','t statistic','color','Method');
%g(1,4).axe_property('YLim',[-10 15],'LineWidth',1,'FontSize',16,'Ygrid','on','XTickLabel',roilabels);
%g(1,4).set_title('left','position',[-.9 1]);
%g(1,4).set_layout_options('legend',false);

% figure('color','w')
% for i = 1:numel(series)
%     subplot(1,4,i)
%     plot(x,mrw_v(:,i),'b.','markersize',20); hold on
%     plot(x,mdn_v(:,i),'r.','markersize',20);
%     title(seriest{i})
%     xlabel('% of total scan time');
%     ylabel('#voxels: z>3')
%     legend('original','denoised','Location','Best'); legend boxoff
%     xlim([40 100])
% end

%% 
% clear 
% close all
% clc
% 
% addpath(genpath('~/Documents/MATLAB'))
% set(0,'defaultAxesFontSize',20)
% set(0,'defaultAxesFontName','Arial')
% subjs={'24'};
% series = {'verbgen','listening','sentcomp','tapping'};
% seriest = {'Verb Generation','Listening Comprehention','Sentence Completion','Finger Tapping'};
% brw_v = zeros(4,4);
% bdn_v = zeros(4,4);
% 
% root=['/Volumes/Research/fieree01lab/labspace/Projects/RMT-fMRI/PROCESSING_T1_v1/NYU-',subjs{1},'/vary_scan_length'];
% inds = floor([.60 .70 .80 .90 1].*160);
% c = 0;
% for i = inds
%     c = c+1;
%     for j = 1:numel(series)
%         try
%             rw = fullfile(root,num2str(i),[series{j},'.feat'],'stats','tstat1.nii');
%             dn = fullfile(root,num2str(i),[series{j},'_dn.feat'],'stats','tstat1.nii');
% 
%             nii = load_untouch_nii(rw); rw_z = nii.img;
%             nii = load_untouch_nii(dn); dn_z = nii.img;
% 
%             brw_v(j,c) = sum(abs(rw_z(:))>3);
%             bdn_v(j,c) = sum(abs(dn_z(:))>3);
%         catch
%             brw_v(j,c) = NaN;
%             bdn_v(j,c) = NaN;
%         end
%     end
% end
% 
% bardata = [brw_v(1,end), flip(bdn_v(1,1:end))];
% bardata = [bardata; brw_v(2,end), flip(bdn_v(2,1:end))];
% bardata = [bardata; brw_v(3,end), flip(bdn_v(3,1:end))];
% bardata = [bardata; brw_v(4,end), bdn_v(4,1:end)];
% bardata(2,2:end) = sort(bardata(2,2:end),'descend');
% bardata(4,2:end) = sort(bardata(4,2:end),'descend');
% 
% 
% figure('color','w')
% bar(bardata); 
% set(gca,'XTickLabels',seriest)
% ylabel('number of voxels: t>3')
% legend('Original','100% denoising','90% denoising','80% denoising','70% denoising','60% denoising','Location','Best'); legend boxoff
% 
