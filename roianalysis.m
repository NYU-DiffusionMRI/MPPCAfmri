addpath(genpath('~/Documents/MATLAB'));
addpath('/Volumes/Research/fieree01lab/labspace/Projects/RMT-fMRI/scripts')
clear 
set(0,'defaultAxesFontSize',20)
set(0,'defaultAxesFontName','Arial')
close all

subjects = [1 2 5 7 8 9 10 11 12 13 14 16 17 18 19 20 21 22 23 24 25 27 28 29 30];
sid = strcat('NYU-',strsplit(num2str(subjects)));

rois = {'LB','LH','LL','LM','LnW','RB','RH','RL','RM','RnW'};
%rois = {'LB','LL','LW','RB'}; %,'LW','RB','RH','RL','RM','RW'};
sphere_r = {'.nii','_sphere10_bin.nii','_sphere15_bin.nii','_sphere20_bin.nii'};
spheret = {'single voxel','sphere r=10v','sphere r=15v','sphere r=20v'};
%sphere_r = {'_sphere20_bin.nii'};
%spheret = {'sphere r=20v'};

root = '/Volumes/Research/fieree01lab/labspace/Projects/RMT-fMRI';
fdir = fullfile(root,'PROCESSING_T1_v1');
rdir = fullfile(root,'ROI_tim');

series = 'VERB';
sertit ='Verb Generation';
comproiraw = []; %zeros(100000,10);
comproiraw_ = {};
comproidn = []; %zeros(100000,10);
count = 1;
countM = 1;

% for each subject, load their relevant zstats and ROIs
c = 0;
for i = 1:numel(sid)
    % zstats
    zdir = dir(fullfile(fdir,sid{i},['*',series,'*.noproc.feat']));
    if isempty(zdir), continue; end
    if ~exist(fullfile(zdir(1).folder,zdir(1).name,'stats','tstat1.nii'),'file'), continue; end
    if ~exist(fullfile(zdir(2).folder,zdir(2).name,'stats','tstat1.nii'),'file'), continue; end
    cd(fullfile(zdir(1).folder,zdir(1).name,'stats'));
    nii = load_untouch_nii('tstat1.nii'); rw = double(nii.img);
    cd(fullfile(zdir(1).folder,zdir(2).name,'stats'));
    nii = load_untouch_nii('tstat1.nii'); dn = double(nii.img);
    cd('..')
    nii = load_untouch_nii('mask.nii'); mask = logical(nii.img);
%     cd('..')
%     nii = load_untouch_nii(fullfile(zdir(1).folder,[zdir(1).name(1:end-12),'_noisemap.nii'])); sigma = double(nii.img);
%     nii = load_untouch_nii(fullfile(zdir(1).folder,[zdir(1).name(1:end-12),'.nii'])); rawbold = double(nii.img);
%     nii = load_untouch_nii(fullfile(zdir(1).folder,[zdir(1).name(1:end-12),'_dn.nii'])); dnbold = double(nii.img);
%     

    compMraw(countM:countM+sum(mask(:))-1) = rw(mask);
    compMdn(countM:countM+sum(mask(:))-1) = dn(mask);
    countM = sum(mask(:)) + 1;
    
    % rois
    roilist = fullfile(rdir,sid{i});
    for j = 1:numel(rois)
        %comproiraw_ = [];
        for k = 1:numel(sphere_r)
            if ~exist(fullfile(roilist,[rois{j},sphere_r{k}]),'file')
                rwdata(i,j,k,:) = [NaN NaN NaN NaN];
                dndata(i,j,k,:) = [NaN NaN NaN NaN];
                pdata(i,j,k) = NaN;
                continue
            end
            nii = load_untouch_nii(fullfile(roilist,[rois{j},sphere_r{k}])); roi = logical(nii.img);
            %roi(abs(dn)<3|abs(rw)<3) = 0;
            roi_ = roi;
            %[roi_, croi_, z_, c1, c2] = thrRoi(roi,mask,3,rw,dn);
                      
%             for m = 1:size(roi,3)
%                 hull = bwconvhull(roi(:,:,m));
%                 roi_(:,:,m) = hull;
%             end
            
%             rwsnr(i) = nanmedian(rw(roi_))/mean(sigma(roi_));
%             dnsnr(i) = nanmedian(dn(roi_))/mean(sigma(roi_));
%             snrdiff(i) = abs(dnsnr-rwsnr)/abs(rwsnr)*100;
% if size(rawbold,4)<150, continue; end
%              count = count+1;
%     snrraw(:,count) = nanmean(vec(rawbold,roi_),2)./nanmean(sigma(roi_));
%     snrdn(:,count) = nanmean(vec(dnbold,roi_),2)./nanmean(sigma(roi_));
% 
%             
            rwdata(i,j,k,1) = nanmean(rw(roi_));
            rwdata(i,j,k,2) = nanmedian(rw(roi_));
            rwdata(i,j,k,3) = nanstd(rw(roi_));
            rwdata(i,j,k,4) = sum(rw(roi_)>3);
            
            dndata(i,j,k,1) = nanmean(dn(roi_));
            dndata(i,j,k,2) = nanmedian(dn(roi_));
            dndata(i,j,k,3) = nanstd(dn(roi_));
            dndata(i,j,k,4) = sum(dn(roi_)>3);
            
            if k == 1 || sum(roi_(:))<5
                P = NaN;
            else
                [H,P,CI] = ttest(rw(roi_),dn(roi_));
            end
            pdata(i,j,k) = P;
            
            c = c+1;
            table{c,1} = sid{i};
            table{c,2} = rois{j};
            table{c,3} = spheret{k};
            table{c,4} = nanmean(rw(roi_));
            table{c,5} = nanmean(dn(roi_));
            table{c,6} = nanmedian(rw(roi_));
            table{c,7} = nanmedian(dn(roi_));
            table{c,8} = nanstd(rw(roi_));
            table{c,9} = nanstd(dn(roi_));
            table{c,10} = sum(abs(rw(roi_))>3);
            table{c,11} = sum(abs(dn(roi_))>3);
            table{c,12} = sum(abs(rw(roi_))>6);
            table{c,13} = sum(abs(dn(roi_))>6);
            table{c,14} = min(rw(roi_));
            table{c,15} = min(dn(roi_));
            table{c,16} = max(rw(roi_));
            table{c,17} = max(dn(roi_));
            %c = c+20;
                
            % vectoize z stats for composite histogram
            if k == 3
%                 comproiraw(count:count+sum(roi_(:))-1,j) = rw(roi_);
%                 comproidn(count:count+sum(roi_(:))-1,j) = dn(roi_);
%                 count = sum(roi_(:)) + 1;
                
                comproiraw_ = cat(2, comproiraw_, rw(roi_));
                
                
            end                        
        end
    end
    comproiraw = cat(1, comproiraw, comproiraw_);
    
end


T = cell2table(table,'VariableNames',{'ID','ROI','ROIsize','meanraw','meanDN','medianraw','medianDN','stdraw','stdDN','rawzgt3','dnzgt3','rawzgt6','dnzgt6','minrawz','mindnz','maxrawz','maxdnz'});
writetable(T,['~/Desktop/',series,'_ROIdata20.csv']);
% generate a table of subject, ROI, ROI size, mean zscore, median
% zscore, z std, count of z>3, count of z>5, min z, max z


figure('color','w');
histogram(compMraw,50,'Normalization','pdf'); hold on
histogram(compMdn,50,'Normalization','pdf')
title(sertit)
ylabel('pdf')
xlabel('Z score')
legend('raw','mppca','location','northwest'); legend boxoff
xlim([-10 10])

%save('~/Desktop/listcomp_dndata.mat','dndata');
%save('~/Desktop/listcomp_rwdata.mat','rwdata');

%%
%comproiraw(comproiraw==0) = NaN;
%comproidn(comproidn==0) = NaN;
% 

% comproiraw = abs(comproiraw);
% comproidn = abs(comproidn);
keyboard

%h = raincloud_plot(comproiraw(:,1), 'box_on', 1);
% figure
% a = num2cell(comproiraw);
% h = rm_raincloud(a ,hot(10));
% keyboard

%%

figure('color','w');
for i = 1:10
    subtightplot(2,5,i,[.1 .1])
    A.raw = comproiraw(:,i);
    A.denoised = comproidn(:,i);
    nhist(A,'smooth','linewidth',1.5,'fsize',15,'nolegend')
    
    
%     
%     subtightplot(2,5,i,[.03 .03])
%     histogram(comproiraw(:,i),50,'Normalization','count'); hold on
%     histogram(comproidn(:,i),50,'Normalization','count')
     title(rois{i})
%     %ylabel('pdf')
%     %xlabel('Z score')
        if i == 10
     legend('raw','denoised','location','northwest'); legend boxoff
        end
     ax = gca;
     ax.XRuler.Exponent = 0;
     ax.LineWidth = 1;
     if ismember(i,6:10)
         xlabel('t-score')
     end
     xlim([0 15])
end
suptitle(sertit)




%%


% 
% rwdata(rwdata == 0) = NaN;
% dndata(dndata == 0) = NaN;


figure('color','w')
for i = 1:4
    clear bardata
    subplot(2,2,i)

    bardata(:,1) = nanmean((rwdata(:,:,i,2)),1);
    bardata(:,2) = nanmean((dndata(:,:,i,2)),1);
   sbardata(:,1) = nanmean((rwdata(:,:,i,3)),1);
   sbardata(:,2) = nanmean((dndata(:,:,i,3)),1);
    
     hBar = bar(1:10,bardata);
    for k1 = 1:size(bardata',1)
        ctr(k1,:) = bsxfun(@plus, hBar(1).XData, [hBar(k1).XOffset]');
        ydt(k1,:) = hBar(k1).YData;
    end
    hold on
    errorbar(ctr, ydt, sbardata', '.k')
    hold off
    set(gca,'XTickLabels',rois)
    legend('raw','denoised','location','northeast')
    title(spheret{i})
    %axis([0 11 0 10])
    grid on
    ylabel('mean Z score')
    
end
suptitle(sertit)










