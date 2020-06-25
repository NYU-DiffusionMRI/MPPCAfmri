function [ mask1, mask2, z, c1, c2] = thrRoi( roi, mask, thr, original, denoised)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    dnroi = roi;
    rwroi = roi;
    thrd = mean(denoised(roi==1)) - thr*std(denoised(roi==1));
    thrr = mean(original(roi==1)) - thr*std(original(roi==1));
    dnroi(denoised>thrd) = 0;
    rwroi(original>thrr) = 0;
    roi_ = logical(dnroi+rwroi);

    box = regionprops(roi_, 'BoundingBox');
    if isempty(box)
        roi_ = roi;
        box = regionprops(roi_, 'BoundingBox');
    end
    for b = 1:numel(box)
        vol(b) = box(b).BoundingBox(4)*box(b).BoundingBox(5)*box(b).BoundingBox(6);
    end
    [~,Ib] = max(vol);
    c1 = box(Ib).BoundingBox;
    cStart = ceil(c1(1:3));
    cLen = floor(c1(4:end)) - 1;   
%     for i = 1:numel(cLen)
%         if cLen(i) == size(roi,i), cLen(i) = cLen(i) - 1; end
%     end
    
    mask1 = zeros(size(roi_));
    mask1(cStart(2):cStart(2)+cLen(2),cStart(1):cStart(1)+cLen(1),cStart(3):cStart(3)+cLen(3)) = 1;
    
    com = regionprops(mask,'Centroid'); 
    com = com(1).Centroid;
    mask2 = zeros(size(roi_));
    xdist = floor(2.*(com(2) - mean([cStart(2),cStart(2)+cLen(2)])));
    c2 = c1;
    c2(2) = c2(2) + xdist;
    mask2(cStart(2)+xdist:cStart(2)+cLen(2)+xdist,cStart(1):cStart(1)+cLen(1),cStart(3):cStart(3)+cLen(3)) = 1;
    
    z = floor(median([cStart(3),cStart(3)+cLen(3)]));

    mask1(mask==0) = 0;
    mask2(mask==0) = 0;
    mask1 = logical(mask1);
    mask2 = logical(mask2);
end

