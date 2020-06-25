#!/bin/sh

#  getsphere.sh
#  
#
#  Created by Benjamin Ades-Aron on 5/1/18.
#
export FSLOUTPUTTYPE=NIFTI
root=/Volumes/Research/fieree01lab/labspace/Projects/RMT-fMRI/ROI_tim
subj=(NYU-1 NYU-2 NYU-5 NYU-7 NYU-8 NYU-9 NYU-10 NYU-11 NYU-12 NYU-13 NYU-14 NYU-16 NYU-17 NYU-18 NYU-19 NYU-20 NYU-21 NYU-22 NYU-23 NYU-24 NYU-25 NYU-27 NYU-28 NYU-29 NYU-30 NYU-31)
rois=(LnW RnW)
# LB_3mm LH_3mm LL_3mm LM_3mm LW_3mm RB_3mm RH_3mm RL_3mm RM_3mm RW_3mm)

for i in ${subj[@]}; do
    cd /Volumes/Research/fieree01lab/labspace/Projects/RMT-fMRI/ROI_tim/$i

    for j in ${rois[@]}; do
        fslmaths ${j}.nii -kernel sphere 5 -fmean ${j}_sphere5
        fslmaths ${j}_sphere5 -bin ${j}_sphere5_bin
        
        fslmaths ${j}.nii -kernel sphere 10 -fmean ${j}_sphere10
        fslmaths ${j}_sphere10 -bin ${j}_sphere10_bin

        fslmaths ${j}.nii -kernel sphere 15 -fmean ${j}_sphere15
        fslmaths ${j}_sphere15 -bin ${j}_sphere15_bin
#
        fslmaths ${j}.nii -kernel sphere 20 -fmean ${j}_sphere20
        mrcalc -force ${j}_sphere20.nii 0.000001 -gt 1 0 -if ${j}_sphere20_bin.nii
#
#        fslmaths ${j}.nii -kernel sphere 25 -fmean ${j}_sphere25
#        fslmaths ${j}_sphere25 -bin ${j}_sphere25_bin
#
#        fslmaths ${j}.nii -kernel sphere 30 -fmean ${j}_sphere30
#        mrcalc -force ${j}_sphere30.nii 0.000001 -gt 1 0 -if ${j}_sphere30_bin.nii

    done
done
