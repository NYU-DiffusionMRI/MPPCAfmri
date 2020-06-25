root=/media/labspace2/Projects/RMT-fMRI
#sid=(5 7 8 13 14 17 19 21 22 24 25 29 30)
#sid=(1 3 5 6 7 8 9 10 11 13 14 17 18 19 22 23 24 25 28 29 30)
sid=(1 3 5 6 7 8 9 10 11 13 14 17 18 19 21 22 23 24 25 28 29 30)
#sid=(1)

subjects=( "${sid[@]/#/NYU-}" )

for i in ${subjects[17]}; do
    echo $i
    proc=$root/PROCESSING_T1_v1/$i

    listening=( `ls $proc/*LISTEN*.nii | grep -v _dn` )
    verbgen=( `ls $proc/*VERB*.nii | grep -v _dn` )
    sentcomp=( `ls $proc/*SENT*.nii | grep -v _dn` )
    tapping=( `ls $proc/*PIANO*.nii | grep -v _dn | grep -v noisemap` )
    mpr=( `ls -S $proc/*brain*` )

    st=(`mrinfo -size $listening`)
    
    cd $proc

   if [ ${#tapping[@]} -eq 0 ]; then
        tapping=( `ls $proc/*INDEX_TAP*.nii | grep -v _dn | grep -v noisemap` )
   fi

   echo $tapping

        sed 's|'$root/PROCESSING/NYU-1/BRAIN_MAPPING_SENTENCECOMP_fMRI_SMS_s15'|'$listening'|g' /media/labspace2/Projects/RMT-fMRI/design_template_raw.fsf > $proc/design_listening_raw.fsf
        sed -i 's|'$root/PROCESSING/NYU-1/BRAINMAPPINGAX3DMPR_s72_brain'|'$mpr'|g' $proc/design_listening_raw.fsf
        sed -i 's|'TEMPDIRECTORY'|'${listening}.noproc.feat'|g' $proc/design_listening_raw.fsf
        sed -i 's|'160'|'${st[3]}'|g' $proc/design_listening_raw.fsf

        sed 's|'$root/PROCESSING/NYU-1/BRAIN_MAPPING_SENTENCECOMP_fMRI_SMS_s15'|'$verbgen'|g' /media/labspace2/Projects/RMT-fMRI/design_template_raw.fsf > $proc/design_verbgen_raw_ica.fsf
        sed -i 's|'$root/PROCESSING/NYU-1/BRAINMAPPINGAX3DMPR_s72_brain'|'$mpr'|g' $proc/design_verbgen_raw_ica.fsf
        sed -i 's|'TEMPDIRECTORY'|'${verbgen}.noprocica.feat'|g' $proc/design_verbgen_raw_ica.fsf
        sed -i 's|'160'|'${st[3]}'|g' $proc/design_verbgen_raw_ica.fsf

        sed 's|'$root/PROCESSING/NYU-1/BRAIN_MAPPING_SENTENCECOMP_fMRI_SMS_s15'|'$sentcomp'|g' /media/labspace2/Projects/RMT-fMRI/design_template_raw.fsf > $proc/design_sentcomp_raw.fsf
        sed -i 's|'$root/PROCESSING/NYU-1/BRAINMAPPINGAX3DMPR_s72_brain'|'$mpr'|g' $proc/design_sentcomp_raw.fsf
        sed -i 's|'TEMPDIRECTORY'|'${sentcomp}.noproc.feat'|g' $proc/design_sentcomp_raw.fsf
        sed -i 's|'160'|'${st[3]}'|g' $proc/design_sentcomp_raw.fsf

        sed 's|'$root/PROCESSING/NYU-1/BRAIN_MAPPING_SENTENCECOMP_fMRI_SMS_s15'|'$tapping'|g' /media/labspace2/Projects/RMT-fMRI/design_template_raw.fsf > $proc/design_piano_raw_ica.fsf
        sed -i 's|'$root/PROCESSING/NYU-1/BRAINMAPPINGAX3DMPR_s72_brain'|'$mpr'|g' $proc/design_piano_raw_ica.fsf
        sed -i 's|'TEMPDIRECTORY'|'${tapping}.noprocica.feat'|g' $proc/design_piano_raw_ica.fsf
        sed -i 's|'160'|'${st[3]}'|g' $proc/design_piano_raw_ica.fsf



          sed 's|'$root/PROCESSING/NYU-1/BRAIN_MAPPING_SENTENCECOMP_fMRI_SMS_s15'|'/media/labspace2/Projects/RMT-fMRI/PROCESSING_T1_v1/NYU-24/BRAIN_MAPPING_RIGHT_HAND_PIANO_s53_dncnn.nii'|g' /media/labspace2/Projects/RMT-fMRI/design_template_raw.fsf > $proc/design_piano_raw_dn3.fsf
        sed -i 's|'$root/PROCESSING/NYU-1/BRAINMAPPINGAX3DMPR_s72_brain'|'$mpr'|g' $proc/design_piano_raw_dn3.fsf
        sed -i 's|'TEMPDIRECTORY'|'${tapping}.noprocdncnn3.feat'|g' $proc/design_piano_raw_dn3.fsf
        sed -i 's|'160'|'${st[3]}'|g' $proc/design_piano_raw_dn3.fsf

        #feat $proc/design_listening_raw.fsf &
        #feat $proc/design_verbgen_raw.fsf &
        #feat $proc/design_sentcomp_raw.fsf &
        feat $proc/design_piano_raw_dn3.fsf &
        #~/Downloads/fix1.067/fix ${tapping}.noprocica.feat ~/Downloads/fix1.067/training_files/Standard.RData 20
        #~/Downloads/fix1.067/fix ${verbgen}.noprocica.feat ~/Downloads/fix1.067/training_files/Standard.RData 20


        wait
    #fi

done


