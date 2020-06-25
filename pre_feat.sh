root=/media/labspace2/Projects/RMT-fMRI
#sid=(1 3 5 6 7 8 9 10 11 13 14 17 18 19 21 22 23 24 25 28 29 30)
sid=(2 12 16 20 27)
subjects=( "${sid[@]/#/NYU-}" )

for i in ${subjects[@]}; do

    echo $i
    proc_t1=$root/PROCESSING_T1_v1/$i
    mkdir $proc_t1
    listening=( `ls $root/NII/$i/*LISTEN*` )
    verbgen=( `ls $root/NII/$i/*VERB*` )
    sentcomp=( `ls $root/NII/$i/*SENT*` )
#listeningdn=( `ls $proc_t1/*LISTEN*dn*` )
#    verbgendn=( `ls $proc_t1/*VERB*dn*` )
#    sentcompdn=( `ls $proc_t1/*SENT*dn*` )
    finger=( `ls $root/NII/$i/*PIANO*` )
    fingert=( `ls $root/NII/$i/*ALT*` )

    echo ${finger[@]}
    cp $finger $proc_t1
    mpr=( `ls -S $root/NII/$i/*MPR* | head -n 1` )

    proc=$root/PROCESSING_T1_v1/$i
    mkdir $proc
    cp $listening $proc
    cp $verbgen $proc
    cp $sentcomp $proc
    cp $finger $proc
    cp $fingert $proc
#   cp $listeningdn $proc
#   cp $verbgendn $proc
#   cp $sentcompdn $proc
    cp $mpr $proc

    cd $proc_t1
    dwidenoise -force -noise $(basename "$listening" .nii)_noisemap.nii $listening $(basename "$listening" .nii)_dn.nii
    dwidenoise -force -noise $(basename "$verbgen" .nii)_noisemap.nii $verbgen $(basename "$verbgen" .nii)_dn.nii
    dwidenoise -force -noise $(basename "$sentcomp" .nii)_noisemap.nii $sentcomp $(basename "$sentcomp" .nii)_dn.nii
    rm $(basename "$mpr" .nii)_brain.nii
    bet $mpr $(basename "$mpr" .nii)_brain -f 0.4
    gunzip $(basename "$mpr" .nii)_brain.nii.gz
    dwidenoise -force -noise $(basename "$finger" .nii)_noisemap.nii $finger $(basename "$finger" .nii)_dn.nii
    dwidenoise -force -noise $(basename "$fingert" .nii)_noisemap.nii $finger $(basename "$finger" .nii)_dn.nii


done


