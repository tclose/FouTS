for f in */images/whole_brain.*.peaks.mif; do 
    d=`dirname $f`; 
    b=`basename $f`; 
    n=${b%.peaks.mif}; 
    n=${n#whole_brain.}; 
    mrmath $f 2 -pow - | mrmathaxis - sum 3 - | mrmath - -sqrt - | mrthreshold - $d/whole_brain.$n.sf_mask.mif --toppercent 0.1 --mask $d/whole_brain.$n.mask.mif -force
done
