for f in */images/whole_brain.*.csd.mif; do 
    d=`dirname $f`; 
    b=`basename $f`; 
    n=${b%.csd.mif}; 
    n=${n#whole_brain.}; 
    find_SH_peaks $f $d/whole_brain.$n.peaks.mif --num 1 --mask $d/whole_brain.$n.mask.mif
done
