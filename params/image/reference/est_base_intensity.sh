for f in */images/whole_brain.*0.no-iso.mif; do 
    d=`dirname $f`; 
    b=`basename $f`; 
    n=${b%.no-iso.mif}; 
    n=${n#whole_brain.}; 
    estimate_intensity $f $d/whole_brain.$n.peaks.mif $d/whole_brain.$n.sf_mask.mif $d/intens.$n.mif 
    image_percentile $d/intens.$n.mif $d/whole_brain.$n.sf_mask.mif --percentile 50 > $d/../intensity.$n.txt;
done
