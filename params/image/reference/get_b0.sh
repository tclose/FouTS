for f in */images/whole_brain.*0.mif; do 
    d=`dirname $f`; 
    b=`basename $f`; 
    n=${b%.mif}; 
    n=${n#whole_brain.}; 
    dwi_extract $f - --bzero | mrmathaxis - sum 3 - | mrmath - 6 -divide $d/whole_brain.$n.b0_avg.mif --force
    image_percentile $d/whole_brain.$n.b0_avg.mif $d/whole_brain.$n.sf_mask.mif --percentile 50 > $d/../b0.$n.txt
done
