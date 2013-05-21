for f in */images/whole_brain.*0.mif; do 
    d=`dirname $f`; 
    b=`basename $f`; 
    n=${b%.mif}; 
    n=${n#whole_brain.}; 
    gen_brain_mask $f $d/whole_brain.$n.mask.mif
done
