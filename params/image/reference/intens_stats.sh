for f in */images/whole_brain.*0.no-iso.mif; do 
    d=`dirname $f`; 
    b=`basename $f`; 
    n=${b%.no-iso.mif}; 
    n=${n#whole_brain.}; 
    mrstats $d/intens.$n.mif --mask $d/whole_brain.$n.sf_mask.mif --dump $d/../intens_stats.$n.txt
done
