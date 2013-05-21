for f in */images/whole_brain.*0.no-iso.mif; do 
    d=`dirname $f`; 
    b=`basename $f`; 
    n=${b%.no-iso.mif}; 
    n=${n#whole_brain.}; 
    estimate_intensity $f $d/intensities.mif > $d/../intens.$n.txt;
    #for i in `seq 5`; do 
        #echo "estimate_intensity $f --reference_tract $d/../ref_tract.$n.$i.tct > $d/../intens.$n.$i.txt;"
    #done;
done
