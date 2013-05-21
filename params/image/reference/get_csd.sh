for f in */images/corpus_callosum.*0.mif; do 
    d=`dirname $f`; 
    b=`basename $f`; 
    n=${b%.mif}; 
    n=${n#corpus_callosum.}; 
    subtract_isotropic $f $d/corpus_callosum.$n.no-iso.mif
    cp $d/corpus_callosum.$n.no-iso.mif $d/fornix.$n.no-iso.mif
done
