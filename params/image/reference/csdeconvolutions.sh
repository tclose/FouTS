for f in */images/whole_brain.*0.mif; do 
    d=`dirname $f`; 
    b=`basename $f`; 
    n=${b%.mif}; 
    n=${n#whole_brain.}; 
    if [ $n == '20' ]; then
        response=tensor_response.1000.txt
    else
        response=tensor_response.3000.txt
    fi
    csdeconv $f $response $d/whole_brain.$n.csd.mif --mask $d/whole_brain.$n.mask.mif
done
