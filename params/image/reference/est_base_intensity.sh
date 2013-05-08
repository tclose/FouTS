for f in */images/corpus_callosum.*0.mif; do d=donald/images; b=corpus_callosum.150.mif; n=corpus_callosum.150; n=150; for i in 1
2
3
4
5; do estimate_intensity donald/images/corpus_callosum.150.mif --reference_tract donald/images/../ref_tract.150.1.tct > donald/images/../intens.150.1.txt; echo done donald/images/corpus_callosum.150.mif; done; done
