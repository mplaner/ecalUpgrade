#!/bin/csh

# was: 
# export ROOTSYS=/usr/bin/root
# by hand:
echo $ROOTSYS
#setenv ROOTSYS /Users/drozdets/ROOT
#setenv PATH $ROOTSYS/bin:$PATH
setenv LD_LIBRARY_PATH $ROOTSYS/lib:$LD_LIBRARY_PATH
setenv DYLD_LIBRARY_PATH $ROOTSYS/lib:$DYLD_LIBRARY_PATH
echo $ROOTSYS

rm -f results.out2
touch results.out2
rm -rf PLOTS
mkdir -p PLOTS/C

foreach item (`cat input_files_list.txt`)
    echo $item | sed 's%\/% %g' | sed 's%\.% %' >> results.out2
	ln -sf $item input.root
	make clean; make test
	set dir=`echo $item | sed 's%\/% %g' | sed 's%\.% %' | awk '{print $4}'`
	cd PLOTS
	mkdir $dir
	mv *.* ../results.root $dir/.
	mv *.* ../*.out $dir/.
	mv C/ $dir/.
	cd $dir/
	mkdir DEBUG
	mv *_Debug_*.png DEBUG/.
	mv debugIndex.html DEBUG/index.html
	cd ..
	mkdir C
	cd ..
end

#cd PLOTS
#l PLOTS/*/results.out | awk '{print $9}' | sed 's%^%cat %' | sed 's%$% >> results.out%'
