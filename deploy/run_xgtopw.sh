#!/bin/bash

XGTOPW_PATH="/home/mateusz/workspace/gtopw/build/Release/"
PHOTO_MISC_PATH="/home/mateusz/workspace/photo_misc/build/Release/"
PHOTO_FIT_PATH="/home/mateusz/workspace/photo_fit/"

HERE=$(dirname $(readlink -f $0))
mkdir ints
INTS="$HERE/ints"
INPS="$HERE/inputs"

cd res
touch kfile.txt
./run_ci.sh -dump
cat log.out >> kfile.txt
./run_hf.sh -dump
cat log.out >> kfile.txt
rm log.out 
cd ..

mv res/kfile.txt .
cp kfile.txt $PHOTO_FIT_PATH
cd $PHOTO_FIT_PATH
./input_gen.py kfile.txt
rm kfile.txt
cd $HERE

$PHOTO_MISC_PATH./main_he -xgtopw $INPS/basis.inp $HERE kfile.txt
rm kfile.txt

mkdir logs_xgtopw

for f in he_k*
do
	echo $f
	sed -i "6 a \$PATH\n$INTS/\n\$END" $f
	cp $f $XGTOPW_PATH
	cd $XGTOPW_PATH
	./xgtopw $f >> log_$f.out
	mv log_$f.out $HERE/logs_xgtopw
	rm $f
	cd $HERE
	mv $f $INPS
done
