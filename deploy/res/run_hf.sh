
#!/bin/bash

rm log.out

PPATH="/home/mateusz/workspace/photo_h2/build/Release/"

for e in 21.5 22.5 
do
	echo $e
	sed -i -e "s/^PHOTON_EN.*/PHOTON_EN $e/g" test_R0_hf.inp
	$PPATH./photo test_R0_hf.inp $1 >> log.out
done


