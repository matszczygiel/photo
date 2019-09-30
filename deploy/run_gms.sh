#!/bin/bash

GAMESS_PATH="/home/mateusz/gamessSTO"
PHOTO_MISC_PATH="/home/mateusz/workspace/photo_misc/build/Release/"


HERE=$(dirname $(readlink -f $0))
mkdir tmp_data
mkdir data
TMP="$HERE/tmp_data"
INPS="$HERE/inputs"
DATA="$HERE/data"

$PHOTO_MISC_PATH./main_h2 -gms $INPS/basis.inp $INPS

cp $INPS/h2p_R${1}.inp $GAMESS_PATH

cd $GAMESS_PATH
./rungms h2p_R${1}.inp 01 > log_i.out
rm h2p_R${1}.inp

mv c.dat h2p_R${1}_c.dat
mv e.dat h2p_R${1}_e.dat

mv h2p_R${1}_c.dat $DATA
mv h2p_R${1}_e.dat $DATA
mv log_i.out $TMP
cd $TMP
grep 'NUMBER OF CARTESIAN GAUSSIAN BASIS FUNCTIONS' log_i.out > basis_l.dat
cd ..
python3 prep/format_length.py
NUM=$(cat basis_l.dat)
mv basis_l.dat basis_R${1}.dat
mv basis_R${1}.dat $DATA

cp $INPS/h2_R${1}.inp $GAMESS_PATH
cd $GAMESS_PATH
./rungms h2_R${1}.inp > log_n.out
rm h2_R${1}.inp

mv c.dat h2_R${1}_c.dat
mv e.dat h2_R${1}_e.dat
mv ci.dat h2_R${1}_ci.dat

mv h2_R${1}_c.dat $DATA
mv h2_R${1}_e.dat $DATA
mv h2_R${1}_ci.dat $DATA
mv log_n.out $TMP

cd $TMP
grep 'STATE   1  ENERGY=' log_n.out > E_CI.dat
grep 'FINAL RHF ENERGY IS' log_n.out > E_HF.dat

cd ..
python3 prep/format_energy.py
ENER_CI=$(cat E_CI.dat)
ENER_HF=$(cat E_HF.dat)
mv E_CI.dat E_CI_R${1}.dat
mv E_HF.dat E_HF_R${1}.dat
mv E_CI_R${1}.dat $DATA
mv E_HF_R${1}.dat $DATA


TARGET_CI="$HERE/res/test_R${1}_ci.inp"
TARGET_HF="$HERE/res/test_R${1}_hf.inp"
sed -i -e "s/^ENERGY_I_STATE .*/ENERGY_I_STATE       ${ENER_CI}/g" $TARGET_CI
sed -i -e "s/^ENERGY_I_STATE .*/ENERGY_I_STATE       ${ENER_HF}/g" $TARGET_HF
sed -i -e "s/^NUMBER_GTO .*/NUMBER_GTO       ${NUM}/g" $TARGET_CI
sed -i -e "s/^NUMBER_GTO .*/NUMBER_GTO       ${NUM}/g" $TARGET_HF

sed -i -e "s|^PATH_IN .*|PATH_IN       ${HERE}/|g" $TARGET_CI
sed -i -e "s|^PATH_IN .*|PATH_IN       ${HERE}/|g" $TARGET_HF
sed -i -e "s|^PATH_OUT .*|PATH_OUT       ${HERE}/res/|g" $TARGET_CI
sed -i -e "s|^PATH_OUT .*|PATH_OUT       ${HERE}/res/|g" $TARGET_HF
