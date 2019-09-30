#!/bin/bash

GAMESS_PATH="/home/mateusz/gamessSTO"
PHOTO_MISC_PATH="/home/mateusz/workspace/photo_misc/build/Release/"


HERE=$(dirname $(readlink -f $0))
mkdir tmp_data
mkdir data
TMP="$HERE/tmp_data"
INPS="$HERE/inputs"
DATA="$HERE/data"

$PHOTO_MISC_PATH./main_he -gms $INPS/basis.inp $INPS

cp $INPS/hep.inp $GAMESS_PATH

cd $GAMESS_PATH
./rungms hep.inp > log_i.out
rm hep.inp

mv c.dat Hep_c.dat
mv e.dat Hep_e.dat

mv Hep_c.dat $DATA
mv Hep_e.dat $DATA
mv log_i.out $TMP

cd $TMP
grep 'NUMBER OF CARTESIAN GAUSSIAN BASIS FUNCTIONS' log_i.out > basis_l.dat
cd ..

python3 prep/format_length.py
NUM=$(cat basis_l.dat)
mv basis_l.dat basis.dat
mv basis.dat $DATA

cp $INPS/he.inp $GAMESS_PATH
cd $GAMESS_PATH
./rungms he.inp > log_n.out
rm he.inp

mv c.dat He_c.dat
mv e.dat He_e.dat
mv ci.dat He_ci.dat

mv He_c.dat  $DATA
mv He_e.dat  $DATA
mv He_ci.dat $DATA
mv log_n.out $TMP

cd $TMP
grep 'STATE   1  ENERGY=' log_n.out > E_CI.dat
grep 'FINAL RHF ENERGY IS' log_n.out > E_HF.dat

cd ..
python3 prep/format_energy.py
ENER_CI=$(cat E_CI.dat)
ENER_HF=$(cat E_HF.dat)
mv E_CI.dat $DATA
mv E_HF.dat $DATA


TARGET_CI="$HERE/res/test_ci.inp"
TARGET_HF="$HERE/res/test_hf.inp"

sed -i -e "s/^ENERGY_I_STATE .*/ENERGY_I_STATE       ${ENER_CI}/g" $TARGET_CI
sed -i -e "s/^ENERGY_I_STATE .*/ENERGY_I_STATE       ${ENER_HF}/g" $TARGET_HF
sed -i -e "s/^NUMBER_GTO .*/NUMBER_GTO       ${NUM}/g" $TARGET_CI
sed -i -e "s/^NUMBER_GTO .*/NUMBER_GTO       ${NUM}/g" $TARGET_HF

sed -i -e "s|^PATH_IN .*|PATH_IN       ${HERE}/|g" $TARGET_CI
sed -i -e "s|^PATH_IN .*|PATH_IN       ${HERE}/|g" $TARGET_HF

mkdir ${HERE}/res/outps
sed -i -e "s|^PATH_OUT .*|PATH_OUT       ${HERE}/res/outps/|g" $TARGET_CI
sed -i -e "s|^PATH_OUT .*|PATH_OUT       ${HERE}/res/outps/|g" $TARGET_HF
