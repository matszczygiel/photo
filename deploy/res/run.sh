#!/bin/bash

sed -i -e "s/^GAUGE.*/GAUGE            dipole/g" test_ci.inp
sed -i -e "s/^GAUGE.*/GAUGE            dipole/g" test_hf.inp

sed -i -e "s/^FILE_OUT.*/FILE_OUT         cNbN_dip_ci.out/g" test_ci.inp
sed -i -e "s/^FILE_OUT.*/FILE_OUT         cNbN_dip_hf.out/g" test_hf.inp

./run_ci.sh
./run_hf.sh

sed -i -e "s/^GAUGE.*/GAUGE            velocity/g" test_ci.inp
sed -i -e "s/^GAUGE.*/GAUGE            velocity/g" test_hf.inp

sed -i -e "s/^FILE_OUT.*/FILE_OUT         cNbN_vel_ci.out/g" test_ci.inp
sed -i -e "s/^FILE_OUT.*/FILE_OUT         cNbN_vel_hf.out/g" test_hf.inp

./run_ci.sh
./run_hf.sh
