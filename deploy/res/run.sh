#!/bin/bash

sed -i -e "s/^GAUGE.*/GAUGE            dipole/g" test_R0_ci.inp
sed -i -e "s/^GAUGE.*/GAUGE            dipole/g" test_R0_hf.inp

sed -i -e "s/^FILE_OUT.*/FILE_OUT         perpendicular_dip_ci_R0.out/g" test_R0_ci.inp
sed -i -e "s/^FILE_OUT.*/FILE_OUT         perpendicular_dip_hf_R0.out/g" test_R0_hf.inp

sed -i -e "s/^POL_THETA.*/POL_THETA            1.570796/g" test_R0_ci.inp
sed -i -e "s/^POL_THETA.*/POL_THETA            1.570796/g" test_R0_hf.inp

./run_ci.sh
./run_hf.sh




sed -i -e "s/^GAUGE.*/GAUGE            velocity/g" test_R0_ci.inp
sed -i -e "s/^GAUGE.*/GAUGE            velocity/g" test_R0_hf.inp

sed -i -e "s/^FILE_OUT.*/FILE_OUT         perpendicular_vel_ci_R0.out/g" test_R0_ci.inp
sed -i -e "s/^FILE_OUT.*/FILE_OUT         perpendicular_vel_hf_R0.out/g" test_R0_hf.inp

sed -i -e "s/^POL_THETA.*/POL_THETA            1.570796/g" test_R0_ci.inp
sed -i -e "s/^POL_THETA.*/POL_THETA            1.570796/g" test_R0_hf.inp

./run_ci.sh
./run_hf.sh




sed -i -e "s/^GAUGE.*/GAUGE            dipole/g" test_R0_ci.inp
sed -i -e "s/^GAUGE.*/GAUGE            dipole/g" test_R0_hf.inp

sed -i -e "s/^FILE_OUT.*/FILE_OUT         parallel_dip_ci_R0.out/g" test_R0_ci.inp
sed -i -e "s/^FILE_OUT.*/FILE_OUT         parallel_dip_hf_R0.out/g" test_R0_hf.inp

sed -i -e "s/^POL_THETA.*/POL_THETA            0/g" test_R0_ci.inp
sed -i -e "s/^POL_THETA.*/POL_THETA            0/g" test_R0_hf.inp

./run_ci.sh
./run_hf.sh




sed -i -e "s/^GAUGE.*/GAUGE            velocity/g" test_R0_ci.inp
sed -i -e "s/^GAUGE.*/GAUGE            velocity/g" test_R0_hf.inp

sed -i -e "s/^FILE_OUT.*/FILE_OUT         parallel_vel_ci_R0.out/g" test_R0_ci.inp
sed -i -e "s/^FILE_OUT.*/FILE_OUT         parallel_vel_hf_R0.out/g" test_R0_hf.inp

sed -i -e "s/^POL_THETA.*/POL_THETA            0/g" test_R0_ci.inp
sed -i -e "s/^POL_THETA.*/POL_THETA            0/g" test_R0_hf.inp

./run_ci.sh
./run_hf.sh