#!/bin/bash

echo START


echo RUN q0p2_e0p0
./q0p2_e0p0 > q0p2_e0p0_N200_100orb_rout2p5.dat
echo PLOT q0p2_e0p0
python ../Analysis/Adaptive_PlotOrbits.py "q0p2_e0p0_N200_100orb_rout4.dat" 2

echo RUN q0p2_e0p1
./q0p2_e0p1 > q0p2_e0p1_N200_100orb_rout2p5.dat
echo PLOT q0p2_e0p1
python ../Analysis/Adaptive_PlotOrbits.py "q0p2_e0p1_N200_100orb_rout4.dat" 2

echo RUN q0p2_e0p2
./q0p2_e0p2 > q0p2_e0p2_N200_100orb_rout2p5.dat
echo PLOT q0p2_e0p2
python ../Analysis/Adaptive_PlotOrbits.py "q0p2_e0p2_N200_100orb_rout4.dat" 2


echo RUN q0p2_e0p3
./q0p2_e0p3 > q0p2_e0p3_N200_100orb_rout2p5.dat
echo PLOT q0p2_e0p3
python ../Analysis/Adaptive_PlotOrbits.py "q0p2_e0p3_N200_100orb_rout4.dat" 2