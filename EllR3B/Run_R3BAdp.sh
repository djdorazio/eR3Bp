#!/bin/bash

echo START


echo RUN qen3_e0p0
#./qen3_e0p0 > qen3_e0p0_N200_100orb_rout4.dat
echo PLOT qen3_e0p0
python ../Analysis/Adaptive_PlotOrbits.py "qen3_e0p1_N200_100orb_rout4.dat" 2




echo RUN q0p01_e0p1
#./q0p01_e0p0 > q0p01_e0p0_N200_100orb_rout4.dat
echo PLOT q0p01_e0p1
python ../Analysis/Adaptive_PlotOrbits.py "q0p01_e0p1_N200_100orb_rout4.dat" 2


echo RUN q0p02_e0p0
#./q0p02_e0p0 > q0p02_e0p0_N200_100orb_rout4.dat
echo PLOT q0p02_e0p0
python ../Analysis/Adaptive_PlotOrbits.py "q0p02_e0p1_N200_100orb_rout4.dat" 2


echo RUN q0p03_e0p0
#./q0p03_e0p0 > q0p03_e0p0_N200_100orb_rout4.dat
echo PLOT q0p03_e0p0
python ../Analysis/Adaptive_PlotOrbits.py "q0p03_e0p1_N200_100orb_rout4.dat" 2

echo RUN q0p036_e0p0
#./q0p036_e0p0 > q0p036_e0p0_N200_100orb_rout4.dat
echo PLOT q0p036_e0p0
python ../Analysis/Adaptive_PlotOrbits.py "q0p036_e0p1_N200_100orb_rout4.dat" 2


echo RUN q0p04_e0p0
#./q0p04_e0p1 > q0p04_e0p1_N200_100orb_rout4.dat
echo PLOT q0p04_e0p0
python ../Analysis/Adaptive_PlotOrbits.py "q0p04_e0p1_N200_100orb_rout4.dat" 2


echo RUN q0p05_e0p0
#./q0p05_e0p0 > q0p05_e0p1_N200_100orb_rout4.dat
echo PLOT q0p05_e0p0
python ../Analysis/Adaptive_PlotOrbits.py "q0p05_e0p1_N200_100orb_rout4.dat" 2



echo RUN q0p1_e0p1
./q0p1_e0p1 > q0p1_e0p1_N200_100orb_rout4.dat
echo PLOT q0p1_e0p1
python ../Analysis/Adaptive_PlotOrbits.py "q0p1_e0p1_N200_100orb_rout4.dat" 2







