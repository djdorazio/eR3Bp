#!/bin/bash

echo START

echo RUN qne4
./q0p01_e0p1 > q0p01_e0p1_N200_100orb_rout4.dat
echo PLOT q0p01_e0p1
python ../Analysis/Adaptive_PlotOrbits.py "q0p01_e0p1_N200_100orb_rout4.dat" 2


echo RUN qne3
#./qne3 > qne3_Om_N1e4_100orb_Adstp.dat
echo PLOT qne3
python Adaptive_PlotOrbits.py "qne3_Om_N1e4_100orb_ADstp.dat"

echo RUN qne2
#./qne2 > qne2_Om_N1e4_100orb_Adstp.dat
echo PLOT qne2
python Adaptive_PlotOrbits.py "qne2_Om_N1e4_100orb_ADstp.dat"

echo RUN qne1
#./qne1 > qne1_Om_N1e4_100orb_Adstp.dat
echo PLOT qne1
python Adaptive_PlotOrbits.py "qne1_Om_N1e4_100orb_ADstp.dat"

echo RUN q0p3
#./q0p3 > q0p3_Om_N1e4_100orb_Adstp.dat
echo PLOT q0p3
python Adaptive_PlotOrbits.py "q0p3_Om_N1e4_100orb_ADstp.dat"

echo RUN q1
#./q1 > q1_Om_N1e4_100orb_Adstp.dat
echo PLOT q1
python Adaptive_PlotOrbits.py "q1_Om_N1e4_100orb_ADstp.dat"
