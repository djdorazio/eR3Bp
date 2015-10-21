import numpy 
from numpy import *
import scipy 
from scipy import *
#import pylab as plt
#from pylab import *
import matplotlib
import matplotlib.pyplot as plt
import sys
import math
from math import pi as M_PI


fname1 = sys.argv[1]
Nstp = int(sys.argv[2])
#fname = "../dat/" + fname1  #OLDin pry 1 then!
fname = "../EllR3B/" + fname1 #OLDin pry 0 then
OLDin = 0  ## =1 uses old input style
  

Panel = 1    #1 Plot individual panels, 2 plot both begin and end together 3 plot all ouputs for movie
PlotICs = False
PlotCJ = False
stf = False




Colors = 0   #0 is old style red blue black green, 1 is red in blue out and inner part split between orange and green 
#2 is all black

###Load in t and x,y positions
if (OLDin==1):
	Vprof   = 4
	q       =  genfromtxt(fname, usecols=[0])[0] 
	ecc 	= 0.0
	np      = int( genfromtxt(fname, usecols=[0])[1] )
	DskSz   =  2.5
	NP      = int(np*np)
	BAD     = int( genfromtxt(fname, usecols=[0])[NP+Nstp] )
	tt    = genfromtxt(fname, usecols=[0], skip_header=2, skip_footer=1)
	xx    = genfromtxt(fname, usecols=[1], skip_header=2, skip_footer=1)
	yy    = genfromtxt(fname, usecols=[2], skip_header=2, skip_footer=1)
	CCJ   = genfromtxt(fname, usecols=[3], skip_header=2, skip_footer=1)
else:
	Vprof = genfromtxt(fname, usecols=[0])[0]    #1 = Kep plus bin quadrapole (MM08), 2 = r*sqrt(Om1 + Om2), 3 = Kep round each hole trans to Kep aroudn binary with vr and vphi comp
	q       =  genfromtxt(fname, usecols=[0])[1] 
	ecc     =  genfromtxt(fname, usecols=[0])[2] 
	np      = int( genfromtxt(fname, usecols=[0])[3] )
	DskSz   =  genfromtxt(fname, usecols=[0])[4]
	nout    =  genfromtxt(fname, usecols=[0])[5]
	NP      = int(np*np)
	BAD     = int( genfromtxt(fname, usecols=[0])[NP+Nstp] )

	tt    = genfromtxt(fname, usecols=[0], skip_header=6, skip_footer=1)
	xx    = genfromtxt(fname, usecols=[1], skip_header=6, skip_footer=1)
	yy    = genfromtxt(fname, usecols=[2], skip_header=6, skip_footer=1)
	CCJ   = genfromtxt(fname, usecols=[3], skip_header=6, skip_footer=1)



M = 1.
a = 1.
Mach=20.

u2 = q/(1.+q) 
u1 = a-u2

xbh1 = -u2
ybh1 = 0.0
xbh2 = u1
ybh2 = 0.0






print "LOADED"

t = zeros(Nstp)
#x = zeros([NP,2])
#y = zeros([NP,2])
#rr = zeros([NP,2])
#CJ = zeros([NP,2])

x = zeros([Nstp,NP])
y = zeros([Nstp,NP])
rr = zeros([Nstp,NP])
CJ = zeros([Nstp,NP])
for i in range(0, NP):
	for j in range(0, Nstp):
		if (i==0): t[j] = tt[j]
		#x[i][0]         = xx[i*Nstp]
		#y[i][0]         = yy[i*Nstp]
		x[j][i]         = xx[i*(Nstp)+j]
		y[j][i]         = yy[i*(Nstp)+j]
		
		#x[i][1]         = xx[i*Nstp+1]
		#y[i][1]         = yy[i*Nstp+1]
		#x[1][i]         = xx[i*Nstp+1]
		#y[1][i]         = yy[i*Nstp+1]
		
		#CJ[i][0]	    = CCJ[i*Nstp]
		#CJ[i][1]	    = CCJ[i*Nstp+1]
		#CJ[0][i]	    = CCJ[i*Nstp]
		CJ[j][i]	    = CCJ[i*(Nstp)+j]
		
		#rr[i][0]        = sqrt(x[i][0]*x[i][0] + y[i][0]*y[i][0]) #+1.e-16
		#rr[i][1]        = sqrt(x[i][1]*x[i][1] + y[i][1]*y[i][1]) #+1.e-16
		
		rr[j][i]        = sqrt(x[j][i]*x[j][i] + y[j][i]*y[j][i]) #+1.e-16
		#rr[1][i]        = sqrt(x[1][i]*x[1][i] + y[1][i]*y[1][i]) #+1.e-16
		
		#x[i][Nstp-1]    = xx[i*Nstp+Nstp-1]
		#y[i][Nstp-1]    = yy[i*Nstp+Nstp-1]


x0x = zeros(np)
y0y = zeros(np)
#vx0x = zeros(np)
#vy0y = zeros(np)

#xfx = zeros(np)
#yfy = zeros(np)
#vxfx = zeros(np)
#vyfy = zeros(np)

#CJ0 = zeros(np)
#CJf = zeros(np)
for i in range(0, np):
	x0x[i]       = xx[i*(np) * Nstp] 
	y0y[i]       = yy[i*Nstp] 
	#vx0x[i]      = vxx[i*(np) *Nstp] 
	#vy0y[i]      = vyy[i*Nstp] 	
	
	#xfx[i]       = xx[i*(np) *Nstp+1] 
	#yfy[i]       = yy[i*Nstp+1] 
	#vxfx[i]      = vxx[i*(np) *Nstp+1] 
	#vyfy[i]      = vyy[i*Nstp+1] 
			
	#CJ0[i]       = CJ[i*Nstp]
	#CJf[i]       = CJ[i*Nstp+1]





FCJ = zeros([np,np])	
CJ0 = zeros([np,np])
CJf = zeros([np,np])
#rr10 = zeros([np,np])
#rr20 = zeros([np,np])
#rr1f = zeros([np,np])
#rr2f = zeros([np,np])


#xy = zeros([np,np])
#rgrd = zeros([np,np])
for i in range(0, np):
	for j in range(0, np):
		#FCJ[i][j] = (CJ[i*np +j*Nstp + 1] - CJ[i*np +j*Nstp])/CJ[i*np +j*Nstp]
		#CJJ[i][j]  = CJ[i*np +j*Nstp]
		#xy[i][j]  =  xx[i*np +j*Nstp]
		#rgrd[j][i] = sqrt(x0x[i]*x0x[i] + y0y[j]*y0y[j])
		#rr10[j][i] = sqrt((x0x[i]-xbh1)*(x0x[i]-xbh1) + (y0y[j]-ybh1)*(y0y[j]-ybh1))
		#rr20[j][i] = sqrt((x0x[i]-xbh2)*(x0x[i]-xbh2) + (y0y[j]-ybh2)*(y0y[j]-ybh2))
		
		#rr1f[j][i] = sqrt((xfx[i]-xbh1)*(xfx[i]-xbh1) + (yfy[j]-ybh1)*(yfy[j]-ybh1))
		#rr2f[j][i] = sqrt((xfx[i]-xbh2)*(xfx[i]-xbh2) + (yfy[j]-ybh2)*(yfy[j]-ybh2))
		
		#VV0[j][i] = rgrd[j][i]*(sqrt(M/rgrd[j][i]/rgrd[j][i]/rgrd[j][i] * ( 1. + 3./4.*(a*a/rgrd[j][i]/rgrd[j][i]) * q/(1.+q)/(1.+q) )**2 - 0.5 * 1./(rgrd[j][i]*rgrd[j][i]*rgrd[j][i]) * 1./Mach/Mach) - sqrt(M/a/a/a))
		#VVf[j][i] = sqrt(Vv2[i* np + j][1])
		
		#CJ0[j][i] = CJ[i* np + j][0]
		#CJf[j][i] = CJ[i* np + j][1]
		
		CJ0[j][i] = CJ[0][i* np + j]
		CJf[j][i] = CJ[Nstp-1][i* np + j]

		#CJ0[j][i] = (x0x[i]*x0x[i] + y0y[j]*y0y[j]) + 2.*(u1/rr10[j][i] + u2/rr20[j][i]) - VV0[j][i]*VV0[j][i]
		#CJf[j][i] = (xfx[i]*xfx[i] + yfy[j]*yfy[j]) + 2.*(u1/rr1f[j][i] + u2/rr2f[j][i]) - VVf[j][i]*VVf[j][i]
		FCJ[j][i] = fabs( (CJf[j][i]-CJ0[j][i]) )/CJ0[j][i]

		#VV0[j][i] = sqrt(vx0x[i]*vx0x[i] + vy0y[j]*vy0y[j])
		#VVf[j][i] = sqrt(vxfx[i]*vxfx[i] + vyfy[j]*vyfy[j])
		#VV0[j][i] = rgrd[j][i]*(sqrt(M/rgrd[j][i]/rgrd[j][i]/rgrd[j][i] * ( 1. + 3./4.*(a*a/rgrd[j][i]/rgrd[j][i]) * q/(1.+q)/(1.+q) )**2 - 0.5 * 1./(rgrd[j][i]*rgrd[j][i]*rgrd[j][i]) * 1./Mach/Mach) - sqrt(M/a/a/a))

	
norb_end = t[len(t)-1]/(2.*pi)
norb = t/(2.*pi)
norbfmt0 = round(norb[0])
norbfmt1 = round(norb[Nstp-1])

norbfmt = zeros(Nstp)
for i in range(Nstp):
	norbfmt[i] = t[i]/(2.*pi)


## IN CASE L points not written in yet
#L3x = 0.
#xL2 = 1.
#CJcrit = 3.
#q=0.01



xL4 = 0.5 - u2
yL4 = 0.866

xL5 = 0.5 - u2
yL5 = -0.866

if (q == 0.0001):
	CJcrit = 3.00886
	xL1 = 0.968066265559007
	xL2 = 1.032424106260663
	xL3 = -1.000041662500362

elif (q == 0.001):
	CJcrit = 3.03859
	xL1 = 0.93130998854097
	xL2 = 1.069892950509135
	xL3 = -1.000416250362032

elif (q == 0.003):
	CJcrit = 3.12563
	xL1 = 0.90047572690107
	xL2 = 1.1001879704990
	xL3 = -1.00124625975948

elif (q == 0.005):
	CJcrit = 3.13882
	xL1 = 0.88137460788153
	xL2 = 1.11800131532167
	xL3 = -1.00207296177854

elif (q == 0.01):
	CJcrit = 3.15345
	xL1 = 0.84862409671786
	xL2 = 1.14631969632798
	xL3 = -1.00412535948298
	
elif (q == 0.02):
	CJcrit = 3.22338
	xL1 = 0.80494927823374
	xL2 = 1.17907038068448
	xL3 = -1.0081695201686

elif (q == 0.025):
	CJcrit = 3.25037894484595
	xL1 = 0.787774496082009
	xL2 = 1.19031198121494
	xL3 = -1.01016180163100

elif (q == 0.03):
	CJcrit = 3.27402
	xL1 = 0.772352667209511
	xL2 = 1.19962909692046
	xL3 = -1.01213455612159

elif (q == 0.04):
	CJcrit = 3.31393
	Xl1 = 0.745098076338014
	xL2 = 1.21435725389812
	xL3 = -1.01602247746864

elif (q == 0.05):
	CJcrit = 3.34669
	xL1 = 0.721126368324148
	xL2 = 1.22557033485370
	xL3 = -1.01983523256860

elif (q == 0.06):
	CJcrit = 3.37426
	xL1 = 0.699433317767269
	xL2 = 1.23440735170966
	xL3 = -1.02357471011260

elif (q == 0.075):
	CJcrit = 3.40837162649910
	xL1 = 0.669969945141183
	xL2 = 1.24457174166952
	xL3 = -1.02905052116936

elif (q == 0.08):
	CJcrit = 3.4182
	xL1 = 0.660795747894620
	xL2 = 1.24733364363021
	xL3 = -1.03084110108191
		
elif (q == 0.1):
	CJcrit = 3.45154
	xL1 = 0.626603496205417
	xL2 = 1.25608290849361
	xL3 = -1.03783564208400

elif (q == 0.3):
	CJcrit = 3.55966
	xL1 = 0.390097485958848
	xL2 = 1.26840220532166
	xL3 = -1.09537841243165

elif (q == 0.5):
	CJcrit = 3.54745813555200
	xL1 = 0.237418238185193
	xL2 = 1.24904738888032
	xL3 = -1.13636129399168

elif (q == 1.0):
	CJcrit = 3.4568
	xL1 = 0.0
	xL2 = 1.1984061445549
	xL3 = -1.1984061445549

else:
	xL2 = 0.0
	xL3 = 0.0
	print "Unknown xL2"
	print "Unknown L3(q)"
		


yL3 = 0.0


## CALCUALTE 2x the rotating frame potential to compare to Jacobi Constants
r1 = sqrt((x-xbh1)*(x-xbh1) + (y-ybh1)*(y-ybh1))
r2 = sqrt((x-xbh2)*(x-xbh2) + (y-ybh2)*(y-ybh2))


if (Vprof == 2):
	#Om = sqrt(M/rr[0]/rr[0]/rr[0] * (1. + 3./4.*(a*a/rr[0]/rr[0]) * q/(1.+q)/(1.+q))**2 - 0.5 * 1./(rr[0]*rr[0]*rr[0]) * 1./Mach/Mach) - sqrt(M/a/a/a)
	Om = sqrt(M/rr[0]/rr[0]/rr[0] * (1. + 3./4.*(a*a/rr[0]/rr[0]) * q/(1.+q)/(1.+q))**2 ) - sqrt(M/a/a/a)
	V = rr[0]*Om
	
if (Vprof == 2):
	Om = sqrt(u1/r1[0]/r1[0]/r1[0] + u2/r2[0]/r2[0]/r2[0]) - sqrt(M/a/a/a)
	V = rr[0]*Om


if (Vprof == 3):
### Vr and Vphi ICs for keplerian aroudn each hole transitioning to keplerian around binary (subtraction of bun ang freq is built in)
	vx0 = (a*(1. + q)*cos(arctan2(y[0],x[0]))*(-(sqrt ((2.*M_PI*q + (1. - 3.*q)*arctan((256.*pow(pow(-(a/(1. + q)) + x[0],2) +(y[0]*y[0]),4))/390625.))/((1. + q)*sqrt(pow(-(a/(1. + q)) + x[0],2) + (y[0]*y[0]))))/(pow(a,2) + pow(1. + q,2)*((x[0]*x[0]) + (y[0]*y[0]))  - 2.*a*(1. + q)*sqrt((x[0]*x[0]) + (y[0]*y[0]))*cos(arctan2(y[0],x[0])))) +  (q*sqrt((2.*M_PI + (-3 + q)*arctan((256.*pow(pow((a*q)/(1. + q) + x[0],2)  + (y[0]*y[0]),4))/390625.))/((1. + q)*sqrt(pow((a*q)/(1. + q) + x[0],2) +  (y[0]*y[0])))))/(pow(a,2)*pow(q,2) + pow(1. + q,2)*((x[0]*x[0]) +  (y[0]*y[0])) + 2.*a*q*(1. + q)*sqrt((x[0]*x[0]) + (y[0]*y[0]))*cos(arctan2(y[0],x[0]))))*sin(arctan2(y[0],x[0])))/sqrt(2.*M_PI) - sqrt((x[0]*x[0]) + (y[0]*y[0]))*(-sqrt(pow(a,-3)) + ((1. + q)*sqrt((2.*M_PI*q + (1. - 3.*q)*arctan((256.*pow(pow(-(a/(1. + q)) +  x[0],2) + (y[0]*y[0]),4))/390625.))/((1. + q)*sqrt(pow(-(a/(1. + q)) +  x[0],2) + (y[0]*y[0]))))*((1. + q)*sqrt((x[0]*x[0]) + (y[0]*y[0])) - a*cos(arctan2(y[0],x[0]))))/(sqrt(2.*M_PI)*(pow(a,2) + pow(1. +q,2)*((x[0]*x[0]) + (y[0]*y[0])) - 2.*a*(1. + q)*sqrt((x[0]*x[0]) +  (y[0]*y[0]))*cos(arctan2(y[0],x[0])))) + ((1. + q)*sqrt((2.*M_PI + (-3 +   q)*arctan((256.*pow(pow((a*q)/(1. + q) + x[0],2) + (y[0]*y[0]),4))/390625.))/((1. + q)*sqrt(pow((a*q)/(1. + q) + x[0],2) +  (y[0]*y[0]))))*((1. + q)*sqrt((x[0]*x[0]) + (y[0]*y[0])) +  a*q*cos(arctan2(y[0],x[0]))))/(sqrt(2.*M_PI)*(pow(a,2)*pow(q,2) + pow(1. + q,2)*((x[0]*x[0]) + (y[0]*y[0])) + 2.*a*q*(1. + q)*sqrt((x[0]*x[0]) + (y[0]*y[0]))*cos(arctan2(y[0],x[0])))))*sin(arctan2(y[0],x[0]))
	vy0 = sqrt((x[0]*x[0]) + (y[0]*y[0]))*cos(arctan2(y[0],x[0]))*(-sqrt(pow(a,-3)) +  ((1. + q)*sqrt((2.*M_PI*q + (1. - 3.*q)*arctan((256.*pow(pow(-(a/(1. +  q)) + x[0],2) + (y[0]*y[0]),4))/390625.))/((1. + q)*sqrt(pow(-(a/(1. + q)) + x[0],2) + (y[0]*y[0]))))*((1. + q)*sqrt((x[0]*x[0]) + (y[0]*y[0])) - a*cos(arctan2(y[0],x[0]))))/(sqrt(2.*M_PI)*(pow(a,2) + pow(1. +  q,2)*((x[0]*x[0]) + (y[0]*y[0])) - 2.*a*(1. + q)*sqrt((x[0]*x[0]) +   (y[0]*y[0]))*cos(arctan2(y[0],x[0])))) + ((1. + q)*sqrt((2.*M_PI + (-3 +  q)*arctan((256.*pow(pow((a*q)/(1. + q) + x[0],2) + (y[0]*y[0]),4))/390625.))/((1. + q)*sqrt(pow((a*q)/(1. + q) + x[0],2) +  (y[0]*y[0]))))*((1. + q)*sqrt((x[0]*x[0]) + (y[0]*y[0])) +  a*q*cos(arctan2(y[0],x[0]))))/(sqrt(2.*M_PI)*(pow(a,2)*pow(q,2) + pow(1. + q,2)*((x[0]*x[0]) + (y[0]*y[0])) + 2.*a*q*(1. + q)*sqrt((x[0]*x[0]) + (y[0]*y[0]))*cos(arctan2(y[0],x[0]))))) + (a*(1. + q)*(-(sqrt((2.*M_PI*q + (1. - 3.*q)*arctan((256.*pow(pow(-(a/(1. + q)) + x[0],2) +  (y[0]*y[0]),4))/390625.))/((1. + q)*sqrt(pow(-(a/(1. + q)) + x[0],2) + (y[0]*y[0]))))/(pow(a,2) + pow(1. + q,2)*((x[0]*x[0]) + (y[0]*y[0]))  - 2.*a*(1. + q)*sqrt((x[0]*x[0]) + (y[0]*y[0]))*cos(arctan2(y[0],x[0])))) +  (q*sqrt((2.*M_PI + (-3 + q)*arctan((256.*pow(pow((a*q)/(1. + q) + x[0],2)  + (y[0]*y[0]),4))/390625.))/((1. + q)*sqrt(pow((a*q)/(1. + q) + x[0],2) +   (y[0]*y[0])))))/(pow(a,2)*pow(q,2) + pow(1. + q,2)*((x[0]*x[0]) +  (y[0]*y[0])) + 2.*a*q*(1. + q)*sqrt((x[0]*x[0]) +  (y[0]*y[0]))*cos(arctan2(y[0],x[0]))))*pow(sin(arctan2(y[0],x[0])),2))/sqrt(2.*M_PI)
	V = sqrt(vx0*vx0 + vy0*vy0)
	
if (Vprof == 4):
	V = sqrt(u1/r1[0] + u2/r2[0]) - rr[0]*sqrt(M/a/a/a)
	
if (Vprof == 5):
	V = sqrt(u1/r1[0] + u2/r2[0])/sqrt(1.+ecc) - rr[0]*sqrt(M/a/a/a)
	
	
	
TU = (x[0]*x[0] + y[0]*y[0]) + 2.*(u1/r1[0] + u2/r2[0])

TUmV2 = TU - V*V



	
	
print "PLOTTING"
###COLOR THE DOTS
reds  = []
blues = []
blacks = []
greens = []
oranges = []

if (ecc!=0.0):
	for i in range(0,NP):
		if (rr[0][i] < 0.8): 
			reds.append(i)
		elif (rr[0][i] > 0.8 and rr[0][i] < 1.2):
			blacks.append(i)
		else:
			blues.append(i)
			
		ff = zeros(Nstp)
		xL45 = zeros(Nstp)
		yL4  = zeros(Nstp)
		for i in range(Nstp):
			ff[i]   = t[i]
			xL45[i] = 0.5*(1.-q)/(1.+q)  # t is f in ER3BpyL4  =
			yL4[i]  = sqrt(1./(1.+ecc*cos(ff[i]))**(4./3.) - 0.25)
		yL5  = -yL4 
else:
	xL45 = ones(Nstp)*0.5
	yL5 = ones(Nstp)*sqrt(3.)/2.
	yL4 = -yL5
	if (Colors==0):
		for i in range(0,NP):
			if (rr[0][i] < xL2 and TUmV2[i] > CJcrit):
				reds.append(i)
			elif (rr[0][i] > xL2 and TUmV2[i] > CJcrit):
				blues.append(i)
			elif(TU[i] <= CJcrit):
				blacks.append(i)
			else:
				greens.append(i)
		
	if (Colors==1):
		for i in range(0,NP):
			if (rr[0][i] < xL2 and TUmV2[i] > CJcrit):
				reds.append(i)
			elif (rr[0][i] > xL2 and TUmV2[i] > CJcrit):
				blues.append(i)
			if (rr[0][i] < xL2 and TUmV2[i] < CJcrit):
				greens.append(i)
			if (rr[0][i] > xL2 and TUmV2[i] < CJcrit):
				oranges.append(i)
	if (Colors==2):
		for i in range(0,NP):
			blacks.append(i)
			
		
##INTERPOLATE SCATTER POINTS TO CONTOUR PLOT

		
		
#zro = list(set(FCJ) & set(zeros(len(FCJ))))		
zro = where(FCJ == 0.0)
		
#----------------------------------------------------#
#----------------------------------------------------#
### PLOTTING ###
#----------------------------------------------------#
#----------------------------------------------------#
if (PlotCJ == True):
	#figure(figsize=[9,3])
	#figure()
	plt.figure(figsize=[4.5,4])

#	subplot(131)
#	title("CJ0 $q=%g$   $t = %i t_{orb}$" %(q, norbfmt1))

#	CT2 = contourf(x0x,y0y,log10(fabs(CJ0)))
	
#	CB2 = colorbar(CT2)
	
#	#scatter(xbh1, ybh1, color='black', s=40)
#	#scatter(xbh2, ybh2, color='blue', s=40)
#	scatter(L3x, L3y, color='blue', marker='x', s=40)
#	scatter(L4x, L4y, color='blue', marker='x', s=40)
#	scatter(L5x, L5y, color='blue', marker='x', s=40)
	
#	xlabel("$x/a$")
#	ylabel("$y/a$")
	
	
		
#	subplot(132)
#	title("CJF $q=%g$   $t = %i t_{orb}$" %(q, norbfmt1))
#	CT3 = contourf(x0x,y0y,log10(fabs(CJf)))
	
#	CB3 = colorbar(CT3)
	
#	xlabel("$x/a$")
#	ylabel("$y/a$")

		
#	subplot(133)
	plt.title("$\log{[\Delta C_J/C_J]}$ $q=%g$   $t = %g t_{orb}$" %(q, norbfmt1))
	CT1 = plt.contourf(x0x, y0y, log10((FCJ)) )
	#CT1 = contourf(x0x,y0y,CJ0 )
 
	

	
	plt.xlim(-(DskSz+2.),DskSz+2.)
	plt.ylim(-(DskSz+2.),DskSz+2.)
	
	CB1 = plt.colorbar(CT1)
	
	plt.xlabel("$x/a$")
	plt.ylabel("$y/a$")
	#clim(-10.,-2.)
	

	if (stf):	
		#savefig("../../Rings_vs_Cavitys/figures/CJRES_q%g_Om_Np1e4_%inorb_500perorb.png" %(q,norbfmt1))
		#savefig("PLOTS/CJRES_q%g_Om_Np1e4_%inorb_Rk5AdptStep.png" %(q,norbfmt1))
		plt.savefig("Visc_Plots/CJRES_q%g_Visc0p00025_Vprof%g_Np1e4_%gnorb_Rk5AdptStep.png" %(q,Vprof,norbfmt1))
	else:
		plt.show() 
	plt.close("all")
#----------------------------------------------------#
#----------------------------------------------------#
#----------------------------------------------------#
		

if (Panel == 1 and PlotICs):
	plt.figure(figsize=[4.5,4])
	plt.title("$q=%g$  $e=%g$  $t = %g t_{orb}$" %(q, ecc, norbfmt0))
	for i in range(0,len(reds)):
		plt.scatter(x[0][reds[i]], y[0][reds[i]], s=1, color='red')#, marker=',')
	for i in range(0,len(blues)):
		plt.scatter(x[0][blues[i]], y[0][blues[i]], s=1, color='blue')#, marker=',')
	for i in range(0,len(blacks)):
		plt.scatter(x[0][blacks[i]], y[0][blacks[i]], s=1, color='black')
		#plt.scatter(x[0][blacks[i]], y[0][blacks[i]], s=1, color='DarkGreen')#, marker=',')
	for i in range(0,len(greens)):
		plt.scatter(x[0][greens[i]], y[0][greens[i]], s=1, color='LightGreen')#, marker=',')
	for i in range(0,len(oranges)):
		plt.scatter(x[0][oranges[i]], y[0][oranges[i]], s=1, color='orange')#, marker=',')
			

	plt.scatter(xbh1, ybh1, color='black', s=40)
	plt.scatter(xbh2, ybh2, color='blue', s=40)
	plt.scatter(xL3, yL3, color='black', marker='o', s=30)
	plt.scatter(xL1, 0.0, color='black', marker='o', s=30)
	plt.scatter(xL2, 0.0, color='black', marker='o', s=30)
	plt.scatter(xL45[1], yL4[1], color='black', marker='o', s=30)
	plt.scatter(xL45[1], yL5[1], color='black', marker='o', s=30)
	


	plt.xlim(-(DskSz+2.),DskSz+2.)
	plt.ylim(-(DskSz+2.),DskSz+2.)
	plt.xlabel("$x/a$")
	plt.ylabel("$y/a$")

	print "SUBPLOT 1 COMPLETE"

	if (stf):	
		#savefig("../../Rings_vs_Cavitys/figures/CJRES_q%g_Om_Np1e4_%inorb_500perorb.png" %(q,norbfmt0))
		#savefig("PLOTS/q%g_Om_Np1e4_%inorb__Rk5AdptStep.png" %(q,norbfmt0))
		plt.savefig("../EccPlots/q%g_ecc%g_norb%g_Vprof%g_Np1e4_Rk5AdptStep.png" %(q,ecc,norbfmt0,Vprof))
	else:
		plt.show() 
	plt.close("all")
else:
	print "Not Plotting ICs"
#----------------------------------------------------#

if (Panel == 1):
	plt.figure(figsize=[4.5,4])
	plt.title("$q=%g$  $e=%g$  $t = %g t_{orb}$" %(q, ecc, norbfmt1))
	for i in range(0,len(reds)):
		plt.scatter(x[Nstp-1][reds[i]], y[Nstp-1][reds[i]], s=2, color='red')#, marker=',')
	for i in range(0,len(blues)):
		plt.scatter(x[Nstp-1][blues[i]], y[Nstp-1][blues[i]], s=2, color='blue')#, marker=',')
	for i in range(0,len(blacks)):
		plt.scatter(x[Nstp-1][blacks[i]], y[Nstp-1][blacks[i]], s=2, color='black')
		#plt.scatter(x[Nstp-1][blacks[i]], y[Nstp-1][blacks[i]], s=1, color='DarkGreen')#, marker=',')
	for i in range(0,len(greens)):
		plt.scatter(x[Nstp-1][greens[i]], y[Nstp-1][greens[i]], s=2, color='LightGreen')#, marker=',')
	for i in range(0,len(oranges)):
		plt.scatter(x[Nstp-1][oranges[i]], y[Nstp-1][oranges[i]], s=2, color='orange')#, marker=',')


	plt.scatter(xbh1, ybh1, color='black', s=40)
	plt.scatter(xbh2, ybh2, color='blue', s=40)
	plt.scatter(xL3, yL3, color='black', marker='x', s=30)
	plt.scatter(xL1, 0.0, color='black', marker='x', s=30)
	plt.scatter(xL2, 0.0, color='black', marker='x', s=30)
	plt.scatter(xL45[1], yL4[1], color='black', marker='x', s=30)
	plt.scatter(xL45[1], yL5[1], color='black', marker='x', s=30)


	
	#plt.xlim(-(DskSz+2.),DskSz+2.)
	#plt.ylim(-(DskSz+2.),DskSz+2.)
	
	plt.xlim(-2.5, 2.5)
	plt.ylim(-2.5, 2.5)

	plt.xlabel("$x/a$")
	plt.ylabel("$y/a$")

	print "SUBPLOT 2 COMPLETE"


	if (stf):	
		#savefig("....//Rings_vs_Cavitys/figures/CJRES_q%g_Om_Np1e4_%inorb_500perorb.png" %(q,norbfmt1))
		#savefig("PLOTS/q%g_Om_Np1e4_%inorb__Rk5AdptStep.png" %(q,norbfmt1))
		plt.savefig("../EccPlots/q%g_ecc%g_norb%g_Vprof%g_Np4e4_r4_Rk5AdptStep_Zoom.png" %(q,ecc,norbfmt1,Vprof))
	else:
		plt.show() 
	plt.close("all")

#----------------------------------------------------#
#----------------------------------------------------#
#----------------------------------------------------#

if (Panel == 2):
	figure(figsize=[4.5,8.])


	title("$q=%g$   $t = %g t_{orb}$" %(q, norbfmt0))
	for i in range(0,len(reds)):
		scatter(x[0][reds[i]], y[0][reds[i]], s=1, color='red')
	for i in range(0,len(blues)):
		scatter(x[0][blues[i]], y[0][blues[i]], s=1, color='blue')
	for i in range(0,len(blacks)):
		scatter(x[0][blacks[i]], y[0][blacks[i]], s=1, color='black')
	for i in range(0,len(greens)):
		scatter(x[0][greens[i]], y[0][greens[i]], s=1, color='green')
	for i in range(0,len(oranges)):
		scatter(x[0][oranges[i]], y[0][oranges[i]], s=1, color='orange')

			

	plt.scatter(xbh1, ybh1, color='black', s=40)
	plt.scatter(xbh2, ybh2, color='blue', s=40)
	plt.scatter(xL3, yL3, color='black', marker='o', s=30)
	plt.scatter(xL1, 0.0, color='black', marker='o', s=30)
	plt.scatter(xL2, 0.0, color='black', marker='o', s=30)
	plt.scatter(xL45[1], yL4[1], color='black', marker='o', s=30)
	plt.scatter(xL45[1], yL5[1], color='black', marker='o', s=30)

	xlim(-(DskSz+2.),DskSz+2.)
	ylim(-(DskSz+2.),DskSz+2.)
	xlabel("$x/a$")
	ylabel("$y/a$")

	print "PLOT COMPLETE"

	

	subplots_adjust(right=0.92)
	subplots_adjust(left=0.16)
	subplots_adjust(top=0.95)
	subplots_adjust(bottom=0.06)
	subplots_adjust(hspace=0.23)

	if (stf):	
		#savefig("../../Rings_vs_Cavitys/figures/CJRES_q%g_Om_Np1e4_%inorb_500perorb.png" %(q,norbfmt1))
		#savefig("PLOTS/q%g_Om_Np1e4_%inorb__Rk5AdptStep.png" %(q,norbfmt1))
		savefig("Visc_Plots/q%g_Visc0p00025_Vprof%g_Np1e4_%inorb__Rk5AdptStep.png" %(q,Vprof,norbfmt1))
	else:
		show() 






if (Panel == 3):
  for ii in range(Nstp):
	plt.figure(figsize=[4.5,4])
	plt.title("$q=%g$  $e=%g$  $t = %g t_{orb}$" %(q, ecc, norbfmt[ii]))
	for i in range(0,len(reds)):
		plt.scatter(x[ii][reds[i]], y[ii][reds[i]], s=1, color='red')#, marker=',')
	for i in range(0,len(blues)):
		plt.scatter(x[ii][blues[i]], y[ii][blues[i]], s=1, color='blue')#, marker=',')
	for i in range(0,len(blacks)):
		#plt.scatter(x[0][blacks[i]], y[0][blacks[i]], s=1, color='black')
		plt.scatter(x[ii][blacks[i]], y[ii][blacks[i]], s=1, color='black')#, marker=',')
	for i in range(0,len(greens)):
		plt.scatter(x[ii][greens[i]], y[ii][greens[i]], s=1, color='DarkGreen')#, marker=',')
	for i in range(0,len(oranges)):
		plt.scatter(x[ii][oranges[i]], y[ii][oranges[i]], s=1, color='orange')#, marker=',')
			

	if (ecc==0.0):
		plt.scatter(xbh1, ybh1, color='black', s=40)
		plt.scatter(xbh2, ybh2, color='blue', s=40)
		plt.scatter(xL3, yL3, color='black', marker='o', s=30)
		plt.scatter(xL1, 0.0, color='black', marker='o', s=30)
		plt.scatter(xL2, 0.0, color='black', marker='o', s=30)
		plt.scatter(xL45[1], yL4[1], color='black', marker='o', s=30)
		plt.scatter(xL45[1], yL5[1], color='black', marker='o', s=30)
	else:
		plt.scatter(xL45[ii], yL4[ii], color='black', marker='o', s=30)
		plt.scatter(xL45[ii], yL5[ii], color='black', marker='o', s=30)
		plt.scatter(xbh1, ybh1, color='black', s=40)
		plt.scatter(xbh2, ybh2, color='black', s=40)
	
	plt.xlim(-(DskSz+2.),DskSz+2.)
	plt.ylim(-(DskSz+2.),DskSz+2.)
	
	plt.xlabel("$x/a$")
	plt.ylabel("$y/a$")

	print "PLOT %i COMPLETE" %ii

	if (stf):	
		#savefig("../../Rings_vs_Cavitys/figures/CJRES_q%g_Om_Np1e4_%inorb_500perorb.png" %(q,norbfmt0))
		#savefig("PLOTS/q%g_Om_Np1e4_%inorb__Rk5AdptStep.png" %(q,norbfmt0))
		plt.savefig("../EccPlots/%3i_q%g_ecc%g_norb%g_Vprof%g_Np1e4_Rk5AdptStep.png" %(ii,q,ecc,norbfmt1,Vprof))
	else:
		plt.show() 
	plt.close("all")








