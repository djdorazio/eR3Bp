/*
 *  R3B_MultP_Adapt.cpp
 *  
 *
 *  Created by Daniel D'Orazio on 6/23/14.
 *  Copyright 2014 Columbia University. All rights reserved. <- LOL
 *
 */

#include "R3B_MultP_Adapt.h"
#include<iostream>
#include<fstream>
#include<math.h> //for sqrt and pow
#include"nr3.h"
#include"odeint.h"
#include"stepper.h"
#include"R3B_Adaptive.h" // the Dormand Price rk5 adaptive step size integrator

int main(){
	// Choose Velocity Profile
	Int Vprof = 4; // 1 = Om, 2 = AddV, 3= Under Constr (to be consistent AddV with vr and vphi)
	Int RandSeed = 1; //randomly seed ICs in a disk or annulus if 1, else uniform square grid
	cout << Vprof << endl;
	const Int NUM_P=100;
	double DskSze = 2.5;
	double DskMin = 0.0;
	
	const Int Norb = 1;
	const Int nvar=6;  // x, y, z, vx, vy, vz
	const Int nout=1;//32*Norb; //number of outputs if NUM_P=1 else output ICs and final step reset nout below - =1 means only records first and final step
	

	const Doub q = 0.1;   //mass ratio
	const Doub ecc = 0.0;  //eccentricity
	const Doub nu  = 0.0;   // coef of kin visc
	
	//unstable q
	//>0.040064205
	//=0.024898786
	//=0.013701202
	//=0.011034092
	
	
	const Doub a = 1.0;
	double u2 = q/(1.+q);
	double u1 = 1.0-u2;
	
	const Doub atol=1.0e-8, rtol=atol, hmin=0.0; //h1=0.01
	const Doub x1 = 0.0;
	const Doub x2 = 100.0*2.*M_PI; // END TIME PASS THIS TO R3B_Adaptive.h
	const Doub h1 = 0.0000000001;//(x2-x1)/(x2*100.);
	VecDoub ystart(nvar);
	

	if (NUM_P == 1){
		ystart[0] = 2.0; //x
		ystart[1] = 1./sqrt(ystart[0]);   //y
		ystart[2] = 0.0;   //z
		ystart[3] = 0.0;   //vx
		ystart[4] = 0.0;   //vy
		ystart[5] = 0.0;   //vz
		Output ouut(nout);  // output nout steps
		
		//rhs_R3B d(q);
		rhs_EllR3B d(q, ecc, nu);
		//rhs_VdP d(1.0e-3);
		
		//Odeint<StepperDopr5<rhs_R3B> > ode(ystart, x1, x2, atol, rtol, h1, hmin, ouut, d);
		Odeint<StepperDopr5<rhs_EllR3B> > ode(ystart, x1, x2, atol, rtol, h1, hmin, ouut, d);
		//Odeint<StepperDopr5<rhs_VdP> > ode(ystart, x1, x2, atol, rtol, h1, hmin, ouut, d);
		
		ode.integrate();
		
		// Print output
		cout << q << endl;
		cout << ecc << endl;
		cout << NUM_P << endl;
		cout << DskSze << endl;
		cout << nout << endl;
		for (Int i=0;i<ouut.count;i++){
			cout.precision(20);
			cout << ouut.xsave[i] << " " << ouut.ysave[0][i] << " " << ouut.ysave[1][i] <<endl;
		}
		
	}else {
		//Int nout=1;
		

		
		double x0;
		double y0;
		double rr;
		double vx0;
		double vy0;
		double Om;
		double Vv;
		double DX = DskSze*2./NUM_P; //0.05; //0.05
		double xL = -(DX*((double)NUM_P+1.))/2.; //+1 offetes from 0
		double yL = -(DX*((double)NUM_P+1.))/2.;
		
		double xbh1 = -u2;  // 1 - u2 - (-u2)= 1 = a/a or r/r, for circ or ecc
		double ybh1 = 0.0;
		double xbh2 = 1. - u2;
		double ybh2 = 0.0;
		
		double r1;
		double r2;
		
		
		
		
		double rr1[nout+1];
		double rr2[nout+1];
		double VVv[nout+1];
		double CJ[nout+1];
		
		cout << q << endl;
		cout << ecc << endl;
		cout << NUM_P << endl;
		cout << DskSze << endl;
		cout << nout << endl;
	//	double nc=10.; //for cored profiles
		//double Mach = 20.;
		//double OM1;
		//double OM2;
		//double OM;
		double rc=0.1;
		for (Int i=0;i<NUM_P;i++){
			for (Int j=0;j<NUM_P;j++){ //i++ or ++i?
				if (RandSeed==1){
					rr = DskSze + 1.;
					while (rr>DskSze or rr<DskMin) {
						//cout << rr << endl;
						x0 = rand()/(double)RAND_MAX * 2.*DskSze - DskSze;
						y0 = rand()/(double)RAND_MAX * 2.*DskSze - DskSze;
						rr  = sqrt(x0*x0 + y0*y0);
					}
					

				}else {
					x0  = xL + DX*(double)i;
					y0  = yL + DX*(double)j;
				}

				
				rr  = sqrt(x0*x0 + y0*y0);
				r1 = sqrt((x0-xbh1)*(x0-xbh1) + (y0-ybh1)*(y0-ybh1));
				r2 = sqrt((x0-xbh2)*(x0-xbh2) + (y0-ybh2)*(y0-ybh2));
				//
				//SET VELOICTY PROFILES
				if (Vprof == 1){
					if (rr > 0.05){
						//Om = sqrt(1./rr/rr/rr * pow( (1. + 3./4.*(a*a/rr/rr) * q/(1.+q)/(1.+q)), 2) - 0.5 * 1./(rr*rr*rr) * 1./Mach/Mach) - sqrt(1./a/a/a); 
						Om = sqrt(1./rr/rr/rr) *  (1. + 3./4.*(a*a/rr/rr) * q/(1.+q)/(1.+q)) -  sqrt(1./a/a/a);
					
						Vv = rr*Om;
						vx0 = -Vv * y0/rr; 
						vy0 = Vv * x0/rr;
						
						ystart[0] = x0; //x
						ystart[1] = y0;   //y
						ystart[2] = 0.0;   //z
						ystart[3] = vx0;   //vx
						ystart[4] = vy0;   //vy
						ystart[5] = 0.0;   //vz
					}else {
						ystart[0] = 2.0; //x
						ystart[1] = 0.0;   //y
						ystart[2] = 0.0;   //z
						ystart[3] = 0.0;   //vx
						ystart[4] = 1.0;   //vy
						ystart[5] = 0.0;   //vz
					}
				}
				if (Vprof == 2){
					if (r1>0.05 && r2>0.05){
						Om = sqrt( 1./(1.+q)/r1/r1/r1  +  1./(1.+1./q)/r2/r2/r2 ) - sqrt(1./a/a/a );
						
						Vv = rr*Om;
						vx0 = -Vv * y0/rr; 
						vy0 = Vv * x0/rr;
						
						ystart[0] = x0; //x
						ystart[1] = y0;   //y
						ystart[2] = 0.0;   //z
						ystart[3] = vx0;   //vx
						ystart[4] = vy0;   //vy
						ystart[5] = 0.0;   //vz
					}else {
						ystart[0] = 2.0; //x
						ystart[1] = 0.0;   //y
						ystart[2] = 0.0;   //z
						ystart[3] = 0.0;   //vx
						ystart[4] = 1.0;   //vy
						ystart[5] = 0.0;   //vz
					}
				//}else{
				//	throw("No velocity profile set!")
				}
				if (Vprof == 3){ //From Mathematica
					//below doesnt quite get KEp at large r
					//vx0 = (a*(1.0 + q)*cos(atan2(y0,x0))*(-(sqrt(q/((1.0 + q)*sqrt(pow(1.0/(1.0 + q) - x0,2) + (y0*y0))))/(pow(a,2) + pow(1.0 + q,2)*((x0*x0) +  (y0*y0)) - 2*a*(1.0 + q)*sqrt((x0*x0) +  (y0*y0))*cos(atan2(y0,x0)))) + (q*sqrt(1.0/((1.0 + q)*sqrt(pow(q/(1.0 + q) + x0,2) +  (y0*y0)))))/(pow(a,2)*pow(q,2) + pow(1.0 + q,2)*((x0*x0)  + (y0*y0)) + 2*a*q*(1.0 + q)*sqrt((x0*x0) +  (y0*y0))*cos(atan2(y0,x0)))) + sqrt((x0*x0) +  (y0*y0))*(sqrt(pow(a,-3)) - ((1.0 + q)*sqrt(q/((1.0 +  q)*sqrt(pow(1.0/(1.0 + q) - x0,2) + (y0*y0))))*((1.0 +  q)*sqrt((x0*x0) + (y0*y0)) -  a*cos(atan2(y0,x0))))/(pow(a,2) + pow(1.0 + q,2)*((x0*x0) + (y0*y0)) - 2*a*(1.0 + q)*sqrt((x0*x0) +  (y0*y0))*cos(atan2(y0,x0))) - ((1.0 + q)*sqrt(1.0/((1.0 +  q)*sqrt(pow(q/(1.0 + q) + x0,2) + (y0*y0))))*((1.0 + q)*sqrt((x0*x0) + (y0*y0)) + a*q*cos(atan2(y0,x0))))/(pow(a,2)*pow(q,2) + pow(1.0 + q,2)*((x0*x0) + (y0*y0)) + 2*a*q*(1.0 + q)*sqrt((x0*x0) + (y0*y0))*cos(atan2(y0,x0)))))*sin(atan2(y0,x0));
					
					//vy0 = sqrt((x0*x0) + (y0*y0))*cos(atan2(y0,x0))*(-sqrt(pow(a,-3)) + ((1.0 + q)*sqrt(q/((1.0 + q)*sqrt(pow(1.0/(1.0 + q) - x0,2) + (y0*y0))))*((1.0 + q)*sqrt((x0*x0) + (y0*y0)) - a*cos(atan2(y0,x0))))/(pow(a,2) + pow(1.0 + q,2)*((x0*x0) +  (y0*y0)) - 2*a*(1.0 + q)*sqrt((x0*x0) + (y0*y0))*cos(atan2(y0,x0))) + ((1.0 + q)*sqrt(1.0/((1.0 + q)*sqrt(pow(q/(1.0 + q) + x0,2) + (y0*y0))))*((1.0 +  q)*sqrt((x0*x0) + (y0*y0)) +  a*q*cos(atan2(y0,x0))))/(pow(a,2)*pow(q,2) + pow(1.0 +  q,2)*((x0*x0) + (y0*y0)) + 2*a*q*(1.0 + q)*sqrt((x0*x0) +   (y0*y0))*cos(atan2(y0,x0)))) + a*(1.0 + q)*(-(sqrt(q/((1.0 +  q)*sqrt(pow(1.0/(1.0 + q) - x0,2) + (y0*y0))))/(pow(a,2) +  pow(1.0 + q,2)*((x0*x0) + (y0*y0)) - 2*a*(1.0 +  q)*sqrt((x0*x0) + (y0*y0))*cos(atan2(y0,x0)))) +  (q*sqrt(1.0/((1.0 + q)*sqrt(pow(q/(1.0 + q) + x0,2) +  (y0*y0)))))/(pow(a,2)*pow(q,2) + pow(1.0 + q,2)*((x0*x0)  + (y0*y0)) + 2*a*q*(1.0 + q)*sqrt((x0*x0) +  (y0*y0))*cos(atan2(y0,x0))))*pow(sin(atan2(y0,x0)),2);
					
					//below approaches KEP at large r and Kep around each BH
					vx0 = (a*(1. + q)*cos(atan2(y0,x0))*(-(sqrt ((2.*M_PI*q + (1. - 3.*q)*atan((256.*pow(pow(-(a/(1. + q)) + x0,2) +(y0*y0),4))/390625.))/((1. + q)*sqrt(pow(-(a/(1. + q)) + x0,2) + (y0*y0))))/(pow(a,2) + pow(1. + q,2)*((x0*x0) + (y0*y0))  - 2.*a*(1. + q)*sqrt((x0*x0) + (y0*y0))*cos(atan2(y0,x0)))) +  (q*sqrt((2.*M_PI + (-3 + q)*atan((256.*pow(pow((a*q)/(1. + q) + x0,2)  + (y0*y0),4))/390625.))/((1. + q)*sqrt(pow((a*q)/(1. + q) + x0,2) +  (y0*y0)))))/(pow(a,2)*pow(q,2) + pow(1. + q,2)*((x0*x0) +  (y0*y0)) + 2.*a*q*(1. + q)*sqrt((x0*x0) + (y0*y0))*cos(atan2(y0,x0))))*sin(atan2(y0,x0)))/sqrt(2.*M_PI) - sqrt((x0*x0) + (y0*y0))*(-sqrt(pow(a,-3)) + ((1. + q)*sqrt((2.*M_PI*q + (1. - 3.*q)*atan((256.*pow(pow(-(a/(1. + q)) +  x0,2) + (y0*y0),4))/390625.))/((1. + q)*sqrt(pow(-(a/(1. + q)) +  x0,2) + (y0*y0))))*((1. + q)*sqrt((x0*x0) + (y0*y0)) - a*cos(atan2(y0,x0))))/(sqrt(2.*M_PI)*(pow(a,2) + pow(1. +q,2)*((x0*x0) + (y0*y0)) - 2.*a*(1. + q)*sqrt((x0*x0) +  (y0*y0))*cos(atan2(y0,x0)))) + ((1. + q)*sqrt((2.*M_PI + (-3 +   q)*atan((256.*pow(pow((a*q)/(1. + q) + x0,2) + (y0*y0),4))/390625.))/((1. + q)*sqrt(pow((a*q)/(1. + q) + x0,2) +  (y0*y0))))*((1. + q)*sqrt((x0*x0) + (y0*y0)) +  a*q*cos(atan2(y0,x0))))/(sqrt(2.*M_PI)*(pow(a,2)*pow(q,2) + pow(1. + q,2)*((x0*x0) + (y0*y0)) + 2.*a*q*(1. + q)*sqrt((x0*x0) + (y0*y0))*cos(atan2(y0,x0)))))*sin(atan2(y0,x0));
						
					
					vy0 = sqrt((x0*x0) + (y0*y0))*cos(atan2(y0,x0))*(-sqrt(pow(a,-3)) +  ((1. + q)*sqrt((2.*M_PI*q + (1. - 3.*q)*atan((256.*pow(pow(-(a/(1. +  q)) + x0,2) + (y0*y0),4))/390625.))/((1. + q)*sqrt(pow(-(a/(1. + q)) + x0,2) + (y0*y0))))*((1. + q)*sqrt((x0*x0) + (y0*y0)) - a*cos(atan2(y0,x0))))/(sqrt(2.*M_PI)*(pow(a,2) + pow(1. +  q,2)*((x0*x0) + (y0*y0)) - 2.*a*(1. + q)*sqrt((x0*x0) +   (y0*y0))*cos(atan2(y0,x0)))) + ((1. + q)*sqrt((2.*M_PI + (-3 +  q)*atan((256.*pow(pow((a*q)/(1. + q) + x0,2) + (y0*y0),4))/390625.))/((1. + q)*sqrt(pow((a*q)/(1. + q) + x0,2) +  (y0*y0))))*((1. + q)*sqrt((x0*x0) + (y0*y0)) +  a*q*cos(atan2(y0,x0))))/(sqrt(2.*M_PI)*(pow(a,2)*pow(q,2) + pow(1. + q,2)*((x0*x0) + (y0*y0)) + 2.*a*q*(1. + q)*sqrt((x0*x0) + (y0*y0))*cos(atan2(y0,x0))))) + (a*(1. + q)*(-(sqrt((2.*M_PI*q + (1. - 3.*q)*atan((256.*pow(pow(-(a/(1. + q)) + x0,2) +  (y0*y0),4))/390625.))/((1. + q)*sqrt(pow(-(a/(1. + q)) + x0,2) + (y0*y0))))/(pow(a,2) + pow(1. + q,2)*((x0*x0) + (y0*y0))  - 2.*a*(1. + q)*sqrt((x0*x0) + (y0*y0))*cos(atan2(y0,x0)))) +  (q*sqrt((2.*M_PI + (-3 + q)*atan((256.*pow(pow((a*q)/(1. + q) + x0,2)  + (y0*y0),4))/390625.))/((1. + q)*sqrt(pow((a*q)/(1. + q) + x0,2) +   (y0*y0)))))/(pow(a,2)*pow(q,2) + pow(1. + q,2)*((x0*x0) +  (y0*y0)) + 2.*a*q*(1. + q)*sqrt((x0*x0) +  (y0*y0))*cos(atan2(y0,x0))))*pow(sin(atan2(y0,x0)),2))/sqrt(2.*M_PI);
					
					ystart[0] = x0; //x
					ystart[1] = y0;   //y
					ystart[2] = 0.0;   //z
					ystart[3] = vx0;   //vx
					ystart[4] = vy0;   //vy
					ystart[5] = 0.0;
				}
//Circular Virial BINARY ICs
				if (Vprof == 4){
					if (r1>0.05 && r2>0.05){
						Vv = sqrt( 1./(1.+q)/r1  +  1./(1.+1./q)/r2 ) - rr*sqrt(1./a/a/a );
						vx0 = -Vv * y0/rr; 
						vy0 = Vv * x0/rr;
						
						ystart[0] = x0;    //x
						ystart[1] = y0;    //y
						ystart[2] = 0.0;   //z
						ystart[3] = vx0;   //vx
						ystart[4] = vy0;   //vy
						ystart[5] = 0.0;   //vz
					}else {
						ystart[0] = 10.0;   //x
						ystart[1] = 0.0;   //y
						ystart[2] = 0.0;   //z
						ystart[3] = 0.0;   //vx
						ystart[4] = 1.0/sqrt(ystart[0]) - rr*sqrt(1./a/a/a );   //vy
						ystart[5] = 0.0;   //vz
					}
				}
//ECCENTRIC Virial BINARY ICs - coordinates are nonuniformly rotating and pulsating!
				if (Vprof == 5){
					double f0 = 0.0;
					double dfdt_p = pow((1. +  ecc*x0/rr ), 2)/ pow((1.-ecc*ecc),1.5);//sqrt(1.-ecc*ecc)/(1.-ecc)/(1.-ecc);
					double VarR_p = (a*(1. - ecc*ecc)) / (1. + ecc*x0/rr);
					double dfdt_B = pow((1. +  ecc*cos(f0) ), 2)/ pow((1.-ecc*ecc),1.5);
					double VarR_B = (a*(1. - ecc*ecc)) / (1. + ecc*cos(f0));
					if (r1>0.05 && r2>0.05){
						double Vp = sqrt( 1./(1.+q)/r1  +  1./(1.+1./q)/r2 ) / sqrt(1. + ecc * x0/rr);
						double Vr = 1./sqrt(1.-ecc*ecc) * ecc*y0/rr;
						//
						//NORMALIZE and change to d/df
						// letting n = 1 
						//
						Vp = Vp - rr;//*dfdt_B/VarR_B;
						Vr = 0.0;//Vr - 1./sqrt(1.-ecc*ecc) * ecc*sin(f0); 
						vx0 = (Vr*x0/rr - Vp * y0/rr);///dfdt_B/VarR_B;
						vy0 = (Vr*y0/rr + Vp * x0/rr);///dfdt_B/VarR_B;
						
						ystart[0] = x0;    //x
						ystart[1] = y0;    //y
						ystart[2] = 0.0;   //z
						ystart[3] = vx0;   //vx
						ystart[4] = vy0;   //vy
						ystart[5] = 0.0;   //vz
					}else {
						ystart[0] = 10.0;   //x
						ystart[1] = 0.0;   //y
						ystart[2] = 0.0;   //z
						ystart[3] = 0.0;   //vx
						ystart[4] = (1.0/sqrt(ystart[0]) - ystart[0]*dfdt_B) / VarR_B / dfdt_B;   //vy
						ystart[5] = 0.0;   //vz
					}
				}
				
				
				//Cored Om
				//OM1 = ( pow(r1, (8. - 3/2)) + pow(rc,(8. - 3/2)) )/(pow(r1, 8.) + pow(rc,8.));
				//OM2 = ( pow(r2, (8. - 3/2)) + pow(rc,(8. - 3/2)) )/(pow(r2, 8.) + pow(rc,8.));
				//Om = ( pow(rr, (nc - 3/2)) + pow(rc,(nc - 3/2)) )/(pow(rr, nc) + pow(rc,nc)) - sqrt(1./a/a/a );

				// VELOCITY PROFILES SET ABOVE^
				

				
				Output ouut(nout);    //ouput length nout long
				
				//// Circ R3B
				//if (ecc = 0.0){
				//	rhs_R3B d(q);  // Give accelerations
				//}else{
				////Ellipt R3B
				rhs_EllR3B d(q, ecc, nu);
				//};

				fop_global = 0;
				//NOE INTEGRATE EQNS
				//Odeint<StepperDopr5<rhs_R3B> > ode(ystart, x1, x2, atol, rtol, h1, hmin, ouut, d);
				Odeint<StepperDopr5<rhs_EllR3B> > ode(ystart, x1, x2, atol, rtol, h1, hmin, ouut, d);
				
				ode.integrate();
				
				// Print output

				for (Int k=0;k<ouut.count;k++){
					rr1[k]=sqrt((ouut.ysave[0][k]-xbh1)*(ouut.ysave[0][k]-xbh1) + (ouut.ysave[1][k]-ybh1)*(ouut.ysave[1][k]-ybh1));
					rr2[k]=sqrt((ouut.ysave[0][k]-xbh2)*(ouut.ysave[0][k]-xbh2) + (ouut.ysave[1][k]-ybh2)*(ouut.ysave[1][k]-ybh2));
					VVv[k]=sqrt(ouut.ysave[3][k]*ouut.ysave[3][k] + ouut.ysave[4][k]*ouut.ysave[4][k]);
					//CJ[k]=2.*((a-u2)/rr1[k] + u2/rr2[k]) + 1./a/a/a * (ouut.ysave[0][k]*ouut.ysave[0][k] + ouut.ysave[1][k]*ouut.ysave[1][k]) - VVv[k]*VVv[k];
					CJ[k] = eCJ( ouut.xsave[k], u2, ecc, ouut.ysave[0][k], ouut.ysave[1][k], rr1[k], rr2[k], VVv[k] );
					// 0-x, 1-y, 2-z;  3-vx, 4-vy, 5-vz
					cout.precision(20);
					cout << ouut.xsave[k] << " " << ouut.ysave[0][k] << " " << ouut.ysave[1][k] << " " << ouut.ysave[3][k] << " " << ouut.ysave[4][k] << " " <<CJ[k] <<endl;
				}
				

			} //for j
		} //for i
		cout << BAD_global << endl;
	} //NUMP else
		
		
} //main
