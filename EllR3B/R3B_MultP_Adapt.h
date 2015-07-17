/*
 *  R3B_MultP_Adapt.h
 *  
 *
 *  Created by Daniel D'Orazio on 6/23/14.
 *  Copyright 2014 Columbia University. All rights reserved.
 *
 */
#include"nr3.h"
Int BAD_global=0;
Int fop_global=0;



double eCJ( double f, double u2, double e, double xx, double yy, double rr1, double rr2, double VVv){
	double a = 1.;
	double PsPot = ( ((a-u2)/rr1 + u2/rr2) + 1./a/a/a * (xx*xx + yy*yy) ) / ( 1.+e*cos(f) );
	//
	//for 
	//
	return 2.*PsPot - VVv*VVv;//  -  2*e*Intf;

}


double root0( double E , double e ){ 
	return( E  - e*sin(E) );
}

double root1( double E , double e ){ 
	return( 1. - e*cos(E) );
}



struct rhs_EllR3B {
	/*
	Doub q;
	rhs_R3B(Doub qq) : q(qq) {}
	void operator() (const Doub x, VecDoub_I &y, VecDoub_O &dydx){
		double u2 = q/(1.+q);
		double u1 = 1.0-u2;
		double a  = 1.0;
		double n  = 1./sqrt(a*a*a); 
		double r1 = sqrt((y[0] + u2)*(y[0] + u2) + y[1]*y[1]);
		double r2 = sqrt((y[0] - u1)*(y[0] - u1) + y[1]*y[1]);
		double dUdx1 = n*n*y[0] - (u2*(y[0] - u1))/(r2*r2*r2) - (u1*(u2 + y[0]))/(r1*r1*r1);
		double dUdx2 = n*n*y[1] - (u2*y[1])/(r2*r2*r2) - (u1*y[1])/(r1*r1*r1);
		
		double rr  = sqrt(y[0]*y[0] + y[1]*y[1]);
		double OmK = pow(rr,-1.5);
		double rho = 1.0;
		double nu = 0.00025;
		//double Fvisc_phi = -9./4.*rho*nu*(OmK - n)/(rr);
		
		double FK = -9./4.*rho*nu;
		//double Fvisc_xi = -Fvisc_phi/rr * (y[0]*sin(x) + y[1]*cos(x));
		//double Fvisc_et = Fvisc_phi/rr *  (y[0]*cos(x) - y[1]*sin(x));
		
		double FvscRx = FK/(rr*rr)*(y[3] - n*y[1]);
		double FvscRy = FK/(rr*rr)*(y[4] + n*y[0]);
		
		//double FvscRx = -Fvisc_phi*y[1]/rr;  //Fvisc_xi *cos(x) + Fvisc_et *sin(x);
		//double FvscRy = Fvisc_phi*y[0]/rr;   //Fvisc_et *cos(x) - Fvisc_xi *sin(x);
		
		double a1 = dUdx1 + 2.*n*y[4] + FvscRx;
		double a2 = dUdx2 - 2.*n*y[3] + FvscRy;
		double a3 = 0.0; ///work in the binary plane
		//acc.setAll(a1,a2,a3);
		
		
		
		// SET DERIVATIVES
		dydx[0] = y[3];        // dx/dt = vx
		dydx[1] = y[4];        // dy/dt = vy
		dydx[2] = y[5];		   // dz/dt = vz
		dydx[3] = a1;         // dvx/dt = ax
		dydx[4] = a2;         // dvy/dt = ay
		dydx[5] = a3;         // dvz/dt = az
	}
	*/
	
	
	
	
// Eliptical R3B
	Doub q;
	Doub ecc;
	rhs_EllR3B(Doub qq, Doub ecce) : q(qq) , ecc(ecce){}
	void operator() (const Doub x, VecDoub_I &y, VecDoub_O &dydx){
		double u2 = q/(1.+q);
		double u1 = 1.0-u2;
		double a  = 1.0;
		double n  = 1./sqrt(a*a*a); 
		
		
		double r1 = sqrt((y[0] + u2)*(y[0] + u2) + y[1]*y[1]);
		double r2 = sqrt((y[0] - u1)*(y[0] - u1) + y[1]*y[1]);
		double dUdx1 = n*n*y[0] - (u2*(y[0] - u1))/(r2*r2*r2) - (u1*(u2 + y[0]))/(r1*r1*r1);
		double dUdx2 = n*n*y[1] - (u2*y[1])/(r2*r2*r2) - (u1*y[1])/(r1*r1*r1);
		
		/*
		// Solve for true anamoly
		
		//minimize for Ecc anam
		//Newton-Rapheson to solve M = E - e*sin(E)
		double M = n * x;
		double E = M;  //Guess value for E is M.
		double TOL = 1e-8;
		double ff = root0( E , ecc ) - M; 
		while( fabs(ff) > TOL ){
			double dfdE = root1( E , ecc ); 
			double dE = -ff/dfdE;
			E += dE;
			ff = root0( E , ecc ) - M; 
		} 
		
		
		//get true anamoly f from EEa
		double fta = 2.* atan2(sqrt(1.+ecc)* tan(E/2.), sqrt(1.- ecc));
		
		
		double rsep = 1./(1. + ecc*cos(fta));
		*/
		
		double rsep = 1./(1. + ecc*cos(x));
		
		
		double a1 = dUdx1*rsep + 2.*n*y[4];
		double a2 = dUdx2*rsep - 2.*n*y[3];
		double a3 = 0.0; ///work in the binary plane
		//acc.setAll(a1,a2,a3);
		
		
		
		// SET DERIVATIVES
		dydx[0] = y[3];        // dx/dt = vx
		dydx[1] = y[4];        // dy/dt = vy
		dydx[2] = y[5];		   // dz/dt = vz
		dydx[3] = a1;         // dvx/dt = ax
		dydx[4] = a2;         // dvy/dt = ay
		dydx[5] = a3;         // dvz/dt = az
	}
	
	
	
	
};
