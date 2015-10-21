/*
 *  odeint.h
 *  
 *
 *  Created by Daniel D'Orazio on 6/23/14.
 *  Copyright 2014 Columbia University. All rights reserved.
 *
 */
#include"nr3.h"

#define GLOBAL_CONST_VAR BAD_global;


//Output Object
struct Output{
	Int kmax;
	Int nvar;
	Int nsave; // number of saves per integration if < 0 save at each tstep
	bool dense;
	Int count;
	Doub x1,x2,xout,dxout;
	VecDoub xsave;
	MatDoub ysave;
	Output() : kmax(-1), dense(false), count(0) {} // default no output
	Output(const Int nsavee) : kmax(500), nsave(nsavee), count(0), xsave(kmax) {
		dense = nsave >0 ? true : false;
	}
	void init(const Int neqn, const Doub xlo, const Doub xhi){
		nvar = neqn;
		if (kmax == -1) return;
		ysave.resize(nvar,kmax);
		if (dense) {
			x1=xlo;
			x2=xhi;
			xout=x1;
			dxout=(x2-x1)/nsave;
		}
	}
	// resize storage arrays by factor of 2 - keep saved data
	void resize() {
		Int kold=kmax;
		kmax *= 2;
		VecDoub tempvec(xsave);
		xsave.resize(kmax);
		for (Int k=0; k<kold; k++)
			xsave[k]=tempvec[k];
		MatDoub tempmat(ysave);
		ysave.resize(nvar,kmax);
		for (Int i=0; i<nvar; i++)
			for (Int k=0; k<kold; k++)
				ysave[i][k]=tempmat[i][k];
	}
	//use dens_out in R3B_Adaptive.h
	template<class Stepper>
	void save_dense(Stepper &s, const Doub xout, const Doub h){
		if (count == kmax) resize();
		for (Int i=0;i<nvar;i++)
			ysave[i][count]=s.dense_out(i,xout,h);
		xsave[count++]=xout;
	}
	// save current values of indep var x and dep vars y
	void save(const Doub x, VecDoub_I &y){
		if (kmax <= 0) return;
		if (count == kmax) resize();
		for (Int i=0;i<nvar;i++)
			ysave[i][count]=y[i];
		xsave[count++]=xout;
	}
	template<class Stepper>
	void ouut(const Int nstp, const Doub x, VecDoub_I &y, Stepper &s, const Doub h){
		if (!dense)
			throw("Dense output is not set!");
		if (nstp == -1){
			save(x,y);
			xout += dxout;
		}else{
			while ( (x-xout)*(x2-x1) > 0.0 ){
				save_dense(s,xout,h);
				xout += dxout;
			}
		}
	}
};




template<class Stepper>
struct Odeint {
	static const Int MAXSTP=100000000;
	Doub EPS;
	Int nok; //good steps
	Int nbad; // bad steps
	Int nvar; // usually 6 variables x,y,z,vx,vy,vz
	Doub x1,x2,hmin;  
	bool dense;  // use dense output?
	VecDoub y, dydx;
	VecDoub &ystart;  //will be replaced by output of integrator
	Output &ouut;  //an object to print step output at desired time
	typename Stepper::Dtype &derivs;  // derivs for specific ODE
	Stepper s;
	Int nstp;
	Doub x,h;
	// x1, x2 = tstart, tend - atol, rtol = abs and rel err tolerance, h1 is guessed starting stepsize, hmin is min stepsize (can be 0)
	Odeint(VecDoub_IO &ystartt, const Doub xx1, const Doub xx2, const Doub atol, const Doub rtol, const Doub h1, const Doub hminn, Output &outt, typename Stepper::Dtype &derivss);
	void integrate();
};

template<class Stepper>
Odeint<Stepper>::Odeint(VecDoub_IO &ystartt, const Doub xx1, const Doub xx2, const Doub atol, const Doub rtol, const Doub h1, const Doub hminn, Output &outt, typename Stepper::Dtype &derivss): nvar(ystartt.size()),y(nvar), dydx(nvar), ystart(ystartt), x(xx1), nok(0), nbad(0), x1(xx1), x2(xx2), hmin(hminn), dense(outt.dense), ouut(outt), derivs(derivss), s(y,dydx,x,atol,rtol,dense){
	EPS=numeric_limits<Doub>::epsilon();
	h=SIGN(h1,x2-x1);
	for (Int i=0;i<nvar;i++) y[i]=ystart[i];
	ouut.init(s.neqn,x1,x2);
}

template<class Stepper>
void Odeint<Stepper>::integrate() {
	derivs(x,y,dydx);
	if (dense)
		ouut.ouut(-1,x,y,s,h);
	else
		ouut.save(x,y);
	
//	int fop = 0;
	for (nstp=0;nstp<MAXSTP;nstp++){
		if ( (x+h*1.0001-x2)*(x2-x1) > 0.0)
			h=x2-x;
		s.step(h,derivs);
		if (s.hdid == h) ++nok; else ++nbad;
		if (dense)
			ouut.ouut(nstp,x,y,s,s.hdid);
		else
			ouut.save(x,y);
		if ((x-x2)*(x2-x1) >= 0.0){
			for (Int i=0;i<nvar;i++) ystart[i]=y[i];
			if (ouut.kmax > 0 && abs(ouut.xsave[ouut.count-1]-x2)>100.0*abs(x2)*EPS)
				ouut.save(x,y); //make sure last step is saved
			return;
		}
		//ADDED BY DAN TO STOP INTEGRATION WHEN PARTICLE GOES FAR AWAY
		//if (sqrt(y[1]*y[1] + y[2]*y[2]) > 100.) {
		//	x=x2;
		//	nstp = MAXSTP;
		//	ouut.save(x,y);	
		//}
		
		
		if (abs(s.hnext) <= hmin) {
			x=x2;
			nstp = MAXSTP;
			// fop_global is just a switch, it starts as 0 and stays 0 as logn as nothing goes wrong
			if (fop_global == 0){
				BAD_global += 1;
				fop_global = 1;
			}
			 
			cout << " step size too SMALL in Odeint "; 
			ouut.save(x,y); // throw("step size too SMALL in Odeint");
			//return;
		}
		h=s.hnext;
	}
	x=x2;
	if (fop_global == 0){
		BAD_global += 1;
		fop_global = 1;
	}
	cout << " Too MANY STEPS in Odeint! ";
	ouut.save(x,y); //throw("Too MANY STEPS in Odeint!");
	//return;
}


			
	
		


