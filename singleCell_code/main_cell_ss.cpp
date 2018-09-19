/*
 *  O'Hara-Rudy Human model with RK method
 *  Created by Mao-Tsuen Jeng on 6/13/11.
 *
 *	Copyright 2011 __Colleen Clancy Lab__. All rights reserved.
 *  Simulated Single cell
 */




#include <fstream>
#include <iostream>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <omp.h>
using namespace std;




const double BCL = 1000;

const int stimuli = 2;


const double base_dt = 0.005;
const int fold = 1.0/base_dt;



typedef struct cell {
	double y0[41];
	double highestV, lowestV, Apd90;
    double highestdvdt, V90;
    double DI, t1, t2;
	double  highestV2, lowestV2, t3, t4, s2APD;
	double dvdt, v_new, v, dvdt2;
	double currents[20];
    int celltype;
    double I_inj;
    double coef_PKA;
    
} Cell;

#include "ord_model_betaStim.h"
#include "integrate_rk2rk1.h"

typedef struct simState {
	double t, tt;
	int tstep, counter, beat;
	Cell cellData;
} SimState;

SimState theState;
SimState *S = &theState;

double Calcu_I_Total(Cell *theCell, double dt );


int main () {
	
	double p = 0;//par;  // Parameter array for passing nondefault conditions
	
	//// Initial conditions
	double v=-87.;
	double nai=7.;
	double nass=nai;
	double ki=145.;
	double kss=ki;
	double cai=1.0e-4;
	double cass=cai;
	double cansr=1.2;
	double cajsr=cansr;
	double m=0.;
	double hf=1.;
	double hs=1.;
	double j=1.;
	double hsp=1.;
	double jp=1.;
	double mL=0.;
	double hL=1.;
	double hLp=1.;
	double a=0.;
	double iF=1.;
	double iS=1.;
	double ap=0.;
	double iFp=1.;
	double iSp=1.;
	double d=0.;
	double ff=1.;
	double fs=1.;
	double fcaf=1.;
	double fcas=1.;
	double jca=1.;
	double nca=0.;
	double ffp=1.;
	double fcafp=1.;
	double xrf=0.;
	double xrs=0.;
	double xs1=0.;
	double xs2=0.;
	double xk1=1.;
	double Jrelnp=0.;
	double Jrelp=0.;
	double CaMKt=0.;
	// %X0 is the vector for initial sconditions for state variables
	//	X0=[v nai nass ki kss cai cass cansr cajsr m hf hs j hsp jp mL hL hLp a iF iS ap iFp iSp d ff fs fcaf fcas jca nca ffp fcafp xrf xrs xs1 xs2 xk1 Jrelnp Jrelp CaMKt]';
    
    int sy0 = 41;
    
    double y0[] = { v, nai, nass, ki, kss,
        cai, cass, cansr, cajsr, m,
        hf, hs, j, hsp, jp,
        mL, hL, hLp, a, iF,
        iS, ap, iFp, iSp, d,
        ff, fs, fcaf, fcas, jca,
        nca, ffp, fcafp, xrf, xrs,
        xs1, xs2, xk1, Jrelnp, Jrelp,
        CaMKt };
	
	
	char name[30];
	FILE *output;
	FILE *output_m;
	FILE *output1;
	FILE *output2;
	FILE *output3;
	int ll, ww;
	
	double V, dvdt, dvdt2;
	Cell *theCell;
	double dt;
	double I_Total;
	double dv;
	double Vm = -86;
	double cycle_length;
	
	double I_inj;
	int done = 0;
	const char *stateFileName = "model_INPUTSteadyStateControl.InitialsENDO"; //FOR initial conditions from steady-state
	const char *fromSingle = "SaveFinal.Cell"; //For final simulated State output name
	
	double waitTime;
	
	time_t startTime;
	time_t previousTime;
	
	const time_t timeSave = 0.1*60*60;
	startTime = time(NULL);
	previousTime = startTime;
    
    //OUTPUT results files
	output_m = fopen("model_OUTPUT_apd90_control1HzENDO.txt", "w");
	output=fopen("model_OUTPUT_ap_gates_control1HzENDO.txt", "w");
    output2=fopen("model_OUTPUT_currents_control1HzENDO.txt", "w");

    //READ in initial conditions from steady-state
	FILE *fp = fopen(stateFileName, "r");
    
	if (fp == NULL) {
		cout << "New start" << endl;
		
		S->counter=0;
		
		
		theCell = &(S->cellData);
		
		for ( int i = 0; i < sy0 ; i++ ) {
			theCell->y0[i] = y0[i];
		}
		
		theCell->Apd90 = 0;
		theCell->DI = 0;
		theCell->s2APD = 0;
		theCell->highestV = -86;
		theCell->lowestV = -86;
		theCell->highestV2 = -86;
		theCell->lowestV2 = -86;
		theCell->t1 = 1E10;
		theCell->t2 = 1E10;
		theCell->t3 = 1E10;
		theCell->t4 = 1E10;
		theCell->dvdt = 0;
		theCell->dvdt2 = 0;
		theCell->v_new = theCell->y0[0];
		theCell->I_inj = 0;
        theCell->celltype = 0; // %endo = 0, epi = 1, M = 2
        theCell->coef_PKA = 0; // beta-stimulation input range from 0 to 1 FOR linear change
        //For Control case, theCell->coef_PKA = 0 (no beta-stimulation)
        
		S->beat = 1;
		S->t = 0.0;
		S->tt = 0.0;
		S->tstep = 0;
	} else {
		fread(S, sizeof(SimState), 1, fp);
		cout << "restarting from the steady-state" << endl;
        theCell = &(S->cellData);

		fclose(fp);
	}

    cout<<"Beta-stimulation: "<< theCell->coef_PKA << endl;
    cout<<"Cell type: "<< theCell->celltype << " (endo = 0, epi = 1, M-cell = 2)" << endl;

	cycle_length = BCL;
    dt = base_dt;
    
	while (!done) {
		
      
		
		theCell->v = theCell->v_new;
		
		if ( S->t < 0.5  ) {
			theCell->I_inj = - 80;
		} else {
			theCell->I_inj=0.0;
		}
		
		V = theCell->v;
		
		I_Total = Calcu_I_Total( theCell, dt );
		
		dv = dt * ( -I_Total );
		
		dvdt = dv / dt;
		    
		
		theCell->v_new = theCell->v + dv;
        S->t += dt;
		S->tt += dt;
		S->counter += 1;

        if ( theCell->highestdvdt < dvdt ) {
           
            theCell->highestdvdt = dvdt;
            theCell->t1 = S->t;
        }
        if ( (S->t < 10) & (theCell->highestV <= theCell->v_new) ) {
            theCell->highestV = theCell->v_new;
            theCell->V90 = theCell->highestV - 0.9 * ( theCell->highestV - theCell->lowestV );
        } else if ( theCell->v_new > theCell->V90 ) {
            theCell->t2 = S->t;
            theCell->Apd90 = S->t - theCell->t1;
        }
        
        
		
		if(  ( S->counter % (fold*1) == 0 ) ) {
            
            int id;
            
			fprintf ( output, "%10.3f\t",
					 S->tt );
            
            for ( id = 0; id < 41; id ++ ) {
                fprintf ( output, "%8.5f\t",
                         S->cellData.y0[id] );
                
                // double y0[] = {
                //    v, nai, nass, ki, kss,
                //    cai, cass, cansr, cajsr, m,
                //    hf, hs, j, hsp, jp,
                //    mL, hL, hLp, a, iF,
                //    iS, ap, iFp, iSp, d,
                //    ff, fs, fcaf, fcas, jca,
                //    nca, ffp, fcafp, xrf, xrs,
                //    xs1, xs2, xk1, Jrelnp, Jrelp,
                //    CaMKt };
            }
            
            fprintf ( output, "\n" );
            
            fprintf ( output2, "%10.3f\t",
                     S->tt );
            
            for ( id = 0; id < 11; id ++ ) {
                fprintf ( output2, "%8.5f\t",
                         S->cellData.currents[id] );
                // Order of currents:
                // INa, ICaL, IKs, IKr, INaK
                // IKb, Ito, INaCa_i, INaCa_ss, IK1
                // INaK
            }
            fprintf ( output2, "\n" );
		}
		
		

		
		
		
		if(S->t >= cycle_length) {
            cout << "APD90 = " << theCell->Apd90 << endl;
            fprintf ( output_m, "%8.5f\t",
                     theCell->Apd90 );
			S->t = 0;
			S->beat++;
            theCell->lowestV = theCell->v_new;
            theCell->highestV = theCell->v_new;
			cout << "now the beat is: " << S->beat << endl;
		}
		
		if (S->beat > stimuli) {
			S->beat = 1;
			done = 1;
		}
		
		
	} // end while
	
	
	
	FILE *fp2 = fopen(fromSingle, "w");
	
	fwrite(S, sizeof(SimState), 1, fp2);
	fclose(fp2);
	
	fclose ( output );
	fclose ( output2 );
	fclose ( output_m );
    
	return 0;
}


double Calcu_I_Total( Cell *theCell, double dt ) {
	
	double tspan[2];
	double tspan_total[] = { 0, dt }; // ms
	double partition = dt;  // ms
	int parts_total = floor( 0.5 + (tspan_total[1] - tspan_total[0]) / partition );
	int part;
	int sy = 1;
	int sy0 = 41;
	
	double y[sy][sy0], t1[sy];
	double * y1[sy];
	for( int i = 0; i < sy; i++ ) {
		y1[i] = y[i];
	}
	
	time_t t_start, t_stop;
	t_start = time(NULL);
	
	
	double V = theCell->y0[0];
	
	for( part = 0; part < parts_total; part++ ) {
		
		
		tspan[0] = partition * part ;
		tspan[1] = tspan[0] + partition ;
        
		integrate_rk2rk1( ord_model_betaStim , tspan, sy0, theCell->y0, dt, sy, y1 ,t1 , theCell );
		
		for( int i = 0; i < sy0; i++ ) {
			theCell->y0[i] = y1[sy-1][i];
		}
		
		
        
	}
	
	return ( - (theCell->y0[0] - V) / dt  );
	
    
}
