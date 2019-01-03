/*
 *  O'Hara-Rudy Human model with RK method
 *  Created by Mao-Tsuen Jeng on 6/13/11.
 *
 *	Copyright 2011 __Colleen Clancy Lab__. All rights reserved.
 *
 *  Two-D simulation
 */




#include <fstream>
#include <iostream>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <omp.h>
using namespace std;




const double BCL = 500;

const int stimuli = 1;

const int tl = 165;  // length of the tissue
const int tw = 165;  // width of the tissue

const double base_dt = 0.005;
const int fold = 1.0/base_dt;

const double dx= 0.01; //(x-direction)
const double dy= 0.01; // (y-direction)


const double gj_x = 2.5; //uS (x-direction)
const double gj_y = 2.5 ; //uS (y-direction)

const double Rmyo = 150.0;       //ohm*cm
const double rad = 0.0011;
const double Rg_x = 3.14159*pow(rad,2)/(gj_x*1.0e-6);
const double Rg_y = 3.14159*pow(rad,2)/(gj_y*1.0e-6);

const double Dfx = ((1e3*rad)/(4.0*(Rmyo+Rg_x/dx)));

const double Dfy = (1e3*rad)/(4.0*(Rmyo+Rg_y/dy));



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
    Cell cellData[tw][tl];
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
    const char *stateFileName = "Save2D.currentState";
    const char *fromCable = "model_INPUT_Control.from1D";
    
    
    
    time_t startTime;
    time_t previousTime;
    
    const time_t timeSave = 0.1*60*60;
    startTime = time(NULL);
    previousTime = startTime;
    
    
    
    FILE *fp = fopen(stateFileName, "r");
    
    if (fp == NULL) {
        cout << "Start from 1D steady-state" << endl;
        
        S->counter=0;
        
        for ( ww = 0; ww < tw; ww++) {
            
            fp = fopen( fromCable, "r" );
            
            for ( ll = 0; ll < tl; ll++) {
                
                theCell = &(S->cellData[ww][ll]);
                fread( theCell, sizeof(Cell), 1, fp );
                
            }
            
            fclose( fp );
        }
        
        S->beat = 1;
        S->t = 0.0;
        S->tt = 0.0;
        S->tstep = 0;
        
    } else {
        
        fread(S, sizeof(SimState), 1, fp);
        cout << "restarting from current state " <<  endl;
        
        fclose(fp);
    }
    
    cycle_length = BCL;
    
    dt = base_dt;
    
    while (!done) {
        
        
        
#pragma omp parallel default(none) private( ww, ll, theCell,  I_Total ) shared( S, dt, cout  )
        {
#pragma omp for schedule(static, 1)
            for ( ww = 0; ww < tw; ww++) {
                
                for ( ll = 0; ll < tl; ll++) {
                    
                    theCell = &(S->cellData[ww][ll]);
                    
                    theCell->v = theCell->v_new;
                    
                    if ( S->t < 0.5 && ( ll < 15 ) ) {
                        theCell->I_inj = - 180;
                    } else {
                        theCell->I_inj=0.0;
                    }
                    
                    
                    I_Total = Calcu_I_Total( theCell, dt );
                    
                    
                    theCell->v = theCell->y0[0];
                    
                }
            }
        }
        
        
#pragma omp parallel default(none) private( ww, ll, theCell, dv ) shared( cout , S, dt )
        {
#pragma omp for // schedule(static, 1)
            for ( ww = 0; ww < tw; ww++) {
                
                for ( ll = 0; ll < tl; ll++ ) {
                    
                    theCell = &(S->cellData[ww][ll]);
                    
                    if(ww>0 && ww<(tw-1) && ll>0 && ll<(tl-1)) {
                        dv=dt*( Dfx*(S->cellData[ww][ll-1].v-2*theCell->v+S->cellData[ww][ll+1].v)/(dx*dx)+Dfy*(S->cellData[ww-1][ll].v-2*theCell->v+S->cellData[ww+1][ll].v)/(dy*dy) ) ;
                    }
                    else if (ww==0 && ll==0) {
                        dv=dt*( Dfx*(-theCell->v+S->cellData[ww][ll+1].v)/(dx*dx)+Dfy*(-theCell->v+S->cellData[ww+1][ll].v)/(dy*dy) ) ;
                    }
                    else if (ww==0 && ll==(tl-1)) {
                        dv=dt*( Dfx*(-theCell->v+S->cellData[ww][ll-1].v)/(dx*dx)+Dfy*(-theCell->v+S->cellData[ww+1][ll].v)/(dy*dy) );
                    }
                    else if (ww==(tw-1) && ll==0) {
                        dv=dt*( Dfx*(-theCell->v+S->cellData[ww][ll+1].v)/(dx*dx)+Dfy*(-theCell->v+S->cellData[ww-1][ll].v)/(dy*dy) );
                    }
                    else if (ww==(tw-1) && ll==(tl-1)) {
                        dv=dt*( Dfx*(-theCell->v+S->cellData[ww][ll-1].v)/(dx*dx)+Dfy*(-theCell->v+S->cellData[ww-1][ll].v)/(dy*dy) );
                    }
                    else if (ww==0 && ll>0 && ll<(tl-1)) {
                        dv=dt*( Dfx*(S->cellData[ww][ll-1].v-2*theCell->v+S->cellData[ww][ll+1].v)/(dx*dx)+Dfy*(-theCell->v+S->cellData[ww+1][ll].v)/(dy*dy) );
                    }
                    else if (ww==(tw-1) && ll>0 && ll<(tl-1)) {
                        dv=dt*( Dfx*(S->cellData[ww][ll-1].v-2*theCell->v+S->cellData[ww][ll+1].v)/(dx*dx)+Dfy*(-theCell->v+S->cellData[ww-1][ll].v)/(dy*dy) );
                    }
                    else if (ww>0 && ww<(tw-1) && ll==0) {
                        dv=dt*( Dfx*(-theCell->v+S->cellData[ww][ll+1].v)/(dx*dx)+Dfy*(S->cellData[ww-1][ll].v-2*theCell->v+S->cellData[ww+1][ll].v)/(dy*dy) );
                    }
                    else if (ww>0 && ww<(tw-1) && ll==(tl-1)) {
                        dv=dt*( Dfx*(-theCell->v+S->cellData[ww][ll-1].v)/(dx*dx)+Dfy*(S->cellData[ww-1][ll].v-2*theCell->v+S->cellData[ww+1][ll].v)/(dy*dy) );
                    }
                    
                    
                    
                    theCell->y0[0] += dv;
                }
            }
            
        } // end openmp
        
        
        
        S->t += dt;
        S->tt += dt;
        S->counter += 1;
        
        
        
        
        //Create voltage (tissue) files for every 1 ms
        if(   ( S->counter % (fold*1) == 0 ) ) {
            
            sprintf(name, "modelOUTPUT_ap%d.dat", S->counter/fold);
            output = fopen(name, "w");
            cout << name << endl;
            
            for (ww=0; ww<=(tw-1); ww+=1) {
                for (ll=0; ll<=(tl-1); ll+=1){
                    theCell = &(S->cellData[ww][ll]);
                    fprintf(output, "%8.2f\t", theCell->y0[0]);
                }
                fprintf(output, "\n");
            }
            fclose (output);
            
            
        }
        
        
        if (S->counter % (fold*100) == 0) {
            cout<< "tt: " << S->tt << ", runtime: " << (time(NULL) - startTime)/60 << " min " << ", tstep: " << S->tstep << ", counter; " << S->counter << endl;
        }
        
        
        
        if(S->t >= cycle_length) {
            S->t = 0;
            S->beat++;
            cout << "now the beat is: " << S->beat << endl;
        }
        
        if (S->beat > stimuli) {
            S->beat = 1;
            done = 1;
        }
        
        if (time(NULL) - previousTime > timeSave) {
            previousTime = time(NULL);
            cout << "Saving current states." << endl;
            FILE *fp = fopen(stateFileName, "w");
            fwrite(S, sizeof(SimState), 1, fp);
            fclose(fp);
            cout << "Saved." << endl;
        }
        
    } // end while
    
    
    
    
    fclose ( output );
    
    
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
