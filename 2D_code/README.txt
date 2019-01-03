Source code contains: main_cell_ss_2D.cpp (main program), ord_model_betaStim.h (O’Hara-Rudy Human model) and integrate_rk2rk1.h (Numerical method)

######################################################################################

model_INPUTSteadyStateControl.1DInitials was generated using one-D cable program and allowing the transmural cable to run for 1000 beats (dt = 0.005, I_inj = -300 at the indicated frequency without beta-stimulation. Cell type and beta-stimulations were set in 1D main program (main_cell_ss_1D.cpp).

Figures of time snapshot Vm were generated using model_INPUTSteadyStateControl.1Dinitials, pacing for 2 beats without beta-stimulation.

######################################################################################
Transmural in one-direction: Cells 1–60 were endo, 61–105 were M, and 106– 165 were epi.

           
######################################################################################

Key Parameters in main_cell_ss_2D.cpp:

1. INPUT initial steady-state condition: Change initials input name at line 162

2. OUTPUT file names: Change OUTPUT file names at line 301 (Create voltage (tissue) files for every 1 ms)

3. Set cycle length (double BCL) at line 24
 
4. Set number of beat (double stimuli) at line 26

5. Set length of the tissue (double tl) at line 28

6. Set width of the tissue (double tw) at line 29

7. Set cell size (double dx: x-direction) at lind 34

8. Set cell size (double dy: y-direction) at lind 35

9. Set gap junction (double gj_x: x-direction ) at line 38

10. Set gap junction (double gj_y: y-direction) at line 39

11. Set inject current ( theCell->I_inj) at line 226



