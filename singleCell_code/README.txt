Source code contains: main_cell_ss.cpp (main program), ord_model_betaStim.h (O’Hara-Rudy Human model) and integrate_rk2rk1.h (Numerical method)

######################################################################################

model_INPUTSteadyStateControl.InitialsENDO were arrived at by starting initial conditions in the main program (line 66 to 107 and line 163 to 191) and allowing the ENDO cell to run for 1000 beats (dt = 0.005, I_inj = -80 at the indicated frequency without beta-stimulation.

Figures of individual AP and concentrations traces were generated using model_INPUTSteadyStateControl.InitialsENDO, pacing for 2 beats without beta-stimulation.

######################################################################################


Key Parameters in main_cell_ss.cpp:

1. INPUT initial steady-state condition: Change initials input name at line 142

2. OUTPUT file names: Change OUTPUT file names at line 154 to line 157

3. Set basic cycle length (double BCL) at line 23
 
4. Set number of beat (double stimuli) at line 25

5. Set cell type (theCell->celltype) at line 189

6. Set beta-stimulation (theCell->coef_PKA) at line 190 

NOTE: Beta-stimulation inputs for linear change: range from 0 to 1 (0: basal level, 1: maximum level)

For example, set theCell->coef_PKA = 0 (no beta-stimulation) for control case.

7. Set inject current ( theCell->I_inj) at line 217



