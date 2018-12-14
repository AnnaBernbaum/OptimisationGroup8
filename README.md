DE4 Opti Coursework 2018/19: Code Submission
============
This repository contains the code and SolidWorks files for Anna Bernbaum and Ina Roll Backe's (Group 8) optimisation project. 

The System
----------
A tumble dryer is an electrical appliance commonly found in homes. It is used to dry clothes by removing moisture through evaporation and centrifugal force. This optimisation problem aims to minimise the cost and maximise the drying efficiency of the tumble dryer, through a multi-objective formulation. 
The problem was split into two subsystems, the first of which aimed to optimise the heating element. Nonlinear solvers found that the optimum transferred heat flux was 3985 W/m2, a 150 % increase on the unoptimized problem and multi-objective methods reduced the cost by 30.6 %.
The second subsystem looked at the belt drive, maximising the angular velocity of the drum. The drum was modelled as a sub-subsystem and a non-linear regression model was fit to data points collected from a finite element analysis. The latter was used to identify optimal material selection and thickness, increasing the speed of the drum by 45%.
The overall, system level optimisation accounted for cost and drying efficiency, finding all combinations of drums and sheets that would cost under £3, and then selecting the best of these for drying efficiency through a weighted equation. This found the optimum heat flux to be 3984 W/m2, drum speed of 61.5 RPM, making a total cost of £2.93, a 67 % improvement on the unoptimized problem.

Directory Structure
-----------
1. System
2. Subsystem_1: The Heating Element
3. Subsystem_2: The Belt Drive

The System directory contains all necessary code and files to reproduce
the results presented for the System optimisation (i.e. Group). 

The Subsystem_x folders contain all necessary code and files such that
the results for the stand-alone subsystem optimisations that you have worked on
individually can be reproduced. 
