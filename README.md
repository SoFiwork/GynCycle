# Menstrual Cycle model published by Fischer-Holzhausen et al. 2022 (https://doi.org/10.1016/j.jtbi.2022.111150)

## Description 
Mathematical model of the female menstrual cycle that couples the dynamics between follicles and hormones.
Expansion of GynCycle (Röbitz et al. (2013)) with follicular growth mechanism from Lange et al. (2019)

![Alt text](https://github.com/SoFiwork/GynCycle_newVersion/blob/main/Flowchart.jpg)

## Installation
The model is written in MATLAB - no additional package installations are necessary.

To use this GynCycle version, please clone the repository:

```
git clone https://github.com/SoFiwork/GynCycle
```

## Usage

run HumanSimulationFollGrowth.m to start a simulation
-needed: 'Parameter.dat', 'InitialValues.txt'
-set specifications for simulation 
  -ShowStuff = 1: show plots
  -runnumb: set number of simulation runs
  -SaveStuff = 1: save data as CSV 
  -DirStuff : choose directory for CSV
  -simulation options
    -NorCycle = 1: run without drug treatment

## Authors
Sophie Fischer-Holzhausen
Susanna Röblitz

## License
This model is licensed under the MIT License.

## Citation

@article{fischer2022hormonal,
  title={Hormonal regulation of ovarian follicle growth in humans: Model-based exploration of cycle variability and parameter sensitivities},
  author={Fischer-Holzhausen, Sophie and R{\"o}blitz, Susanna},
  journal={Journal of Theoretical Biology},
  pages={111150},
  year={2022},
  publisher={Elsevier}
}

@article{fischer2021mathematical,
  title={Mathematical modeling and simulation provides evidence for new strategies of ovarian stimulation},
  author={Fischer, Sophie and Ehrig, Rainald and Sch{\"a}fer, Stefan and Tronci, Enrico and Mancini, Toni and Egli, Marcel and Ille, Fabian and Kr{\"u}ger, Tillmann HC and Leeners, Brigitte and R{\"o}blitz, Susanna},
  journal={Frontiers in endocrinology},
  volume={12},
  pages={613048},
  year={2021},
  publisher={Frontiers Media SA}
}










ODE model of the female menstrual cycle

How to run an simulation

-make sure that all files are available
    -InititalValues.txt
    -Parameter.dat
    -StartTimesPoiss.txt
    -FSH.txt
-set specifications in HumanSimulationFollGrowth.m
    -number of simulation runs 
    -specify save options
    -selected type of simulation 
      -normal cycle
      -treatment protocol
         -stimulation during the late folicular phase
         -stimulation during the luteal phase
-run!

HumanSimulationFollGrowth.m
-needed: 'Parameter.dat', 'InitialValues.txt'
-set specifications for simulation 
  -ShowStuff = 1: show plots
  -runnumb: set number of simulation runs
  -SaveStuff = 1: save data as CSV 
  -DirStuff : choose directory for CSV
  -simulation options
    -NorCycle = 1: run without drug treatment
    LutStim = 1: luteal stimulation protocol from Kuang et al. (2014b)
    LateFollPhase = 1: stimulation in the late follicular phase protocol from 
    Zhu and Fu (2019)

parameter vectors
-> para
para(1): ODE function called to test(0) or not (1)
para(2): number of equations not related to the follicular maturation -> number is higher than the number of used equations
para(7): mean FSH sensitivity
para(8): standard deviation FSH sensitivity
para(9): threshold LH concentration for ovulation
para(10): lifetime of follicle after reaching the maximal size

-> paraOde: parameter for growth equation
v: fractal dimensio
gamma: growth rate
xi: maximal follicular diameter
mu: proportion of self-harm
k: strength of competition
rho: rate of decline
Folmax: min. follicular size for ovulation

-> paraPoi: follicle start points follow Poissan distribution 
lambda: number of follicles per days
intervallPoi: intervall per day in which follicles appear

CreateLoadFolTimes.m: creates time points for follicle appearance and FSH sensitivity for each follicle
total number of follicles set to 1000 

Simulation_X.m: 
-calls ODE solve
-checks the state and growth behaviour of all follicles and evaluate the state
-there are several versions of the function - one for each simulation type (normal Cycle or with drug administration)
-Tovu: time point of the last ovulation -> 14 at the beginning of the simulation, in order to start the simulation -> will be set after each ovulation 
-Follicles: class in which all follicles and their properties (start point, destiny, growth behaviour) are saved
  -destiny of follicles
  1: ovulation
  -1: growing
  -2: decreasing
  3: bin enough but not ovulated due to low LH level
  4: will ovulate within the next 12 h 
-integration time is from one follicle appearance to the next one
-options: integration stops whenever a follicle ovulates 
-dosing_timeIdx: find the count of the current dosing event
-dosing_events: vector that includes the time point of the dosing event and the number of the administration 
->if more than one drug is administered than there is more than one dosing_events vector 
-output
  -FollOvulInfo: [number of follicle; starttime; time of ovulation; lifetime]
  -solutions
  -CycleInfo: [[0 Cyclelength]; [rest FollperCycle]; OvuT]

ODE_Model_X.m: function includes all ODEs

testfun_X:
-calculation of the follicles’ growth and decline
-calculation of E2 and P4
-includes also the calculation of drug concentrations
