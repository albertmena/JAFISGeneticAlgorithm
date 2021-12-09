# JAFISGeneticAlgorithm
This repository calculates minimum hysteresis paths from JAFIS .txt file

# Run

Edit the parameter configuration and run the batch.py or batch_2.py 
- batch: initial version of the code.  Step through the holes and calculate the minimum hysteresis path.
- batch_2: Paths start with "INITIAL_GEN_NUMBER" and grow to "CHANGES_ALLOWED_PER_LENGTH". Terminate if "NGEN" or find zero hysteresis in the "MAXSPOT" path.

# Parameters

MAX_VALUE_ALLOWED = 65535

MAXSPOT = 25 #maximum number of spots in the system

MAXDISTANCE = 50000 #max distance allowed between holes

FILE_0 = 'data/3-9-2021---17-22-21_Llacer335_LR2xyshots.txt'

FILE_1 = 'data/11-10-2021----18-31-38_grid345_xy_shots.txt'

FILE_2 = 'data/11-11-2021----15-5-9_Skpyke369_shot_run.txt.txt'

FILE_3 = 'data/29-10-2021----11-50-58_3030_2nd_XY_shots_run.txt.txt'

RADIATION = 0.2 #0-no mutations. 1-all chromosomes mute

ELITE = 0.2 #percentage individuals  in next generation

NGEN = 500 # number of generatios alloweds

N_POB = 100 #size of the population

INITIAL_GEN_NUMBER = 1 #ges in chomosomes in first generation

CHANGES_ALLOWED_PER_LENGTH = 0.2 #changes / length-chromosome allowed
