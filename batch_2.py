# -*- coding: utf-8 -*-
"""
Created on Sat Oct  9 22:44:40 2021
@author: coss
        agarcia-mena
"""

import numpy as np
import copy

MAX_VALUE_ALLOWED = 65535
MAXSPOT = 25
MAXDISTANCE = 50000
FILE_0 = 'data/3-9-2021---17-22-21_Llacer335_LR2xyshots.txt'
FILE_1 = 'data/11-10-2021----18-31-38_grid345_xy_shots.txt'
FILE_2 = 'data/11-11-2021----15-5-9_Skpyke369_shot_run.txt.txt'
FILE_3 = 'data/29-10-2021----11-50-58_3030_2nd_XY_shots_run.txt.txt'
RADIATION = 0.2 #0-no mutations. 1-all chromosomes mute
ELITE = 0.2 #percentage individuals  in next generation
NGEN = 500
N_POB = 100
INITIAL_GEN_NUMBER = 1 #ges in chomosomes in first generation
CHANGES_ALLOWED_PER_LENGTH = 0.2 #changes / length-chromosome allowed
DEBUG = False

class Spot:
    def __init__(self, pos, cla1, is1):
        self.pos = pos
        self.cla1 = cla1
        self.is1 = is1

    def __add__(self, other):
        return Spot(self.pos, 
                    (self.cla1[0]+other.cla1[0], self.cla1[1]+other.cla1[1]),
                    (self.is1[0]+other.is1[0], self.is1[1]+other.is1[1]))
    
    def same_location(self, other):
        return self.cla1[0] == other.cla1[0] and \
               self.cla1[1] == other.cla1[1] and \
               self.is1[0] == other.is1[0] and \
               self.is1[1] == other.is1[1]
    
    def is_within_distance(self, maxDist):
        return abs(self.cla1[0]) < maxDist and \
               abs(self.cla1[1]) < maxDist and \
               abs(self.is1[0]) < maxDist and \
               abs(self.is1[1]) < maxDist
    
    def isValid(self):
        return 0 < self.cla1[0] < MAX_VALUE_ALLOWED and \
               0 < self.cla1[1] < MAX_VALUE_ALLOWED and \
               0 < self.is1[0] < MAX_VALUE_ALLOWED and \
               0 < self.is1[1] < MAX_VALUE_ALLOWED
    
    def printSpot(self, change=False):
        changeStr = ""
        if change:
            changeStr = " === "
        print("Pos=[%d,%d] CLA1=[%d,%d] IS1=[%d,%d] %s"%\
              (self.pos[0], self.pos[1],
               self.cla1[0], self.cla1[1],
               self.is1[0], self.is1[1], changeStr))


def read_spots(fn, maxDistance, maxSpot):
    # Read all spots
    fh = open(fn, 'r')
    allSpots = []
    refSpot = None
    for line in fh.readlines():
        if line.startswith('[['):
            tokens = [x.strip() for x in line.replace('[', ' ').replace(']', ' ').
                split(',')]
            pos = (int(tokens[0]), int(tokens[1]))
            if pos[1] != 0 or pos[0] > maxSpot:
                continue
            cla1x = int(tokens[8])
            cla1y = int(tokens[9])
            
            is1x = int(tokens[16])
            is1y = int(tokens[17])
            newSpot = Spot(pos, (cla1x, cla1y), (is1x, is1y))
            allSpots.append(newSpot)
            #print(str(pos[0]), ', ' + str(pos[1]) + '\n')
            if pos == (0, 0):
                refSpot = Spot(pos, (cla1x, cla1y), (is1x, is1y))
    fh.close()
    if refSpot is None:
        raise Exception("Cannot find reference spot")
    
    # Add the refSpot to all spots to get the valid ones
    # Check they are within a distance
    allValidSpots=[]
    def already_in_list(newSpot, allValidSpots):
        for spot in allValidSpots:
            if spot.same_location(newSpot):
                return True
        return False
    
    for spot in allSpots:
        valid = False
        if spot.pos != (0, 0):
            if spot.is_within_distance(maxDistance):
                newSpot=spot+refSpot
                valid = newSpot.isValid()
        else:
            newSpot = copy.copy(spot)
            valid = True
        if valid and not already_in_list(newSpot, allValidSpots):
            allValidSpots.append(newSpot)
    return allValidSpots

# Individual
class Individual():
    def __init__(self, chromosome=None, n=None):
        if not n is None:
            self.initRandom(n)
        else:
            self.chromosome = chromosome
    
    def initRandom(self, n):
        self.chromosome = np.random.permutation(n)[:INITIAL_GEN_NUMBER].tolist() #TO_CHANGE

    def mutation(self, p, allValidSpots):
        if np.random.uniform() < p:
            i = int(np.random.uniform(0, len(self.chromosome)))
            randomVaule = np.random.uniform()
            if randomVaule < 0.25:
                j = int(np.random.uniform(0, len(self.chromosome)))
                aux = self.chromosome[i]
                self.chromosome[i] = self.chromosome[j]
                self.chromosome[j] = aux
            elif 0.25 <= randomVaule < 0.5:
                self.chromosome = self.chromosome[i:] + self.chromosome[0:i]
            elif 0.5 <= randomVaule < 0.75:
                self.chromosome.append(
                    int(np.random.uniform(0, len(allValidSpots) - 1)))
            else:
                self.chromosome.insert(0,
                          int(np.random.uniform(0, len(allValidSpots) - 1)))

    def crossover(self, other):
        N = len(self.chromosome)
        i = int(np.random.uniform(0, N))
        return Individual(self.chromosome[0:i]+other.chromosome[i:N])

    def splice(self, other):
        i = int(np.random.uniform(0, len(other.chromosome)))
        np_random = np.random.uniform()
        if np_random < 0.2:
            return Individual(self.chromosome + other.chromosome)
        elif 0.2 < np_random <= 0.4:
            return Individual(other.chromosome + self.chromosome)
        elif 0.4 < np_random <= 0.6:
            return Individual(other.chromosome[:i] + self.chromosome)
        elif 0.6 < np_random <= 0.8:
            return Individual(other.chromosome[i:] + self.chromosome)
        else:
            return Individual(self.chromosome + other.chromosome[:i])

    def isValidChromosome(self, allValidSpots):
        set_chr = set(self.chromosome)
        if len(self.chromosome) != len(set_chr) or \
                len(self.chromosome) > len(allValidSpots):
                return False
        else:
            return True

    def changei(self, allValidSpots, i):
        if i >= len(self.chromosome)-2:
            return 0
        icla1x = allValidSpots[self.chromosome[i]].cla1[0]
        icla1y = allValidSpots[self.chromosome[i]].cla1[1]
        jcla1x = allValidSpots[self.chromosome[i+1]].cla1[0]
        jcla1y = allValidSpots[self.chromosome[i+1]].cla1[1]
        kcla1x = allValidSpots[self.chromosome[i+2]].cla1[0]
        kcla1y = allValidSpots[self.chromosome[i+2]].cla1[1]
        
        sign1x = np.sign(jcla1x-icla1x)
        sign1y = np.sign(jcla1y-icla1y)
        sign2x = np.sign(kcla1x-jcla1x)
        sign2y = np.sign(kcla1y-jcla1y)
        
        change = False
        if sign1x != sign2x:
            change = True
        if sign1y != sign2y:
            change = True
            
        iis1x = allValidSpots[self.chromosome[i]].is1[0]
        iis1y = allValidSpots[self.chromosome[i]].is1[1]
        jis1x = allValidSpots[self.chromosome[i+1]].is1[0]
        jis1y = allValidSpots[self.chromosome[i+1]].is1[1]
        kis1x = allValidSpots[self.chromosome[i+2]].is1[0]
        kis1y = allValidSpots[self.chromosome[i+2]].is1[1]
        
        sign1x = np.sign(jis1x-iis1x)
        sign1y = np.sign(jis1y-iis1y)
        sign2x = np.sign(kis1x-jis1x)
        sign2y = np.sign(kis1y-jis1y)
        
        if sign1x != sign2x:
            change = True
        if sign1y != sign2y:
            change = True
        return change
    
    def fitnessAllChanges(self, allValidSpots):
        retval = 0
        N = len(self.chromosome)
        for i in range(0, N-2):
            if self.changei(allValidSpots, i):
                retval += 1
        return retval

    def fitnessChanges(self, allValidSpots):
        retval = 0
        N = len(self.chromosome)
        for i in range(0, N-2):
            if self.changei(allValidSpots, i):
                return 1
        return 0


    def print(self, allValidSpots):
        for i in range(len(self.chromosome)):
            allValidSpots[self.chromosome[i]].printSpot(
                self.changei(allValidSpots, i))



class Population():
    def __init__(self, allValidSpots, p=0, debug=False):
        self.allValidSpots = allValidSpots
        self.debug = debug
        self.allIndividuals = []
        self.lenAllValidSpots = len(allValidSpots)
        for i in range(p):
            self.addIndividual(Individual(n=self.lenAllValidSpots))
        if p > 0:
            self.evaluateFitness()
            self.evaluateLength()


    def evaluateLength(self):
        N = len(self.allIndividuals)
        self.individualLength = np.zeros(N)
        for n in range(N):
            self.individualLength[n] = len(self.allIndividuals[n].chromosome)
        self.individualLength = [int(x) for x in self.individualLength]
        self.print(ev='length')

    def evaluateFitness(self):
        N = len(self.allIndividuals)
        self.individualFitness = np.zeros(N)
        for n in range(N):
            self.individualFitness[n] = self.allIndividuals[n].fitnessAllChanges(
                self.allValidSpots)
        self.individualFitness = [int(x) for x in self.individualFitness]
        if self.debug:
            self.print(ev='fitnes')
    
    def addIndividual(self, individual):
        self.allIndividuals.append(individual)
    
    def print(self, ev, N = int(ELITE * 100)):
        s = np.sort(self.individualFitness)
        fitness_sort = np.argsort(self.individualFitness)
        if ev == 'fitnes':
            print('Elite fitness: ', list(s[0:N]))
        else:
            listLe = []
            for x in fitness_sort[:N]:
                listLe.append(len(self.allIndividuals[x].chromosome))
            if self.debug:
                print('Elite length: ', np.sort(listLe))
    
    def newGeneration(self):
        newPopulation = Population(self.allValidSpots, debug=self.debug)
        sortedIdx = list(reversed(np.argsort(self.individualLength)))
        # Elite 20% longest and zero change
        elite_idx = 0
        for i in sortedIdx:
            if len(newPopulation.allIndividuals) < int(ELITE * len(sortedIdx))\
                and self.individualFitness[i] / len(self.allIndividuals[i].chromosome) <= \
                    CHANGES_ALLOWED_PER_LENGTH:
                    newPopulation.addIndividual(self.allIndividuals[i])
                    # if self.debug:
                    #     print(newPopulation.allIndividuals[elite_idx].chromosome)
                    elite_idx += 1

        # Sexual reproduction
        n = len(self.allIndividuals)
        p = self.individualLength
        if np.sum(p) == 0:
            p = p + (1 / len(p))
        else:
            p = p/np.sum(p)
        allIdx = np.arange(0, n)
        while len(newPopulation.allIndividuals) != n:
            i = int(np.random.choice(allIdx, p=p))
            j = int(np.random.choice(allIdx, p=p))
            newIndividual = self.allIndividuals[i].splice(self.allIndividuals[j])
            newIndividual.mutation(RADIATION, self.allValidSpots)
            if newIndividual.isValidChromosome(self.allValidSpots) and\
                    self.is_not_repeated(newPopulation, newIndividual):
                newPopulation.addIndividual(newIndividual)
        newPopulation.evaluateFitness()
        newPopulation.evaluateLength()
        return newPopulation

    def is_not_repeated(self, newPopulation, newIndividual):
        for p in newPopulation.allIndividuals:
            if p.chromosome == newIndividual.chromosome:
                return False
        return True

    def finished(self, N = int(ELITE * 100)):
        fitness_sort = np.argsort(self.individualFitness)
        listLe = []
        for x in fitness_sort[:N]:
            len_x = len(self.allIndividuals[x].chromosome)
            listLe.append(len_x)
            if len_x == self.lenAllValidSpots:
                print('Elite length: ', np.sort(listLe))
                return True

    def getBest(self):
        list_elems = []
        len_init = 0
        for i, n in enumerate(self.individualFitness):
            coef = self.individualFitness[i] / len(self.allIndividuals[i].chromosome)
            if coef <= CHANGES_ALLOWED_PER_LENGTH:
                list_elems.append((self.allIndividuals[i].chromosome,
                                   self.individualFitness[i]))
                if int(len(self.allIndividuals[i].chromosome)) > len_init:
                    len_init = int(len(self.allIndividuals[i].chromosome))
        for x in list_elems:
            if len(x[0]) == len_init:
                print('\nChromosome: {}\nchanges: {} '
                      'length: {} coef: {} <= {}'.
                      format(x[0], x[1], len(x[0]), round(x[1] / len(x[0]), 2),
                             CHANGES_ALLOWED_PER_LENGTH))
                return x[0]



def optimize_JAFIS(fn, maxSpot, maxDistance=MAXDISTANCE, Npop=100, Ngen=30, debug=False):
    allValidSpots = read_spots(fn, maxDistance, maxSpot)
    if debug:
        print("all valid spots")
        for spot in allValidSpots:
            spot.printSpot()


    population = Population(allValidSpots, Npop, debug)
    for n in range(Ngen):
        if debug:
            print("\nGeneration %d" % (n+1))
        population = population.newGeneration()
        if population.finished():
            break
    return (population.getBest(), allValidSpots)


def get_num_combs(numHoles=35, numGens=3):
    import itertools

    stuff = range(numHoles)
    listSubset = []
    Numcombs = 0
    for subset in itertools.combinations(stuff, numGens):
        permutations = list(itertools.permutations(subset))
        Numcombs += len(permutations)
        listSubset.append(permutations)
    return Numcombs


if __name__ == "__main__":
    print('Number of 3-gene combinations to go through {} holes: {}'
          .format(21, get_num_combs(numHoles=21, numGens=3)))
    print('Number of 4-gene combinations to go through {} holes: {}'
          .format(21, get_num_combs(numHoles=21, numGens=4)))
    print('Number of 5-gene combinations to go through {} holes: {}\n'
          .format(21, get_num_combs(numHoles=21, numGens=5)))

    if not DEBUG:
        print('debug = ', DEBUG)
    FILES = [FILE_0]
    for F in FILES:
        bestPath, allValidSpots = optimize_JAFIS(
            F, maxSpot=MAXSPOT, Npop=N_POB, debug=DEBUG, Ngen=NGEN)
