# -*- coding: utf-8 -*-
"""
Created on Sat Oct  9 22:44:40 2021
Edited on Fry Nov 26 08:08
@author: coss
        agarcia-mena
"""

import numpy as np
import copy

MAX_VALUE_ALLOWED = 65535
MAXSPOT = 25
MAXDISTANCE = 50000
DEBUG = True
NGEN = 300
RADIATION = 0.2
ELITE = 0.2
FILE = '3-9-2021---17-22-21_Llacer335_LR2xyshots.txt'

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
        self.chromosome = np.random.permutation(n).tolist() #TO_CHANGE
    
    def mutation(self, p):
        if np.random.uniform() < p:
            i = int(np.random.uniform(0, len(self.chromosome)))
            if np.random.uniform() < 0.5:
                j = int(np.random.uniform(0, len(self.chromosome)))
                aux = self.chromosome[i]
                self.chromosome[i] = self.chromosome[j]
                self.chromosome[j] = aux
            else:
                self.chromosome = self.chromosome[i:] + self.chromosome[0:i]
                
    def combine(self, other):
        N = len(self.chromosome)
        i = int(np.random.uniform(0, N))
        return Individual(self.chromosome[0:i]+other.chromosome[i:N])
    
    def isValidChromosome(self):
        N = len(self.chromosome)
        for n in range(N):
            if not n in self.chromosome:
                return False
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
    
    def fitnessLargestSeq(self, allValidSpots):
        N = len(self.chromosome)
        retval = N
        for i in range(0, N-2):
            if not self.changei(allValidSpots, i):
                retval -= 1
            else:
                break
        return retval

    def print(self, allValidSpots):
        for i in range(len(self.chromosome)):
            allValidSpots[self.chromosome[i]].printSpot(
                self.changei(allValidSpots, i))

class Population():
    def __init__(self, allValidSpots, p=0, debug=False):
        self.allValidSpots = allValidSpots
        self.debug = debug
        self.allIndividuals = []
        n = len(allValidSpots)
        for i in range(p):
            self.addIndividual(Individual(n=n))
        if p > 0:
            self.evaluateFitness()
    
    def evaluateFitness(self):
        N = len(self.allIndividuals)
        self.individualFitness = np.zeros(N)
        for n in range(N):
            self.individualFitness[n] = self.allIndividuals[n].fitnessAllChanges(
                self.allValidSpots)
            # self.individualFitness[n]=self.allIndividuals[n].fitnessLargestSeq(self.allValidSpots)
        if self.debug:
            self.print()
    
    def addIndividual(self, individual):
        self.allIndividuals.append(individual)
    
    def print(self, N = 30):
        s = np.sort(self.individualFitness)
        print(s[0:N])
    
    def newGeneration(self):
        newPopulation = Population(self.allValidSpots, debug=self.debug)
        sortedIdx = np.argsort(self.individualFitness)
        
        # Elite 20%
        for i in range(int(ELITE*len(sortedIdx))):
            newPopulation.addIndividual(self.allIndividuals[sortedIdx[i]])
        
        # Sexual reproduction
        n = len(self.allIndividuals)
        p = np.max(self.individualFitness)-self.individualFitness+1
        p = p/np.sum(p)
        allIdx = np.arange(0, n)
        while len(newPopulation.allIndividuals) != n:
            i = int(np.random.choice(allIdx, p=p))
            j = int(np.random.choice(allIdx, p=p))
            newIndividual = self.allIndividuals[i].combine(self.allIndividuals[j])
            newIndividual.mutation(RADIATION)
            if newIndividual.isValidChromosome():
                newPopulation.addIndividual(newIndividual)
        newPopulation.evaluateFitness()
        
        return newPopulation

    def finished(self):
        return np.min(self.individualFitness) == 0
    
    def getBest(self):
        sortedIdx = np.argsort(self.individualFitness)
        return self.allIndividuals[sortedIdx[0]]        
    
    def printBest(self):
        bestIndividual = self.getBest()
        bestIndividual.print(self.allValidSpots)
        

def optimize_JAFIS(fn, maxSpot, maxDistance=MAXDISTANCE, Npop=100, Ngen=3000, debug=False):
    allValidSpots = read_spots(fn, maxDistance, maxSpot)
    if debug:
        print("all valid spots")
        for spot in allValidSpots:
            spot.printSpot()

    print('Number of 3-gene combinations to go through {} holes: {}\n'
          .format(len(allValidSpots), get_num_combs(numHoles=len(allValidSpots), numGens=3)))

    population = Population(allValidSpots, 100, debug)
    for n in range(Ngen):
        if debug:
            print("Generation %d"%(n+1))
        population = population.newGeneration()
        if population.finished():
            break
    if debug:
        population.printBest()
    return (population.getBest(), allValidSpots)


def get_num_combs(numHoles=25, numGens=3):
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
    bestPath, allValidSpots = optimize_JAFIS(FILE, maxSpot=MAXSPOT, debug=DEBUG,
                                             Ngen=NGEN)
    bestPath.print(allValidSpots)
