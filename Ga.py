import main
import random
from matplotlib.animation import FuncAnimation
import matplotlib.pyplot as plt
import sys


class Pop:

    def __init__(self,pep,k,c,t=None):
        """
        :param pep: list of peptides in the population
        :param k: int, selection pressure: number of individuals for tournament selection
        :param t: float, mutation rate, [0,1]
        :param c: float, crossover rate, [0,1]
        """

        self.cross = c
        self.pop = pep
        self.k = k
        self.n = len(self.pop)
        if t is None:
            self.t = 1/self.n
        else:
            self.t = t

    @classmethod
    def random_init(cls,seq,n,k,c,t=None):
        """
        Random initiation of the population
        :param seq: string sequence of 0's and 1's representing if a aa is H or P respectably
        :param n: number of individuals in the population
        :param k: int, selection pressure: number of individuals for tournament selection
        :param t: float, mutation rate, [0,1]
        :param c: float, crossover rate, [0,1]
        :return: Pop object
        """
        temp = []
        for pept in range(n):
            temp.append(main.Pep([], seq))

        return cls(temp,k,c,t)

    def selection(self):
        """
        :return: individual chosen by k tournament
        """
        selected = []
        fittest = 0

        # select two parents
        for ind in range(2):

            # each from a k tournament
            for i in range(self.k-1):
                pept = random.choice(self.pop)

                if type(fittest) == int:
                    fittest = pept

                elif pept.fitness() > fittest.fitness():
                    fittest = pept

            selected.append(fittest)

        return selected

    def replacement(self):
        random_n = random.uniform(0, 1)
        new_pop = []

        # select two parents from population until new_pop has the same number of individuals as old_pop
        while len(new_pop) < self.n:
            parent1, parent2 = self.selection()

            # crossover the two selected parents with probability of self.cross
            if random_n < self.cross:
                child1, child2 = parent1.crossover(parent2)
            else:
                child1, child2 = parent1, parent2

            # mutate each new individual with a probability of self.t
            for child in [child1,child2]:
                random_n = random.uniform(0, 1)
                if random_n < self.t:
                    child = child.mutate()
                new_pop.append(child)

        return Pop(new_pop,self.k,self.cross,self.t)

    def meanfitness(self):
        fitnessavg = 0
        highst_fit = 0
        for pep in self.pop:
            fit_pep = pep
            fit = pep.fitness()
            fitnessavg += fit
            if fit > highst_fit:
                highst_fit = fit
                fit_pep = pep

        return (fitnessavg/self.n),highst_fit,fit_pep


class GA:

    def __init__(self,population,gen,iterations,fitness,n,k,c,t=None):
        """
        :param population: Population of this generation
        :param gen: Number of the present generation
        :param iterations: Num
        :param fitness: Fitness maximum the algorithm can reach
        :param n: Number of individuals in the population
        :param k: Number of individuals chosen for k tournament selection
        :param c: Rate of crossover
        :param t: Rate of mutation
        """
        self.gen = gen
        self.iterations = iterations
        self.fitness = fitness
        self.n = n
        self.k = k
        self.c = c
        self.t = t
        self.population = population
        self.fitness = self.population.meanfitness()

    @classmethod
    def random_init(cls,iterations,fitness,seqHP,n,k,c,t):
        """
        :param iterations: Number of maxmn iterations/generations the algorithm can reach
        :param fitness: Fitness maximum the algorithm can reach
        :param seqHP: Sequence of aa's in format eg. "HPHPHHHP"
        :param n: Number of individuals in the population
        :param k: Number of individuals chosen for k tournament selection
        :param c: Rate of crossover
        :param t: Rate of mutation
        :return: Ga object with gen = 0
        """
        seq = ""
        for i in seqHP:
            if i == "H":
                seq = seq + "0"
            else:
                seq = seq + "1"

        population = Pop.random_init(seq, n, k, c, t)
        print("The mean fitness is:", population.meanfitness())
        return cls(population,0,iterations,fitness,n,k,c,t=None)

    def oneiteration(self):
        gen = self.gen + 1
        new_c = self.population.cross
        new_t = self.population.t

        new_pop = self.population.replacement()
        print("The mean fitness is:",new_pop.meanfitness()[0])
        print("The highest fitness is:",new_pop.meanfitness()[1])

        return GA(new_pop,gen,self.iterations,self.fitness,self.n,self.k,new_c,new_t)

    def plotbest(self):
        best = self.fitness[2]
        best.plot(1,1)



def runGa(iterations,fitness,seqHP,n,k,c,t=None):
    pop_avg_fitness = []
    pop_highest_fitness = []
    pop = GA.random_init(iterations,fitness,seqHP,n,k,c,t)
    pop_avg_fitness.append(pop.fitness[0])
    pop_highest_fitness.append(pop.fitness[1])

    while pop.gen < iterations and pop.fitness[0] < fitness:
        print("Generation:",pop.gen)
        pop = pop.oneiteration()
        pop_avg_fitness.append(pop.fitness[0])
        pop_highest_fitness.append(pop.fitness[1])

    if pop.gen == iterations:
        print("Reached maximum number of iterations")
    if pop.fitness == fitness:
        print("Reached fitness maximum")


    pop.plotbest()
    gens = [i for i in range(pop.gen+1)]
    plt.plot(gens, pop_avg_fitness)
    plt.plot(gens, pop_highest_fitness)
    plt.legend(["Mean fitness", "Highest fitness"], loc='upper left')

    plt.xlabel("Generations")
    plt.ylabel("Fitness")


    plt.show()


    pop.population.meanfitness()

    return pop_avg_fitness,pop_highest_fitness,pop.fitness[2]


# ll = runGa(400,50,"HPHPPHHPHPPHPHHPPHPH",100,4,0.5,0.8)
iterations, fitness, sequence, size, k, cross, mutation = sys.argv[1:]

runGa(int(iterations),int(fitness),sequence,int(size),int(k),float(cross),float(mutation))

