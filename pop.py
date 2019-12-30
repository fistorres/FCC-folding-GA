import main
import random

class Pop:

    def __init__(self,seq,n,k,c,t=None):
        """
        :param n: int, number of individuals in the population
        :param k: int, selection pressure: number of individuals for tournament selection
        :param t: float, mutation rate, [0,1]
        :param c: float, crossover rate, [0,1]
        :param seq: list of aa's of the peptide sequence in binary
        """

        self.n = n
        self.cross = c
        self.pop = []
        for pept in range(n):
            self.pop.append(main.Pep([],seq))

        self.k = k
        if t is None:
            self.t = 1/n
        else:
            self.t = t

    def mutation(self):
        random_n = random.uniform(0, 1)
        new_pop = []
        for pept in self.pop:
            if random_n < self.cross:
                pept = pept.mutate()

            new_pop.append(pept)

        self.pop = new_pop

    def selection(self):
        selected = []
        fittest = 0

        # select two parents
        for ind in range(2):

            # each from a k tournament
            for i in range(self.k-1):
                pept = random.choice(self.pop)

                if fittest == 0:
                    fittest = pept

                elif pept.fitness() > fittest.fitness():
                    fittest = pept

            selected.append(fittest)

        return selected

    def replacement(self):
        random_n = random.uniform(0, 1)
        new_pop = []

        while len(new_pop) < self.n:
            parent1, parent2 = self.selection()
            if random_n < self.cross:
                child1, child2 = parent1.crossover(parent2)
            else:
                child1, child2 = parent1, parent2

            new_pop.extend([child1,child2])

        self.pop = new_pop



class GA:

    def __init__(self):
        pass

    def meanfitness(self):
        pass

    def plotbest(self):
        pass



