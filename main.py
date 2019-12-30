from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt
import numpy as np
import random
from copy import deepcopy


class Aa:

    def __init__(self,bin,n,pos):
        """
        :param bin: 0 if H. 1 if P.
        :param pos: from 1 to pept len
        """
        self.n = n
        self.bin = bin
        self.pos = pos
        self.vect = [[1, 1, 0], [1, -1, 0], [-1, -1, 0], [-1, 1, 0],
                     # up
                     [1, 0, 1], [0, -1, 1], [-1, 0, 1], [0, 1, 1],
                     # down
                     [1, 0, -1], [0, -1, -1], [-1, 0, -1], [0, 1, -1]]

        self.neigh = []
        for vect in self.vect:
            self.neigh.append(list(np.add(np.array(self.pos),np.array(vect))))

    def consecutive(self,otheraa):
        """
        :param self:
        :param otheraa: other aa object
        :return: are two aa neighbourss
        """
        if abs(otheraa.n - self.n) == 1:
            return True

    # |xi − xj|+|yi − yj|+|zi − zj| = 2 and |i - j| > 1
    def neighbours(self,otheraa):
        """
        :param otheraa:
        :return: true if two aa are neighbours
        """

        distance = np.linalg.norm(np.array(self.pos)-np.array(otheraa.pos))
        return distance == 2**0.5

    def setPos(self,pos):
        self.pos = pos

    def poss_neigh(self,aalist):
        """
        returns all free neighboring positions given a list of aa's
        """
        availablespots = deepcopy(self.neigh)
        for aa in aalist:
            if aa.pos in availablespots:
                availablespots.remove(aa.pos)

        return availablespots

    def __str__(self):
        return "(({},{}))".format(self.pos, self.n)

class Pep:

    def __init__(self,aas,seq=None):
        """
        :param seq: string of sequence of aas with bin value...
        :param aas: list od aa's object
        """
        self.aas = aas
        self.seq = seq
        # vector to neighbourss
                     # same plane
        self.vect = [[1, 1, 0], [1, -1, 0], [-1, -1, 0], [-1, 1, 0],
                     # up
                     [1, 0, 1], [0, -1, 1], [-1, 0, 1], [0, 1, 1],
                     # down
                     [1, 0, -1], [0, -1, -1], [-1, 0, -1], [0, 1, -1]]

        if len(aas) == 0:
            #self.taken = {}
            self.aaspos = []

            ## initialize with random positioning
            anum = 1
            aminoa = Aa(self.seq[0],0,[0,0,0])

            while anum < len(self.seq):
                self.aas.append(aminoa)
                newpos = list(np.add(np.array(aminoa.pos),
                                 np.array(random.choice(self.vect))))

                # while the position is ocupied generate new one
                while newpos in self.aaspos:
                    newpos = list(np.add(np.array(aminoa.pos),
                                     np.array(random.choice(self.vect))))

                aminoa = Aa(self.seq[anum],anum,newpos)

                # self.taken is "ordered" by aa.n
                #self.taken[tuple(newpos)] = aminoa
                self.aaspos.append(newpos)

                anum += 1

        else:
            #self.taken = {}
            self.aaspos = []
            for aa in self.aas:
                # do i need this dictionary ?
                #self.taken[tuple(aa.pos)] = aa
                self.aaspos.append(aa.pos)

    def setAas(self,aaslist):
        self.aas = aaslist

    def setAasPos(self,aasplist):
        self.aaspos = aasplist

    def fitness(self):
        HHbonds = 0

        for aa in self.aas:
            # if aa is H then
            if aa.bin == "0":
                index_aa = self.aas.index(aa)

                # compute every possible position for bonding
                for v in self.vect:
                    possibleHH = list(np.add(np.array(aa.pos),np.array(v)))

                    # if there is aa in that position
                    # search only furthers aa's
                    if possibleHH in self.aaspos:
                        index_HH = self.aaspos.index(possibleHH)
                        other_aa = self.aas[index_HH]

                        # if that aa is also a H and aa+other_aa are not neighbourss:
                        if other_aa.bin == "0" and not aa.consecutive(other_aa):
                            HHbonds += 1

            return HHbonds

    # by pull move
    def mutate(self):
        """
        Mutates a peptide in a random aa by pull move
        """
        self.plot()

        # choose a random aa

        random_n = random.randrange(0, len(self.aas))
        new_aa = self.aas[random_n]

        # give the chosen aa a new free neighboring position
        possible_position = new_aa.poss_neigh(self.aas)
        random_pos = random.choice(possible_position)
        new_aa.setPos(random_pos)

        # makes it impossible to check the same position twice
        possible_position.remove(random_pos)
        print("N chosen, pos chosen",random_n,random_pos)
        print("Aa chosen:",new_aa)

        # check only left aa
        if random_n == 0:

            print(new_aa.pos)
            print(new_aa.neighbours(self.aas[random_n + 1]))
            print(self.aas[random_n + 1])

            new_aa.setPos(random_pos)

            if new_aa.neighbours(self.aas[random_n + 1]):

                self.aas[random_n] = new_aa
                self.aaspos[random_n] = new_aa.pos

                self.plot()

            # pull left
            else:
                dummy_aas = []
                dummy_pos = []

                dummy_pos.append(new_aa)

                for a in range(1,len(self.aas)):

                    aa = self.aas[a]
                    aa.setPos(self.aaspos[a-1])
                    dummy_aas.append(aa)
                    dummy_pos.append(aa.pos)

                self.setAas(dummy_aas)
                self.setAasPos(dummy_pos)

                self.plot()

        # check only right aa
        elif random_n == len(self.aas)-1:

            print("New position:",new_aa.pos)
            print("Is right aa neighbour?:",new_aa.neighbours(self.aas[random_n - 1]))
            print("Right neighbour:",self.aas[random_n - 1])

            new_aa.setPos(random_pos)

            if new_aa.neighbours(self.aas[random_n - 1]):

                self.aas[random_n] = new_aa
                self.aaspos[random_n] = random_pos

                self.plot()

            # pull right
            else:
                dummy_aas = []
                dummy_pos = []

                for a in range(0, len(self.aas)-1):
                    aa = self.aas[a]
                    aa.setPos(self.aaspos[a + 1])
                    dummy_aas.append(aa)
                    dummy_pos.append(aa.pos)

                dummy_aas.append(new_aa)
                dummy_pos.append(new_aa.pos)

                self.setAas(dummy_aas)
                self.setAasPos(dummy_pos)

                self.plot()

        # check both
        else:
            # while the random_pos is taken and there is no aa to be neighbours, generate a new position

            while not any([new_aa.neighbours(self.aas[random_n+1]),new_aa.neighbours(self.aas[random_n-1])]):
                random_pos = random.choice(possible_position)
                possible_position.remove(random_pos)
                new_aa.setPos(random_pos)
                print("wtf")

            print("New pos:",new_aa.pos)
            print(new_aa.neighbours(self.aas[random_n+1]),new_aa.neighbours(self.aas[random_n-1]))
            print("Left and right neighbour:",self.aas[random_n+1],self.aas[random_n-1])

            # if True then the move is already legal (self avoiding and connected to both)
            if all([new_aa.neighbours(self.aas[random_n+1]),new_aa.neighbours(self.aas[random_n-1])]):

                self.aas[random_n] = new_aa
                self.aaspos[random_n] = new_aa.pos

                self.plot()

            else:
                # if the aa to the left is still neighbours then pull all the right aas to the left
                if new_aa.neighbours(self.aas[random_n - 1]):
                    dummy_aas = []
                    dummy_pos = []

                    dummy_aas.extend(self.aas[:random_n])
                    dummy_pos.extend(self.aaspos[:random_n])

                    dummy_aas.append(new_aa)
                    dummy_pos.append(new_aa.pos)

                    for a in range(random_n+1, len(self.aas)):
                        aa = self.aas[a]
                        aa.setPos(self.aaspos[a - 1])
                        dummy_aas.append(aa)
                        dummy_pos.append(aa.pos)

                    self.setAas(dummy_aas)
                    self.setAasPos(dummy_pos)

                    self.plot()


                # if the aa to the right is still neighbours then pull all the left aas to the right
                else:
                    dummy_aas = []
                    dummy_pos = []

                    for a in range(0, random_n):
                        aa = self.aas[a]
                        aa.setPos(self.aaspos[a + 1])
                        dummy_aas.append(aa)
                        dummy_pos.append(aa.pos)

                    dummy_aas.append(new_aa)
                    dummy_pos.append(new_aa.pos)

                    dummy_aas.extend(self.aas[random_n+1:])
                    dummy_pos.extend(self.aaspos[random_n+1:])

                    self.setAas(dummy_aas)
                    self.setAasPos(dummy_pos)

                    self.plot()

                return 1

    def crossover(self,other_pep):
        pass

    def plot(self):
        ax = plt.axes(projection='3d')

        zline, xline, yline = [],[],[]


        for point in self.aaspos:
            xline.append(point[0])
            yline.append(point[1])
            zline.append(point[2])

        # lines
        ax.plot3D(xline, yline, zline)
        # points
        ax.scatter3D(xline, yline, zline)
        ax.scatter3D(xline[0], yline[0], zline[0], color="k")
        ax.scatter3D(xline[-1], yline[-1], zline[-1], color="r")

        axes = plt.gca()
        axes.set_xlim([-10, 10])
        axes.set_ylim([-10, 10])
        axes.set_zlim([-10, 10])

        # ax.set_axis_off()

        plt.show()
