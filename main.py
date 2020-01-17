from mpl_toolkits import mplot3d
from mpl_toolkits.mplot3d.axes3d import get_test_data
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

        # all neighboring positions
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
        return otheraa.pos in self.neigh

    def setPos(self,pos):

        return Aa(self.bin,self.n,pos)

    def poss_neigh(self,aalist):
        """
        returns all free neighboring positions given a list of aa's
        """
        availablespots = deepcopy(self.neigh)
        # print("!!",self)
        for aa in aalist:
            # print(aa)
            if aa.pos in availablespots:
                # print(aa)
                availablespots.remove(aa.pos)

        return availablespots

    def __str__(self):
        return "(({},{}))".format(self.pos, self.n)

    def __eq__(self, other):
        return self.n == other.n and self.bin == other.bin and self.pos == other.pos

    def __ne__(self,other):
        return self.n != other.n or self.bin != other.bin or self.pos != other.pos


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
            self.aaspos = []

            ## initialize with random positioning
            anum = 1
            aminoa = Aa(self.seq[0], 0, [0, 0, 0])
            self.aas.append(aminoa)
            self.aaspos.append(aminoa.pos)

            while anum < len(self.seq):

                newpos = list(np.add(np.array(aminoa.pos),
                                np.array(random.choice(self.vect))))

                # while the position is occupied generate new one
                while newpos in self.aaspos:
                    newpos = list(np.add(np.array(aminoa.pos),
                                     np.array(random.choice(self.vect))))

                aminoa = Aa(self.seq[anum], anum, newpos)

                self.aas.append(aminoa)
                self.aaspos.append(aminoa.pos)

                anum += 1

        else:
            self.aaspos = []
            for aa in self.aas:
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

        a_copy = deepcopy(self.aas)
        a_pos_copy = deepcopy(self.aaspos)
        # choose a random aa
        random_n = random.randrange(0, len(a_copy))
        new_aa = a_copy[random_n]

        # print(random_n)
        # give the chosen aa a new free neighboring position
        possible_position = new_aa.poss_neigh(a_copy)
        while len(possible_position) == 0:
            random_n = random.randrange(0, len(a_copy))
            new_aa = a_copy[random_n]

        random_pos = random.choice(possible_position)
        new_aa = new_aa.setPos(random_pos)

        # makes it impossible to check the same position twice
        possible_position.remove(random_pos)

        # check only left aa
        if random_n == 0:

            a_copy[random_n] = new_aa

            dummy_aas = deepcopy(a_copy)
            dummy_pos = deepcopy(a_pos_copy)

            a_pos_copy[random_n] = new_aa.pos

            place = random_n  # equal to 0
            while place < len(self.aas)-1 and not a_copy[place].neighbours(a_copy[place+1]):

                aa = dummy_aas[place+1]
                aa = aa.setPos(dummy_pos[place])

                a_pos_copy[place+1] = aa.pos
                a_copy[place+1] = aa

                place += 1

            return Pep(a_copy)

        # check only right aa
        elif random_n == len(self.aas)-1:

            new_aa = new_aa.setPos(random_pos)

            a_copy[random_n] = new_aa

            dummy_aas = deepcopy(a_copy)
            dummy_pos = deepcopy(a_pos_copy)

            a_pos_copy[random_n] = new_aa.pos

            place = random_n  # equal to len(self.aas)-1
            while place > 0 and not a_copy[place].neighbours(a_copy[place - 1]):
                aa = dummy_aas[place - 1]
                aa = aa.setPos(dummy_pos[place])

                a_pos_copy[place - 1] = aa.pos
                a_copy[place - 1] = aa

                place -= 1

            return Pep(a_copy)

        # check both

        else:
            # while there is no aa to be neighbour, generate a new position
            try:
                chosen_n = [a for a in range(1,len(a_copy)-1)]
                while not any([new_aa.neighbours(a_copy[random_n+1]), new_aa.neighbours(a_copy[random_n-1])]):
                    if len(chosen_n) == 0:
                        print("break loop")
                        return self
                    if len(possible_position) == 0:
                        print("happend")
                        random_n = random.choice(chosen_n)
                        chosen_n.remove(random_n)

                        new_aa = a_copy[random_n]
                        possible_position = new_aa.poss_neigh(a_copy)

                    random_pos = random.choice(possible_position)
                    possible_position.remove(random_pos)
                    new_aa = new_aa.setPos(random_pos)


            except IndexError:
                print(random_n)

            # print("right,left")
            # print(new_aa.neighbours(self.aas[random_n+1]), new_aa.neighbours(self.aas[random_n-1]))
            # print("!!!", self.aas[random_n + 1], self.aas[random_n - 1])
            # print("!",new_aa)

            # if the aa to the left is still neighbours then pull all the right aas to the left
            if new_aa.neighbours(a_copy[random_n - 1]):

                new_aa = new_aa.setPos(random_pos)

                a_copy[random_n] = new_aa

                dummy_aas = deepcopy(a_copy)
                dummy_pos = deepcopy(a_pos_copy)

                place = random_n
                while place < len(a_copy)-1 and not a_copy[place].neighbours(a_copy[place + 1]):
                    aa = dummy_aas[place + 1]
                    aa = aa.setPos(dummy_pos[place])

                    a_pos_copy[place + 1] = aa.pos
                    a_copy[place + 1] = aa

                    place += 1

                a_pos_copy[random_n] = new_aa.pos

                return Pep(a_copy)

            # if the aa to the right is still neighbours then pull all the left aas to the right
            elif new_aa.neighbours(a_copy[random_n + 1]):
                new_aa.setPos(random_pos)
                a_copy[random_n] = new_aa

                dummy_aas = deepcopy(a_copy)
                dummy_pos = deepcopy(a_pos_copy)

                place = random_n
                while place > 0 and not a_copy[place].neighbours(a_copy[place - 1]):
                    aa = dummy_aas[place - 1]
                    aa.setPos(dummy_pos[place])

                    a_pos_copy[place - 1] = aa.pos
                    a_copy[place - 1] = aa

                    place -= 1

                a_pos_copy[random_n] = new_aa.pos

                return Pep(a_copy)

            # if both consecutive aas are already neighbours then no need to change all others aas
            else:
                a_copy[random_n] = new_aa
                return Pep(a_copy)

    def crossover(self, other_pep):
        random_aa = int(random.uniform(1, len(self.aas)-1))
        parent1 = deepcopy(self.aas)
        parent2 = deepcopy(other_pep.aas)

        child1 = parent1[:random_aa] + parent2[random_aa:]
        child2 = parent2[:random_aa] + parent1[random_aa:]

        # basic chance with no rotation
        distance = np.subtract(np.array(parent2[random_aa-1].pos), np.array(parent1[random_aa-1].pos))
        leng = random_aa
        # child 1
        done = False

        while leng < len(parent1):

            child1_aa = child1[leng]
            child2_aa = child2[leng]

            new_pos_c1 = list(np.subtract(np.array(child1_aa.pos),distance))
            new_pos_c2 = list(np.add(np.array(child2_aa.pos),distance))

            child1_aa = child1_aa.setPos(new_pos_c1)
            child2_aa = child2_aa.setPos(new_pos_c2)

            if child1_aa in child1[:random_aa+1] or child2_aa in child2[:random_aa+1]:
                return self, other_pep


            child1[leng] = child1_aa
            child2[leng] = child2_aa

            leng += 1

        return Pep(child1), Pep(child2)


    def plot(self,x,y):


        fig = plt.figure(figsize=plt.figaspect(1/y))
        ax = fig.add_subplot(1, y, x, projection='3d')
        # ax = plt.axes(projection='3d')

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
        axes.set_xlim([-5, 5])
        axes.set_ylim([-5, 5])
        axes.set_zlim([-5, 5])

        #ax.set_axis_off()
        if x == y:
            plt.show()

    def __eq__(self, other):
        return self.aas == other.aas

    # def __str__(self):
    #     return ""
    #