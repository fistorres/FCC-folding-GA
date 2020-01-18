from Ga import *
from main import *

### Benchmarking sequences
# length of 20
S1 = "HPHPPHHPHPPHPHHPPHPH"
# length of 48
H1 = "HPHHPPHHHHPHHHPPHHPPHPHHHPHPHHPPHHPPPHPPPPPPPPHH"

# length of 90
F_90_1 = "PPHHHPPPHHPPPPHHPHHHHHHPHPHPHHPHHHHHPHHHPHPHHHHPHHPPPPHHHPHPHPPHHHPHHPHPHPPHHHPPPPHHPPHPPP"
# length of 135
S1_2 = "HHHHPHHHHHHPPHHPHHHHHHHHPHHPHHHHHHHHHHPPHHPPPPPHHPHHHHHHHHPHPPHH" \
       "PPPHHHHHHHHPHHHHHHPPHHHHHHHPHPPHHHHHHHHHPPHHPPPHHHHHHHPHHPHHHHHHHPPHHHH"
# length of 200
R1 = "PPPHPHHPHHPPPHPHPPPPHPHHPPHPHHHHHPPHH" \
     "PPHHHHHHPPHPPHHPPHPHPHHHHHPHHPHHHPPPHHHPHHPPHPHPPHPPP" \
     "HPPHPPHPPHHHPHHHPHPPHPHHPHHHHPHPHHHPHHHPPPPPPHHHHHHPPPPP" \
     "PPPHHHPPHPHPPPHPHPHPHHPPHHPPPPHHHHHHPPPHHPPPPPHPPPHHPP"


# demonstrating neighbours
def plotneighbours():
    vect = [[1, 1, 0], [1, -1, 0], [-1, -1, 0], [-1, 1, 0],
            # up
            [1, 0, 1], [0, -1, 1], [-1, 0, 1], [0, 1, 1],
            # down
            [1, 0, -1], [0, -1, -1], [-1, 0, -1], [0, 1, -1]]
    ax = plt.axes(projection='3d')

    zline, xline, yline = [], [], []

    for point in vect:
        print(point)

        xline.append(point[0])
        yline.append(point[1])
        zline.append(point[2])

    xline.append(0)
    yline.append(0)
    zline.append(0)

    # lines
    # ax.plot3D(xline, yline, zline)
    # points
    ax.scatter3D(xline, yline, zline)
    # ax.scatter3D(xline[0], yline[0], zline[0])
    ax.scatter3D(xline[-1], yline[-1], zline[-1])

    axes = plt.gca()
    axes.set_xlim([-2, 2])
    axes.set_ylim([-2, 2])
    axes.set_zlim([-2, 2])

    plt.show()

# demostrating mutation
def plotmutation():
    Pep1 = Pep([],"0101101001")

    Pep2 = Pep1.mutate()

    Pep1.plot(1,2)
    Pep2.plot(2,2)

    for a,b in zip(Pep1.aas,Pep2.aas):
        print(a,b)

# demonstrating crossover
def plotcrossover():
    Pep1 = Pep([],"010111110101001")
    Pep2 = Pep([],"010110111001001")

    c1,c2 = Pep1.crossover(Pep2)

    Pep1.plot(1,4)
    Pep2.plot(2,4)
    c1.plot(3,4)
    c2.plot(4,4)

    for a,b,c,d in zip(Pep1.aas,Pep2.aas,c1.aas,c2.aas):
        print(a,b,c,d)
