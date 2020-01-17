from pop import *
from main import *

a = Aa("0",0,[0,0,0])
b = Aa("1",1,[1,1,0])
c = Aa("1",2,[1,2,1])
d = Aa("0",3,[1,3,1])
new = [a,b,c,d]
#
# newPep = Pep(new)
# newPep.aaspos
Pep1 = Pep([],"0011101000")
Pep2 = Pep([],"0011101000")
# Pep.plot()
# Pep.plot()
#
#
# # print(type(newPep),type(Pep),Pep.fitness)
#
# random_fit = Pep.fitness()
# fit = newPep.fitness()
# print(random_fit,fit)

pop = Pop.random_init("0110110000101101111000001110000",100,3,0.4)

# Pep2 = Pep.mutate()

for a,b in zip(Pep2.aas,Pep1.aas):
    print(a,b)



pep1, pep2 = Pep1.crossover(Pep2)

for a,b in zip(pep2.aas,pep1.aas):
    print("child2:",a,"child1:",b)


for a in range(len(pep2.aas)-1):
    print(pep2.aas[a].neighbours(pep2.aas[a+1]))
    print(pep1.aas[a].neighbours(pep1.aas[a+1]))


Pep1.plot(1,4)
Pep2.plot(2,4)
pep1.plot(3,4)
pep2.plot(4,4)