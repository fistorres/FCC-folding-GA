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
Pep = Pep([],"0011101000")
# Pep.plot()
Pep.mutate()
# Pep.plot()
#
#
# # print(type(newPep),type(Pep),Pep.fitness)
#
# random_fit = Pep.fitness()
# fit = newPep.fitness()
# print(random_fit,fit)

pop = Pop("0110110000101101111000001110000",100,2,0.4)

# ll = pop.selection()
# print(ll[0].fitness(),ll[1].fitness())
# for i in pop.pop:
#     print(i.fitness())


for i in Pep.aas:
    print(i)