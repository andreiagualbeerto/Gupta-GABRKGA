import numpy as np
#from scen import genDumbScenario1
from scen import genDumbScenario_wsn, genDumbScenario_toy
from plotScen import plotScenario, plotSolutionOnScenario,plot_hist
from utils import *
import matplotlib.pyplot as plt

from deap import algorithms
from deap import base
from deap import creator
from deap import tools
import random
import networkx as nx
scenario = genDumbScenario_toy()

plotScenario(scenario)

sim = 0

def modFit(ind, sim):
    '''
    Esta versão de função fitness possui quatro pesos. O último peso agora diz respeito ao grafo de connectividade, ou seja,
    retorna o número de componentes conexas do grafo de proximidade induzido pela posição dos pontos selecionados.
    Como a função é de maximização, então queremos maximizar o termo 1/(# de componentes conexas), ou seja, o grafo que tiver
    apenas 1 componente é o de máximo valor do termo.
    '''

    numCC = getNumOfCC(ind, sim['CommMatrix'])
    F4 = 1.0/numCC
    w1, w2, w3, w4 = sim['weights']
    FOther, = fit2(ind, sim)
    return FOther+w4*F4,

def fit2(ind, sim):
    f1 = sum(ind) / len(ind)
    N, k, SenMatrix = sim['scenario']['quantTargets'], sim['k'], sim['SenMatrix']
    cov = [sum(np.multiply(SenMatrix[i, :], ind)) for i in range(N)]
    covCost = [min(k, cov[i])/k for i in range(N)]
    f2 = (1.0/N) * sum(covCost)
    m, CommMatrix = sim['m'], sim['CommMatrix']
    com = [sum(np.multiply(CommMatrix[i, :], ind)) for i in range(len(ind)) if ind[i] == 1]
    commCost = [min(m, com[i])/m for i in range(len(com))]
    f3 = (1.0/sum(ind)) * sum(commCost)
    # print('F3: {}'.format(f3))
    w1, w2, w3, w4 = sim['weights']
    #plotSolutionOnScenario(ind, sim['scenario'], plotConn=False)
    return w1*(1.0 - f1) + w2*(f2) + w3*f3,

def fit(ind, sim):
    f1 = sum(ind) / len(ind)
    N, k, SenMatrix = sim['scenario']['quantTargets'], sim['k'], sim['SenMatrix']
    cov = [sum(np.multiply(SenMatrix[i, :], ind)) for i in range(N)]
    covCost = [k if cov[i] >= k else k - cov[i] for i in range(N)]
    f2 = (1.0 / (N * k)) * sum(covCost)
    m, CommMatrix = sim['m'], sim['CommMatrix']
    com = [sum(np.multiply(CommMatrix[i, :], ind)) for i in range(len(ind)) if ind[i] == 1]
    # as proposed by the paper
    commCost = [m if com[i] >= m else m - com[i] for i in range(len(com))]
    # as proposed by the paper
    f3 = (1.0 / (sum(ind) * m)) * sum(commCost)
    w1, w2, w3, w4 = sim['weights']
    #plotSolutionOnScenario(ind, sim['scenario'], plotConn=False)
    return w1*(1.0 - f1) + w2*(f2) + w3*f3,

def mutZeroBit(individual, indpb):
    for i in range(len(individual)):
        if random.random() < indpb:
            individual[i] = type(individual[i])(0)
            #print('Mutate on bit %d' % i)

    return individual,

# creator.create("FitnessMax", base.Fitness, weights=(0.3, 0.4, 0.3))

k, m = 2, 2
sim = {
    'GAParams': {'ngen': 50, 'popSize': 60, 'cxProb': 0.9, 'mutProb': 0.03,
                 'fitFunction':fit2},
    'k': k,
    'm': m,
    'scenario': scenario,
    'SenMatrix': initSenMatrix(scenario),
    'CommMatrix': initCommMatrix(scenario),
    #'weights': (0.2, 0.2, 0.2, 0.4)
    'weights': (0.4, 0.3, 0.3, 0.0),#fit2
    #'weights':(0.25,0.25,0.25,0.25), #modfit
}
creator.create("FitnessMax", base.Fitness, weights=(1.0,))
# creator.create("FitnessMax", base.Fitness, weights=(0.02, 0.02, 0.02, 0.94))
# creator.create("Individual", array.array, typecode='b', fitness=creator.FitnessMax)
creator.create("Individual", list, typecode='b', fitness=creator.FitnessMax)
#def getGfromIndividual(ind, scen):

countConex = {}
countSensorUtilizados= {}
countSensorDescobertos= {}
countTargetDescobertos= {}
countTarget={}
for exec in range(20):

    toolbox = base.Toolbox()
    # Attribute generator
    toolbox.register("attr_bool", random.randint, 0, 1)
    # Structure initializers
    toolbox.register("individual", tools.initRepeat, creator.Individual, toolbox.attr_bool, n=scenario['quantPot'])
    toolbox.register("population", tools.initRepeat, list, toolbox.individual)
    toolbox.register("mate", tools.cxTwoPoint)
    toolbox.register("mutate", mutZeroBit, indpb=sim['GAParams']['mutProb'])
    toolbox.register("select", tools.selRoulette)
    #toolbox.register("select", tools.selTournament, tournsize=3)
    toolbox.register("evaluate", sim['GAParams']['fitFunction'], sim=sim)
    #print("Running for (k={},m={})".format(k, m))
    pop = toolbox.population(n=sim['GAParams']['popSize'])
    hof = tools.HallOfFame(3)
    stats = tools.Statistics(lambda ind: ind.fitness.values)
    stats.register("avg", np.mean)
    stats.register("std", np.std)
    stats.register("min", np.min)
    stats.register("max", np.max)

    pop, log = algorithms.eaSimple(pop, toolbox, cxpb=sim['GAParams']['cxProb'], mutpb=sim['GAParams']['mutProb'],
                                   ngen=sim['GAParams']['ngen'],
                                   stats=stats, halloffame=hof, verbose=False)
    #print(hof[0])
    plotSolutionOnScenario(hof[0], sim['scenario'], plotSen=True, connMatrix=sim['CommMatrix'], plotConn=True)
    #ha = [1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    #print([(sim['GAParams']['fitFunction'](i, sim)) for i in hof])
    # Cálculo do número de componentes conexas
    numCC = getNumOfCC(hof[0], sim['CommMatrix'])
    numSU=sum(hof[0])#ok
numSD = getNumSensorDescobertos(hof[0],sim)
numTD = getNumTargetDescobertos(hof[0],sim)

#print("num cc no fora", numCC)
#x, y = scenario['potential'][0]
#nx.draw(G, [(x1, y1) for x1, y1, i in zip(x, y, hof[0]) if i])
#plt.show()
if numCC in countConex.keys():
    countConex[numCC] += 1
else:
    countConex[numCC] = 1

if numSU in countSensorUtilizados.keys():
    countSensorUtilizados[numSU] += 1
else:
    countSensorUtilizados[numSU] = 1

if numSD in countSensorDescobertos.keys():
    countSensorDescobertos[numSD] += 1
else:
    countSensorDescobertos[numSD] = 1

if numTD in countTargetDescobertos.keys():
    countTargetDescobertos[numTD] += 1
else:
    countTargetDescobertos[numTD] = 1

print("Qualidade da solução - BRKGA - TOY\n",
      "Número de  número de componentes conexas:",numCC,"\n",
      "Número de Sensores Utilizados:",numSU,"\n",
      "Número de alvos descobertos (total ou parcial):",numTD,"\n",
      "Número de Sensores descobertos:",numSD,"\n")

print("countConex",countConex)
print("countSensorUtilizados",countSensorUtilizados)
print("countSensorDescobertos",countSensorDescobertos)
print("countTargetDescobertos",countTargetDescobertos)


plot_hist(countConex, countSensorUtilizados,countSensorDescobertos, countTargetDescobertos)


print(countConex)
#print(fitValue(fit2(hof[0], sim), sim['weights']))
