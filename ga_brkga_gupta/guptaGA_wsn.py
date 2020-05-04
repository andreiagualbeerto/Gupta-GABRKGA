import numpy as np
#from scen import genDumbScenario1
from scen import genDumbScenario_wsn, genDumbScenario_toy
from plotScen import plotScenario, plotSolutionOnScenario,plot_hist
from utils import *
import matplotlib.pyplot as pl

from deap import algorithms
from deap import base
from deap import creator
from deap import tools
import random
import networkx as nx

WSN1, WSN2 = genDumbScenario_wsn()
sim = 0


def modFit(ind, sim):
    '''
    Esta versão de função fitness possui quatro pesos. O último peso agora diz respeito ao grafo de connectividade, ou seja,
    retorna o número de componentes conexas do grafo de proximidade induzido pela posição dos pontos selecionados.
    Como a função é de maximização, então queremos maximizar o termo 1/(# de componentes conexas), ou seja, o grafo que tiver
    apenas 1 componente é o de máximo valor do termo.
    '''
    CommMatrix = sim['CommMatrix']
    # CommMatrix after multiplied each row or column by the ind
    A = np.multiply(np.array(ind), CommMatrix)
    G = nx.from_numpy_matrix(A)
    F4 = 1.0/nx.number_connected_components(G)
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

creator.create("FitnessMax", base.Fitness, weights=(1.0,))
creator.create("Individual", list, typecode='b', fitness=creator.FitnessMax)

"""#Simulations"""
def main(WSN1, WSN2):
    kmList = []
    ppd = []
    pps = []

    random.seed(32)

    kmSim1 = [(2, 2)]
    #print("atual 1", WSN1)
    #WSN = [WSN1, WSN2]
    #WSN = [WSN2]
    WSN = [WSN1]
    summaryResults = {}
    ix=0
    for wsn in WSN:
        ix=+1
        for k, m in kmSim1:
            for scen in wsn:
                sim = {
                    'GAParams': {'ngen': 50, 'popSize': 60, 'cxProb': 0.9, 'mutProb': 0.03,'fitFunction':modFit},
                    'k': k,
                    'm': m,
                    'scenario': scen,
                    'SenMatrix': initSenMatrix(scen),
                    'CommMatrix': initCommMatrix(scen),
                    #'weights': [0.4, 0.3, 0.3, 0.0]
                    'weights':(0.25,0.25,0.25,0.25), #modfit
                }


                countConex = {}
                countSensorUtilizados= {}
                countSensorDescobertos= {}
                countTargetDescobertos= {}
                quantPot = scen['quantPot']
                for exec in range(20):
                    toolbox = base.Toolbox()
                    # Attribute generator
                    toolbox.register("attr_bool", random.randint, 0, 1)
                    # Structure initializers
                    toolbox.register("individual", tools.initRepeat, creator.Individual, toolbox.attr_bool,
                                     scen['quantPot'])
                    toolbox.register("population", tools.initRepeat, list, toolbox.individual)
                    toolbox.register("mate", tools.cxTwoPoint)
                    toolbox.register("mutate", mutZeroBit, indpb=sim['GAParams']['mutProb'])
                    toolbox.register("select", tools.selRoulette)
                    toolbox.register("evaluate", sim['GAParams']['fitFunction'], sim=sim)
                    pop = toolbox.population(n=sim['GAParams']['popSize'])

                    hof = tools.HallOfFame(3)
                    stats = tools.Statistics(lambda ind: ind.fitness.values)
                    stats.register("avg", np.mean)
                    stats.register("std", np.std)
                    stats.register("min", np.min)
                    stats.register("max", np.max)

                    pop, log = algorithms.eaSimple(pop, toolbox, cxpb=sim['GAParams']['cxProb'],
                                                   mutpb=sim['GAParams']['mutProb'],
                                                   ngen=sim['GAParams']['ngen'],
                                                   stats=stats, halloffame=hof, verbose=False)
                    print("Melhor Indivíduo", hof[0])

                    summaryResults[(k, m)] = sum(hof[0])
                    print("km={},{} best: {}, Pot:{} Targets:{}".format(k, m, sum(hof[0]), scen['quantPot'],
                                                                        scen['quantTargets']))

                    kmList.append(k)
                    kmList.append(m)
                    d,s = plotSolutionOnScenario(hof[0], scen, plotConn=True,connMatrix=sim['CommMatrix'], plotSen=True)
                    ppd.append(d)
                    pps.append(s)

                    # Cálculo do número de componentes conexas
                    numCC = getNumOfCC(hof[0], sim['CommMatrix'])  # ok
                    numSU = sum(hof[0])  # ok
                    numSD = getNumSensorDescobertos(hof[0], sim)
                    numTD = getNumTargetDescobertos(hof[0], sim)

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
                          "Número de  número de componentes conexas:", numCC, "\n",
                          "Número de Sensores Utilizados:", numSU, "\n",
                          "Número de alvos descobertos (total ou parcial):", numTD, "\n",
                          "Número de Sensores descobertos:", numSD, "\n")

                print("countConex", countConex)
                print("countSensorUtilizados", countSensorUtilizados)
                print("countSensorDescobertos", countSensorDescobertos)
                print("countTargetDescobertos", countTargetDescobertos)

                plot_hist(countConex, countSensorUtilizados,countSensorDescobertos, countTargetDescobertos, ix,quantPot)

if __name__ == "__main__":
    WSN1,WSN2 = genDumbScenario_wsn()
    main(WSN1, WSN2)