import numpy as np
#from scen import genDumbScenario1
from scen import genDumbScenario_wsn, genDumbScenario_toy
from plotScen import plotScenario, plotSolutionOnScenario,plot_hist
from utils import *
import matplotlib.pyplot as pl
import sys
np.set_printoptions(threshold=sys.maxsize)
from deap import algorithms
from deap import base
from deap import creator
from deap import tools
import random
import networkx as nx

WSN1, WSN2 = genDumbScenario_wsn()
sim = 0
global quantPot

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
    #print("covCost", covCost)
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

def chaveIndividual(sim):
    global quantPot
    chave = []
    for i in range(quantPot):
        chave.append(random.uniform(0, 1))
    return chave

# 1. Decoder
def decode(chave, sim):
    pindice = 0.3 #resolver
    # 2. Order in ascending order the initial vector, represented by chromosome
    chave_a_ordenados = sorted(chave)
    #3.Check in the initial vector the indexes of the first 3 values of the ordered vector;
    indices = chave_a_ordenados[:(int(pindice* sim['scenario']['quantPot']))]
    #4. Define a new vector and assign the value 1 in the key figures verified in the previous step
    # and value 0 for the other positions of the new vector;
    vetor_decodificado = [0] * sim['scenario']['quantPot']
    for i in range(len(indices)):
        x = indices[i]
        posicao = chave.index(x)
        vetor_decodificado[posicao] = 1
    #5.return the chromosome and calculate the fitness
    return vetor_decodificado

creator.create("FitnessMax", base.Fitness,weights=(1.0,))
creator.create("Individual", list, typecode='b', fitness=creator.FitnessMax)

"""#Simulations"""
def main(WSN1, WSN2):
    kmList = []
    allKmList = []
    allPpd = []
    allPps = []
    ppd = []
    pps = []

    random.seed(32)

    kmSim1 = [(2, 2)]
    #print("atual 1", WSN1)
   # WSN = [WSN1, WSN2]
    WSN = [WSN2]
    summaryResults = {}

    for wsn in WSN:
        for k, m in kmSim1:
            for scen in wsn:
                sim = {
                    'GAParams': {'ngen': 50, 'popSize': 60, 'cxProb': 0.9, 'mutProb': 0.03,
                                 'fitFunction':modFit,'peSize':0.25,'pmSize':0.25,'pENEsize':0.5,
                                 'porIndice':0.3, 'decodeFunction':decode},
                    'k': k,
                    'm': m,
                    'scenario': scen,
                    'SenMatrix': initSenMatrix(scen),
                    'CommMatrix': initCommMatrix(scen),
                    #'weights': [0.4, 0.3, 0.3, 0.0]#fit2
                    'weights': (0.25, 0.25, 0.25, 0.25),  # modfit
                }
                countConex = {}
                countSensorUtilizados= {}
                countSensorDescobertos= {}
                countTargetDescobertos= {}
                quantPot =scen['quantPot']
                for exec in range(20):
                    toolbox = base.Toolbox()
                    # Attribute generator
                    toolbox.register("attr_bool", random.randint, 0, 1)
                    toolbox.register("attr_key", random.uniform, 0, 1)

                    # Structure initializers
                    toolbox.register("individual", tools.initRepeat, creator.Individual, toolbox.attr_bool,
                                     scen['quantPot'])
                    #toolbox.register("chave_individual", chaveIndividual,creator.Individual)
                    toolbox.register("chave_individual",  tools.initRepeat, creator.Individual,toolbox.attr_bool,
                                     scen['quantPot'])
                    toolbox.register("population", tools.initRepeat, list, toolbox.individual)
                    toolbox.register("pop_decoder", tools.initRepeat, list, toolbox.individual)
                    toolbox.register("population_chave", tools.initRepeat, list, toolbox.chave_individual)
                    # toolbox.register("mate", tools.cxTwoPoint)
                    # toolbox.register("mate", mate_uniforme_parametrizado,cxpb=sim['GAParams']['cxProb'] )
                    toolbox.register("mutate", mutZeroBit, indpb=sim['GAParams']['mutProb'])
                    # toolbox.register("select", tools.selRoulette)
                    # toolbox.register("select",selectParents_brkga,peSize=sim['GAParams']['peSize'],  sim=sim)
                    toolbox.register("evaluate", sim['GAParams']['fitFunction'], sim=sim)
                    # print("Running for (k={},m={})".format(k, m))
                    # Decodificador - BRKGA
                    # toolbox.register("decode",sim['GAParams']['decodeFunction'],ipb = sim['GAParams']['porIndice'],sim=sim)
                    toolbox.register("decode", sim['GAParams']['decodeFunction'], sim=sim)
                    pop = toolbox.population(n=sim['GAParams']['popSize'])
                    popChave = toolbox.population_chave(n=sim['GAParams']['popSize'])

                    hof = tools.HallOfFame(3)
                    stats = tools.Statistics(lambda ind: ind.fitness.values)
                    stats.register("avg", np.mean)
                    stats.register("std", np.std)
                    stats.register("min", np.min)
                    stats.register("max", np.max)

                    pop, log = algorithms.eaSimple_BRKGA(pop, popChave, toolbox, cxpb=sim['GAParams']['cxProb'],
                                                   mutpb=sim['GAParams']['mutProb'],
                                                   ngen=sim['GAParams']['ngen'],
                                                   peSize=sim['GAParams']['peSize'],
                                                   pmSize=sim['GAParams']['pmSize'],
                                                   quantPot=scen['quantPot'],
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
                    numCC = getNumOfCC(hof[0], sim['CommMatrix'])#ok
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


                plot_hist(countConex, countSensorUtilizados,countSensorDescobertos, countTargetDescobertos, quantPot)

if __name__ == "__main__":
    WSN1,WSN2 = genDumbScenario_wsn()
    main(WSN1, WSN2)