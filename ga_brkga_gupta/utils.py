import numpy as np
import networkx as nx
# Implementation of Gupta 2016
def dist(x1, y1, x2, y2):
    return np.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2)

def getNumOfCC(ind, CommMatrix):#ok
    index = np.array(ind, np.bool)
    indexes = np.array([i for i, bit in zip(range(len(index)), index) if bit])
    A = CommMatrix[:, indexes]
    A = A[indexes, :]
    G = nx.from_numpy_matrix(A)
    return nx.number_connected_components(G)

def getNumSensorDescobertos(ind, sim):
    todos=[]
    k, teste = sim['k'], Teste(sim['scenario'])
    ativos = np.array(ind, np.bool)
    indexes_ativos = np.array([i for i, bit in zip(range(len(ativos)), ativos) if bit])
    for i in range(len(ind)):
        if ind[i] == 1:
            index_teste = np.array(teste[i,:], np.bool)
            indexes_teste = np.array([i for i, bit in zip(range(len(index_teste)), index_teste) if bit])
            for x in range(len(indexes_teste)):
                 for indice2 in indexes_teste:
                     if indice2 not in todos:
                         todos.append(indice2)
    cobertos = np.setdiff1d(todos, indexes_ativos)
    qtdSDescobertos =len(ind) - len(cobertos)-len(indexes_ativos)
    return qtdSDescobertos



def getNumTargetDescobertos(ind, sim):#0k
    N, k, SenMatrix = sim['scenario']['quantTargets'], sim['k'], sim['SenMatrix']
    cov = [sum(np.multiply(SenMatrix[i, :], ind))for i in range(N)]
    index = np.array(cov, np.bool)
    covDescobertos = np.array([i for i, bit in zip(range(len(index)), index) if bit==False])
    return len(covDescobertos)

# retorna uma matriz binaria de sensoriamento entre ntargets x potential
def initSenMatrix(scenario):
    Xt, Yt = scenario['targets'][0]
    Xpt, Ypt = scenario['potential'][0]
    M = np.zeros((len(Xt), len(Xpt)), dtype=np.int32)
    for t in range(scenario['quantTargets']):
        xt, yt = Xt[t], Yt[t]
        for p in range(scenario['quantPot']):
            xpt, ypt = Xpt[p], Ypt[p]
            d = dist(xt, yt, xpt, ypt)
            # print("Dist: target {} to potential: {} is: {}".format(t, p, d))
            if dist(xt, yt, xpt, ypt) <= scenario['Rsen']:
                M[t, p] = 1.0
    #print("M",M)
    return M

# retorna uma matriz de booleanos onde i,j = 1, se i se conecta com j
def initCommMatrix(scenario):
    Xpt, Ypt = scenario['potential'][0]
    M = np.zeros((len(Xpt), len(Xpt)), dtype=np.int32)
    for p1 in range(len(Xpt)):
        x1, y1 = Xpt[p1], Ypt[p1]
        for p2 in range(len(Xpt)):
            x2, y2 = Xpt[p2], Ypt[p2]
            # p1 should be different from p2 in order to avoid misinformation about p1 being closed to p1
            M[p1, p2] = (dist(x1, y1, x2, y2) <= scenario['Rcomm']) and p2 != p1
    return M

def Teste(scenario):
    todos = []
    Xpt, Ypt = scenario['potential'][0]
    M = np.zeros((len(Xpt), len(Xpt)), dtype=np.int32)
    for p1 in range(len(Xpt)):
        x1, y1 = Xpt[p1], Ypt[p1]
        for p2 in range(len(Xpt)):
            x2, y2 = Xpt[p2], Ypt[p2]
            # p1 should be different from p2 in order to avoid misinformation about p1 being closed to p1
            M[p1, p2] = (dist(x1, y1, x2, y2) <= scenario['Rsen']) and p2 != p1
    return M