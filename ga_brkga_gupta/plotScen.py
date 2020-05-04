import matplotlib.pyplot as pl


# Plotting
def plotBarra(scenario):
    vkm = []
    pos = 0
    qtd = int(len(kmList) / 2)

    print("vsol", vsol)
    print("qtd", qtd)

    for i in range(qtd):
        # print(i)
        vkm.append('(' + str(kmList[pos]) + ',' + str(kmList[pos + 1]) + ')')
        pos += 2
        # print(i)

    # print("vkm", vkm)
    pl.bar(vkm, vsol, color="#ff9f43")

    pl.xticks(vkm)
    pl.ylabel('Número de Posições Potenciais selecionadas')

    pl.xlabel('valores(k,m)')
    pl.title('WSN#1')
    pl.show()


def plotLinha(allPps, allKmList, allPpd):
    vkm = []
    vpps = []
    vpps_final = []
    pos = 0

    v = []
    print('allPps', allPps)
    print('len(allKmList)', len(allKmList))

    for i in range(len(allKmList)):
        for j in range(len(allPps)):
            v += allPps[j]

    v.sort()
    set(v)
    print('v.sort()', v)

    for i in range(len(allKmList)):
        print('XX', allPpd[i], allPps[i])
        pl.plot(allPpd[i], allPps[i])
        # pl.plot(allPpd[0],allPps[i], 'o')

    pl.title('WSN#1')
    pl.xticks(allPpd[0])
    pl.yticks(v)

    pl.ylabel('Número de Posições Potenciais Selecionadas')
    pl.xlabel('Número de Posições Potenciais Dadas')

    pl.legend(['k=1, m=1', 'k=2, m=1', 'k=2, m=2', 'k=3, m=2'], loc=4)
    pl.show()


def plotScenario(scenario1, indexes=True):
    # only plot the scenario  with no solution selected
    x, y = scenario1['targets'][0]
    pl.plot(x, y, 'b*', label='Target point')
    if indexes:
        index = 0
        for xx, yy in zip(x, y):
            pl.text(xx-0.2, yy-0.2, "{}".format(index), fontsize=8)
            index += 1
    # aleatorio
    # x,y = genPotential([0, 300, 0, 300], 300)
    # grid
    x, y = scenario1['potential'][0]
    pl.plot(x, y, 'ro', fillstyle='none', label='Potential Position')

    if indexes:
        index = 0
        for xx, yy in zip(x, y):
            pl.text(xx-0.2, yy-0.2, "{}".format(index), fontsize=8)
            index += 1

    pl.legend(loc='upper center', bbox_to_anchor=(0.5, -0.08),
              fancybox=True, shadow=True, ncol=2)
    pl.show()


def plotSolutionOnScenario(binarySol, scenario, **kargs):
    # only plot the scenario  with no solution selected
    plotConn = False if 'plotConn' not in kargs.keys() else kargs['plotConn']
    plotSen = False if 'plotSen' not in kargs.keys() else kargs['plotSen']
    connMatrix = [] if 'connMatrix' not in kargs.keys() else kargs['connMatrix']


    x, y = scenario['targets'][0]
    pl.plot(x, y, 'b*', label='Target point')
    x, y = scenario['potential'][0]

    if len(connMatrix) > 0:
        printed = {}
        for i in range(len(binarySol)):
            if binarySol[i]:
                for j in range(i+1, len(binarySol)):
                    if binarySol[j] and connMatrix[i, j]:
                        # yeah yeah, I know, I know, too expensive.
                        if (i, j) not in printed.keys():
                            printed[(i, j)] = True
                            pl.plot([x[i], x[j]], [y[i], y[j]], 'k-', lineWidth=1, alpha=0.8)

    xun = [x[i] for i in range(len(binarySol)) if binarySol[i] == 0]
    yun = [y[i] for i in range(len(binarySol)) if binarySol[i] == 0]
    pl.plot(xun, yun, 'ro', fillstyle='none', alpha=0.5, label='Potential Position')

    xsel = [x[i] for i in range(len(binarySol)) if binarySol[i] == 1]
    ysel = [y[i] for i in range(len(binarySol)) if binarySol[i] == 1]
    pl.plot(xsel, ysel, 'o', color='tab:red', alpha=1, markeredgecolor='r', label='Placed Sensor')

    ppd = len(xun) + len(xsel)
    # print('tam PD', ppd)

    pps = len(xsel)
    # print('tam PS', pps)


    if plotConn:
        Rcomm = scenario['Rcomm']
        for xsel_x, ysel_y in zip(xsel, ysel):
            pl.gca().add_artist(pl.Circle((xsel_x, ysel_y), Rcomm, color='g', fill=True, alpha=0.08))
        #pl.plot(xsel, ysel, 'o', color='none', alpha=0.8, markersize=Rcomm, markeredgecolor='g')

    if plotSen:
        Rsen = scenario['Rsen']
        #pl.plot(xsel, ysel, 'o', color='tab:red', alpha=0.15, markersize=Rsen, markeredgecolor='r')
        for xsel_x, ysel_y in zip(xsel, ysel):
            pl.gca().add_artist(pl.Circle((xsel_x, ysel_y), Rsen, color='tab:green', fill=True, alpha=0.08))

    pl.legend(loc='upper center', bbox_to_anchor=(0.5, -0.08),
              fancybox=True, shadow=True, ncol=3)

    pl.show()
    #print("plotou")
    return ppd, pps

def plot_hist(countConex,countSensorUtilizados,countSensorDescobertos,countTargetDescobertos,quantPot):
    import matplotlib.pyplot as pl
    pl.bar(list(countConex.keys()), countConex.values(), color='#607c8e')
    pl.bar(list(countSensorUtilizados.keys()), countSensorUtilizados.values(), color='#FFB6C1')
    pl.bar(list(countSensorDescobertos.keys()), countSensorDescobertos.values(), color='#D2691E')
    pl.bar(list(countTargetDescobertos.keys()), countTargetDescobertos.values(), color='#3CB371')
    #pl.title('Qualidade de Soluções - WSN#'+str(ix)+'-'+str(quantPot)+' - GA - fit')
    pl.title('Qualidade de Soluções - WSN#2' + '-' + str(quantPot) + ' - BRKGA - modFIT')
    pl.xlabel('Qtd. pontos')
    pl.ylabel('Qtd.de Soluções')
    pl.grid(axis='y', alpha=0.75)
    pl.legend(['Componentes Conexas', 'Sensores Utilizados', 'Sensor Descobertos','Target Descobertos'], loc=9)

    pl.show()

# Plotar Gráfico
def plot_log(logbook):
    gen = logbook.select("gen")
    min = logbook.select("min")
    avg = logbook.select("avg")
    max = logbook.select("max")

    import matplotlib.pyplot as plt

    fig, ax1 = plt.subplots()
    line1 = ax1.plot(gen, min, "b-", label="Minimum Fitness")
    ax1.set_xlabel("Generation")
    ax1.set_ylabel("Fitness", color="b")
    for tl in ax1.get_yticklabels():
        tl.set_color("b")

    ax2 = ax1.twinx()
    line2 = ax2.plot(gen, avg, "g-", label="Average Fitness")
    for tl in ax2.get_yticklabels():
        tl.set_color("g")

    ax3 = ax1.twinx()
    line3 = ax3.plot(gen, max, "y-", label="Maximum Fitness")
    ax3.set_ylabel("Size")
    for tl in ax3.get_yticklabels():
        tl.set_color("y")

    lns = line1 + line2 + line3
    labs = [l.get_label() for l in lns]
    ax1.legend(lns, labs, loc="center right")

    plt.show()

def plot_teste():
    import pandas as pd

    # Generate data on commute times.
    dicionario = {1: 30, 2: 28, 3: 3}
    size, scale = 1000, 10
    commutes = pd.Series(dicionario.keys(), dicionario.values())


    commutes.plot.hist(grid=True, bins=20, rwidth=0.9,
                       color='#607c8e')
    pl.title('Commute Times for 1,000 Commuters')
    pl.xlabel('Counts')
    pl.ylabel('Commute Time')
    pl.grid(axis='y', alpha=0.75)

    pl.show()

def plotLinha2(allKmList,allPps,allPpd):
    print("no plot linhaa")
    v = []
    print('allPps', allPps)
    print('len(allKmList)', len(allKmList))

    #for i in range(len(allKmList)):
    for j in range(len(allPps)):
        v += allPps[j]

    print('v', v)
    v.sort()
    set(v)
    print('v.sort()', v)

    for i in range(len(allKmList)):
        print('XX', allPpd[i], allPps[i])
        pl.plot(allPpd[i], allPps[i])
        # pl.plot(allPpd[0],allPps[i], 'o')

    pl.title('WSN#1')
    pl.xticks(allPpd[0])
    pl.yticks(v)

    pl.ylabel('Número de Posições Potenciais Selecionadas')
    pl.xlabel('Número de Posições Potenciais Dadas')

   # pl.legend(['k=1, m=1', 'k=2, m=1', 'k=2, m=2', 'k=3, m=2'], loc=4)
    pl.legend(['k=2, m=2'], loc=1)
    pl.show()
