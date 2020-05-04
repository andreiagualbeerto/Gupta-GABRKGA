# generating the scenarios
import numpy as np
import copy


def genDumbScenario_wsn():
    # WSN1 Grid scenario
    scenario1 = genScenario(limits=[0, 300, 0, 300], grid=True, quantPot=100, quantTargets=100, Rsen=50, Rcomm=100)
    WSN1 = [scenario1]
    for qtPt in range(200, 501, 100):
        scen = copy.deepcopy(scenario1)
        qt = int(np.sqrt(qtPt))
        qt = qt ** 2
        scen['quantPot'] = qt
        scen['potential'] = [genPotential(limits=[0, 300, 0, 300], grid=True, quantPot=qt)]
        WSN1.append(scen)

    #print("WSN1", WSN1)


    # WSN2
    scenario2 = genScenario(limits=[0, 300, 0, 300], quantPot=100, quantTargets=100, Rsen=50, Rcomm=100)
    WSN2 = [scenario2]
    for qtPt in range(200, 501, 100):
        scen = copy.deepcopy(scenario2)
        scen['quantPot'] = qtPt
        scen['potential'] = [genPotential(limits=[0, 300, 0, 300], grid=False, quantPot=qtPt)]
        WSN2.append(scen)

    #print("WSN2",WSN2)
    #print("WSN2type", type(WSN2))

    return WSN1,WSN2

def genDumbScenario_toy():
    quantPot = 62
    quantTargets = 4
    Rsen = 1.0
    Rcomm = 1.8
    limits = [0, 16, 0, 16]
    X_target = [3, 3, 11, 11]
    Y_target = [11, 14, 11, 14]
    x_pot = []
    y_pot = []
    for x in range(2, 13):
        for y in range(10, 16):
            app = True
            # This is to avoid the generation
            #I know, too expensive...
            for x_t, y_t in zip(X_target, Y_target):
                if x_t == x and y_t == y:
                    app = False
                    break
            if app:
                x_pot.append(x)
                y_pot.append(y)

    X = {
        'quantPot': quantPot,
        'quantTargets': quantTargets,
        'Rsen': Rsen,
        'Rcomm': Rcomm,
        'limits': limits,
        'targets': [(X_target, Y_target)],
        'potential': [(x_pot, y_pot)],
    }

    return X


def genTargets(**kargs):
    limits = [0, 300, 0, 300] if 'limits' not in kargs.keys() else kargs['limits']
    quantTargets = 10 if 'quantTargets' not in kargs.keys() else kargs['quantTargets']
    xinf, xsup, yinf, ysup = limits

    X = (xsup - xinf) * np.random.random_sample((quantTargets,)) + xinf
    Y = (ysup - yinf) * np.random.random_sample((quantTargets,)) + yinf
    return (X, Y)


def genPotential(**kargs):
    # default params
    limits = [0, 300, 0, 300] if 'limits' not in kargs.keys() else kargs['limits']
    quantPot = 10 if 'quantPot' not in kargs.keys() else kargs['quantPot']
    grid = False if 'grid' not in kargs.keys() else kargs['grid']
    space = 25 if 'space' not in kargs.keys() else kargs['space']

    xinf, xsup, yinf, ysup = limits

    X, Y = np.array([]), np.array([])
    if grid:
        x, y = xinf + space, yinf + space
        xvector = np.linspace(xinf, xsup, np.sqrt(quantPot))
        yvector = np.linspace(yinf, ysup, np.sqrt(quantPot))
        X, Y = np.meshgrid(xvector, yvector)
        X = X.flatten()
        Y = Y.flatten()
    else:
        kargs['quantTargets'] = quantPot
        return genTargets(**kargs)

    return X, Y

def genScenario(**kargs):
    quantPot = 10 if 'quantPot' not in kargs.keys() else kargs['quantPot']
    quantTargets = 10 if 'quantTargets' not in kargs.keys() else kargs['quantTargets']
    Rsen = 30 if 'Rsen' not in kargs.keys() else kargs['Rsen']
    Rcomm = 10 if 'Rcomm' not in kargs.keys() else kargs['Rcomm']
    limits = [0, 300, 0, 300] if 'limits' not in kargs.keys() else kargs['limits']
    if 'grid' in kargs and kargs['grid']:
        quantPot = int(np.sqrt(quantPot)) ** 2
    X = {
        'quantPot': quantPot,
        'quantTargets': quantTargets,
        'Rsen': Rsen,
        'Rcomm': Rcomm,
        'limits': limits,
        'targets': [genTargets(**kargs)],
        'potential': [genPotential(**kargs)],
    }

    return X

