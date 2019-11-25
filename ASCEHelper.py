# Notes about the helper
from __future__ import absolute_import

from six.moves import range
def ComputeCuTa(h_n, Ct=0.02, x=0.75):
    return Ct*h_n**x

def ComputeELFForces(FloorWeights, FloorHeights, Tn ):
    import numpy as np
    k = min(2.,max(1. ,1. + 1. / 2. * (Tn - 0.5)))
    Cvx = FloorWeights*np.array(FloorHeights)**k/np.sum(FloorWeights*np.array(FloorHeights)**k)
    return Cvx

def ComputeELFShearsMomentsLoads(FloorWeights, FloorHeights, CuTa, Cs):
    import numpy as np
    Cvx = ComputeELFForces(FloorWeights, np.array(FloorHeights)/12., CuTa)

    Loads = Cvx * np.sum(FloorWeights) * Cs
    StoryShear = []
    StoryMoment = []

    import numpy as np
    for i in range(len(Loads)):
        StoryShear.append(np.sum(Loads[i:]))
        if i == 0:
            StoryMoment.append(np.sum(Loads[i:]*(np.array(FloorHeights[i:]))))
        else:
            StoryMoment.append(np.sum(Loads[i:] * (np.array(FloorHeights[i:]) - FloorHeights[i-1])))

    class Output():
        def __init__(self):
            self.StoryShear = StoryShear
            self.StoryMoment = StoryMoment
            self.PushOverLoads = Loads

    return Output()

def ComputeMRSAShearsMomentsLoads(FloorWeights, FloorHeights, MRSALoads, Cs):
    import numpy as np

    Loads = MRSALoads
    StoryShear = []
    StoryMoment = []

    import numpy as np
    for i in range(len(Loads)):
        StoryShear.append(np.sum(Loads[i:]))
        if i == 0:
            StoryMoment.append(np.sum(Loads[i:]*(np.array(FloorHeights[i:]))))
        else:
            StoryMoment.append(np.sum(Loads[i:] * (np.array(FloorHeights[i:]) - FloorHeights[i-1])))

    class Output():
        def __init__(self):
            self.StoryShear = StoryShear
            self.StoryMoment = StoryMoment
            self.PushOverLoads = Loads

    return Output()

def GetDesignSa(t, S1, Sds, Sd1, TL, R, I, IgnoreMinBaseShear = False):
    Cs=min(Sds/R*I,Sd1/R*I/t,Sd1*TL/t**2/R*I)
    if not IgnoreMinBaseShear:
        Csmin = max(0.044*Sds*I,0.01)
    else:
        Csmin = 0.01

    if S1 > 0.6:
        Csmin = max(Csmin,0.5*S1/R*I)
    return max(Csmin,Cs)

def ComputeWallDesignAxialLoads(DeadLoad, LiveLoad, Sds, LLLessThan100psf = True):
    import numpy as np

    DeadLoad = np.array(DeadLoad)
    LiveLoad = np.array(LiveLoad)

    if LLLessThan100psf:
        Comb1 = (1.2+0.2*Sds)*DeadLoad + 0.5*LiveLoad
        Comb2 = (0.9-0.2*Sds)*DeadLoad
    else:
        Comb1 = (1.2+0.2*Sds)*DeadLoad + 1.0*LiveLoad
        Comb2 = (0.9-0.2*Sds)*DeadLoad

    def TakeSum(data):
        temp = []
        for i in range(len(data)):
            temp.append(np.sum(data[i:]))
        return np.array(temp)

    import numpy as np
    if np.sum(Comb1) < np.sum(Comb2): #take the smaller of the two, thats more conservative
        return TakeSum(Comb1)
    else:
        return TakeSum(Comb2)

def GetTBIDamping(Height):
    return min(0.05, max(0.36/(Height)**.5, 0.025))