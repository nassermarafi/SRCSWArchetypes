from __future__ import absolute_import
from __future__ import print_function
from six.moves import range
__author__ = 'marafi'
import os
# from numba import jit

def ExtractGroundMotion(Folder_Location):
    GMids = []
    GMFiles = {}
    GMData = {}
    Dt = {}
    NumPoints = {}
    # for subdir, dirs, files in os.walk(Folder_Location, True):

    files = [f for f in os.listdir(Folder_Location) if os.path.isfile(os.path.join(Folder_Location, f))]
    for file in files:

        def GetGMid(file):
            index = file.find('(')
            indexend = file.find(')')
            return file[index+1:indexend]

        def ExtractFirstLine(file):
            line = open(Folder_Location+'//'+file,'r').readline()
            return line

        def ExtractTable(file):
            lines = open(Folder_Location+'//'+file,'r').readlines()
            GMData = []
            for line in lines:
                for point in line.split():
                    GMData.append(float(point))
            import numpy as np
            return np.array(GMData)

        if file.startswith('SortedEQFile_'):
            id = GetGMid(file)
            GMids.append(id)
            GMFiles[id]=file
            GMData[id]=ExtractTable(file)

        if file.startswith('DtFile_('):
            id = GetGMid(file)
            Dt[id]=float(ExtractFirstLine(file))

        if file.startswith('NumPointsFile_('):
            id = GetGMid(file)
            NumPoints[id]=int(ExtractFirstLine(file))

    return GMids, GMFiles, Dt, NumPoints, GMData

def ExtractGroundMotionIds(Folder_Location):
    GMids = []

    files = [f for f in os.listdir(Folder_Location) if os.path.isfile(os.path.join(Folder_Location, f))]
    for file in files:

        def GetGMid(file):
            index = file.find('(')
            indexend = file.find(')')
            return file[index+1:indexend]

        if file.startswith('NumPointsFile_('):
            id = GetGMid(file)

            GMids.append(id)

    return GMids

def ExtractOneGroundMotion(Folder_Location, GMid):
    import numpy as np
    try:
        GMData = np.loadtxt(Folder_Location+'SortedEQFile_(%s).dat'%GMid)
    except:
        GMData = np.loadtxt(Folder_Location+'SortedEQFile_(%s).txt'%GMid)
    try:
        Dt = np.loadtxt(Folder_Location+'DtFile_(%s).dat'%GMid)
    except:
        Dt = np.loadtxt(Folder_Location+'DtFile_(%s).txt'%GMid)
    try:
        NumPoints =  np.loadtxt(Folder_Location+'NumPointsFile_(%s).dat'%GMid)
    except:
        NumPoints =  np.loadtxt(Folder_Location+'NumPointsFile_(%s).txt'%GMid)

    return Dt, NumPoints, GMData

def ExtractGroundMotion2(Folder_Location):
    GMids = []
    GMFiles = {}
    GMData = {}
    Dt = {}
    NumPoints = {}
    # for subdir, dirs, files in os.walk(Folder_Location, True):

    files = [f for f in os.listdir(Folder_Location) if os.path.isfile(os.path.join(Folder_Location, f))]
    for file in files:

        def GetGMid(file):
            index = file.find('(')
            indexend = file.find(')')
            return file[index+1:indexend]

        def ExtractFirstLine(file):
            line = open(Folder_Location+'//'+file,'r').readline()
            return line

        def ExtractTable(file):
            lines = open(Folder_Location+'//'+file,'r').readlines()
            GMData = []
            for line in lines:
                for point in line.split():
                    GMData.append(float(point))
            import numpy as np
            return np.array(GMData)

        if file.startswith('SortedEQFile_'):
            id = GetGMid(file)
            GMids.append(id)
            GMFiles[id]=file
            GMData[id]=ExtractTable(file)

        if file.startswith('DtFile_('):
            id = GetGMid(file)
            Dt[id]=float(ExtractFirstLine(file))

        if file.startswith('NumPointsFile_('):
            id = GetGMid(file)
            NumPoints[id]=int(ExtractFirstLine(file))
    try:
    #Get Mad and R from GM Parameter File
        Output = GetGMParameters2(Folder_Location)
    except:
        class output():
            pass
        Output = output()

        Output.M = None
        Output.R = None

    class Ouput:
        def __init__(self):
            self.GMids = GMids
            self.GMFiles = GMFiles
            self.Dt = Dt
            self.NumPoints = NumPoints
            self.GMData = GMData
            self.M = Output.M
            self.R = Output.R

    O = Ouput()

    return O

# @jit
# def FindSa(GMData, Dt, T, Zeta):
#     import numpy as np
#
#     #Using Newmarks: Linear System
#     #Using Linear Acceleration Method
#     #From Chopra, Page 177
#
#     #Mass
#     m = 1
#     wn = 2*np.pi/T
#
#     #Stiffness
#     k = wn**2.0*m
#
#     #Damping 2%
#     c=Zeta*2*m*wn
#
#     gamma = 0.5
#     beta = 1./4.
#
#     u = np.zeros(len(GMData))
#     v = np.zeros(len(GMData))
#     a = np.zeros(len(GMData))
#     p = np.array(GMData)
#
#     a[0]=p[0]/m
#
#     a1 = 1/beta/Dt**2*m+gamma/beta/Dt*c
#     a2 = 1/beta/Dt*m+(gamma/beta-1)*c
#     a3 = (1./2./beta-1.)*m+Dt*(gamma/2./beta-1.)*c
#
#     khat = k + a1
#
#     for i in range(1,len(GMData)):
#         phat = p[i] + a1*u[i-1] + a2*v[i-1] + a3*a[i-1]
#         u[i] = phat/khat
#         v[i] = gamma/beta/Dt*(u[i]-u[i-1]) + (1-gamma/beta)*v[i-1]+Dt*(1.-gamma/2./beta)*a[i-1]
#         a[i] = 1./beta/Dt**2.0*(u[i]-u[i-1])-1/beta/Dt*v[i-1]-(1/beta/2-1)*a[i-1]
#
#     a = a-p
#
#     return max(abs(max(u)),abs(min(u))),max(abs(max(v)),abs(min(v))),max(abs(max(a)),abs(min(a)))

import numba
import numpy as np
# nopython=True means an error will be raised
# if fast compilation is not possible.
@numba.jit(nopython=True)
def FindSa(GMData, Dt, T, Zeta):
    # Using Newmarks: Linear System
    # Using Linear Acceleration Method
    # From Chopra, Page 177

    #Mass
    m = 1
    wn = 2*np.pi/T

    #Stiffness
    k = wn**2.0*m

    #Damping 2%
    c=Zeta*2*m*wn

    gamma = 0.5
    beta = 1./4.

    u = np.zeros(len(GMData))
    v = np.zeros(len(GMData))
    a = np.zeros(len(GMData))
    p = GMData #np.array(GMData)

    a[0]=p[0]/m

    a1 = 1/beta/Dt**2*m+gamma/beta/Dt*c
    a2 = 1/beta/Dt*m+(gamma/beta-1)*c
    a3 = (1./2./beta-1.)*m+Dt*(gamma/2./beta-1.)*c

    khat = k + a1

    for i in range(1,len(GMData)):
        phat = p[i] + a1*u[i-1] + a2*v[i-1] + a3*a[i-1]
        u[i] = phat/khat
        v[i] = gamma/beta/Dt*(u[i]-u[i-1]) + (1-gamma/beta)*v[i-1]+Dt*(1.-gamma/2./beta)*a[i-1]
        a[i] = 1./beta/Dt**2.0*(u[i]-u[i-1])-1/beta/Dt*v[i-1]-(1/beta/2-1)*a[i-1]

    a = a-p
    a = np.abs(a)
    u = np.abs(u)
    v = np.abs(v)
    return np.max(u), np.max(v), np.max(a)

# nopython=True means an error will be raised
# if fast compilation is not possible.
@numba.jit
def FindSaQuick(GMData, Dt, T, Zeta):
    # Using Newmarks: Linear System
    # Using Linear Acceleration Method
    # From Chopra, Page 177

    #Mass
    m = 1
    wn = 2*np.pi/T

    #Stiffness
    k = wn**2.0*m

    #Damping 2%
    c=Zeta*2*m*wn

    gamma = 0.5
    beta = 1./4.

    u = np.zeros(len(GMData))
    v = np.zeros(len(GMData))
    a = np.zeros(len(GMData))
    p = GMData #np.array(GMData)

    a[0]=p[0]/m

    a1 = 1/beta/Dt**2*m+gamma/beta/Dt*c
    a2 = 1/beta/Dt*m+(gamma/beta-1)*c
    a3 = (1./2./beta-1.)*m+Dt*(gamma/2./beta-1.)*c

    khat = k + a1

    for i in range(1,len(GMData)):
        phat = p[i] + a1*u[i-1] + a2*v[i-1] + a3*a[i-1]
        u[i] = phat/khat
        v[i] = gamma/beta/Dt*(u[i]-u[i-1]) + (1-gamma/beta)*v[i-1]+Dt*(1.-gamma/2./beta)*a[i-1]
        a[i] = 1./beta/Dt**2.0*(u[i]-u[i-1])-1/beta/Dt*v[i-1]-(1/beta/2-1)*a[i-1]

    a = a-p
    a = np.abs(a)
    u = np.abs(u)
    v = np.abs(v)
    return np.max(a)

@numba.jit
def FindSaHistory(GMData, Dt, T, Zeta):
    # Using Newmarks: Linear System
    # Using Linear Acceleration Method
    # From Chopra, Page 177

    #Mass
    m = 1
    wn = 2*np.pi/T

    #Stiffness
    k = wn**2.0*m

    #Damping 2%
    c=Zeta*2*m*wn

    gamma = 0.5
    beta = 1./4.

    u = np.zeros(len(GMData))
    v = np.zeros(len(GMData))
    a = np.zeros(len(GMData))
    p = GMData #np.array(GMData)

    a[0]=p[0]/m

    a1 = 1/beta/Dt**2*m+gamma/beta/Dt*c
    a2 = 1/beta/Dt*m+(gamma/beta-1)*c
    a3 = (1./2./beta-1.)*m+Dt*(gamma/2./beta-1.)*c

    khat = k + a1

    for i in range(1,len(GMData)):
        phat = p[i] + a1*u[i-1] + a2*v[i-1] + a3*a[i-1]
        u[i] = phat/khat
        v[i] = gamma/beta/Dt*(u[i]-u[i-1]) + (1-gamma/beta)*v[i-1]+Dt*(1.-gamma/2./beta)*a[i-1]
        a[i] = 1./beta/Dt**2.0*(u[i]-u[i-1])-1/beta/Dt*v[i-1]-(1/beta/2-1)*a[i-1]

    a = a-p
    return np.array(a)

def ComputeResponseSpectrum(GMData, Dt, Periods, Zeta=0.05):
    Sa = []
    for T in Periods:
        u,v,a = FindSa(GMData, Dt, T, Zeta)
        Sa.append(a)

    return Sa

def FindSaOutput(GMData, Dt, T, Zeta):
    import numpy as np

    #Using Newmarks: Linear System
    #Using Linear Acceleration Method
    #From Chopra, Page 177

    #Mass
    m = 1
    wn = 2*np.pi/T

    #Stiffness
    k = wn**2.0*m

    #Damping 2%
    c=Zeta*2*m*wn

    gamma = 0.5
    beta = 1./4.

    u = np.zeros(len(GMData))
    v = np.zeros(len(GMData))
    a = np.zeros(len(GMData))
    p = np.array(GMData)

    a[0]=p[0]/m

    a1 = 1/beta/Dt**2*m+gamma/beta/Dt*c
    a2 = 1/beta/Dt*m+(gamma/beta-1)*c
    a3 = (1./2./beta-1.)*m+Dt*(gamma/2./beta-1.)*c

    khat = k + a1

    for i in range(1,len(GMData)):
        phat = p[i] + a1*u[i-1] + a2*v[i-1] + a3*a[i-1]
        u[i] = phat/khat
        v[i] = gamma/beta/Dt*(u[i]-u[i-1]) + (1-gamma/beta)*v[i-1]+Dt*(1.-gamma/2./beta)*a[i-1]
        a[i] = 1./beta/Dt**2.0*(u[i]-u[i-1])-1/beta/Dt*v[i-1]-(1/beta/2-1)*a[i-1]

    a = a-p

    class Output():
        def __init__(self):
            self.Sd = max(abs(u))
            self.Sv = max(abs(v))
            self.Sa = max(abs(a))
            self.TimeAtMaxSa = float(list(abs(a)).index(self.Sa))/len(a)

    return Output()

def RunElasticSDOF(GMData, Dt, T, Zeta):
    import numpy as np

    #Using Newmarks: Linear System
    #Using Linear Acceleration Method
    #From Chopra, Page 177

    #Mass
    m = 1
    wn = 2*np.pi/T

    #Stiffness
    k = wn**2.0*m

    #Damping 2%
    c=Zeta*2*m*wn

    gamma = 0.5
    beta = 1./4.

    u = np.zeros(len(GMData))
    v = np.zeros(len(GMData))
    a = np.zeros(len(GMData))
    p = np.array(GMData)

    a[0]=p[0]/m

    a1 = 1/beta/Dt**2*m+gamma/beta/Dt*c
    a2 = 1/beta/Dt*m+(gamma/beta-1)*c
    a3 = (1./2./beta-1.)*m+Dt*(gamma/2./beta-1.)*c

    khat = k + a1

    for i in range(1,len(GMData)):
        phat = p[i] + a1*u[i-1] + a2*v[i-1] + a3*a[i-1]
        u[i] = phat/khat
        v[i] = gamma/beta/Dt*(u[i]-u[i-1]) + (1-gamma/beta)*v[i-1]+Dt*(1.-gamma/2./beta)*a[i-1]
        a[i] = 1./beta/Dt**2.0*(u[i]-u[i-1])-1/beta/Dt*v[i-1]-(1/beta/2-1)*a[i-1]

    a = a-p/m

    class Output:
        def __init__(self):
            self.Displacement = u
            self.Velocity = v
            self.Accelecation = a
            self.GMData = GMData
            self.TimeStep = Dt
            self.Period = T
            self.Zeta = Zeta
            self.MaxU = max(abs(u))
            self.MaxV = max(abs(v))
            self.MaxA = max(abs(a))
    O = Output()

    return O

def RunElastoPlasticSDOF(GMData, Dt, T, Dy, Zeta=0.05):
    def getFs(k, delta_y, u, fOLD, uOLD):
        Fmax = k*delta_y
        fnew = fOLD + k*(u-uOLD)

        if abs(fnew) > Fmax:
            return Fmax*fnew/abs(fnew), 0.0
        else:
            return fnew, k

    import numpy as np

    #Mass
    g = 386.4
    m = 1
    wn = 2*np.pi/T

    #Stiffness
    k = wn**2.0*m

    #Damping 2%
    c=Zeta*2*m*wn

    #Using Average Acceleration Method
    gamma = 0.5
    beta = 1./4.

    #How many iterations in cycle
    maxj = 1000

    #Arrays
    u = np.zeros(len(GMData))
    v = np.zeros(len(GMData))
    a = np.zeros(len(GMData))
    p = np.array(GMData)
    phat = np.zeros(len(GMData))
    fs = np.zeros(len(GMData))
    kt = np.zeros(len(GMData))
    kthat = np.zeros(len(GMData))
    Rhat = np.zeros((len(GMData),maxj))

    #Initial Calculations
    a[0]=(p[0])/m
    fs[0]=0.0
    kt[0]=k

    a1 = 1/beta/Dt**2*m+gamma/beta/Dt*c
    a2 = 1/beta/Dt*m+(gamma/beta-1)*c
    a3 = (1./2./beta-1.)*m+Dt*(gamma/2./beta-1.)*c

    #Convergence
    tol = 1*10**-5

    for i in range(1,len(GMData)):
        u[i] = u[i-1]
        fs[i]=fs[i-1]
        kt[i]=kt[i-1]
        phat[i] = p[i] + a1*u[i-1]+ a2*v[i-1] + a3*a[i-1]
        for j in range(0,maxj):
            Rhat[i,j] = phat[i] - fs[i] - a1*u[i]
            if abs(Rhat[i,j]) < tol:
                break
            kthat[i] = kt[i] + a1
            deltau = Rhat[i,j]/kthat[i]
            u[i] = u[i] + deltau
            fs[i], kt[i] = getFs(k, Dy, u[i], fs[i-1], u[i-1])
        v[i]=gamma/beta/Dt*(u[i]-u[i-1])+(1.0-gamma/beta)*v[i-1]+Dt*(1-gamma/2.0/beta)*a[i-1]
        a[i]=1/beta/(Dt**2.0)*(u[i]-u[i-1])-1.0/(beta*Dt)*v[i-1]-(1.0/2.0/beta-1.0)*a[i-1]

    ########################## Plot Push Over ##########################

    t = np.arange(0,len(GMData),1.0)*Dt
    a = a-p

    class Output:
        def __init__(self):
            self.Displacement = u
            self.Velocity = v
            self.Accelecation = a
            self.Reaction = fs
            self.GMData = GMData
            self.TimeStep = Dt
            self.time = t
            self.Period = T
            self.Zeta = Zeta
            self.MaxU = max(abs(u))
            self.MaxV = max(abs(v))
            self.MaxA = max(abs(a))

    O = Output()

    return O

def FindSDI_NumericalIntegration(GMData, Dt, T, Dy):
    import numpy as np

    # def getFs(x,ke,dy):
    #     if x<dy:
    #         return ke*x
    #     else:
    #         return ke*dy+(x-dy)*0.05*ke
    #
    # def getKt(x,ke,dy):
    #     if x<dy:
    #         return ke
    #     else:
    #         return 0.05*ke

    def getFs(next_u,E,H,dy,alphaKIN,uP):
        Sy = E*dy
        TrialStress = E*(next_u-uP)
        TrialqKIN = -H*alphaKIN
        TrialEta = TrialStress+TrialqKIN
        TrialF = abs(TrialEta)-Sy

        if TrialF < 0:
            next_alphaKIN = alphaKIN
            next_Stress = TrialStress
            if TrialEta != 0:
                Kt = E*TrialEta/abs(TrialEta)
            else:
                Kt = E
            return next_Stress, uP, next_alphaKIN, Kt
        else:
            dGamma = TrialF/(E+H)
            next_uP = uP + dGamma*TrialEta/abs(TrialEta)
            next_alphaKIN = alphaKIN + dGamma*TrialEta/abs(TrialEta)
            next_Stress = TrialStress - dGamma*E*TrialEta/abs(TrialEta)
            Kt = H*TrialEta/abs(TrialEta)
            return next_Stress, next_uP, next_alphaKIN, Kt

    #Using Newmarks: Linear System
    #Using Linear Acceleration Method
    #From Chopra, Page 177

    #Mass
    m = 1
    wn = 2*np.pi/T

    #Stiffness
    k = wn**2.0*m

    #Damping 2%
    c=0.05*2*m*wn

    gamma = 0.5
    beta = 1./4.

    u = np.zeros(len(GMData))
    v = np.zeros(len(GMData))
    a = np.zeros(len(GMData))
    p = np.array(GMData)
    fs = np.zeros(len(GMData))
    alphaKIN = np.zeros(len(GMData))
    uP = np.zeros(len(GMData))

    a[0]=(p[0])/m

    a1 = 1/beta/Dt**2*m+gamma/beta/Dt*c
    a2 = 1/beta/Dt*m+(gamma/beta-1)*c
    a3 = (1./2./beta-1.)*m+Dt*(gamma/2./beta-1.)*c

    for i in range(1,len(GMData)):
        maxj = 1000
        tol = 1e-6
        j = 0

        utrial = np.zeros(maxj)
        fstrial = np.zeros(maxj)
        kttrial = np.zeros(maxj)
        kthattrial = np.zeros(maxj)
        alphaKINtrial = np.zeros(maxj)
        uPtrial = np.zeros(maxj)

        # for i in range(0,999):
        #     utrial[i]=.1*np.sin(2*np.pi*i*0.001)
        #     fstrial[i+1],uPtrial[i+1],alphaKINtrial[i+1],kttrial[i+1] = getFs(utrial[i], k, 0.05*k, Dy, alphaKINtrial[i], uPtrial[i])
        #
        # import matplotlib.pyplot as plt
        # plt.plot(utrial,fstrial,'.-')
        # plt.xlabel(r'$T_{n}, seconds$')
        # plt.ylabel(r'$S_{DI}, in$')
        # plt.title('Inelastic Response Spectrum, Dy=0.01')
        # ax = plt.gca()
        # ax.grid(True)
        # plt.show()


        utrial[j] = u[i]
        fstrial[j],uPtrial[j],alphaKINtrial[j],kttrial[j] = getFs(u[i], k, 0.05*k, Dy, alphaKIN[i], uP[i])
        phat = p[i] + a1*u[i-1] + a2*v[i-1] + a3*a[i-1]

        R = 1
        while abs(R) > tol:
            fstrial[j],uPtrial[j],alphaKINtrial[j],kttrial[j] = getFs(utrial[j], k, 0.05*k, Dy, alphaKIN[i], uP[i])
            R = phat - fstrial[j] - a1*utrial[j]
            kthattrial[j] = kttrial[j] + a1
            deltau = R / kthattrial[j]
            utrial[j+1] = utrial[j]+deltau
            fs[i] = fstrial[j]

            j=j+1

        alphaKIN[i] = alphaKINtrial[j-1]
        uP[i] = uPtrial[j-1]
        u[i] = utrial[j]
        v[i] = gamma/beta/Dt*(u[i]-u[i-1]) + (1-gamma/beta)*v[i-1]+Dt*(1.-gamma/2./beta)*a[i-1]
        a[i] = 1./beta/Dt**2.0*(u[i]-u[i-1])-1/beta/Dt*v[i-1]-(1/beta/2-1)*a[i-1]

    import matplotlib.pyplot as plt
    plt.plot(u,fs,'.-')
    plt.xlabel(r'$T_{n}, seconds$')
    plt.ylabel(r'$S_{DI}, in$')
    plt.title('Inelastic Response Spectrum, Dy=0.01')
    ax = plt.gca()
    ax.grid(True)
    plt.show()

    return max(abs(max(u)),abs(min(u))),max(abs(max(v)),abs(min(v))),max(abs(max(a)),abs(min(a)))

def ProcessDataFile(File_Location,File_Name):
    for subdir, dirs, files in os.walk(File_Location):

        def GetGMid(file):
                text = file.split('_')
                for i in range(0,len(text)):
                    if text[i-1] == 'no':
                        return text[i]

        for file in files:
            #Skip file names not matching
            if file != File_Name:
                continue

            lines = open(File_Location+'//'+file,'r').readlines()
            StartReading = False
            NGANo = GetGMid(file)
            Dt = None
            Num = None
            GMData = []
            for line in lines:
                if StartReading == False:
                    temp = line.split()
                    if len(temp)==4:
                        if temp[3]=='DT':
                            Dt = float(temp[1])
                            Num = int(temp[0])
                            StartReading =True
                else:
                    temp = line.split()
                    for num in temp:
                        GMData.append('%f'%(float(num)))

            Dt = 0.01
            Num = len(GMData)

            #Now Write Data Files
            outFile = open(File_Location+'//'+'DtFile_(%s).dat'%NGANo, 'w')
            outFile.write('%f'%Dt)
            outFile.close()
            outFile = open(File_Location+'//'+'NumPointsFile_(%s).dat'%NGANo, 'w')
            outFile.write('%d'%Num)
            outFile.close()
            outFile = open(File_Location+'//'+'SortedEQFile_(%s).dat'%NGANo, 'w')
            for point in GMData:
                outFile.write(point)
                outFile.write('\n')
            outFile.close()
    return None

def ProcessUSGSDataFile(File_Location,File_Name):
    for subdir, dirs, files in os.walk(File_Location):

        for file in files:
            #Skip file names not matching
            if file != File_Name:
                continue

            if not(File_Name.endswith('.acc')):
                continue

            lines = open(File_Location+'//'+file,'r').readlines()
            StartReading = False
            NGANo = File_Name[:-4]
            Dt = None
            Num = None

            import numpy as np
            data = np.loadtxt(File_Location+File_Name,skiprows=2)

            GMData = data[:,1]/np.mean(data[:,1])*0.784/100./9.81

            Dt = data[1,0]-data[0,0]
            Num = len(GMData)

            #Now Write Data Files
            outFile = open(File_Location+'//'+'DtFile_(%s).dat'%NGANo, 'w')
            outFile.write('%f'%Dt)
            outFile.close()
            outFile = open(File_Location+'//'+'NumPointsFile_(%s).dat'%NGANo, 'w')
            outFile.write('%d'%Num)
            outFile.close()
            outFile = open(File_Location+'//'+'SortedEQFile_(%s).dat'%NGANo, 'w')
            for point in GMData.tolist():
                outFile.write('%f'%point)
                outFile.write('\n')
            outFile.close()
    return None

def ProcessNGAWest2DataFiles(File_Location):

    def GetGMid(file):
        text = file.split('_')
        for i in range(0,len(text)):
            if text[i-1] == 'no':
                return text[i]
    import pandas as pd
    NGAWest2FF = pd.read_csv(os.getcwd() + '/GroundMotions/NGA-West-2-FlatFiles/ngawest2flatfile_filenamesonly.csv')

    for subdir, dirs, files in os.walk(File_Location):
        for file in files:
            #Make Sure Acceleration File
            if not(file.endswith('.AT2')):
                continue
            #Find if H or V
            if file.endswith('Z.AT2'):
                continue #Ignoring Vertical GMs
            if file.endswith('UP.AT2'):
                continue #Ignoring Vertical GMs
            #Find Out Directions
            if file in set(NGAWest2FF['File Name (Horizontal 1)']):
                Direction = '1'
            elif file in set(NGAWest2FF['File Name (Horizontal 2)']):
                Direction = '2'
            else:
                RSN = int(file[3:file.find('_')])
                record = NGAWest2FF[NGAWest2FF['Record Sequence Number'] == RSN]
                if file[-7:] == record['File Name (Horizontal 1)'].tolist()[0][-7:]:
                    Direction = '1'
                elif file[-7:] == record['File Name (Horizontal 2)'].tolist()[0][-7:]:
                    Direction = '2'
                else:
                    continue
            #Make Find Out GMId
            ind = file.index('_')
            GMid = file[3:ind]+Direction

            # Check if File Already Created
            if os.path.isfile(File_Location+'//'+'SortedEQFile_(%s).dat'%GMid):
                continue

            lines = open(File_Location+'//'+file,'r').readlines()
            StartReading = False
            Dt = None
            Num = None
            GMData = []
            for line in lines:
                if StartReading == False:
                    temp = line.split()
                    for i in range(len(temp)):
                        if temp[i]=='DT=' or temp[i]=='dt=':
                            try:
                                Dt = float(temp[i+1])
                            except:
                                Dt = float(temp[i+1][:-4])
                            StartReading =True
                else:
                    temp = line.split()
                    for num in temp:
                        GMData.append('%f'%(float(num)))

            Num = len(GMData)

            #Now Write Data Files
            try:
                outFile = open(File_Location+'//'+'DtFile_(%s).dat'%GMid, 'w')
                outFile.write('%f'%Dt)
                outFile.close()
                outFile = open(File_Location+'//'+'NumPointsFile_(%s).dat'%GMid, 'w')
                outFile.write('%d'%Num)
                outFile.close()
                outFile = open(File_Location+'//'+'SortedEQFile_(%s).dat'%GMid, 'w')
                for point in GMData:
                    outFile.write(point)
                    outFile.write('\n')
                outFile.close()
            except:
                print(('Error in GMid: ', GMid))
    return None

def ProcessNGASubductionDataFiles(File_Location):
    def GetGMid(file):
        text = file.split('_')
        for i in range(0,len(text)):
            if text[i-1] == 'no':
                return text[i]

    import pandas as pd
    NGAWest2FF = pd.read_csv(File_Location + 'NGASubductionIntraslabMetadata.csv')

    for subdir, dirs, files in os.walk(File_Location):
        for file in files:
            print(file)

            #Make Sure Acceleration File
            if not(file.endswith('.AT2')):
                continue
            #Find if H or V
            if file.endswith('Z.AT2'):
                continue #Ignoring Vertical GMs
            if file.endswith('UP.AT2'):
                continue #Ignoring Vertical GMs
            #Find Out Directions
            if file in set(NGAWest2FF['Filename_H1']):
                Direction = '1'
            elif file in set(NGAWest2FF['Filename_H2']):
                Direction = '2'
            else:
                continue
            #Make Find Out GMId
            ind = file.index('_')
            GMid = file[3:ind]+Direction

            # Check if File Already Created
            if os.path.isfile(File_Location+'//'+'SortedEQFile_(%s).dat'%GMid):
                continue

            lines = open(File_Location+'//'+file,'r').readlines()
            StartReading = False
            Dt = None
            Num = None
            GMData = []
            for line in lines:
                if StartReading == False:
                    temp = line.split()
                    for i in range(len(temp)):
                        if temp[i]=='DT=' or temp[i]=='dt=':
                            try:
                                Dt = float(temp[i+1])
                            except:
                                Dt = float(temp[i+1][:-4])
                            StartReading =True
                else:
                    temp = line.split()
                    for num in temp:
                        GMData.append('%f'%(float(num)))

            Num = len(GMData)

            #Now Write Data Files
            try:
                outFile = open(File_Location+'//'+'DtFile_(%s).dat'%GMid, 'w')
                outFile.write('%f'%Dt)
                outFile.close()
                outFile = open(File_Location+'//'+'NumPointsFile_(%s).dat'%GMid, 'w')
                outFile.write('%d'%Num)
                outFile.close()
                outFile = open(File_Location+'//'+'SortedEQFile_(%s).dat'%GMid, 'w')
                for point in GMData:
                    outFile.write(point)
                    outFile.write('\n')
                outFile.close()
            except:
                print(('Error in GMid: ', GMid))
    return None

def ProcessNGAEastDataFiles(File_Location):

    def GetGMid(file):
        text = file.split('_')
        for i in range(0,len(text)):
            if text[i-1] == 'no':
                return text[i]
    import pandas as pd
    NGAWest2FF = pd.read_csv(File_Location + 'FileNames.csv')

    for subdir, dirs, files in os.walk(File_Location):
        for file in files:
            #Make Sure Acceleration File
            if not(file.endswith('.AT2')):
                continue
            #Find if H or V
            if file.endswith('Z.AT2'):
                continue #Ignoring Vertical GMs
            if file.endswith('UP.AT2'):
                continue #Ignoring Vertical GMs
            #Find Out Directions
            if file in set(NGAWest2FF['Horizontal-1Acc.Filename']):
                Direction = '1'
            elif file in set(NGAWest2FF['Horizontal-2Acc.Filename']):
                Direction = '2'
            else:
                RSN = int(file[3:file.find('_')])
                record = NGAWest2FF[NGAWest2FF['RecordSequenceNumber'] == RSN]
                if file[-7:] == record['Horizontal-1Acc.Filename'].tolist()[0][-7:]:
                    Direction = '1'
                elif file[-7:] == record['Horizontal-2Acc.Filename'].tolist()[0][-7:]:
                    Direction = '2'
                else:
                    continue
            #Make Find Out GMId
            ind = file.index('_')
            GMid = file[3:ind]+Direction

            # Check if File Already Created
            if os.path.isfile(File_Location+'//'+'SortedEQFile_(%s).dat'%GMid):
                continue

            lines = open(File_Location+'//'+file,'r').readlines()
            StartReading = False
            Dt = None
            Num = None
            GMData = []
            for line in lines:
                if StartReading == False:
                    temp = line.split()
                    for i in range(len(temp)):
                        if temp[i]=='DT=' or temp[i]=='dt=':
                            try:
                                Dt = float(temp[i+1])
                            except:
                                Dt = float(temp[i+1][:-4])
                            StartReading =True
                else:
                    temp = line.split()
                    for num in temp:
                        GMData.append('%f'%(float(num)))

            Num = len(GMData)

            #Now Write Data Files
            try:
                outFile = open(File_Location+'//'+'DtFile_(%s).dat'%GMid, 'w')
                outFile.write('%f'%Dt)
                outFile.close()
                outFile = open(File_Location+'//'+'NumPointsFile_(%s).dat'%GMid, 'w')
                outFile.write('%d'%Num)
                outFile.close()
                outFile = open(File_Location+'//'+'SortedEQFile_(%s).dat'%GMid, 'w')
                for point in GMData:
                    outFile.write(point)
                    outFile.write('\n')
                outFile.close()
            except:
                print(('Error in GMid: ', GMid))
    return None

def ProcessDataFile_KNET(File_Location):
    import pandas as pd
    GMDatabase = pd.DataFrame()

    for subdir, dirs, files in os.walk(File_Location):
        for file in files:
            if not (file.endswith('EW1') or file.endswith('EW') or file.endswith('EW2')
                    or file.endswith('NS1') or file.endswith('NS') or file.endswith('NS2')):
                continue
                # EndNo = '3'

            print(('Processing File Name: %s'%file))

            lines = open(File_Location+'/'+file, 'r').readlines()

            StartReading = False
            file = file.replace('.', '')

            NGANo = file

            OriginalTime = lines[0][18:].replace('\r','')
            Latitude = float(lines[1][18:])
            Longitude = float(lines[2][18:])
            Depth = float(lines[3][18:])
            Magnitude = float(lines[4][18:])
            StationCode = lines[5][18:].replace('\r','')
            StationLat = float(lines[6][18:])
            StationLong = float(lines[7][18:])
            StationHeight = float(lines[8][18:])
            RecordTime = lines[9][18:].replace('\r','')
            SamplingFrequency = float(lines[10][18:-3])
            DurationTime = float(lines[11][18:])
            Dir = lines[12][18:].replace('\r','')
            LastCorrection = lines[15][18:].replace('\r','')

            #Compute R_JB Distance
            from geopy.distance import vincenty
            R_JB = vincenty((Latitude, Longitude),(StationLat, StationLong)).kilometers

            Dt = None
            Num = None
            ScaleFactor = None

            GMData = []
            for line in lines:
                if StartReading == False:
                    temp = line.split()
                    if temp[0] == 'Sampling':
                        Dt = 1.0/float(temp[2].replace('Hz',''))
                    if temp[0] == 'Scale':
                        text = temp[2].split('(gal)/')
                        ScaleFactor = float(text[0])/float(text[1])
                    if temp[0]=='Memo.':
                            StartReading =True
                else:
                    temp = line.split()
                    for num in temp:
                        GMData.append(float(num)/386.4/2.54*ScaleFactor)

            Num = len(GMData)

            import numpy as np
            GMData = np.array(GMData) - float(GMData[0])

            #Now Write Data Files
            outFile = open(File_Location+'/'+'DtFile_(%s).dat'%NGANo, 'w')
            outFile.write('%f'%Dt)
            outFile.close()
            outFile = open(File_Location+'/'+'NumPointsFile_(%s).dat'%NGANo, 'w')
            outFile.write('%d'%Num)
            outFile.close()
            outFile = open(File_Location+'/'+'SortedEQFile_(%s).dat'%NGANo, 'w')
            for point in GMData:
                outFile.write('%f'%point)
                outFile.write('\n')
            outFile.close()

            GMDatabase = GMDatabase.append({'GMid': NGANo, 'Station Code': StationCode, 'OriginalTime' : OriginalTime, 'Latitude' : Latitude, 'Longitude' : Longitude, 'Depth' : Depth, 'Magnitude' : Magnitude, 'StationLat' : StationLat, 'StationLong' : StationLong, 'StationHeight': StationHeight,  'RecordTime' : RecordTime, 'SamplingFrequency' : SamplingFrequency, 'Dt': Dt, 'DurationTime' : DurationTime, 'Dir':Dir, 'LastCorrection':LastCorrection, 'NumPoints': Num, 'R_JB':R_JB}, ignore_index=True)

    GMDatabase.to_csv(File_Location+'/GMDatabase.csv')

    return None

def ProcessDataFile_KNET_GMDataBaseFileOnly(File_Location):
    import pandas as pd
    GMDatabase = pd.DataFrame()

    for subdir, dirs, files in os.walk(File_Location):
        for file in files:
            if not (file.endswith('EW1') or file.endswith('EW') or file.endswith('EW2')
                    or file.endswith('NS1') or file.endswith('NS') or file.endswith('NS2')):
                continue
                # EndNo = '3'

            print(('Processing File Name: %s'%file))

            lines = open(File_Location+'/'+file, 'r').readlines()

            StartReading = False
            file = file.replace('.', '')

            NGANo = file

            OriginalTime = lines[0][18:].replace('\r','').replace('\n','')
            Latitude = float(lines[1][18:])
            Longitude = float(lines[2][18:])
            Depth = float(lines[3][18:])
            Magnitude = float(lines[4][18:])
            StationCode = lines[5][18:].replace('\r','').replace('\n','')
            StationLat = float(lines[6][18:])
            StationLong = float(lines[7][18:])
            StationHeight = float(lines[8][18:])
            RecordTime = lines[9][18:].replace('\r','').replace('\n','')
            SamplingFrequency = float(lines[10][18:-3])
            DurationTime = float(lines[11][18:])
            Dir = lines[12][18:].replace('\r','')
            LastCorrection = lines[15][18:].replace('\r','').replace('\n','')

            #Compute R_JB Distance
            from geopy.distance import vincenty
            R_JB = vincenty((Latitude, Longitude),(StationLat, StationLong)).kilometers

            Dt = None
            Num = None
            ScaleFactor = None

            GMData = []
            for line in lines:
                if StartReading == False:
                    temp = line.split()
                    if temp[0] == 'Sampling':
                        Dt = 1.0/float(temp[2].replace('Hz',''))
                    if temp[0] == 'Scale':
                        text = temp[2].split('(gal)/')
                        ScaleFactor = float(text[0])/float(text[1])
                    if temp[0]=='Memo.':
                            StartReading =True
                else:
                    temp = line.split()
                    for num in temp:
                        GMData.append(float(num)/386.4/2.54*ScaleFactor)

            Num = len(GMData)

            GMDatabase = GMDatabase.append({'GMid': NGANo, 'Station Code': StationCode, 'OriginalTime' : OriginalTime, 'Latitude' : Latitude, 'Longitude' : Longitude, 'Depth' : Depth, 'Magnitude' : Magnitude, 'StationLat' : StationLat, 'StationLong' : StationLong, 'StationHeight': StationHeight,  'RecordTime' : RecordTime, 'SamplingFrequency' : SamplingFrequency, 'Dt': Dt, 'DurationTime' : DurationTime, 'Dir':Dir, 'LastCorrection':LastCorrection, 'NumPoints': Num, 'R_JB':R_JB}, ignore_index=True)

    GMDatabase.to_csv(File_Location+'/GMDatabase.csv')

    return None

def CreateResponseSpectrumFile_OLD(File_Location, NGANo, GMData, Dt, Zeta, Ry, PostYieldPer, Comments=''):
    import numpy as np

    g = 386.4

    accuracy = 0.1
    Tmax = 20

    T = np.zeros(len(list(range(0,int(Tmax/accuracy)))))

    uElastic = np.zeros(len(list(range(0,int(Tmax/accuracy)))))
    vElastic = np.zeros(len(list(range(0,int(Tmax/accuracy)))))
    aElastic = np.zeros(len(list(range(0,int(Tmax/accuracy)))))

    u = np.zeros(len(list(range(0,int(Tmax/accuracy)))))
    v = np.zeros(len(list(range(0,int(Tmax/accuracy)))))
    a = np.zeros(len(list(range(0,int(Tmax/accuracy)))))

    for i in range(0,int(Tmax/accuracy)):
        T[i] = accuracy+i*accuracy
        uElastic[i], vElastic[i], aElastic[i] = FindSa(GMData*g, Dt, T[i], Zeta)
        # uElastic[i], vElastic[i], aElastic[i] = FindSDI_OpenSees(GMData*g, Dt, T[i], 1000, Zeta, 1.0)
        u[i], v[i], a[i] = FindSDI_OpenSees(GMData*g, Dt, T[i], uElastic[i]/Ry, Zeta, PostYieldPer)

    #Find Ds
    Ds = FindDs(GMData, Dt)

    outFile = open(File_Location+'//'+'ResponseSpectrum_(%s).dat'%NGANo, 'w')
    outFile.write('NGANo: %s \n'%NGANo)
    outFile.write('Zeta: %f \n'%Zeta)
    outFile.write('Ry: %f \n'%Ry)
    outFile.write('Post Yield Percentage of Elastic Stiffness: %f \n'%PostYieldPer)
    outFile.write('Ds5-95: %f \n'%Ds)
    outFile.write('Comments: %s \n'%Comments)
    outFile.write('Period, SD, SV, SA, SDI, SVI, SAI \n')
    for i in range(0,int(Tmax/accuracy)):
        if i != int(Tmax/accuracy)-1:
            outFile.write('%f %f %f %f %f %f %f \n'%(T[i],uElastic[i], vElastic[i], aElastic[i],u[i],v[i],a[i]))
        else:
            outFile.write('%f %f %f %f %f %f %f'%(T[i],uElastic[i], vElastic[i], aElastic[i],u[i],v[i],a[i]))
    outFile.close()

def CreateResponseSpectrumFile(File_Location, NGANo, GMData, Dt, Zeta, Comments='', Tmax = 40):
    import numpy as np

    g = 386.4

    accuracy = 0.01
    # Tmax = 40

    T = np.zeros(len(list(range(0,int(Tmax/accuracy)))))

    uElastic = np.zeros(len(list(range(0,int(Tmax/accuracy)))))
    vElastic = np.zeros(len(list(range(0,int(Tmax/accuracy)))))
    aElastic = np.zeros(len(list(range(0,int(Tmax/accuracy)))))

    u = np.zeros(len(list(range(0,int(Tmax/accuracy)))))
    v = np.zeros(len(list(range(0,int(Tmax/accuracy)))))
    a = np.zeros(len(list(range(0,int(Tmax/accuracy)))))

    for i in range(0,int(Tmax/accuracy)):
        T[i] = accuracy+i*accuracy
        uElastic[i], vElastic[i], aElastic[i] = FindSa(GMData*g, Dt, T[i], Zeta)

    #Find Ds
    Ds = FindDs(GMData, Dt)

    outFile = open(File_Location+'//'+'ResponseSpectrum_(%s).dat'%NGANo, 'w')
    outFile.write('NGANo: %s \n'%NGANo)
    outFile.write('Zeta: %f \n'%Zeta)
    try:
        outFile.write('Ds5-95: %f \n'%Ds)
    except:
        outFile.write('Ds5-95: FAILED \n')
    outFile.write('Comments: %s \n'%Comments)
    outFile.write('Period, SD, SV, SA\n')
    for i in range(0,int(Tmax/accuracy)):
        if i != int(Tmax/accuracy)-1:
            outFile.write('%f %f %f %f  \n'%(T[i],uElastic[i], vElastic[i], aElastic[i]))
        else:
            outFile.write('%f %f %f %f '%(T[i],uElastic[i], vElastic[i], aElastic[i]))
    outFile.close()

def ReadResponseSpectrumFile(File_Location, File_Name):
    f = open(File_Location+'//'+File_Name)
    lines = f.readlines()
    ResponseSpectrum = []
    #Read Zeta
    Zeta = float(lines[1].split(': ')[1])
    # Read Ds
    Ds = float(lines[2].split(': ')[1])
    #Read Data
    for i in range(7,len(lines)):
        temp = lines[i].split()
        fixed = []
        for j in range(0,len(temp)):
            fixed.append(float(temp[j]))
        ResponseSpectrum.append(fixed)
    import numpy as np
    ResponseSpectrum = np.array(ResponseSpectrum)
    f.close()
    return ResponseSpectrum, Zeta, Ds

def ProcessSynFiles(File_Location, File_Name, NGA_Number, NGA_Name):
    f = open(File_Location+'//'+File_Name)
    lines = f.readlines()
    ResponseSpectrum = []
    for i in range(0,len(lines)):
        ResponseSpectrum.append(lines[i].split())
    import numpy as np
    ResponseSpectrum = np.array(ResponseSpectrum)
    CheckDt = float(ResponseSpectrum[1][0])-float(ResponseSpectrum[0][0])
    # Find Dt and Create File
    for i in range(1,len(ResponseSpectrum)):
        Dt = float(ResponseSpectrum[i][0])-float(ResponseSpectrum[i-1][0])
        if round(CheckDt,3) != round(Dt,3):
            print('Error with DT')
            return None
    for j in range(1,4):
        #Create Dt File
            f = open(File_Location+'//'+'DtFile_(%s%d).dat'%(NGA_Number,j),'w')
            f.write('%f'%Dt)
            f.close()
        # Create GMData File
            f = open(File_Location+'//'+'SortedEQFile_(%s%d).dat'%(NGA_Number,j),'w')
            for i in range(0,len(ResponseSpectrum)):
                f.write('%f\n'%(float(ResponseSpectrum[i][j])/386.4/2.54))
            f.close()
        # Create Num Point File
            f = open(File_Location+'//'+'NumPointsFile_(%s%d).dat'%(NGA_Number,j),'w')
            f.write('%d'%len(ResponseSpectrum))
            f.close()

def ProcessGMFileFromGMData(File_Location, GMData, Dt, GMid):
    # Find Dt and Create File
    #Create Dt File
    f = open(File_Location+'//'+'DtFile_(%s).dat'%(GMid),'w')
    f.write('%f'%Dt)
    f.close()
    # Create GMData File
    f = open(File_Location+'//'+'SortedEQFile_(%s).dat'%(GMid),'w')
    for i in range(0,len(GMData)):
        f.write('%f\n'%(float(GMData[i])))
    f.close()
    # Create Num Point File
    f = open(File_Location+'//'+'NumPointsFile_(%s).dat'%(GMid),'w')
    f.write('%d'%len(GMData))
    f.close()

def PlotResponseSpecrum(File_Location, File_Name, axarr, htmlcolor, Scale_Factor=1.0, Label=None):
    g = 386.4
    Data, Zeta, Ds = ReadResponseSpectrumFile(File_Location, File_Name)

    Plot = axarr[0].plot(Data[:,0], Data[:,1]*Scale_Factor,'-',linewidth=2.0,color=htmlcolor,label=Label)
    Plot = axarr[1].plot(Data[:,0], Data[:,2]*Scale_Factor,'-',linewidth=2.0,color=htmlcolor,label=Label)
    Plot = axarr[2].plot(Data[:,0], Data[:,3]*Scale_Factor/g,'-',linewidth=2.0,color=htmlcolor,label=Label)

    axarr[0].grid(True)
    axarr[1].grid(True)
    axarr[2].grid(True)

    axarr[0].set_ylabel(r'$S_{D},\; in$')
    axarr[1].set_ylabel(r'$S_{V},\; {in}/{sec}$')
    axarr[2].set_ylabel(r'$S_{A},\; g$')
    axarr[2].set_xlabel(r'$T_{n},\; seconds$')

    axarr[0].set_xlim(0,10)
    axarr[1].set_xlim(0,10)
    axarr[2].set_xlim(0,10)

    return Plot

def PlotPseudoResponseSpecrum(File_Location, File_Name, axarr, htmlcolor, Scale_Factor=1.0, PlotInelastic=False, Label=None):
    import numpy as np
    g = 386.4
    Data, Zeta, Ry, PostYieldStiffness, Ds = ReadResponseSpectrumFile(File_Location, File_Name)

    if PlotInelastic == False:
        Wn = 2*np.pi/Data[:,0]
        Plot = axarr[0].plot(Data[:,0], Data[:,1]*Scale_Factor,'-',color=htmlcolor,label=Label, linewidth=2.0)
        Plot = axarr[1].plot(Data[:,0], Data[:,1]*Scale_Factor*Wn,'-',color=htmlcolor,label=Label, linewidth=2.0)
        Plot = axarr[2].plot(Data[:,0], Data[:,1]*Scale_Factor/g*Wn*Wn,'-',color=htmlcolor,label=Label, linewidth=2.0)
    else:
        Wn = 2*np.pi/Data[:,0]
        # Plot = axarr[0].plot(Data[:,0], Data[:,4]*Scale_Factor,'.-',color=htmlcolor,label=Label)
        # Plot = axarr[1].plot(Data[:,0], Data[:,4]*Scale_Factor*Wn,'.-',color=htmlcolor,label=Label)
        # Plot = axarr[2].plot(Data[:,0], Data[:,4]*Scale_Factor/g*Wn*Wn,'.-',color=htmlcolor,label=Label)

        PsuedoAcc = np.zeros(len(Data[:,0]))
        PsuedoVel = np.zeros(len(Data[:,0]))
        for i in range(0,len(Data[:,0])):
            wn = 2*np.pi/Data[i,0]
            k = wn**2
            if Data[i,4] > Data[i,1]/Ry:
                PsuedoAcc[i] = ((Data[i,4]-Data[i,1]/Ry)*k*PostYieldStiffness + (Data[i,1]/Ry)*k)/g
                PsuedoVel[i] = PsuedoAcc[i]/wn*g
            else:
                PsuedoAcc[i] = (Data[i,4])*k/g
                PsuedoVel[i] = PsuedoAcc[i]/wn*g

        Plot = axarr[0].plot(Data[:,0], Data[:,4]*Scale_Factor,'.-',color=htmlcolor,label=Label)
        Plot = axarr[1].plot(Data[:,0], PsuedoVel,'.-',color=htmlcolor,label=Label)
        Plot = axarr[2].plot(Data[:,0], PsuedoAcc,'.-',color=htmlcolor,label=Label)

    axarr[0].grid(True)
    axarr[1].grid(True)
    axarr[2].grid(True)

    axarr[0].set_ylabel(r'$S_\mathrm{d}$, in', fontsize=20)
    axarr[1].set_ylabel(r'$S_\mathrm{v}$, in/sec', fontsize=20)
    axarr[2].set_ylabel(r'$S_\mathrm{a}$, g', fontsize=20)
    axarr[2].set_xlabel(r'$T_\mathrm{n}$, seconds', fontsize=20)

    return Plot

def FindDs(GMData, Dt, min=.05, max=0.95):
    import numpy as np
    from scipy.integrate import trapz
    T = np.array(list(range(0,len(GMData))))*Dt
    Total = trapz(np.array(GMData)**2,dx=Dt)
    Ds = []
    mini = None
    maxi = None
    for i in range(0,len(GMData)):
        Ds.append(trapz(GMData[0:i]**2,dx=Dt)/Total)
        if min < Ds[i] and min > Ds[i-1]:
            mini = i
        if max < Ds[i] and max > Ds[i-1]:
            maxi = i
    DsMinMax = T[maxi]-T[mini]
    Ds = np.array(Ds)

    class Output:
        def __init__(self):
            self.Ds = DsMinMax
            self.MaxIndex = maxi
            self.MinIndex = mini

    O = Output()

    return DsMinMax

def FindDsFromGMid(GMFolder, GMid, min=.05, max=0.95):
    Dt, NumPoints, GMData = ExtractOneGroundMotion(GMFolder, GMid)

    import numpy as np
    from scipy.integrate import trapz
    T = np.array(list(range(0,len(GMData))))*Dt
    Total = trapz(np.array(GMData)**2,dx=Dt)
    Ds = []
    mini = None
    maxi = None
    for i in range(0,len(GMData)):
        Ds.append(trapz(GMData[0:i]**2,dx=Dt)/Total)
        if min < Ds[i] and min > Ds[i-1]:
            mini = i
        if max < Ds[i] and max > Ds[i-1]:
            maxi = i
    DsMinMax = T[maxi]-T[mini]
    Ds = np.array(Ds)

    class Output:
        def __init__(self):
            self.Ds = DsMinMax
            self.MaxIndex = maxi
            self.MinIndex = mini

    O = Output()

    return DsMinMax

@numba.jit(nopython=True)
def FindDsQuick(GMData, Dt, min=.05, max=0.95):
    GMData = np.array(GMData)**2.
    GMData = GMData/GMData.sum()
    GMData = list(np.cumsum(GMData))
    DsMinMax = (len([x for x in GMData if x>min and x < max])-1)*Dt
    return DsMinMax

def FindBracketedDuration(GMData, Dt, Threshold=0.05):
    import numpy as np
    GMData = np.abs(np.array(GMData))
    Index = np.argwhere(GMData>Threshold)
    return (max(Index)-min(Index))*Dt

def FindCAV(GMData, Dt, Threshold=0.05):
    import numpy as np
    from scipy.integrate import trapz
    GMData = np.abs(np.array(GMData))
    GMData[GMData < Threshold] = 0
    return np.trapz(GMData, dx = Dt)

def FindDsQuickFromGMid(GMFolder, GMid, min=.05, max=0.95):
    Dt, NumPoints, GMData = ExtractOneGroundMotion(GMFolder, GMid)

    import numpy as np
    from scipy.integrate import trapz
    GMData = np.array(GMData)**2.
    GMData = GMData/GMData.sum()
    GMData = np.cumsum(GMData).tolist()
    DsMinMax = len([x for x in GMData if x>min and x < max])*Dt
    return DsMinMax

def FindDs595FromRS(GMFolder, GMid,):
    f = open(GMFolder+'/'+'ResponseSpectrum_(%s).dat'%GMid, 'r')
    for i in range(10):
        line = f.readline()
        if line.startswith('Ds5-95'):
            return float(line.split()[1])

def GetDs(GMData, Dt, min=.05, max=0.95):
    import numpy as np
    from scipy.integrate import trapz
    T = np.array(list(range(0,len(GMData))))*Dt
    Total = trapz(np.array(GMData)**2,dx=Dt)
    Ds = []
    mini = None
    maxi = None
    for i in range(0,len(GMData)):
        Ds.append(trapz(GMData[0:i]**2,dx=Dt)/Total)
        if min < Ds[i] and min > Ds[i-1]:
            mini = i
        if max < Ds[i] and max > Ds[i-1]:
            maxi = i
    DsMinMax = T[maxi]-T[mini]
    Ds = np.array(Ds)

    class Output:
        def __init__(self):
            self.Ds = DsMinMax
            self.MaxIndex = maxi
            self.MinIndex = mini

    O = Output()

    return O

def FindScaleFactorToFit(Response_Spectrum, HazardFit, Tmin=None, Tmax=None):
    T = Response_Spectrum[:,0]
    Sa = Response_Spectrum[:,3]
    HT = HazardFit[:,0]
    HSa = HazardFit[:,1]

    import numpy as np
    HSa = np.interp(T,HT,HSa)
    HT = T

    def AdjustMin(Tadjust,Saadjust):
        if Tmin != None:
            for i in range(0,len(Tadjust)):
                if Tmin > Tadjust[i]:
                    continue
                else:
                    Tadjust = Tadjust[i:]
                    Saadjust = Saadjust[i:]
                    return Tadjust, Saadjust

    def AdjustMax(Tadjust,Saadjust):
        if Tmax != None:
            for i in range(0,len(Tadjust)):
                if Tmax > Tadjust[i]:
                    continue
                else:
                    Tadjust = Tadjust[:i]
                    Saadjust = Saadjust[:i]
                    return Tadjust, Saadjust

    T, Sa = AdjustMin(T,Sa)
    HT, HSa = AdjustMin(HT,HSa)
    T, Sa = AdjustMax(T,Sa)
    HT, HSa = AdjustMax(HT,HSa)

    def func(x,a):
        Results=np.zeros(len(x))
        for j in range(len(x)):
            for i in range(1,len(T)):
                if T[i] >= x[j] and T[i-1] < x[j]:
                     Results[j]=a*Sa[i]
        return Results

    g = 386.4

    xData = HT
    yData = HSa*g

    x0 = np.array([1])
    import scipy.optimize as optimization
    print('Optimizing Curve')
    SF = optimization.curve_fit(func, xData, yData)
    return SF

def HazardCurve_OLD(Lat, Long):
    def round_to(n, precision):
        correction = 0.5 if n >= 0 else -0.5
        return int(n/precision+correction ) * precision
    Lat = round_to(Lat, 0.05)
    Long = round_to(Long, 0.05)
    import os
    import numpy as np
    Folder_Location = os.getcwd()+'//GroundMotions//UniformHazard'
    T = []
    Sa = []
    for subdir, dirs, files in os.walk(Folder_Location):
        for file in files:
            print(('Opening File %s'%file))
            lines = open(Folder_Location+'//'+file,'r').readlines()
            #Get Period
            T.append(float(lines[1].split()[6]))

            for i in range(2,len(lines)):
                templine = lines[i].split()
                LatRow = round_to(float(templine[1]),0.05)
                LongRow = round_to(float(templine[0]),0.05)

                if LatRow == Lat and LongRow == Long:
                    Sa.append(float(templine[2]))

    Hazard = np.transpose(np.array([T, Sa]))
    Hazard = Hazard[np.argsort(Hazard[:,0])]
    Hazard = np.transpose(Hazard)
    return Hazard

def HazardCurve(Lat, Long):
    import os
    import numpy as np

    Folder_Location = os.getcwd()+'//GroundMotions//UniformHazard'
    T = []
    Sa = []

    for subdir, dirs, files in os.walk(Folder_Location):
        for file in files:
            print(('Opening File %s'%file))
            lines = open(Folder_Location+'//'+file,'r').readlines()
            #Get Period
            T.append(float(lines[1].split()[6]))

            import OpenSees.OpenSeesPostProcessor as OPP
            Data = OPP.get_array_from_string_vector(lines[2:])
            Data = np.double(Data)
            LatData = Data[:,1]
            LongData = Data[:,0]
            HazardData = Data[:,2]

            ReducedData = []
            rowsCopied = []
            for i in range(0,len(LatData)-1):
                if abs(LatData[i]-Lat) <= 0.05 and abs(LatData[i+1]-Lat) <= 0.05:
                    if abs(LongData[i]-Long) <= 0.05 and abs(LongData[i+1]-Long) <= 0.05:
                        if not(i in rowsCopied):
                            ReducedData.append([LatData[i],LongData[i],HazardData[i]])
                            rowsCopied.append(i)
                        if not(i+1 in rowsCopied):
                            ReducedData.append([LatData[i+1],LongData[i+1],HazardData[i+1]])
                            rowsCopied.append(i+1)
                # if LongData[i]<Long and LongData[i+1]>Long:
                #     if abs(LongData[i]-Long) < 0.05:
                #         if not(i in rowsCopied):
                #             ReducedData.append([LatData[i],LongData[i],HazardData[i]])
                #             rowsCopied.append(i)
                #         if not(i+1 in rowsCopied):
                #             ReducedData.append([LatData[i+1],LongData[i+1],HazardData[i+1]])
                #             rowsCopied.append(i+1)

            ReducedData = np.array(ReducedData)

            import scipy.interpolate as spint
            f = spint.interp2d(ReducedData[:,0],ReducedData[:,1],ReducedData[:,2])
            Sa.append(f(Lat,Long))

    Hazard = np.transpose(np.array([T, Sa]))
    Hazard = Hazard[np.argsort(Hazard[:,0])]
    Hazard = np.transpose(Hazard)
    return Hazard

def GetSa2in50(Lat, Long, Period):
    # Return 2 Percent in 50 Hazard SA
    Hazard = HazardCurve(Lat, Long)
    import numpy as np
    return np.interp(Period, Hazard[0],Hazard[1])

def PlotGMHysteresis_OpenSees(GMData, Dt, T, Dy, Zeta, PostYieldStrengthPercentage, Plot_Name, Plot=False):
    import numpy as np
    from OpenSees import OpenSeesHelper as oh
    from OpenSees import OpenSeesAnimation as osa
    from OpenSees import OpenSeesCollector as oc
    from OpenSees import OpenSeesPostProcessor

    g=386.4

    #Mass
    m = 1
    wn = 2*np.pi/T

    #Stiffness
    k = wn**2.0*m

    #Damping 2%
    c=Zeta*2*m*wn

    import numpy as np
    import geometryhelper as gh
    import StickModelHelper
    import tclhelper as th

    Height = 1

    ########################## Initializing ##########################
    import time
    timestamp = time.strftime("%y%m%d-%H%M-B")
    ModelName = 'FindHysteresis'
    tcl = open(os.getcwd()+'//tcl//%s-%s.tcl'%(ModelName,timestamp), 'w') # Saving TCL File Location
    OpenSees = oh.OpenSeesHelper() # Opening the OpenSees Helper Class
    OpenSeesCollector = oc.OpenSeesCollector() #opening the OpenSees Node/Element/Material Collector
    OpenSees.Clock().write_to_tcl_start(tcl)

    ########################## Setup and Source Definition ##########################
    OpenSees.TclTitle('File Setup').write_to_tcl(tcl)
    OpenSees.Model(2, 3).open_write_to_tcl(tcl) # Number of Dim, Number of DOFs

    ########################## Define Building Geometry, Nodes and Constraints ##########################
    OpenSees.TclTitle('Geometry Setup').write_to_tcl(tcl)

    #Defining Grids
    XGrids = np.array([0])
    YGrids = np.array([0])
    ZGrids = np.array(list(range(0, 2)))*Height

    #Defining Floor Weights
    FloorMass = [0,m]

    #Creating OpenSees Nodal Elements
    ONodes = gh.get_node_points(XGrids, YGrids, ZGrids, True)
    OpenSeesCollector.assign_ONodes(ONodes)

    ########################## Define Geometric Transformations ##########################
    OpenSees.TclTitle('Define Geometric Transformations').write_to_tcl(tcl)

    #Define Geometry Transformations for Beams and Column
    LinearTransformation = OpenSees.LinearGeometryTransformation('1', 0, 0, 0, True)
    LinearTransformation.write_to_tcl(tcl)

    ########################## Define Materials and Sections ##########################
    OpenSees.TclTitle('Define Materials and Sections').write_to_tcl(tcl)

    ########################## Define Rotational Springs for Plastic Hinge ##########################

    #Define Rotational Spring
    # Spring = StickModelHelper.get_rotational_spring_elastic_plastic_kinematic_hardening(OpenSeesCollector, k, Dy, PostYieldStrengthPercentage*k, Height)
    Spring = StickModelHelper.get_rotational_spring_elastic_hardening(OpenSeesCollector, k, Dy, PostYieldStrengthPercentage*k, Height)


    ########################## Define Elements ##########################
    OpenSees.TclTitle('Define Elements').write_to_tcl(tcl)

    StickModelHelper.get_stick_element(OpenSeesCollector, Spring, Height, LinearTransformation)

    ########################## Writing to TCL File ##########################

    # Setting Nodes as Used
    for object in set(OpenSeesCollector._OElements):
        object._ONode_I.setAsUsed()
        object._ONode_J.setAsUsed()

    #Writing Nodes to File
    OpenSees.TclTitle('Defining Nodes').write_to_tcl(tcl)
    for i in range(0, len(ONodes)):
        ONodes[i].write_to_tcl(tcl)

    #Defining Fixity
    OpenSees.TclTitle('Defining Fixity').write_to_tcl(tcl)
    OpenSees.Supports(OpenSeesCollector.get_ONodes(gh.get_node_name(0,0,0)),1,0,1,0,0,1,True).write_to_tcl(tcl)

    #Defining Masses
    OpenSees.TclTitle('Defining Masses').write_to_tcl(tcl)
    for i in range(1, len(ZGrids)):
        Nodes = OpenSeesCollector.get_ONodes_at_Elevation(ZGrids[i])
        for node in Nodes:
            if node._Type == 'Main':
                OpenSees.Mass(node, FloorMass[i], 1e-6, 1e-6, 1e-6, 1e-6, 1e-6, True).write_to_tcl(tcl)

    #Write Element from OpenSees Collector
    OpenSees.TclTitle('Writing Materials').write_to_tcl(tcl)
    for i in range(0,len(OpenSeesCollector._OMaterials)):
        OpenSeesCollector._OMaterials[i].write_to_tcl(tcl)

    #Write Sections from OpenSees Collector
    OpenSees.TclTitle('Writing Sections').write_to_tcl(tcl)
    for i in range(0,len(OpenSeesCollector._OSections)):
        OpenSeesCollector._OSections[i].write_to_tcl(tcl)

    #Write Elements from OpenSees Collector
    OpenSees.TclTitle('Writing Elements').write_to_tcl(tcl)
    for i in range(0,len(OpenSeesCollector._OElements)):
        OpenSeesCollector._OElements[i].write_to_tcl(tcl)

    #Write Constraints from OpenSees Collector
    OpenSees.TclTitle('Writing Constraints').write_to_tcl(tcl)
    for i in range(0,len(OpenSeesCollector._OConstraints)):
        OpenSeesCollector._OConstraints[i].write_to_tcl(tcl)

    #Write Shells from OpenSees Collector
    OpenSees.TclTitle('Writing Shells').write_to_tcl(tcl)
    for i in range(0,len(OpenSeesCollector._OQuadrilateral)):
        OpenSeesCollector._OQuadrilateral[i].write_to_tcl(tcl)

    ########################## Eigenvalue Analysis ##########################
    OpenSees.TclTitle('Eigenvalue Analysis').write_to_tcl(tcl)

    OpenSees.EigenValueAnalysis(1).write_to_tcl(tcl)

    ########################## Rayleigh Damping ##########################
    OpenSees.TclTitle('Rayleigh Damping').write_to_tcl(tcl)

    OpenSees.RayleighDamping(Zeta, 2*np.pi/T).write_to_tcl(tcl)

    ########################## Loads ##########################
    OpenSees.TclTitle('Loads').write_to_tcl(tcl)

    ########################## Recorders ##########################
    OpenSees.TclTitle('Recorder Setup').write_to_tcl(tcl)

    Nodes = OpenSeesCollector.get_ONodes_at_Elevation(ZGrids[len(ZGrids)-1])
    OutputFolder = 'Results'
    Displacement_File_Name = '%s-NodeD-%s.dat'%(ModelName,timestamp)
    OpenSees.NodeRecorder('Node', Displacement_File_Name, OutputFolder, Nodes, '1 2 3', 'disp').write_to_tcl(tcl)
    Velocity_File_Name = '%s-NodeV-%s.dat'%(ModelName,timestamp)
    OpenSees.NodeRecorder('Node', Velocity_File_Name, OutputFolder, Nodes, '1 2 3', 'vel').write_to_tcl(tcl)

    #Define TimeSeries
    OpenSees.TimeSeries(1,GMData,Dt, 1.0).write_to_tcl_start(tcl)
    Acceleration_File_Name = '%s-NodeA-%s.dat'%(ModelName,timestamp)
    OpenSees.NodeRecorder('Node', Acceleration_File_Name, OutputFolder, Nodes, '1', 'accel', extra='-timeSeries 1').write_to_tcl(tcl)

    #Reactions Recorder

    Reaction_File_Name = '%s-NodeReact-%s.dat'%(ModelName,timestamp)
    OpenSees.NodeRecorder('Node',Reaction_File_Name,OutputFolder,OpenSeesCollector.get_Main_ONodes_at_Elevation(ZGrids[len(ZGrids)-2]),'3','reaction').write_to_tcl(tcl)

    ########################## Display Results ##########################
    OpenSees.TclTitle('Display Results').write_to_tcl(tcl)

    ########################## Gravity Analysis ##########################
    OpenSees.TclTitle('Gravity Analysis').write_to_tcl(tcl)

    ########################## Pushover Analysis ##########################
    OpenSees.TclTitle('Pushover Analysis').write_to_tcl(tcl)

    ########################## Time History Analysis ##########################
    OpenSees.TclTitle('Time History Analysis').write_to_tcl(tcl)

    OpenSees.DynamicAnalysis_GMData(GMData,Dt,1.0,0).write_to_tcl(tcl)

    ########################## Close File ##########################
    OpenSees.TclTitle('Closing File').write_to_tcl(tcl)

    OpenSees.Clock().write_to_tcl_end(tcl)
    th.print_line(tcl, 'wipe al b l;')
    th.print_line(tcl, 'puts "Models Run Complete";')
    tcl.close()

    ########################## Plot Geometry ##########################

    ########################## Run OpenSees Script ##########################
    import subprocess
    #Single Processor
    run = subprocess.Popen([os.getcwd()+'//tcl//OpenSees', os.getcwd()+'//tcl//%s-%s.tcl'%(ModelName,timestamp)],cwd=os.getcwd()+'//tcl')
    run.wait()

    #Multiple Processors
    # run = subprocess.Popen(['mpiexec','-np','1',os.getcwd()+'//tcl//OpenSeesMP', os.getcwd()+'//tcl//%s-%s.tcl'%(ModelName,timestamp)],cwd=os.getcwd()+'//tcl')
    # run.wait()

    ########################## Plot Push Over ##########################
    Displ = OpenSeesPostProcessor.get_table_from_data_file(os.getcwd()+'//tcl//'+OutputFolder, Displacement_File_Name)
    Vel = OpenSeesPostProcessor.get_table_from_data_file(os.getcwd()+'//tcl//'+OutputFolder, Velocity_File_Name)
    Acc = OpenSeesPostProcessor.get_table_from_data_file(os.getcwd()+'//tcl//'+OutputFolder, Acceleration_File_Name)
    Reac = OpenSeesPostProcessor.get_table_from_data_file(os.getcwd()+'//tcl//'+OutputFolder, Reaction_File_Name)

    import matplotlib.pylab as plt

    t = Acc[:,0]
    ag = GMData[:]/g
    a = Acc[:,1]
    fs = Reac[:,1]/(k*Dy)
    u = Displ[:,1]/Dy

    if Plot == True:
        ani = osa.HysteresisAnimation(t, ag, a, fs, u)
        # ani.save('Figures/Hysteresis.mp4', fps=100)
        plt.show()

    # import matplotlib.pyplot as plt
    #
    # fig, ((ax1, ax2), (ax3, ax4), (ax5, ax6)) = plt.subplots(nrows=3, ncols=2)
    # fig = plt.figure(1, figsize=(12,12))
    # ax1.set_title(Plot_Name)
    #
    # ax1.plot(Acc[:,0],GMData/g)
    # ax1.set_xlabel('Time, seconds')
    # ax1.set_ylabel('Ground Acc, g')
    # ax1.grid(True)
    #
    # ax2.axis('off')
    #
    # ax3.plot(Acc[:,0],Reac[:,1])
    # ax3.set_xlabel('Time, seconds')
    # ax3.set_ylabel('Base Shear/W')
    # ax3.grid(True)
    #
    # ax4.plot(Displ[:,1],Reac[:,1])
    # ax4.set_xlabel('Displacement, in')
    # ax4.set_ylabel('Base Shear/W')
    # ax4.grid(True)
    #
    # ax5.axis('off')
    #
    # ax6.plot(Displ[:,1],Acc[:,0])
    # ax6.set_xlabel('Displacement, in')
    # ax6.set_ylabel('Time, seconds')
    # ax6.grid(True)
    #
    # plt.tight_layout()
    #
    # plt.show()
    # # plt.savefig('HysteresisAndAll.png')

    os.remove(os.getcwd()+'//tcl//%s-%s.tcl'%(ModelName,timestamp))
    os.remove(os.getcwd()+'//tcl//'+OutputFolder+'//'+Displacement_File_Name)
    os.remove(os.getcwd()+'//tcl//'+OutputFolder+'//'+Velocity_File_Name)
    os.remove(os.getcwd()+'//tcl//'+OutputFolder+'//'+Acceleration_File_Name)
    os.remove(os.getcwd()+'//tcl//'+OutputFolder+'//'+Reaction_File_Name)

    return

def PlotGMHysteresis_OpenSees_Deterioration(GMData, Dt, T, Dy, Zeta, PostYieldStrengthPercentage, LamdaSCA, LamdaK, rateDeterioration, R, Mu, Kappa, PostCappingStiffness, Plot_Name, Plot=False, OpenSeesCommand = os.getcwd()+'//tcl//OpenSees', WriteAnimation=False):
    import numpy as np
    from OpenSees import OpenSeesHelper as oh
    from OpenSees import OpenSeesCollector as oc
    from OpenSees import OpenSeesPostProcessor

    g=386.4

    #Mass
    m = 1
    wn = 2*np.pi/T

    #Stiffness
    k = wn**2.0*m

    #Damping 2%
    c=Zeta*2*m*wn

    import numpy as np
    import geometryhelper as gh
    import StickModelHelper
    import tclhelper as th
    import os

    ########################## Initializing ##########################
    import time
    # randomnumber = random.randint(1, 10000000000)
    import uuid
    randomnumber = str(uuid.uuid4().hex.upper()[0:12])
    timestamp = time.strftime("%y%m%d-%H%M-")+'%s'%randomnumber
    ModelName = 'FindHysteresis'
    tcl = open(os.getcwd()+'//tcl//%s-%s.tcl'%(ModelName,timestamp), 'w') # Saving TCL File Location
    logfilename = '%s-%s-LOG.dat'%(ModelName,timestamp)
    OpenSees = oh.OpenSeesHelper() # Opening the OpenSees Helper Class
    OpenSeesCollector = oc.OpenSeesCollector() # Opening the OpenSees Node/Element/Material Collector
    OpenSees.Clock().write_to_tcl_start(tcl)

    ########################## Setup and Source Definition ##########################
    OpenSees.TclTitle('File Setup').write_to_tcl(tcl)
    OpenSees.Model(2, 3).open_write_to_tcl(tcl,logfilename) # Number of Dim, Number of DOFs

    ########################## Define Building Geometry, Nodes and Constraints ##########################
    OpenSees.TclTitle('Geometry Setup').write_to_tcl(tcl)

    #Defining Grids
    XGrids = np.array([0])
    YGrids = np.array([0])
    ZGrids = np.array([0])

    #Defining Floor Weights
    FloorMass = [m]

    #Creating OpenSees Nodal Elements
    ONodes = gh.get_node_points(XGrids, YGrids, ZGrids, True)
    OpenSeesCollector.assign_ONodes(ONodes)

    ########################## Define Geometric Transformations ##########################
    OpenSees.TclTitle('Define Geometric Transformations').write_to_tcl(tcl)

    ########################## Define Materials and Sections ##########################
    OpenSees.TclTitle('Define Materials and Sections').write_to_tcl(tcl)

    ########################## Define Rotational Springs for Plastic Hinge ##########################

    #Define Rotational Spring
    Spring = StickModelHelper.get_rotational_spring_elastic_plastic_kinematic_hardening(OpenSeesCollector, k, Dy, PostYieldStrengthPercentage, LamdaSCA, LamdaK, rateDeterioration, R, Mu, Kappa, PostCappingStiffness)

    ########################## Define Elements ##########################
    OpenSees.TclTitle('Define Elements').write_to_tcl(tcl)

    #Define P-Delta Column
    SupportNode = OpenSeesCollector.get_ONodes(gh.get_node_name(0,0,0))
    MassNode = OpenSeesCollector.create_ONode(0,0,0,0,0,0,'Constraint',True)

    SpringElement = OpenSees.ElementZeroLength(OpenSeesCollector.get_ZeroLength_Name(0),SupportNode,MassNode,[Spring],'1','Non-Linear Spring')

    OpenSeesCollector.add_OElement(SpringElement)

    ########################## Writing to TCL File ##########################

    # Setting Nodes as Used
    for object in set(OpenSeesCollector._OElements):
        object._ONode_I.setAsUsed()
        object._ONode_J.setAsUsed()

    #Writing Nodes to File
    OpenSees.TclTitle('Defining Nodes').write_to_tcl(tcl)
    for i in range(0, len(ONodes)):
        ONodes[i].write_to_tcl(tcl)

    #Defining Fixity
    OpenSees.TclTitle('Defining Fixity').write_to_tcl(tcl)
    OpenSees.Supports(SupportNode,1,0,1,0,0,1,True).write_to_tcl(tcl)
    OpenSees.Supports(MassNode,0,0,1,0,0,1,True).write_to_tcl(tcl)

    #Defining Masses
    OpenSees.TclTitle('Defining Masses').write_to_tcl(tcl)
    OpenSees.Mass(MassNode, FloorMass[0], 1e-6, 1e-6, 1e-6, 1e-6, 1e-6, True).write_to_tcl(tcl)

    #Write Element from OpenSees Collector
    OpenSees.TclTitle('Writing Materials').write_to_tcl(tcl)
    for i in range(0,len(OpenSeesCollector._OMaterials)):
        OpenSeesCollector._OMaterials[i].write_to_tcl(tcl)

    #Write Sections from OpenSees Collector
    OpenSees.TclTitle('Writing Sections').write_to_tcl(tcl)
    for i in range(0,len(OpenSeesCollector._OSections)):
        OpenSeesCollector._OSections[i].write_to_tcl(tcl)

    #Write Elements from OpenSees Collector
    OpenSees.TclTitle('Writing Elements').write_to_tcl(tcl)
    for i in range(0,len(OpenSeesCollector._OElements)):
        OpenSeesCollector._OElements[i].write_to_tcl(tcl)

    #Write Constraints from OpenSees Collector
    OpenSees.TclTitle('Writing Constraints').write_to_tcl(tcl)
    for i in range(0,len(OpenSeesCollector._OConstraints)):
        OpenSeesCollector._OConstraints[i].write_to_tcl(tcl)

    #Write Shells from OpenSees Collector
    OpenSees.TclTitle('Writing Shells').write_to_tcl(tcl)
    for i in range(0,len(OpenSeesCollector._OQuadrilateral)):
        OpenSeesCollector._OQuadrilateral[i].write_to_tcl(tcl)

    ########################## Eigenvalue Analysis ##########################
    OpenSees.TclTitle('Eigenvalue Analysis').write_to_tcl(tcl)

    # OpenSees.EigenValueAnalysis(1).write_to_tcl(tcl)

    ########################## Rayleigh Damping ##########################
    OpenSees.TclTitle('Rayleigh Damping').write_to_tcl(tcl)

    OpenSees.RayleighDamping(Zeta, wn).write_to_tcl(tcl)

    ########################## Loads ##########################
    OpenSees.TclTitle('Loads').write_to_tcl(tcl)

    ########################## Recorders ##########################
    OpenSees.TclTitle('Recorder Setup').write_to_tcl(tcl)

    Nodes = [MassNode]
    OutputFolder = 'Results'
    Displacement_File_Name = '%s-NodeD-%s.dat'%(ModelName,timestamp)
    OpenSees.NodeRecorder('Node', Displacement_File_Name, OutputFolder, Nodes, '1 2 3', 'disp').write_to_tcl(tcl)
    Velocity_File_Name = '%s-NodeV-%s.dat'%(ModelName,timestamp)
    OpenSees.NodeRecorder('Node', Velocity_File_Name, OutputFolder, Nodes, '1 2 3', 'vel').write_to_tcl(tcl)

    #Define TimeSeries
    GMData = np.array(np.array(GMData).tolist()+np.zeros(int(20/Dt)).tolist()) #Add 20 Seconds of Zeros
    OpenSees.TimeSeries(1,GMData,Dt, 1.0).write_to_tcl_start(tcl)
    Acceleration_File_Name = '%s-NodeA-%s.dat'%(ModelName,timestamp)
    OpenSees.NodeRecorder('Node', Acceleration_File_Name, OutputFolder, Nodes, '1', 'accel', extra='-timeSeries 1').write_to_tcl(tcl)

    #Reactions Recorder

    Nodes = [SupportNode]
    Reaction_File_Name = '%s-NodeReact-%s.dat'%(ModelName,timestamp)
    OpenSees.NodeRecorder('Node',Reaction_File_Name,OutputFolder,Nodes,'1 2 3','reaction').write_to_tcl(tcl)

    ########################## Display Results ##########################
    OpenSees.TclTitle('Display Results').write_to_tcl(tcl)

    ########################## Gravity Analysis ##########################
    OpenSees.TclTitle('Gravity Analysis').write_to_tcl(tcl)

    # OpenSees.Display().write_to_tcl(tcl)

    ########################## Time History Analysis ##########################
    OpenSees.TclTitle('Time History Analysis').write_to_tcl(tcl)

    OpenSees.DynamicAnalysis_GMData(GMData, Dt, 1.0, 0).write_to_tcl(tcl)

    # ########################## Pushover Analysis ##########################
    OpenSees.TclTitle('Pushover Analysis').write_to_tcl(tcl)

    NodeLoads = []
    Nodes = [MassNode]
    NodeLoads.append(OpenSees.NodeLoadsPoint(Nodes, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, True))
    PushLoadPattern = OpenSees.LoadPatternPlain('300', NodeLoads)
    PushLoadPattern.write_to_tcl(tcl)

    TestDrift = Dy*0.01

    OpenSees.Pushover_Simple(Nodes[0], '1', TestDrift, TestDrift, SupportNode).write_to_tcl(tcl)

    ########################## Eigenvalue Analysis ##########################
    OpenSees.TclTitle('Eigenvalue Analysis').write_to_tcl(tcl)

    OpenSees.EigenValueAnalysis(1,fullGenLapack=True).write_to_tcl(tcl)

    ########################## Close File ##########################
    OpenSees.TclTitle('Closing File').write_to_tcl(tcl)

    OpenSees.Clock().write_to_tcl_end(tcl)
    th.print_line(tcl, 'wipe al b l;')
    th.print_line(tcl, 'puts "Models Run Complete";')
    tcl.close()

    ########################## Plot Geometry ##########################

    ########################## Run OpenSees Script ##########################
    import subprocess
    #Single Processor
    OutputFile = open(os.getcwd()+'//tcl//%s-%s.out'%(ModelName,timestamp),'w')
    print(('Running OpenSees '+os.getcwd()+'//tcl//%s-%s.tcl'%(ModelName,timestamp)))
    run = subprocess.Popen([OpenSeesCommand, os.getcwd()+'//tcl//%s-%s.tcl'%(ModelName,timestamp)], cwd=os.getcwd()+'//tcl', stdout=OutputFile, stderr=subprocess.STDOUT)

    run.wait()
    OutputFile.close()
    os.remove(os.getcwd()+'//tcl//%s-%s.out'%(ModelName,timestamp))

    #Open Log File and Get NL Period
    Log = OpenSeesPostProcessor.get_string_table_from_data_file(os.getcwd()+'//tcl//',logfilename)
    ElongatedPeriod = None
    PushOverStiffness = []
    for i in range(len(Log)):
        for j in range(len(Log[i])):
            if Log[i][j]=="T1":
                ElongatedPeriod = float(Log[i][j+2])
            if Log[i][j]=='SDOFStiffness:':
                PushOverStiffness.append(float(Log[i][j+1]))

    #Find Stiffness From PushOver
    if len(PushOverStiffness)>0:
        AveragePushOverStiffness = np.mean(PushOverStiffness)
        ElongatedPeriodFromPushOver = (2.0*np.pi/AveragePushOverStiffness**.5)
    else:
        AveragePushOverStiffness = 0
        ElongatedPeriodFromPushOver = 1000000.0

    #Display Some Results for Back Checking
    print(('Elastic Period: %.2f'%T))
    print(('Elongated Period: %.2f'%ElongatedPeriod))
    print(('Elastic Stiffness: %.2f'%k))
    print(('NL Stiffness: %.2f'%AveragePushOverStiffness))
    if AveragePushOverStiffness != 0:
        print(('Calculated NL Period: %.2f'%(2.0*np.pi/AveragePushOverStiffness**.5)))
    else:
        print('Calculated NL Period: infinity')

    ########################## Plot Push Over ##########################
    Displ = OpenSeesPostProcessor.get_table_from_data_file(os.getcwd()+'//tcl//'+OutputFolder, Displacement_File_Name)
    Vel = OpenSeesPostProcessor.get_table_from_data_file(os.getcwd()+'//tcl//'+OutputFolder, Velocity_File_Name)
    Acc = OpenSeesPostProcessor.get_table_from_data_file(os.getcwd()+'//tcl//'+OutputFolder, Acceleration_File_Name)
    Reac = OpenSeesPostProcessor.get_table_from_data_file(os.getcwd()+'//tcl//'+OutputFolder, Reaction_File_Name)

    DataPoints = len(GMData)

    t = Acc[0:DataPoints,0]
    ag = GMData[:]/g
    a = Acc[0:DataPoints,1]
    fs = -1*Reac[0:DataPoints,1]/(k*Dy)
    u = Displ[0:DataPoints,1]/Dy

    MaxU = max(abs(Displ[0:DataPoints,1]))
    MaxV = max(abs(Vel[0:DataPoints,1]))
    MaxA = max(abs(Acc[0:DataPoints,1]))
    MaxF = max(abs(Reac[0:DataPoints,1]))

    DArray = Displ[:,1]
    for i in range(0,len(DArray)):
        if abs(DArray[i]) == max(abs(DArray)):
            TimeAtMax = float(i)/float(len(DArray))
            break

    Stiffness = []
    Collapse = 0

    for i in range(1,len(Displ[0:DataPoints,1])):
        if (Displ[i,1]-Displ[i-1,1]) != 0.0:
            Stiffness.append(abs(Reac[i,1]-Reac[i-1,1])/abs(Displ[i,1]-Displ[i-1,1]))
        else:
            Stiffness.append(k)


    # Collapse = []
    # for i in range(1,len(Displ[:,1])):
    #     if Stiffness[i-1] < 0.01*k and abs(fs[i]) < 0.05:
    #         Collapse.append(1)
    # if len(Collapse) > 10:
    #     Collapse = 1
    # else:
    #     Collapse = 0

    PeriodThreshholdForCollapse = 4
    if ElongatedPeriodFromPushOver > PeriodThreshholdForCollapse*T:
        Collapse = 1
    else:
        Collapse = 0

    if Plot == True:
        import matplotlib.pylab as plt
        from OpenSees import OpenSeesAnimation as osa

        #Another Plot
        fig, ((ax1),(ax2)) = plt.subplots(nrows=2, ncols=1, figsize=(3.5,6))
        ax1.plot(t[1:],Stiffness,'--',linewidth=2.0,label='Stiffness')
        ax2.plot(t[1:],abs(fs[1:]),'--',linewidth=2.0,label='Stiffness')
        ax1.set_title('Collaspe = %d'%Collapse)
        plt.show(block=False)

        PlotText = 'Ground Motion ID: %s\nT = %.2f \n$\zeta$ = %.2f \nR = %.0f \n$K_y$ = %.2f $K_e$ \n$K_c$ = %.2f $K_e$ \n$\lambda_{s,c,a}$ = %.2f \n$\lambda_k$ = %.2f '%(Plot_Name, T, Zeta, R, PostYieldStrengthPercentage, PostCappingStiffness, LamdaSCA, LamdaK)
        ani = osa.HysteresisAnimation(t, ag, a, fs, u, PlotText)
        if WriteAnimation:
            ani.save('Figures/Hysteresis.mp4', fps=30, bitrate=20000) #The Script already divides the frame rate by 10
            ani.ConvertToGIF(os.getcwd()+'/Figures/','Hysteresis.mp4')
        plt.show()

    os.remove(os.getcwd()+'//tcl//%s-%s.tcl'%(ModelName,timestamp))
    os.remove(os.getcwd()+'//tcl//'+OutputFolder+'//'+Displacement_File_Name)
    os.remove(os.getcwd()+'//tcl//'+OutputFolder+'//'+Velocity_File_Name)
    os.remove(os.getcwd()+'//tcl//'+OutputFolder+'//'+Acceleration_File_Name)
    os.remove(os.getcwd()+'//tcl//'+OutputFolder+'//'+Reaction_File_Name)
    os.remove(os.getcwd()+'//tcl//'+logfilename)

    return MaxU, MaxV, MaxA, MaxF, TimeAtMax, Collapse, ElongatedPeriodFromPushOver

def PlotGMHysteresis_OpenSees_Deterioration_Peak_Oriented(GMData, Dt, T, Dy, Zeta, PostYieldStrengthPercentage, LamdaSCA, LamdaK, rateDeterioration, R, Mu, Kappa, PostCappingStiffness, Plot_Name, Plot=False, OpenSeesCommand = 'OpenSeesSP', WriteAnimation=False):
    import numpy as np
    from OpenSees import OpenSeesHelper as oh
    from OpenSees import OpenSeesCollector as oc
    from OpenSees import OpenSeesPostProcessor

    g=386.4

    #Mass
    m = 1
    wn = 2*np.pi/T

    #Stiffness
    k = wn**2.0*m

    #Damping 2%
    c=Zeta*2*m*wn

    import numpy as np
    import geometryhelper as gh
    import StickModelHelper
    import tclhelper as th
    import os

    ########################## Initializing ##########################
    import time
    # randomnumber = random.randint(1, 10000000000)
    import uuid
    randomnumber = str(uuid.uuid4().hex.upper()[0:12])
    timestamp = time.strftime("%y%m%d-%H%M-")+'%s'%randomnumber
    ModelName = 'FindHysteresis'
    tcl = open(os.getcwd()+'//tcl//%s-%s.tcl'%(ModelName,timestamp), 'w') # Saving TCL File Location
    logfilename = '%s-%s-LOG.dat'%(ModelName,timestamp)
    OpenSees = oh.OpenSeesHelper() # Opening the OpenSees Helper Class
    OpenSeesCollector = oc.OpenSeesCollector() # Opening the OpenSees Node/Element/Material Collector
    OpenSees.Clock().write_to_tcl_start(tcl)

    ########################## Setup and Source Definition ##########################
    OpenSees.TclTitle('File Setup').write_to_tcl(tcl)
    OpenSees.Model(2, 3).open_write_to_tcl(tcl,logfilename) # Number of Dim, Number of DOFs

    ########################## Define Building Geometry, Nodes and Constraints ##########################
    OpenSees.TclTitle('Geometry Setup').write_to_tcl(tcl)

    #Defining Grids
    XGrids = np.array([0])
    YGrids = np.array([0])
    ZGrids = np.array([0])

    #Defining Floor Weights
    FloorMass = [m]

    #Creating OpenSees Nodal Elements
    ONodes = gh.get_node_points(XGrids, YGrids, ZGrids, True)
    OpenSeesCollector.assign_ONodes(ONodes)

    ########################## Define Geometric Transformations ##########################
    OpenSees.TclTitle('Define Geometric Transformations').write_to_tcl(tcl)

    ########################## Define Materials and Sections ##########################
    OpenSees.TclTitle('Define Materials and Sections').write_to_tcl(tcl)

    ########################## Define Rotational Springs for Plastic Hinge ##########################
    #Define Rotational Spring
    Spring = StickModelHelper.get_rotational_spring_elastic_plastic_ibarra_peak_oriented(OpenSeesCollector, k, Dy, PostYieldStrengthPercentage, LamdaSCA, LamdaK, rateDeterioration, R, Mu, Kappa, PostCappingStiffness)
    #get_rotational_spring_elastic_plastic_ibarra_peak_oriented
    #get_rotational_spring_elastic_plastic_kinematic_hardening

    ########################## Define Elements ##########################
    OpenSees.TclTitle('Define Elements').write_to_tcl(tcl)

    #Define P-Delta Column
    SupportNode = OpenSeesCollector.get_ONodes(gh.get_node_name(0,0,0))
    MassNode = OpenSeesCollector.create_ONode(0,0,0,0,0,0,'Constraint',True)

    SpringElement = OpenSees.ElementZeroLength(OpenSeesCollector.get_ZeroLength_Name(0),SupportNode,MassNode,[Spring],'1','Non-Linear Spring')

    OpenSeesCollector.add_OElement(SpringElement)


    ########################## Writing to TCL File ##########################

    # Setting Nodes as Used
    for object in set(OpenSeesCollector._OElements):
        object._ONode_I.setAsUsed()
        object._ONode_J.setAsUsed()

    #Writing Nodes to File
    OpenSees.TclTitle('Defining Nodes').write_to_tcl(tcl)
    for i in range(0, len(ONodes)):
        ONodes[i].write_to_tcl(tcl)

    #Defining Fixity
    OpenSees.TclTitle('Defining Fixity').write_to_tcl(tcl)
    OpenSees.Supports(SupportNode,1,0,1,0,0,1,True).write_to_tcl(tcl)
    OpenSees.Supports(MassNode,0,0,1,0,0,1,True).write_to_tcl(tcl)

    #Defining Masses
    OpenSees.TclTitle('Defining Masses').write_to_tcl(tcl)
    OpenSees.Mass(MassNode, FloorMass[0], 1e-6, 1e-6, 1e-6, 1e-6, 1e-6, True).write_to_tcl(tcl)

    #Write Element from OpenSees Collector
    OpenSees.TclTitle('Writing Materials').write_to_tcl(tcl)
    for i in range(0,len(OpenSeesCollector._OMaterials)):
        OpenSeesCollector._OMaterials[i].write_to_tcl(tcl)

    #Write Sections from OpenSees Collector
    OpenSees.TclTitle('Writing Sections').write_to_tcl(tcl)
    for i in range(0,len(OpenSeesCollector._OSections)):
        OpenSeesCollector._OSections[i].write_to_tcl(tcl)

    #Write Elements from OpenSees Collector
    OpenSees.TclTitle('Writing Elements').write_to_tcl(tcl)
    for i in range(0,len(OpenSeesCollector._OElements)):
        OpenSeesCollector._OElements[i].write_to_tcl(tcl)

    #Write Constraints from OpenSees Collector
    OpenSees.TclTitle('Writing Constraints').write_to_tcl(tcl)
    for i in range(0,len(OpenSeesCollector._OConstraints)):
        OpenSeesCollector._OConstraints[i].write_to_tcl(tcl)

    #Write Shells from OpenSees Collector
    OpenSees.TclTitle('Writing Shells').write_to_tcl(tcl)
    for i in range(0,len(OpenSeesCollector._OQuadrilateral)):
        OpenSeesCollector._OQuadrilateral[i].write_to_tcl(tcl)

    ########################## Eigenvalue Analysis ##########################
    OpenSees.TclTitle('Eigenvalue Analysis').write_to_tcl(tcl)

    # OpenSees.EigenValueAnalysis(1).write_to_tcl(tcl)

    ########################## Rayleigh Damping ##########################
    OpenSees.TclTitle('Rayleigh Damping').write_to_tcl(tcl)

    OpenSees.RayleighDamping(Zeta, wn).write_to_tcl(tcl)

    ########################## Loads ##########################
    OpenSees.TclTitle('Loads').write_to_tcl(tcl)

    ########################## Recorders ##########################
    OpenSees.TclTitle('Recorder Setup').write_to_tcl(tcl)

    Nodes = [MassNode]
    OutputFolder = 'Results'
    Displacement_File_Name = '%s-NodeD-%s.dat'%(ModelName,timestamp)
    OpenSees.NodeRecorder('Node', Displacement_File_Name, OutputFolder, Nodes, '1 2 3', 'disp').write_to_tcl(tcl)
    Velocity_File_Name = '%s-NodeV-%s.dat'%(ModelName,timestamp)
    OpenSees.NodeRecorder('Node', Velocity_File_Name, OutputFolder, Nodes, '1 2 3', 'vel').write_to_tcl(tcl)

    #Define TimeSeries
    GMData = np.array(np.array(GMData).tolist()+np.zeros(int(20/Dt)).tolist()) #Add 20 Seconds of Zeros for Free Vibration at the end
    OpenSees.TimeSeries(1,GMData,Dt, 1.0).write_to_tcl_start(tcl)
    Acceleration_File_Name = '%s-NodeA-%s.dat'%(ModelName,timestamp)
    OpenSees.NodeRecorder('Node', Acceleration_File_Name, OutputFolder, Nodes, '1', 'accel', extra='-timeSeries 1').write_to_tcl(tcl)

    #Reactions Recorder

    Nodes = [SupportNode]
    Reaction_File_Name = '%s-NodeReact-%s.dat'%(ModelName,timestamp)
    OpenSees.NodeRecorder('Node',Reaction_File_Name,OutputFolder,Nodes,'1 2 3','reaction').write_to_tcl(tcl)

    ########################## Display Results ##########################
    OpenSees.TclTitle('Display Results').write_to_tcl(tcl)

    ########################## Gravity Analysis ##########################
    OpenSees.TclTitle('Gravity Analysis').write_to_tcl(tcl)

    # OpenSees.Display().write_to_tcl(tcl)

    ########################## Time History Analysis ##########################
    OpenSees.TclTitle('Time History Analysis').write_to_tcl(tcl)

    OpenSees.DynamicAnalysis_GMData(GMData, Dt, 1.0, 0).write_to_tcl(tcl)

    ########################## Eigenvalue Analysis ##########################
    OpenSees.TclTitle('Eigenvalue Analysis').write_to_tcl(tcl)

    OpenSees.EigenValueAnalysis(1,fullGenLapack=True).write_to_tcl(tcl)

    ########################### Pushover Analysis ##########################
    OpenSees.TclTitle('Pushover Analysis').write_to_tcl(tcl)

    NodeLoads = []
    Nodes = [MassNode]
    NodeLoads.append(OpenSees.NodeLoadsPoint(Nodes, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, True))
    PushLoadPattern = OpenSees.LoadPatternPlain('300', NodeLoads)
    PushLoadPattern.write_to_tcl(tcl)

    TestDrift = Dy*0.01

    OpenSees.Pushover_Simple(Nodes[0], '1', TestDrift, TestDrift, SupportNode).write_to_tcl(tcl)

    ########################## Close File ##########################
    OpenSees.TclTitle('Closing File').write_to_tcl(tcl)

    OpenSees.Clock().write_to_tcl_end(tcl)
    th.print_line(tcl, 'wipe al b l;')
    th.print_line(tcl, 'puts "Models Run Complete";')
    tcl.close()

    ########################## Plot Geometry ##########################

    ########################## Run OpenSees Script ##########################
    import subprocess
    #Single Processor
    OutputFile = open(os.getcwd()+'//tcl//%s-%s.out'%(ModelName,timestamp),'w')
    print(('Running OpenSees '+os.getcwd()+'/tcl/%s-%s.tcl'%(ModelName,timestamp)))
    run = subprocess.Popen([OpenSeesCommand, os.getcwd()+'/tcl/%s-%s.tcl'%(ModelName,timestamp)], cwd=os.getcwd()+'/tcl/', stdout=OutputFile, stderr=subprocess.STDOUT)
    run.wait()
    OutputFile.close()
    os.remove(os.getcwd()+'//tcl//%s-%s.out'%(ModelName,timestamp))

    #Open Log File and Get NL Period
    Log = OpenSeesPostProcessor.get_string_table_from_data_file(os.getcwd()+'/tcl/',logfilename)
    ElongatedPeriod = None
    PushOverStiffness = []
    for i in range(len(Log)):
        for j in range(len(Log[i])):
            if Log[i][j]=="T1":
                ElongatedPeriod = float(Log[i][j+2])
            if Log[i][j]=='SDOFStiffness:':
                PushOverStiffness.append(float(Log[i][j+1]))

    #Find Stiffness From PushOver
    if len(PushOverStiffness)>0:
        AveragePushOverStiffness = np.mean(PushOverStiffness)
        ElongatedPeriodFromPushOver = (2.0*np.pi/AveragePushOverStiffness**.5)
    else:
        AveragePushOverStiffness = 0
        ElongatedPeriodFromPushOver = 1000000.0

    #Displau Some Results for Back Checking
    # print 'Elastic Period: %.2f'%T
    # print 'Elongated Period: %.2f'%ElongatedPeriod
    # print 'Elastic Stiffness: %.2f'%k
    # print 'NL Stiffness: %.2f'%AveragePushOverStiffness
    # if AveragePushOverStiffness != 0:
    #     print 'Calculated NL Period: %.2f'%(2.0*np.pi/AveragePushOverStiffness**.5)
    # else:
    #     print 'Calculated NL Period: infinity'

    ########################## Plot Push Over ##########################
    Displ = OpenSeesPostProcessor.get_table_from_data_file(os.getcwd()+'//tcl//'+OutputFolder, Displacement_File_Name)
    Vel = OpenSeesPostProcessor.get_table_from_data_file(os.getcwd()+'//tcl//'+OutputFolder, Velocity_File_Name)
    Acc = OpenSeesPostProcessor.get_table_from_data_file(os.getcwd()+'//tcl//'+OutputFolder, Acceleration_File_Name)
    Reac = OpenSeesPostProcessor.get_table_from_data_file(os.getcwd()+'//tcl//'+OutputFolder, Reaction_File_Name)

    DataPoints = len(GMData)

    t = Acc[0:DataPoints,0]
    ag = GMData[:]/g
    a = Acc[0:DataPoints,1]
    fs = -1*Reac[0:DataPoints,1]/(k*Dy)
    u = Displ[0:DataPoints,1]/Dy

    MaxU = max(abs(Displ[0:DataPoints,1]))
    MaxV = max(abs(Vel[0:DataPoints,1]))
    MaxA = max(abs(Acc[0:DataPoints,1]))
    MaxF = max(abs(Reac[0:DataPoints,1]))

    DArray = Displ[:,1]
    for i in range(0,len(DArray)):
        if abs(DArray[i]) == max(abs(DArray)):
            TimeAtMax = float(i)/float(len(DArray))
            break

    Stiffness = []
    Collapse = 0

    for i in range(1,len(Displ[0:DataPoints,1])):
        if (Displ[i,1]-Displ[i-1,1]) != 0.0:
            Stiffness.append(abs(Reac[i,1]-Reac[i-1,1])/abs(Displ[i,1]-Displ[i-1,1]))
        else:
            Stiffness.append(k)

    # Collapse = []
    # for i in range(1,len(Displ[:,1])):
    #     if Stiffness[i-1] < 0.01*k and abs(fs[i]) < 0.05:
    #         Collapse.append(1)
    # if len(Collapse) > 10:
    #     Collapse = 1
    # else:
    #     Collapse = 0

    PeriodThreshholdForCollapse = 4
    if ElongatedPeriodFromPushOver > PeriodThreshholdForCollapse*T:
        Collapse = 1
    else:
        Collapse = 0

    # if AveragePushOverStiffness < 0:
    #     Collapse = 0
    # else:
    #     PeriodEstimateFromPushOver = 2*np.pi/(AveragePushOverStiffness)**.5
    #     if PeriodEstimateFromPushOver > 100*T:
    #         Collapse = 1
    #     else:
    #         Collapse = 0

    if Plot == True:
        import matplotlib.pylab as plt
        from OpenSees import OpenSeesAnimation as osa

        #Another Plot
        fig, ((ax1),(ax2)) = plt.subplots(nrows=2, ncols=1, figsize=(3.5,6))
        ax1.plot(t[1:],Stiffness,'--',linewidth=2.0,label='Stiffness')
        ax2.plot(t[1:],abs(fs[1:]),'--',linewidth=2.0,label='Stiffness')
        ax1.set_title('Collaspe = %d'%Collapse)
        plt.show(block=False)

        PlotText = 'Ground Motion ID: %s\nT = %.2f \n$\zeta$ = %.2f \nR = %.0f \n$K_y$ = %.2f $K_e$ \n$K_c$ = %.2f $K_e$ \n$\lambda_{s,c,a}$ = %.2f \n$\lambda_k$ = %.2f '%(Plot_Name, T, Zeta, R, PostYieldStrengthPercentage, PostCappingStiffness, LamdaSCA, LamdaK)
        ani = osa.HysteresisAnimation(t, ag, a, fs, u, PlotText)
        if WriteAnimation:
            ani.save('Figures/Hysteresis.mp4', fps=30, bitrate=20000) #The Script already divides the frame rate by 10
            ani.ConvertToGIF(os.getcwd()+'/Figures/','Hysteresis.mp4')
        plt.show()

    os.remove(os.getcwd()+'//tcl//%s-%s.tcl'%(ModelName,timestamp))
    os.remove(os.getcwd()+'//tcl//'+OutputFolder+'//'+Displacement_File_Name)
    os.remove(os.getcwd()+'//tcl//'+OutputFolder+'//'+Velocity_File_Name)
    os.remove(os.getcwd()+'//tcl//'+OutputFolder+'//'+Acceleration_File_Name)
    os.remove(os.getcwd()+'//tcl//'+OutputFolder+'//'+Reaction_File_Name)
    os.remove(os.getcwd()+'//tcl//'+logfilename)

    # print 'EigT/PushT: %.2f'%(ElongatedPeriod/ElongatedPeriodFromPushOver)

    return MaxU, MaxV, MaxA, MaxF, TimeAtMax, Collapse, ElongatedPeriodFromPushOver

def PlotGMHysteresis_OpenSees_Deterioration_Peak_Oriented_Class_Output(GMData, Dt, T, Dy, Zeta, PostYieldStrengthPercentage, LamdaSCA, LamdaK, rateDeterioration, R, Mu, Kappa, PostCappingStiffness, Plot_Name, Plot=False, OpenSeesCommand = os.getcwd()+'//tcl//OpenSees', WriteAnimation=False, FreeVibrationTime = 20):
    import numpy as np
    from OpenSees import OpenSeesHelper as oh
    from OpenSees import OpenSeesCollector as oc
    from OpenSees import OpenSeesPostProcessor

    g=386.4

    #Mass
    m = 1
    wn = 2*np.pi/T

    #Stiffness
    k = wn**2.0*m

    #Damping 2%
    c=Zeta*2*m*wn

    import numpy as np
    import geometryhelper as gh
    import StickModelHelper
    import tclhelper as th
    import os

    ########################## Initializing ##########################
    import time
    # randomnumber = random.randint(1, 10000000000)
    import uuid
    randomnumber = str(uuid.uuid4().hex.upper()[0:12])
    timestamp = time.strftime("%y%m%d-%H%M-")+'%s'%randomnumber
    ModelName = 'FindHysteresis'
    tcl = open(os.getcwd()+'//tcl//%s-%s.tcl'%(ModelName,timestamp), 'w') # Saving TCL File Location
    logfilename = '%s-%s-LOG.dat'%(ModelName,timestamp)
    OpenSees = oh.OpenSeesHelper() # Opening the OpenSees Helper Class
    OpenSeesCollector = oc.OpenSeesCollector() # Opening the OpenSees Node/Element/Material Collector
    OpenSees.Clock().write_to_tcl_start(tcl)

    ########################## Setup and Source Definition ##########################
    OpenSees.TclTitle('File Setup').write_to_tcl(tcl)
    OpenSees.Model(2, 3).open_write_to_tcl(tcl,logfilename) # Number of Dim, Number of DOFs

    ########################## Define Building Geometry, Nodes and Constraints ##########################
    OpenSees.TclTitle('Geometry Setup').write_to_tcl(tcl)

    #Defining Grids
    XGrids = np.array([0])
    YGrids = np.array([0])
    ZGrids = np.array([0])

    #Defining Floor Weights
    FloorMass = [m]

    #Creating OpenSees Nodal Elements
    ONodes = gh.get_node_points(XGrids, YGrids, ZGrids, True)
    OpenSeesCollector.assign_ONodes(ONodes)

    ########################## Define Geometric Transformations ##########################
    OpenSees.TclTitle('Define Geometric Transformations').write_to_tcl(tcl)

    ########################## Define Materials and Sections ##########################
    OpenSees.TclTitle('Define Materials and Sections').write_to_tcl(tcl)

    ########################## Define Rotational Springs for Plastic Hinge ##########################
    #Define Rotational Spring
    Spring = StickModelHelper.get_rotational_spring_elastic_plastic_ibarra_peak_oriented(OpenSeesCollector, k, Dy, PostYieldStrengthPercentage, LamdaSCA, LamdaK, rateDeterioration, R, Mu, Kappa, PostCappingStiffness)
    #get_rotational_spring_elastic_plastic_ibarra_peak_oriented
    #get_rotational_spring_elastic_plastic_kinematic_hardening

    ########################## Define Elements ##########################
    OpenSees.TclTitle('Define Elements').write_to_tcl(tcl)

    #Define P-Delta Column
    SupportNode = OpenSeesCollector.get_ONodes(gh.get_node_name(0,0,0))
    MassNode = OpenSeesCollector.create_ONode(0,0,0,0,0,0,'Constraint',True)

    SpringElement = OpenSees.ElementZeroLength(OpenSeesCollector.get_ZeroLength_Name(0),SupportNode,MassNode,[Spring],'1','Non-Linear Spring')

    OpenSeesCollector.add_OElement(SpringElement)


    ########################## Writing to TCL File ##########################

    # Setting Nodes as Used
    for object in set(OpenSeesCollector._OElements):
        object._ONode_I.setAsUsed()
        object._ONode_J.setAsUsed()

    #Writing Nodes to File
    OpenSees.TclTitle('Defining Nodes').write_to_tcl(tcl)
    for i in range(0, len(ONodes)):
        ONodes[i].write_to_tcl(tcl)

    #Defining Fixity
    OpenSees.TclTitle('Defining Fixity').write_to_tcl(tcl)
    OpenSees.Supports(SupportNode,1,0,1,0,0,1,True).write_to_tcl(tcl)
    OpenSees.Supports(MassNode,0,0,1,0,0,1,True).write_to_tcl(tcl)

    #Defining Masses
    OpenSees.TclTitle('Defining Masses').write_to_tcl(tcl)
    OpenSees.Mass(MassNode, FloorMass[0], 1e-6, 1e-6, 1e-6, 1e-6, 1e-6, True).write_to_tcl(tcl)

    #Write Element from OpenSees Collector
    OpenSees.TclTitle('Writing Materials').write_to_tcl(tcl)
    for i in range(0,len(OpenSeesCollector._OMaterials)):
        OpenSeesCollector._OMaterials[i].write_to_tcl(tcl)

    #Write Sections from OpenSees Collector
    OpenSees.TclTitle('Writing Sections').write_to_tcl(tcl)
    for i in range(0,len(OpenSeesCollector._OSections)):
        OpenSeesCollector._OSections[i].write_to_tcl(tcl)

    #Write Elements from OpenSees Collector
    OpenSees.TclTitle('Writing Elements').write_to_tcl(tcl)
    for i in range(0,len(OpenSeesCollector._OElements)):
        OpenSeesCollector._OElements[i].write_to_tcl(tcl)

    #Write Constraints from OpenSees Collector
    OpenSees.TclTitle('Writing Constraints').write_to_tcl(tcl)
    for i in range(0,len(OpenSeesCollector._OConstraints)):
        OpenSeesCollector._OConstraints[i].write_to_tcl(tcl)

    #Write Shells from OpenSees Collector
    OpenSees.TclTitle('Writing Shells').write_to_tcl(tcl)
    for i in range(0,len(OpenSeesCollector._OQuadrilateral)):
        OpenSeesCollector._OQuadrilateral[i].write_to_tcl(tcl)

    ########################## Eigenvalue Analysis ##########################
    OpenSees.TclTitle('Eigenvalue Analysis').write_to_tcl(tcl)

    # OpenSees.EigenValueAnalysis(1).write_to_tcl(tcl)

    ########################## Rayleigh Damping ##########################
    OpenSees.TclTitle('Rayleigh Damping').write_to_tcl(tcl)

    OpenSees.RayleighDamping(Zeta, wn).write_to_tcl(tcl)

    ########################## Loads ##########################
    OpenSees.TclTitle('Loads').write_to_tcl(tcl)

    ########################## Recorders ##########################
    OpenSees.TclTitle('Recorder Setup').write_to_tcl(tcl)

    Nodes = [MassNode]
    OutputFolder = 'Results'
    Displacement_File_Name = '%s-NodeD-%s.dat'%(ModelName,timestamp)
    OpenSees.NodeRecorder('Node', Displacement_File_Name, OutputFolder, Nodes, '1 2 3', 'disp').write_to_tcl(tcl)
    Velocity_File_Name = '%s-NodeV-%s.dat'%(ModelName,timestamp)
    OpenSees.NodeRecorder('Node', Velocity_File_Name, OutputFolder, Nodes, '1 2 3', 'vel').write_to_tcl(tcl)

    #Define TimeSeries
    GMData = np.array(np.array(GMData).tolist()+np.zeros(int(FreeVibrationTime/Dt)).tolist()) #Add 20 Seconds of Zeros for Free Vibration at the end
    OpenSees.TimeSeries(1,GMData,Dt, 1.0).write_to_tcl_start(tcl)
    Acceleration_File_Name = '%s-NodeA-%s.dat'%(ModelName,timestamp)
    OpenSees.NodeRecorder('Node', Acceleration_File_Name, OutputFolder, Nodes, '1', 'accel', extra='-timeSeries 1').write_to_tcl(tcl)

    #Reactions Recorder

    Nodes = [SupportNode]
    Reaction_File_Name = '%s-NodeReact-%s.dat'%(ModelName,timestamp)
    OpenSees.NodeRecorder('Node',Reaction_File_Name,OutputFolder,Nodes,'1 2 3','reaction').write_to_tcl(tcl)

    ########################## Display Results ##########################
    OpenSees.TclTitle('Display Results').write_to_tcl(tcl)

    ########################## Gravity Analysis ##########################
    OpenSees.TclTitle('Gravity Analysis').write_to_tcl(tcl)

    # OpenSees.Display().write_to_tcl(tcl)

    ########################## Time History Analysis ##########################
    OpenSees.TclTitle('Time History Analysis').write_to_tcl(tcl)

    OpenSees.DynamicAnalysis_GMData(GMData, Dt, 1.0, 0).write_to_tcl(tcl)

    ########################## Eigenvalue Analysis ##########################
    OpenSees.TclTitle('Eigenvalue Analysis').write_to_tcl(tcl)

    OpenSees.EigenValueAnalysis(1,fullGenLapack=True).write_to_tcl(tcl)

    # ########################## Pushover Analysis ##########################
    OpenSees.TclTitle('Pushover Analysis').write_to_tcl(tcl)

    NodeLoads = []
    Nodes = [MassNode]
    NodeLoads.append(OpenSees.NodeLoadsPoint(Nodes, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, True))
    PushLoadPattern = OpenSees.LoadPatternPlain('300', NodeLoads)
    PushLoadPattern.write_to_tcl(tcl)

    TestDrift = Dy*0.01

    OpenSees.Pushover_Simple(Nodes[0], '1', TestDrift, TestDrift, SupportNode).write_to_tcl(tcl)

    ########################## Close File ##########################
    OpenSees.TclTitle('Closing File').write_to_tcl(tcl)

    OpenSees.Clock().write_to_tcl_end(tcl)
    th.print_line(tcl, 'wipe al b l;')
    th.print_line(tcl, 'puts "Models Run Complete";')
    tcl.close()

    ########################## Plot Geometry ##########################

    ########################## Run OpenSees Script ##########################
    import subprocess
    #Single Processor
    OutputFile = open(os.getcwd()+'//tcl//%s-%s.out'%(ModelName,timestamp),'w')
    print(('Running OpenSees '+os.getcwd()+'//tcl//%s-%s.tcl'%(ModelName,timestamp)))
    run = subprocess.Popen([OpenSeesCommand, os.getcwd()+'//tcl//%s-%s.tcl'%(ModelName,timestamp)], cwd=os.getcwd()+'//tcl', stdout=OutputFile, stderr=subprocess.STDOUT)

    run.wait()
    OutputFile.close()
    os.remove(os.getcwd()+'//tcl//%s-%s.out'%(ModelName,timestamp))

    #Open Log File and Get NL Period
    Log = OpenSeesPostProcessor.get_string_table_from_data_file(os.getcwd()+'//tcl//',logfilename)
    ElongatedPeriod = None
    PushOverStiffness = []
    for i in range(len(Log)):
        for j in range(len(Log[i])):
            if Log[i][j]=="T1":
                ElongatedPeriod = float(Log[i][j+2])
            if Log[i][j]=='SDOFStiffness:':
                PushOverStiffness.append(float(Log[i][j+1]))

    #Find Stiffness From PushOver
    if len(PushOverStiffness)>0:
        AveragePushOverStiffness = np.mean(PushOverStiffness)
        ElongatedPeriodFromPushOver = (2.0*np.pi/AveragePushOverStiffness**.5)
    else:
        AveragePushOverStiffness = 0
        ElongatedPeriodFromPushOver = 1000000.0

    #Displau Some Results for Back Checking
    # print 'Elastic Period: %.2f'%T
    # print 'Elongated Period: %.2f'%ElongatedPeriod
    # print 'Elastic Stiffness: %.2f'%k
    # print 'NL Stiffness: %.2f'%AveragePushOverStiffness
    # if AveragePushOverStiffness != 0:
    #     print 'Calculated NL Period: %.2f'%(2.0*np.pi/AveragePushOverStiffness**.5)
    # else:
    #     print 'Calculated NL Period: infinity'

    ########################## Plot Push Over ##########################
    Displ = OpenSeesPostProcessor.get_table_from_data_file(os.getcwd()+'//tcl//'+OutputFolder, Displacement_File_Name)
    Vel = OpenSeesPostProcessor.get_table_from_data_file(os.getcwd()+'//tcl//'+OutputFolder, Velocity_File_Name)
    Acc = OpenSeesPostProcessor.get_table_from_data_file(os.getcwd()+'//tcl//'+OutputFolder, Acceleration_File_Name)
    Reac = OpenSeesPostProcessor.get_table_from_data_file(os.getcwd()+'//tcl//'+OutputFolder, Reaction_File_Name)

    DataPoints = len(GMData)

    t = Acc[0:DataPoints,0]
    ag = GMData[:]/g
    a = Acc[0:DataPoints,1]
    fs = -1*Reac[0:DataPoints,1]/(k*Dy)
    u = Displ[0:DataPoints,1]/Dy

    MaxU = max(abs(Displ[0:DataPoints,1]))
    MaxV = max(abs(Vel[0:DataPoints,1]))
    MaxA = max(abs(Acc[0:DataPoints,1]))
    MaxF = max(abs(Reac[0:DataPoints,1]))

    DArray = Displ[:,1]
    for i in range(0,len(DArray)):
        if abs(DArray[i]) == max(abs(DArray)):
            TimeAtMax = float(i)/float(len(DArray))
            break

    Stiffness = []
    Collapse = 0

    for i in range(1,len(Displ[0:DataPoints,1])):
        if (Displ[i,1]-Displ[i-1,1]) != 0.0:
            Stiffness.append(abs(Reac[i,1]-Reac[i-1,1])/abs(Displ[i,1]-Displ[i-1,1]))
        else:
            Stiffness.append(k)

    # Collapse = []
    # for i in range(1,len(Displ[:,1])):
    #     if Stiffness[i-1] < 0.01*k and abs(fs[i]) < 0.05:
    #         Collapse.append(1)
    # if len(Collapse) > 10:
    #     Collapse = 1
    # else:
    #     Collapse = 0

    PeriodThreshholdForCollapse = 4
    if ElongatedPeriodFromPushOver > PeriodThreshholdForCollapse*T:
        Collapse = 1
    else:
        Collapse = 0

    # if AveragePushOverStiffness < 0:
    #     Collapse = 0
    # else:
    #     PeriodEstimateFromPushOver = 2*np.pi/(AveragePushOverStiffness)**.5
    #     if PeriodEstimateFromPushOver > 100*T:
    #         Collapse = 1
    #     else:
    #         Collapse = 0

    if Plot == True:
        import matplotlib.pylab as plt
        from OpenSees import OpenSeesAnimation as osa

        #Another Plot
        fig, ((ax1),(ax2)) = plt.subplots(nrows=2, ncols=1, figsize=(3.5,6))
        ax1.plot(t[1:],Stiffness,'--',linewidth=2.0,label='Stiffness')
        ax2.plot(t[1:],abs(fs[1:]),'--',linewidth=2.0,label='Stiffness')
        ax1.set_title('Collaspe = %d'%Collapse)
        plt.show(block=False)

        PlotText = 'Ground Motion ID: %s\nT = %.2f \n$\zeta$ = %.2f \nR = %.0f \n$K_y$ = %.2f $K_e$ \n$K_c$ = %.2f $K_e$ \n$\lambda_{s,c,a}$ = %.2f \n$\lambda_k$ = %.2f '%(Plot_Name, T, Zeta, R, PostYieldStrengthPercentage, PostCappingStiffness, LamdaSCA, LamdaK)
        ani = osa.HysteresisAnimation(t, ag, a, fs, u, PlotText)
        if WriteAnimation:
            ani.save('Figures/Hysteresis.mp4', fps=30, bitrate=20000) #The Script already divides the frame rate by 10
            ani.ConvertToGIF(os.getcwd()+'/Figures/','Hysteresis.mp4')
        plt.show()

    os.remove(os.getcwd()+'//tcl//%s-%s.tcl'%(ModelName,timestamp))
    os.remove(os.getcwd()+'//tcl//'+OutputFolder+'//'+Displacement_File_Name)
    os.remove(os.getcwd()+'//tcl//'+OutputFolder+'//'+Velocity_File_Name)
    os.remove(os.getcwd()+'//tcl//'+OutputFolder+'//'+Acceleration_File_Name)
    os.remove(os.getcwd()+'//tcl//'+OutputFolder+'//'+Reaction_File_Name)
    os.remove(os.getcwd()+'//tcl//'+logfilename)

    # print 'EigT/PushT: %.2f'%(ElongatedPeriod/ElongatedPeriodFromPushOver)
    class Output:
        def __init__(self):
            self.MaxU = MaxU
            self.MaxV = MaxV
            self.MaxA = MaxA
            self.MaxF = MaxF
            self.TimeAtMax = TimeAtMax
            self.Collapse = Collapse
            self.ElongatedPeriodFromPushOver = ElongatedPeriodFromPushOver
            self.Acceleration = a
            self.GroundAcc = ag
            self.Displacement = u
            self.SpringForce = fs
            self.Time = t

    o = Output()

    return o

def PlotGMHysteresis_OpenSees_ElastoPlastic(GMData, Dt, T, Dy, Zeta, R, Plot_Name, Plot=False, OpenSeesCommand = os.getcwd()+'//tcl//OpenSees'):
    import numpy as np
    from OpenSees import OpenSeesHelper as oh
    from OpenSees import OpenSeesCollector as oc
    from OpenSees import OpenSeesPostProcessor

    g=386.4

    #Mass
    m = 1
    wn = 2*np.pi/T

    #Stiffness
    k = wn**2.0*m

    #Damping 2%
    c=Zeta*2*m*wn

    import numpy as np
    import geometryhelper as gh
    import StickModelHelper
    import tclhelper as th
    import os

    ########################## Initializing ##########################
    import time
    # randomnumber = random.randint(1, 10000000000)
    import uuid
    randomnumber = str(uuid.uuid4().hex.upper()[0:12])
    timestamp = time.strftime("%y%m%d-%H%M-")+'%s'%randomnumber
    ModelName = 'FindHysteresis'
    tcl = open(os.getcwd()+'//tcl//%s-%s.tcl'%(ModelName,timestamp), 'w') # Saving TCL File Location
    OpenSees = oh.OpenSeesHelper() # Opening the OpenSees Helper Class
    OpenSeesCollector = oc.OpenSeesCollector() # Opening the OpenSees Node/Element/Material Collector
    OpenSees.Clock().write_to_tcl_start(tcl)

    ########################## Setup and Source Definition ##########################
    OpenSees.TclTitle('File Setup').write_to_tcl(tcl)
    OpenSees.Model(2, 3).open_write_to_tcl(tcl) # Number of Dim, Number of DOFs

    ########################## Define Building Geometry, Nodes and Constraints ##########################
    OpenSees.TclTitle('Geometry Setup').write_to_tcl(tcl)

    #Defining Grids
    XGrids = np.array([0])
    YGrids = np.array([0])
    ZGrids = np.array([0])

    #Defining Floor Weights
    FloorMass = [m]

    #Creating OpenSees Nodal Elements
    ONodes = gh.get_node_points(XGrids, YGrids, ZGrids, True)
    OpenSeesCollector.assign_ONodes(ONodes)

    ########################## Define Geometric Transformations ##########################
    OpenSees.TclTitle('Define Geometric Transformations').write_to_tcl(tcl)

    ########################## Define Materials and Sections ##########################
    OpenSees.TclTitle('Define Materials and Sections').write_to_tcl(tcl)

    ########################## Define Rotational Springs for Plastic Hinge ##########################

    #Define Rotational Spring
    Spring = StickModelHelper.get_rotational_spring_elasto_plastic(OpenSeesCollector, k, Dy)

    ########################## Define Elements ##########################
    OpenSees.TclTitle('Define Elements').write_to_tcl(tcl)

    #Define P-Delta Column
    SupportNode = OpenSeesCollector.get_ONodes(gh.get_node_name(0,0,0))
    MassNode = OpenSeesCollector.create_ONode(0,0,0,0,0,0,'Constraint',True)

    SpringElement = OpenSees.ElementZeroLength(OpenSeesCollector.get_ZeroLength_Name(0),SupportNode,MassNode,[Spring],'1','Non-Linear Spring')

    OpenSeesCollector.add_OElement(SpringElement)


    ########################## Writing to TCL File ##########################

    # Setting Nodes as Used
    for object in set(OpenSeesCollector._OElements):
        object._ONode_I.setAsUsed()
        object._ONode_J.setAsUsed()

    #Writing Nodes to File
    OpenSees.TclTitle('Defining Nodes').write_to_tcl(tcl)
    for i in range(0, len(ONodes)):
        ONodes[i].write_to_tcl(tcl)

    #Defining Fixity
    OpenSees.TclTitle('Defining Fixity').write_to_tcl(tcl)
    OpenSees.Supports(SupportNode,1,0,1,0,0,1,True).write_to_tcl(tcl)
    OpenSees.Supports(MassNode,0,0,1,0,0,1,True).write_to_tcl(tcl)

    #Defining Masses
    OpenSees.TclTitle('Defining Masses').write_to_tcl(tcl)
    OpenSees.Mass(MassNode, FloorMass[0], 1e-6, 1e-6, 1e-6, 1e-6, 1e-6, True).write_to_tcl(tcl)

    #Write Element from OpenSees Collector
    OpenSees.TclTitle('Writing Materials').write_to_tcl(tcl)
    for i in range(0,len(OpenSeesCollector._OMaterials)):
        OpenSeesCollector._OMaterials[i].write_to_tcl(tcl)

    #Write Sections from OpenSees Collector
    OpenSees.TclTitle('Writing Sections').write_to_tcl(tcl)
    for i in range(0,len(OpenSeesCollector._OSections)):
        OpenSeesCollector._OSections[i].write_to_tcl(tcl)

    #Write Elements from OpenSees Collector
    OpenSees.TclTitle('Writing Elements').write_to_tcl(tcl)
    for i in range(0,len(OpenSeesCollector._OElements)):
        OpenSeesCollector._OElements[i].write_to_tcl(tcl)

    #Write Constraints from OpenSees Collector
    OpenSees.TclTitle('Writing Constraints').write_to_tcl(tcl)
    for i in range(0,len(OpenSeesCollector._OConstraints)):
        OpenSeesCollector._OConstraints[i].write_to_tcl(tcl)

    #Write Shells from OpenSees Collector
    OpenSees.TclTitle('Writing Shells').write_to_tcl(tcl)
    for i in range(0,len(OpenSeesCollector._OQuadrilateral)):
        OpenSeesCollector._OQuadrilateral[i].write_to_tcl(tcl)

    ########################## Eigenvalue Analysis ##########################
    OpenSees.TclTitle('Eigenvalue Analysis').write_to_tcl(tcl)

    # OpenSees.EigenValueAnalysis(1).write_to_tcl(tcl)

    ########################## Rayleigh Damping ##########################
    OpenSees.TclTitle('Rayleigh Damping').write_to_tcl(tcl)

    OpenSees.RayleighDamping(Zeta, wn).write_to_tcl(tcl)

    ########################## Loads ##########################
    OpenSees.TclTitle('Loads').write_to_tcl(tcl)

    ########################## Recorders ##########################
    OpenSees.TclTitle('Recorder Setup').write_to_tcl(tcl)

    Nodes = [MassNode]
    OutputFolder = 'Results'
    Displacement_File_Name = '%s-NodeD-%s.dat'%(ModelName,timestamp)
    OpenSees.NodeRecorder('Node', Displacement_File_Name, OutputFolder, Nodes, '1 2 3', 'disp').write_to_tcl(tcl)
    Velocity_File_Name = '%s-NodeV-%s.dat'%(ModelName,timestamp)
    OpenSees.NodeRecorder('Node', Velocity_File_Name, OutputFolder, Nodes, '1 2 3', 'vel').write_to_tcl(tcl)

    #Define TimeSeries
    OpenSees.TimeSeries(1,GMData,Dt, 1.0).write_to_tcl_start(tcl)
    Acceleration_File_Name = '%s-NodeA-%s.dat'%(ModelName,timestamp)
    OpenSees.NodeRecorder('Node', Acceleration_File_Name, OutputFolder, Nodes, '1', 'accel', extra='-timeSeries 1').write_to_tcl(tcl)

    #Reactions Recorder

    Nodes = [SupportNode]
    Reaction_File_Name = '%s-NodeReact-%s.dat'%(ModelName,timestamp)
    OpenSees.NodeRecorder('Node',Reaction_File_Name,OutputFolder,Nodes,'1 2 3','reaction').write_to_tcl(tcl)

    ########################## Display Results ##########################
    OpenSees.TclTitle('Display Results').write_to_tcl(tcl)

    ########################## Gravity Analysis ##########################
    OpenSees.TclTitle('Gravity Analysis').write_to_tcl(tcl)

    # OpenSees.Display().write_to_tcl(tcl)

    ########################## Pushover Analysis ##########################
    OpenSees.TclTitle('Pushover Analysis').write_to_tcl(tcl)

    # NodeLoads = []
    # Nodes = [MassNode]
    # NodeLoads.append(OpenSees.NodeLoadsPoint(Nodes, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, True))
    # PushLoadPattern = OpenSees.LoadPatternPlain('300', NodeLoads)
    # PushLoadPattern.write_to_tcl(tcl)
    #
    # Mc = k*Dy+k*PostYieldStrengthPercentage*(Dy*Ductility-Dy)
    #
    # MaxDrift =  Dy*Ductility+Mc/k/PostCappingStiffness
    # OpenSees.Pushover(Nodes[0], '1', MaxDrift, MaxDrift/1000 ).write_to_tcl(tcl)

    ########################## Time History Analysis ##########################
    OpenSees.TclTitle('Time History Analysis').write_to_tcl(tcl)

    OpenSees.DynamicAnalysis_GMData(GMData, Dt, 1.0, 0).write_to_tcl(tcl)

    ########################## Close File ##########################
    OpenSees.TclTitle('Closing File').write_to_tcl(tcl)

    OpenSees.Clock().write_to_tcl_end(tcl)
    th.print_line(tcl, 'wipe al b l;')
    th.print_line(tcl, 'puts "Models Run Complete";')
    tcl.close()

    ########################## Plot Geometry ##########################

    ########################## Run OpenSees Script ##########################
    import subprocess
    #Single Processor
    OutputFile = open(os.getcwd()+'//tcl//%s-%s.out'%(ModelName,timestamp),'w')
    print(('Running OpenSees '+os.getcwd()+'//tcl//%s-%s.tcl'%(ModelName,timestamp)))
    run = subprocess.Popen([OpenSeesCommand, os.getcwd()+'//tcl//%s-%s.tcl'%(ModelName,timestamp)],cwd=os.getcwd()+'//tcl', stdout=OutputFile, stderr=subprocess.STDOUT)

    run.wait()
    OutputFile.close()
    os.remove(os.getcwd()+'//tcl//%s-%s.out'%(ModelName,timestamp))

    ########################## Plot Push Over ##########################
    Displ = OpenSeesPostProcessor.get_table_from_data_file(os.getcwd()+'//tcl//'+OutputFolder, Displacement_File_Name)
    Vel = OpenSeesPostProcessor.get_table_from_data_file(os.getcwd()+'//tcl//'+OutputFolder, Velocity_File_Name)
    Acc = OpenSeesPostProcessor.get_table_from_data_file(os.getcwd()+'//tcl//'+OutputFolder, Acceleration_File_Name)
    Reac = OpenSeesPostProcessor.get_table_from_data_file(os.getcwd()+'//tcl//'+OutputFolder, Reaction_File_Name)

    t = Acc[:,0]
    ag = GMData[:]/g
    a = Acc[:,1]
    fs = -1*Reac[:,1]/(k*Dy)
    u = Displ[:,1]/Dy

    if Plot == True:
        import matplotlib.pylab as plt
        from OpenSees import OpenSeesAnimation as osa

        PlotText = 'Ground Motion ID: %s\nT = %.2f \n$\zeta$ = %.2f \nR = %.0f '%(Plot_Name, T, Zeta, R)
        ani = osa.HysteresisAnimation(t, ag, a, fs, u, PlotText)
        ani.save('Figures/Hysteresis.mp4', fps=30, bitrate=20000) #The Script already divides the frame rate by 10
        ani.ConvertToGIF(os.getcwd()+'/Figures/','Hysteresis.mp4')
        plt.show()

        # # import matplotlib.pyplot as plt
        #
        # fig, ((ax1, ax2), (ax3, ax4), (ax5, ax6)) = plt.subplots(nrows=3, ncols=2)
        # fig = plt.figure(1, figsize=(12,12))
        # ax1.set_title('Time History GM: %s'%Plot_Name)
        #
        # ax1.plot(Acc[:,0],GMData/g)
        # ax1.set_xlabel('Time, seconds')
        # ax1.set_ylabel('Ground Acc, g')
        # ax1.grid(True)
        #
        # # ax2.axis('off')
        # Periods = np.arange(0.1,5,0.1)
        # Sa = []
        # for period in Periods:
        #     u,v,a = FindSa(GMData,Dt,period,Zeta)
        #     Sa.append(a/g)
        # ax2.plot(Periods,Sa)
        # ax2.set_xlabel('Period, sec')
        # ax2.set_ylabel('$S_a, g$')
        # ax2.set_title('Elastic Response Spectrum GM: %s'%Plot_Name)
        #
        # ax3.set_title('Period: %f'%T)
        # ax3.plot(Acc[:,0],Reac[:,1]/g)
        # ax3.set_xlabel('Time, seconds')
        # ax3.set_ylabel('Base Shear/W, g')
        # ax3.grid(True)
        #
        # ax4.set_title('K: %f, Dy: %f'%(k,Dy))
        # ax4.plot(Displ[:,1]/Dy,Reac[:,1]/g)
        # ax4.set_xlabel('Displacement/Dy, in')
        # ax4.set_ylabel('Base Shear/W, g')
        # ax4.grid(True)
        #
        # ax5.axis('off')
        #
        # ax6.set_title('S_design = %f'%(Dy*k/g/m))
        # ax6.plot(Displ[:,1]/Dy,Acc[:,0])
        # ax6.set_xlabel('Displacement/Dy, in')
        # ax6.set_ylabel('Time, seconds')
        # ax6.grid(True)
        #
        # plt.tight_layout()
        # plt.savefig('Figures/HysteresisAndAll.png')
        # plt.show()

    # Temp for plotting Transparent Graphs
    # fig = plt.figure(1, figsize=(8,8))
    # plt.plot(u,fs,linewidth=2.0,color='#4F82C9')
    # plt.xlim(-15,15)
    # plt.grid()
    # plt.tick_params(axis='x',labelbottom='off')
    # plt.tick_params(axis='y',labelleft='off')
    # plt.savefig('Figures/HysteresisTrans.png', transparent=True)
    # plt.show()
    #
    # fig = plt.figure(1, figsize=(8,8))
    # plt.plot(t,ag,linewidth=2.0,color='#4F82C9')
    # # plt.xlim(-15,15)
    # plt.axis('off')
    # plt.grid(False)
    # plt.tick_params(axis='x',labelbottom='off')
    # plt.tick_params(axis='y',labelleft='off')
    # plt.savefig('Figures/HysteresisTrans.png', transparent=True)
    # plt.show()

    os.remove(os.getcwd()+'//tcl//%s-%s.tcl'%(ModelName,timestamp))
    os.remove(os.getcwd()+'//tcl//'+OutputFolder+'//'+Displacement_File_Name)
    os.remove(os.getcwd()+'//tcl//'+OutputFolder+'//'+Velocity_File_Name)
    os.remove(os.getcwd()+'//tcl//'+OutputFolder+'//'+Acceleration_File_Name)
    os.remove(os.getcwd()+'//tcl//'+OutputFolder+'//'+Reaction_File_Name)

    MaxU = max(abs(Displ[:,1]))
    MaxV = max(abs(Vel[:,1]))
    MaxA = max(abs(Acc[:,1]))
    MaxF = max(abs(Reac[:,1]))

    DArray = Displ[:,1]
    for i in range(0,len(DArray)):
        if abs(DArray[i]) == max(abs(DArray)):
            TimeAtMax = float(i)/float(len(DArray))
            break

    Stiffness = []
    Collapse = 0

    for i in range(1,len(Displ[:,1])):
        if (Displ[i,1]-Displ[i-1,1]) != 0.0:
            Stiffness.append(abs(Reac[i,1]-Reac[i-1,1])/abs(Displ[i,1]-Displ[i-1,1]))
        else:
            Stiffness.append(100000)

    for i in range(1,len(Displ[:,1])):
        if Stiffness[i-1] < 0.001 and abs(Reac[i,1]) < 0.0001:
            Collapse = 1

    return MaxU, MaxV, MaxA, MaxF, TimeAtMax, Collapse

def PlotGMHysteresis_OpenSees_Deterioration_PushOver(GMData, Dt, T, Dy, Zeta, PostYieldStrengthPercentage, LamdaSCA, LamdaK, rateDeterioration, Ductility, PostCappingStiffness, PDelta,  Plot_Name):
    import numpy as np
    from OpenSees import OpenSeesHelper as oh
    from OpenSees import OpenSeesCollector as oc
    from OpenSees import OpenSeesPostProcessor

    g=386.4

    #Mass
    m = 1
    wn = 2*np.pi/T

    #Stiffness
    k = wn**2.0*m

    #Damping 2%
    c=Zeta*2*m*wn

    import numpy as np
    import geometryhelper as gh
    import StickModelHelper
    import tclhelper as th

    ########################## Initializing ##########################
    import time
    timestamp = time.strftime("%y%m%d-%H%M-B")
    ModelName = 'FindHysteresis'
    tcl = open(os.getcwd()+'//tcl//%s-%s.tcl'%(ModelName,timestamp), 'w') # Saving TCL File Location
    OpenSees = oh.OpenSeesHelper() # Opening the OpenSees Helper Class
    OpenSeesCollector = oc.OpenSeesCollector() # Opening the OpenSees Node/Element/Material Collector
    OpenSees.Clock().write_to_tcl_start(tcl)

    ########################## Setup and Source Definition ##########################
    OpenSees.TclTitle('File Setup').write_to_tcl(tcl)
    OpenSees.Model(2, 3).open_write_to_tcl(tcl) # Number of Dim, Number of DOFs

    ########################## Define Building Geometry, Nodes and Constraints ##########################
    OpenSees.TclTitle('Geometry Setup').write_to_tcl(tcl)

    #Defining Grids
    XGrids = np.array([0])
    YGrids = np.array([0])
    ZGrids = np.array([0])

    #Defining Floor Weights
    FloorMass = [m]

    #Creating OpenSees Nodal Elements
    ONodes = gh.get_node_points(XGrids, YGrids, ZGrids, True)
    OpenSeesCollector.assign_ONodes(ONodes)

    ########################## Define Geometric Transformations ##########################
    OpenSees.TclTitle('Define Geometric Transformations').write_to_tcl(tcl)

    ########################## Define Materials and Sections ##########################
    OpenSees.TclTitle('Define Materials and Sections').write_to_tcl(tcl)

    ########################## Define Rotational Springs for Plastic Hinge ##########################

    #Define Rotational Spring
    Spring = StickModelHelper.get_rotational_spring_elastic_plastic_kinematic_hardening(OpenSeesCollector, k, Dy, PostYieldStrengthPercentage, LamdaSCA, LamdaK, rateDeterioration, Ductility, PostCappingStiffness)

    ########################## Define Elements ##########################
    OpenSees.TclTitle('Define Elements').write_to_tcl(tcl)

    #Define P-Delta Column
    SupportNode = OpenSeesCollector.get_ONodes(gh.get_node_name(0,0,0))
    MassNode = OpenSeesCollector.create_ONode(0,0,0,0,0,0,'Constraint',True)

    SpringElement = OpenSees.ElementZeroLength(OpenSeesCollector.get_ZeroLength_Name(0),SupportNode,MassNode,[Spring],'1','Non-Linear Spring')

    OpenSeesCollector.add_OElement(SpringElement)


    ########################## Writing to TCL File ##########################

    # Setting Nodes as Used
    for object in set(OpenSeesCollector._OElements):
        object._ONode_I.setAsUsed()
        object._ONode_J.setAsUsed()

    #Writing Nodes to File
    OpenSees.TclTitle('Defining Nodes').write_to_tcl(tcl)
    for i in range(0, len(ONodes)):
        ONodes[i].write_to_tcl(tcl)

    #Defining Fixity
    OpenSees.TclTitle('Defining Fixity').write_to_tcl(tcl)
    OpenSees.Supports(SupportNode,1,0,1,0,0,1,True).write_to_tcl(tcl)
    OpenSees.Supports(MassNode,0,0,1,0,0,1,True).write_to_tcl(tcl)

    #Defining Masses
    OpenSees.TclTitle('Defining Masses').write_to_tcl(tcl)
    OpenSees.Mass(MassNode, FloorMass[0], 1e-6, 1e-6, 1e-6, 1e-6, 1e-6, True).write_to_tcl(tcl)

    #Write Element from OpenSees Collector
    OpenSees.TclTitle('Writing Materials').write_to_tcl(tcl)
    for i in range(0,len(OpenSeesCollector._OMaterials)):
        OpenSeesCollector._OMaterials[i].write_to_tcl(tcl)

    #Write Sections from OpenSees Collector
    OpenSees.TclTitle('Writing Sections').write_to_tcl(tcl)
    for i in range(0,len(OpenSeesCollector._OSections)):
        OpenSeesCollector._OSections[i].write_to_tcl(tcl)

    #Write Elements from OpenSees Collector
    OpenSees.TclTitle('Writing Elements').write_to_tcl(tcl)
    for i in range(0,len(OpenSeesCollector._OElements)):
        OpenSeesCollector._OElements[i].write_to_tcl(tcl)

    #Write Constraints from OpenSees Collector
    OpenSees.TclTitle('Writing Constraints').write_to_tcl(tcl)
    for i in range(0,len(OpenSeesCollector._OConstraints)):
        OpenSeesCollector._OConstraints[i].write_to_tcl(tcl)

    #Write Shells from OpenSees Collector
    OpenSees.TclTitle('Writing Shells').write_to_tcl(tcl)
    for i in range(0,len(OpenSeesCollector._OQuadrilateral)):
        OpenSeesCollector._OQuadrilateral[i].write_to_tcl(tcl)

    ########################## Eigenvalue Analysis ##########################
    OpenSees.TclTitle('Eigenvalue Analysis').write_to_tcl(tcl)

    # OpenSees.EigenValueAnalysis(1).write_to_tcl(tcl)

    ########################## Rayleigh Damping ##########################
    OpenSees.TclTitle('Rayleigh Damping').write_to_tcl(tcl)

    OpenSees.RayleighDamping(Zeta,wn).write_to_tcl(tcl)

    ########################## Loads ##########################
    OpenSees.TclTitle('Loads').write_to_tcl(tcl)

    ########################## Recorders ##########################
    OpenSees.TclTitle('Recorder Setup').write_to_tcl(tcl)

    Nodes = [MassNode]
    OutputFolder = 'Results'
    Displacement_File_Name = '%s-NodeD-%s.dat'%(ModelName,timestamp)
    OpenSees.NodeRecorder('Node', Displacement_File_Name, OutputFolder, Nodes, '1 2 3', 'disp').write_to_tcl(tcl)
    Velocity_File_Name = '%s-NodeV-%s.dat'%(ModelName,timestamp)
    OpenSees.NodeRecorder('Node', Velocity_File_Name, OutputFolder, Nodes, '1 2 3', 'vel').write_to_tcl(tcl)

    #Define TimeSeries
    OpenSees.TimeSeries(1,GMData,Dt, 1.0).write_to_tcl_start(tcl)
    Acceleration_File_Name = '%s-NodeA-%s.dat'%(ModelName,timestamp)
    OpenSees.NodeRecorder('Node', Acceleration_File_Name, OutputFolder, Nodes, '1', 'accel', extra='-timeSeries 1').write_to_tcl(tcl)

    #Reactions Recorder

    Nodes = [SupportNode]
    Reaction_File_Name = '%s-NodeReact-%s.dat'%(ModelName,timestamp)
    OpenSees.NodeRecorder('Node',Reaction_File_Name,OutputFolder,Nodes,'1 2 3','reaction').write_to_tcl(tcl)

    ########################## Display Results ##########################
    OpenSees.TclTitle('Display Results').write_to_tcl(tcl)

    ########################## Gravity Analysis ##########################
    OpenSees.TclTitle('Gravity Analysis').write_to_tcl(tcl)

    # OpenSees.Display().write_to_tcl(tcl)

    ########################## Pushover Analysis ##########################
    OpenSees.TclTitle('Pushover Analysis').write_to_tcl(tcl)

    NodeLoads = []
    Nodes = [MassNode]
    NodeLoads.append(OpenSees.NodeLoadsPoint(Nodes, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, True))
    PushLoadPattern = OpenSees.LoadPatternPlain('300', NodeLoads)
    PushLoadPattern.write_to_tcl(tcl)

    Mc = k*Dy+k*PostYieldStrengthPercentage*(Dy*Ductility-Dy)

    MaxDrift =  Dy*Ductility+Mc/k/PostCappingStiffness
    OpenSees.Pushover(Nodes[0], '1', MaxDrift, MaxDrift/1000 ).write_to_tcl(tcl)

    ########################## Time History Analysis ##########################
    OpenSees.TclTitle('Time History Analysis').write_to_tcl(tcl)

    # OpenSees.DynamicAnalysis_GMData(GMData,Dt,1.0,0).write_to_tcl(tcl)

    ########################## Close File ##########################
    OpenSees.TclTitle('Closing File').write_to_tcl(tcl)

    OpenSees.Clock().write_to_tcl_end(tcl)
    th.print_line(tcl, 'wipe al b l;')
    th.print_line(tcl, 'puts "Models Run Complete";')
    tcl.close()

    ########################## Plot Geometry ##########################

    ########################## Run OpenSees Script ##########################
    import subprocess
    #Single Processor
    run = subprocess.Popen([os.getcwd()+'//tcl//OpenSees', os.getcwd()+'//tcl//%s-%s.tcl'%(ModelName,timestamp)],cwd=os.getcwd()+'//tcl')
    run.wait()

    #Multiple Processors
    # run = subprocess.Popen(['mpiexec','-np','1',os.getcwd()+'//tcl//OpenSeesMP', os.getcwd()+'//tcl//%s-%s.tcl'%(ModelName,timestamp)],cwd=os.getcwd()+'//tcl')
    # run.wait()

    ########################## Plot Push Over ##########################
    Displ = OpenSeesPostProcessor.get_table_from_data_file(os.getcwd()+'//tcl//'+OutputFolder, Displacement_File_Name)
    Vel = OpenSeesPostProcessor.get_table_from_data_file(os.getcwd()+'//tcl//'+OutputFolder, Velocity_File_Name)
    # Acc = OpenSeesPostProcessor.get_table_from_data_file(os.getcwd()+'//tcl//'+OutputFolder, Acceleration_File_Name)
    Reac = OpenSeesPostProcessor.get_table_from_data_file(os.getcwd()+'//tcl//'+OutputFolder, Reaction_File_Name)

    import matplotlib.pyplot as plt

    plt.plot(Displ[:,1],-Reac[:,1])
    plt.xlabel('Displacement, in')
    plt.ylabel('Fz, kips')
    plt.title('Dy: %.2f, k: %.2f, Mc: %.2f, ThetaU: %.2f'%(Dy, k, Mc, MaxDrift))

    plt.show()
    # plt.savefig('HysteresisAndAll.png')

    os.remove(os.getcwd()+'//tcl//%s-%s.tcl'%(ModelName,timestamp))
    os.remove(os.getcwd()+'//tcl//'+OutputFolder+'//'+Displacement_File_Name)
    os.remove(os.getcwd()+'//tcl//'+OutputFolder+'//'+Velocity_File_Name)
    # os.remove(os.getcwd()+'//tcl//'+OutputFolder+'//'+Acceleration_File_Name)
    os.remove(os.getcwd()+'//tcl//'+OutputFolder+'//'+Reaction_File_Name)

    MaxDuoDy = max(abs(Displ[:,1]/Dy))
    MaxSpringForce = max(abs(Reac[:,1]/g))

    DArray = Displ[:,1]
    for i in range(0,len(DArray)):
        if abs(DArray[i]) == max(abs(DArray)):
            TimeAtMax = float(i)/float(len(DArray))
            continue

    return MaxDuoDy, MaxSpringForce, TimeAtMax

def PlotGMHysteresis_Numerical_ElastoPlastic(GMData, Dt, T, Dy, Zeta, R, Plot_Name, Plot=False, SaveAnimation = False):

    def getFs(k, delta_y, u, fOLD, uOLD):
        Fmax = k*delta_y
        fnew = fOLD + k*(u-uOLD)

        if abs(fnew) > Fmax:
            return Fmax*fnew/abs(fnew), 0.0
        else:
            return fnew, k

    import numpy as np
    import os

    #Mass
    g = 386.4
    m = 1
    wn = 2*np.pi/T

    #Stiffness
    k = wn**2.0*m

    #Damping 2%
    c=Zeta*2*m*wn

    #Using Average Acceleration Method
    gamma = 0.5
    beta = 1./4.

    #How many iterations in cycle
    maxj = 1000

    #Arrays
    u = np.zeros(len(GMData))
    v = np.zeros(len(GMData))
    a = np.zeros(len(GMData))
    p = np.array(GMData)
    phat = np.zeros(len(GMData))
    fs = np.zeros(len(GMData))
    kt = np.zeros(len(GMData))
    kthat = np.zeros(len(GMData))
    Rhat = np.zeros((len(GMData),maxj))

    #Initial Calculations
    a[0]=(p[0])/m
    fs[0]=0.0
    kt[0]=k

    a1 = 1/beta/Dt**2*m+gamma/beta/Dt*c
    a2 = 1/beta/Dt*m+(gamma/beta-1)*c
    a3 = (1./2./beta-1.)*m+Dt*(gamma/2./beta-1.)*c

    #Convergence
    tol = 1*10**-5

    for i in range(1,len(GMData)):
        u[i] = u[i-1]
        fs[i]=fs[i-1]
        kt[i]=kt[i-1]
        phat[i] = p[i] + a1*u[i-1]+ a2*v[i-1] + a3*a[i-1]
        for j in range(0,maxj):
            Rhat[i,j] = phat[i] - fs[i] - a1*u[i]
            if abs(Rhat[i,j]) < tol:
                break
            kthat[i] = kt[i] + a1
            deltau = Rhat[i,j]/kthat[i]
            u[i] = u[i] + deltau
            fs[i], kt[i] = getFs(k, Dy, u[i], fs[i-1], u[i-1])
        v[i]=gamma/beta/Dt*(u[i]-u[i-1])+(1.0-gamma/beta)*v[i-1]+Dt*(1-gamma/2.0/beta)*a[i-1]
        a[i]=1/beta/(Dt**2.0)*(u[i]-u[i-1])-1.0/(beta*Dt)*v[i-1]-(1.0/2.0/beta-1.0)*a[i-1]

    ########################## Plot Push Over ##########################
    Displ = u
    Vel = v
    Acc = a
    Reac = fs

    t = np.arange(0,len(GMData),1.0)*Dt
    ag = GMData/g
    a = Acc-p
    fs = Reac/(k*Dy)
    u = Displ/Dy

    if Plot == True:
        import matplotlib.pylab as plt
        from OpenSees import OpenSeesAnimation as osa

        PlotText = 'Ground Motion ID: %s\nT = %.2f \n$\zeta$ = %.2f \nR = %.2f '%(Plot_Name, T, Zeta, R)
        ani = osa.HysteresisAnimation(t, ag, a, fs, u, PlotText)
        if SaveAnimation:
            ani.save('Figures/Hysteresis.mp4', fps=30, bitrate=20000) #The Script already divides the frame rate by 10
            ani.ConvertToGIF(os.getcwd()+'/Figures/','Hysteresis.mp4')
        plt.show()

    MaxU = max(abs(Displ))
    MaxV = max(abs(Vel))
    MaxA = max(abs(a))
    MaxF = max(abs(Reac))

    DArray = Displ
    for i in range(0,len(DArray)):
        if abs(DArray[i]) == max(abs(DArray)):
            TimeAtMax = float(i)/float(len(DArray))
            break

    Stiffness = []
    Collapse = 0

    for i in range(1,len(Displ[:])):
        if (Displ[i]-Displ[i-1]) != 0.0:
            Stiffness.append(abs(Reac[i]-Reac[i-1])/abs(Displ[i]-Displ[i-1]))
        else:
            Stiffness.append(100000)

    for i in range(1,len(Displ[:])):
        if Stiffness[i-1] < 0.001 and abs(Reac[i]) < 0.0001:
            Collapse = 1

    return MaxU, MaxV, MaxA, MaxF, TimeAtMax, Collapse

def GetSaDesignASCE(t, S1, Sds, Sd1, TL, R, I):
    Cs=min(Sds/R*I,Sd1/R*I/t,Sd1*TL/t**2/R*I)
    Csmin = max(0.044*Sds*I,0.01)
    if S1 > 0.6:
        Csmin = max(Csmin,0.5*S1/R*I)

    return max(Csmin,Cs)

def GetASCEMCE(t, S1, Sds, Sd1, TL):
    R = 1
    I = 1
    Cs=min(Sds/R*I,Sd1/R*I/t,Sd1*TL/t**2/R*I)
    Csmin = max(0.044*Sds*I,0.01)
    if S1 > 0.6:
        Csmin = max(Csmin,0.5*S1/R*I)

    return max(Csmin,Cs)

def GetMCE(t, PGA, SMS, SM1, TL):
    Ts = SM1/SMS
    T0 = 0.2*Ts

    if t < T0:
        return PGA + (SMS-PGA)/T0*(t)
    elif t < Ts:
        return SMS
    elif t < TL:
        return SM1/t
    else:
        return SM1*TL/t**2.

def ASCEResponseSpectrum(S1, Sds, Sd1, TL):
    import numpy as np
    Samin = Sds*0.4
    Tmin = 0.0
    T0 = 0.2*Sd1/Sds
    Ts = Sd1/Sds
    PeriodRange = np.arange(Ts,TL*2,0.01)
    Periods = [Tmin,T0]
    Periods.extend(PeriodRange)

    Sa = []
    Sa.append(Samin)
    Sa.append(Sds)
    for t in PeriodRange:
        Sa.append(GetSaDesignASCE(t,S1,Sds,Sd1,TL,1,1))

    return [Periods, Sa]

def ASHTOResponseSpectrum(As, T0, SDS, SD1):
    import numpy as np
    Period = np.arange(0.0, 5.01, 0.01)

    def ComputeSa(t):
        if t <= T0:
            return As + (SDS - As)/T0*t
        else:
            return min(SDS,SD1/t)

    Sa = []
    for t in Period:
        Sa.append(ComputeSa(t))

    return [Period, Sa]

def FindCollapseSa(Drift, Sa, StiffnessDifference, R):
    Gradient = []
    for i in range(1,len(Sa)-1):
        Gradient.append((Sa[i]-Sa[i-1])/(Drift[i]-Drift[0]))

    for i in range(1,len(Gradient)):
        PerDiff = abs((Gradient[i]-Gradient[i-1])/Gradient[i-1])
        if PerDiff > StiffnessDifference and Drift[i]>R:
            return Sa[i]
    return Sa[len(Drift)-1]

def FindCollapseSa_Method2(Drift, Sa, CollapseGradient, R, TimeAtMax):
    import numpy as np

    Gradient = []
    for i in range(1,len(Sa)-1):
        Gradient.append((Sa[i]-Sa[i-1])/(Drift[i]-Drift[i-1]))
    AverageGradient = np.mean(Gradient[0],Gradient[1],Gradient[2])

    SaTemp = 1000

    for i in range(2,len(Gradient)):
        if (Gradient[i]-Gradient[i-1])/Gradient[i-1] < CollapseGradient:
            SaTemp = Sa[i]
            break

    for i in range(1,len(TimeAtMax)):
        if TimeAtMax[i] > 0.75:
            SaTemp = min(Sa[i],SaTemp)
            break

    for i in range(1,len(Drift)):
        if Drift[i] > 2*R:
            SaTemp = min(Sa[i],SaTemp)
            break

    if SaTemp != 1000:
        return SaTemp

    return Sa[len(Drift)-1]

def GetFragility(SaCollapse):
    import numpy as np

    Mean = np.mean(SaCollapse)
    STD = np.std(SaCollapse)
    COV = STD/Mean

    Mu = np.log(Mean**2/(STD**2+Mean**2)**.5)
    Sigma = (np.log(1+COV**2))**.5

    Median = np.exp(Mu)

    Beta = np.std(np.log(SaCollapse))

    Sas = np.arange(0.01,max(SaCollapse),0.01)

    def Prob(x):
        return 1.0/x/(2.0*np.pi)**.5/Sigma*np.exp(-1.0*(np.log(x)-Mu)**2/2/Sigma**2)

    PDF = []
    for i in range(len(Sas)):
        PDF.append(Prob(Sas[i]))

    import scipy.integrate

    CDF = scipy.integrate.cumtrapz(PDF, Sas, initial=0)

    return Sas, CDF, Median, Beta

def GetSS(GMData, Dt, Zeta, T1, T2, dt=0.01):
    import numpy as np
    Sas = []
    T = np.arange(T1,T2+dt,dt)
    for t in T:
        SD, SV, SA = FindSa(GMData, Dt, t, Zeta)
        Sas.append(SA)

    import scipy.integrate as scipyint
    SS = scipyint.trapz(Sas, T)
    SS = SS/(Sas[0]*(T2-T1))
    return SS

def ComputeSSa(ResponseSpectrumPeriods, ResponseSpectrumSa, T1, T2):
    Sas = []
    T = []
    for i in range(len(ResponseSpectrumPeriods)):
        if T1 <= ResponseSpectrumPeriods[i] and T2 >= ResponseSpectrumPeriods[i]:
            T.append(ResponseSpectrumPeriods[i])
            Sas.append(ResponseSpectrumSa[i])

    import scipy.integrate as scipyint
    SS = scipyint.trapz(Sas, T)
    if Sas[0]*(T2-T1) == 0.:
        SS = 1.
    else:
        SS = SS/(Sas[0]*(T2-T1))

    return SS

def ComputeSaSSa(ResponseSpectrumPeriods, ResponseSpectrumSa, T1, T2):
    Sas = []
    T = []
    for i in range(len(ResponseSpectrumPeriods)):
        if T1 <= ResponseSpectrumPeriods[i] and T2 >= ResponseSpectrumPeriods[i]:
            T.append(ResponseSpectrumPeriods[i])
            Sas.append(ResponseSpectrumSa[i])

    import scipy.integrate as scipyint
    SS = scipyint.trapz(Sas, T)
    if Sas[0]*(T2-T1) == 0.:
        SS = 1.
    else:
        SS = SS/(Sas[0]*(T2-T1))

    return SS*Sas[0]

def ComputeSSaUsingGeoMean(ResponseSpectrumPeriods, ResponseSpectrumSa, T1, T2):
    Sas = []
    T = []
    for i in range(len(ResponseSpectrumPeriods)):
        if T1 <= ResponseSpectrumPeriods[i] and T2 >= ResponseSpectrumPeriods[i]:
            T.append(ResponseSpectrumPeriods[i])
            Sas.append(ResponseSpectrumSa[i])

    import scipy.stats.mstats as sm

    SS = sm.gmean(Sas)

    SS = SS/Sas[0]

    return SS

def ComputeSSaVer2(ResponseSpectrumPeriods, ResponseSpectrumSa, T1, T2):

    import numpy as np
    dt = ResponseSpectrumPeriods[1]-ResponseSpectrumPeriods[0]
    T2 = min(np.max(ResponseSpectrumPeriods),T2)

    T = np.arange(T1,T2+dt,dt)
    Sas = np.interp(T, ResponseSpectrumPeriods, ResponseSpectrumSa)

    import scipy.integrate as scipyint
    SS = scipyint.trapz(Sas, T)
    SS = SS/(Sas[0]*(T2-T1))

    return SS

def ComputeSSaFromGMid(GMFolder, GMid, T1, T2):
    O = GetResponseSpectrumDataPoints(GMFolder, GMid)
    ResponseSpectrumPeriods = O.Period
    ResponseSpectrumSa = O.Sa
    Sas = []
    T = []
    for i in range(len(ResponseSpectrumPeriods)):
        if T1 <= ResponseSpectrumPeriods[i] and T2 >= ResponseSpectrumPeriods[i]:
            T.append(ResponseSpectrumPeriods[i])
            Sas.append(ResponseSpectrumSa[i])

    import scipy.integrate as scipyint
    SS = scipyint.trapz(Sas, T)
    SS = SS/(Sas[0]*(T2-T1))

    return SS

def ComputeSSaFromGMidQuick(GMFolder, GMid, T1, T2):
    O = GetResponseSpectrumDataPointsQuick(GMFolder, GMid)
    ResponseSpectrumPeriods = O.Period
    ResponseSpectrumSa = O.Sa

    import numpy as np
    dt = ResponseSpectrumPeriods[1] - ResponseSpectrumPeriods[0]
    T2 = min(np.max(ResponseSpectrumPeriods), T2)

    T = np.arange(T1, T2 + dt, dt)
    Sas = np.interp(T, ResponseSpectrumPeriods, ResponseSpectrumSa)

    import scipy.integrate as scipyint
    SS = scipyint.trapz(Sas, T)
    SS = SS / (Sas[0] * (T2 - T1))

    return SS

def ComputeSSaFromGMDataQuick(GMData, Dt, T1, T2, SaDt = 0.01):
    import numpy as np
    ResponseSpectrumPeriods = np.arange(T1, T2+SaDt, SaDt)

    ResponseSpectrumSa = []
    Zeta = 0.05
    for t1 in ResponseSpectrumPeriods:
        SD, SV, SA = FindSa(GMData, Dt, t1, Zeta)
        ResponseSpectrumSa.append(SA)

    T = ResponseSpectrumPeriods
    Sas = ResponseSpectrumSa

    import scipy.integrate as scipyint
    SS = scipyint.trapz(Sas, T)
    SS = SS / (Sas[0] * (T2 - T1))

    return SS

def GetSSThroughSprectrum(Folder_Location, GMid, T1, T2):
    O = GetResponseSpectrumDataPoints(Folder_Location,GMid)
    Period = O.Period
    Acc = O.Sa
    Vel = O.Sv
    Disp = O.Sd

    PeriodTemp = []
    AccTemp = []

    if T1 == T2:
        return 0.0

    for i in range(len(Period)):
        if Period[i] <= T2 and Period[i] >= T1:
            PeriodTemp.append(Period[i])
            AccTemp.append(Acc[i])

    import scipy.integrate as scipyint
    SS = scipyint.trapz(AccTemp, PeriodTemp)
    SS = SS/(AccTemp[0]*(T2-T1))

    return SS

def MethodNameGotDeleted(Folder_Location, GMid, T1, T2, Tnorm):
    O = GetResponseSpectrumDataPoints(Folder_Location,GMid)
    Period = O.Period
    Acc = O.Sa
    Vel = O.Sv
    Disp = O.Sd

    PeriodTemp = []
    AccTemp = []

    if T1 == T2:
        return 0.0

    import numpy as np

    SaNorm = np.interp(Tnorm,Period,Acc)

    for i in range(len(Period)):
        if Period[i] <= T2 and Period[i] >= T1:
            PeriodTemp.append(Period[i])
            AccTemp.append(Acc[i])



    import scipy.integrate as scipyint
    SS = scipyint.trapz(AccTemp, PeriodTemp)
    SS = SS/(SaNorm*(T2-T1))

    return SS

def GetCOVofSpectrumThroughSprectrum(Folder_Location, GMid, T1, T2):
    import numpy as np
    O = GetResponseSpectrumDataPoints(Folder_Location,GMid)
    Period = O.Period
    Acc = O.Sa
    Vel = O.Sv
    Disp = O.Sd

    PeriodTemp = []
    AccTemp = []
    VelTemp = []
    DispTemp = []

    if T1 == T2:
        return 0.0

    for i in range(len(Period)):
        if Period[i] <= T2 and Period[i] >= T1:
            PeriodTemp.append(Period[i])
            AccTemp.append(Acc[i])
            VelTemp.append(Vel[i])
            DispTemp.append(Disp[i])

    class Output:
        def __init__(self):
            self.COV_Sa = np.std(AccTemp)/np.mean(AccTemp)
            self.COV_Sv = np.std(VelTemp)/np.mean(VelTemp)
            self.COV_Sd = np.std(DispTemp)/np.mean(DispTemp)

    Out = Output()

    return Out

def GetSSLogLogThroughSprectrum(Folder_Location, GMid, T1, T2):
    import numpy as np
    O = GetResponseSpectrumDataPoints(Folder_Location,GMid)
    Period = O.Period
    Acc = O.Sa
    Vel = O.Sv
    Disp = O.Sd

    PeriodTemp = []
    AccTemp = []

    if T1 == T2:
        return 0.0

    for i in range(len(Period)):
        if Period[i] <= T2 and Period[i] >= T1:
            PeriodTemp.append(Period[i])
            AccTemp.append(Acc[i])

    PeriodTemp = np.log(PeriodTemp)
    AccTemp = np.log(AccTemp)

    import scipy.integrate as scipyint
    SS = scipyint.trapz(AccTemp, PeriodTemp)
    SS = SS/(AccTemp[0]*(T2-T1))

    return SS

def GetSSDSemiLogThroughSprectrum(Folder_Location, GMid, T1, T2):
    import numpy as np
    O = GetResponseSpectrumDataPoints(Folder_Location,GMid)
    Period = O.Period
    Acc = O.Sa
    Vel = O.Sv
    Disp = O.Sd

    PeriodTemp = []
    DispTemp = []

    if T1 == T2:
        return 0.0

    for i in range(len(Period)):
        if Period[i] <= T2 and Period[i] >= T1:
            PeriodTemp.append(Period[i])
            DispTemp.append(Disp[i])

    PeriodTemp = np.log(PeriodTemp)

    import scipy.integrate as scipyint
    SS = scipyint.trapz(DispTemp, PeriodTemp)
    SS = SS/(DispTemp[0]*(T2-T1))

    return SS

def GetSSVThroughVelocitySprectrum(Folder_Location, GMid, T1, T2):
    O = GetResponseSpectrumDataPoints(Folder_Location,GMid)
    Period = O.Period
    Acc = O.Sa
    Vel = O.Sv
    Disp = O.Sd

    PeriodTemp = []
    VelTemp = []

    if T1 == T2:
        return 0.0

    for i in range(len(Period)):
        if Period[i] <= T2 and Period[i] >= T1:
            PeriodTemp.append(Period[i])
            VelTemp.append(Vel[i])

    import scipy.integrate as scipyint
    SS = scipyint.trapz(VelTemp, PeriodTemp)
    SS = SS/(VelTemp[0]*(T2-T1))

    return SS

def GetSSDThroughDisplacementSprectrum(Folder_Location, GMid, T1, T2):
    O = GetResponseSpectrumDataPoints(Folder_Location,GMid)
    Period = O.Period
    Acc = O.Sa
    Vel = O.Sv
    Disp = O.Sd

    PeriodTemp = []
    DispTemp = []

    if T1 == T2:
        return 0.0

    for i in range(len(Period)):
        if Period[i] <= T2 and Period[i] >= T1:
            PeriodTemp.append(Period[i])
            DispTemp.append(Disp[i])

    import scipy.integrate as scipyint
    SS = scipyint.trapz(DispTemp, PeriodTemp)
    SS = SS/(DispTemp[0]*(T2-T1))

    return SS

def GetSaThroughSpectrumFile(Folder_Location, GMid, T1):
    import numpy as np
    Period = GetResponseSpectrumDataPoints(Folder_Location,GMid).Period
    Acc = GetResponseSpectrumDataPoints(Folder_Location,GMid).Sa

    Sa = np.interp(T1, Period, Acc)

    return Sa

def GetMedian(Data):
    import numpy as np

    Mean = np.mean(Data)
    STD = np.std(Data)
    COV = STD/Mean

    Mu = np.log(Mean**2/(STD**2+Mean**2)**.5)
    Sigma = (np.log(1+COV**2))**.5

    Median = np.exp(Mu)

    Beta = np.std(np.log(Data))

    return Median, Beta

def GetSaRAlpha(GMData, Dt, Zeta, T1, Alpha=0.5):
    Sas = []
    T2 = 2*T1
    T = [T1,T2]
    for t in T:
        SD, SV, SA = FindSa(GMData, Dt, t, Zeta)
        Sas.append(SA)

    R = (Sas[1]/Sas[0])**Alpha

    return R*Sas[0]

def GetSstar(GM_Folder, GMid, Zeta, T1, Alpha = 0.5):
    Sas = []
    T2 = 2*T1
    T = [T1,T2]
    for t in T:
        SA = GetSaThroughSpectrumFile(GM_Folder, GMid, t)
        Sas.append(SA)

    R = (Sas[1]/Sas[0])**Alpha

    return R*Sas[0]

def GetEpsilon(GM_Folder, GMid, Zeta, T1):
    import numpy as np
    import GMPEs

    import GMHelper
    GMids, GMFiles, Dt, NumPoints, GMData = GMHelper.ExtractGroundMotion(GM_Folder)

    # Parameters = GetGMParameters(GM_Folder,GMids)
    # GMidsPara = Parameters[0]
    # M = Parameters[1]
    # R = Parameters[2]
    #
    # index = GMidsPara.index(GMid)

    O = GetGMParameters2(GM_Folder)

    SD, SV, SA = FindSa(GMData[GMid], Dt[GMid], T1, Zeta)

    SaGMPE, SigmaGMPE = GMPEs.AbrahamsonSilva_1997_horiz(O.M[GMid],O.R[GMid],T1)

    epsilon = (np.log(SA)-np.log(SaGMPE))/SigmaGMPE

    return epsilon

def GetGMParameters(GM_Location, GMids):
    from OpenSees import OpenSeesPostProcessor

    Parameters = OpenSeesPostProcessor.get_string_table_from_data_file(GM_Location,'/GMParameters.txt')
    # ID, M, R

    M = []
    R = []

    for GMid in GMids:
        for i in range(len(Parameters)):
            if GMid[0:5] == Parameters[i][0]:
                M.append(float(Parameters[i][1]))
                R.append(float(Parameters[i][2]))

    return [GMids,M,R]

def GetGMParameters2(GM_Location):
    from OpenSees import OpenSeesPostProcessor

    Parameters = OpenSeesPostProcessor.get_string_table_from_data_file(GM_Location,'/GMParameters.txt')
    # ID, M, R

    M = {}
    R = {}

    for i in range(len(Parameters)):
        M[Parameters[i][0]]=float(Parameters[i][1])
        R[Parameters[i][0]]=float(Parameters[i][2])

    class Output:
        def __init__(self):
            self.M = M
            self.R = R

    O = Output()
    return O

def GetEpsilons(GM_Folder, Period, Zeta):
    import GMHelper
    GMids, GMFiles, Dt, NumPoints, GMData = GMHelper.ExtractGroundMotion(GM_Folder)

    Parameters = GetGMParameters(GM_Folder,GMids)
    GMidsPara = Parameters[0]
    M = Parameters[1]
    R = Parameters[2]

    Epsilon = {}
    for GMid in GMids:
        ParaIndex = GMidsPara.index(GMid)
        Epsilon[GMid]=GetEpsilon(GM_Folder,GMid,Zeta,Period)

    return Epsilon

def GetSaRAlphas(GM_Folder, GMid, Period, Zeta, Alpha=0.5):


    import GMHelper
    GMids, GMFiles, Dt, NumPoints, GMData = GMHelper.ExtractGroundMotion(GM_Folder)

    SaRAlphas = {}
    SaRAlphas=GetSaRAlpha(GMData[GMid], Dt[GMid], Zeta, Period, Alpha)

    return SaRAlphas

def CreateGroundMotionFiles(Time, AccInG, GMName, Folder_Location, CreateResponseSpectrum=True):
    import OpenSees.OpenSeesPostProcessor as opp
    #Find Dt
    Dt = Time[1]-Time[0]
    for i in range(1,len(Time)):
        tempDt = Time[i]-Time[i-1]
        if round(Dt,3) != round(tempDt,3):
            print(('Error in creating Dt, uneven time steps at point %d'%i))
            return None
    opp.save_table_to_file(Folder_Location, 'DtFile_(%s).dat'%GMName, [[round(Dt,3)]])
    #FindNumber of Points
    if len(Time) == len(AccInG):
        opp.save_table_to_file(Folder_Location, 'NumPointsFile_(%s).dat'%GMName, [[len(Time)]])
    #Write
    opp.save_table_to_file(Folder_Location, 'SortedEQFile_(%s).dat'%GMName, [AccInG])

    #Create Response Spectrum
    if CreateResponseSpectrum:
        CreateResponseSpectrumFile(Folder_Location,GMName, AccInG, round(Dt,3), 0.05, Tmax=20)

# @jit
def CreateResponseSpectrumFiles(Folder_Location):
    GMids, GMFiles, Dt, NumPoints, GMData = ExtractGroundMotion(Folder_Location)

    for GMid in GMids:
        print(('Running GMid: %s'%GMid))
        import os
        if not(os.path.isfile(Folder_Location+'/ResponseSpectrum_(%s).dat'%GMid)):
            CreateResponseSpectrumFile(Folder_Location, GMid, GMData[GMid], round(Dt[GMid],3), 0.05)

def CreateResponseSpectrumFromGMid(Folder_Location, GMid):
    Dt, NumPoints, GMData = ExtractOneGroundMotion(Folder_Location, GMid)

    CreateResponseSpectrumFile(Folder_Location, GMid, GMData, np.round(Dt,3), 0.05)

def GetResponseSpectrumDataPoints(Folder_Location, GMid):
    Table = open(Folder_Location+'ResponseSpectrum_(%s).dat'%GMid).readlines()
    Period = []
    Acceleration = []
    Displacement = []
    Velocity = []
    g = 386.4

    for i in range(5,len(Table)):
        row = Table[i].split()
        Period.append(float(row[0]))
        Displacement.append(float(row[1]))
        Velocity.append(float(row[2]))
        Acceleration.append(float(row[3])/g)

    class Output:
        def __init__(self):
            self.Sa = Acceleration
            self.Sv = Velocity
            self.Sd = Displacement
            self.Period = Period

    O = Output()

    return O

def GetResponseSpectrumDataPointsQuick(Folder_Location, GMid):
    import numpy as np
    lines = open(Folder_Location+'ResponseSpectrum_(%s).dat'%GMid, 'r').readlines()[5:]
    lines = [x.split() for x in lines]
    Table = np.array(lines, dtype=float)
    g = 386.4

    class Output:
        def __init__(self):
            self.Sa = Table[:,3]/g
            self.Sv = Table[:,2]
            self.Sd = Table[:,1]
            self.Period = Table[:,0]

    O = Output()

    return O

def FindRyForGivenDuctilityDemand(GMid, GMData, Dt, T, Zeta, DuctilityDemand, StartingR=1.01, ShowPlot=False):
    import numpy as np

    # SDOF Properties
    m = 1
    wn = 2*np.pi/T
    #Stiffness
    k = wn**2.0*m

    #Find Elastic Sa
    g = 386.4
    u, v, Sa = FindSa(GMData*g, Dt, T, Zeta)

    #Define Function
    def function(_R):
        Fy = Sa*m/_R
        Dy = Fy/k

        uI, vI, aI, FI, TatMax, Collapse = PlotGMHysteresis_Numerical_ElastoPlastic(GMData*g, Dt, T, Dy, Zeta, _R, '%s'%GMid)
        # uI, vI, aI, FI, TatMax, Collapse = PlotGMHysteresis_OpenSees_ElastoPlastic(GMData[GMid]*g*SF, Dt[GMid], T, Dy, Zeta, R,'%s'%GMid, OpenSeesCommand=OpenSeesCommand)

        DD = uI/Dy
        return DD

    #Try to Estimate A good Starting Point for R
    MuEstimate = []
    REstimate = np.arange(1,int(DuctilityDemand)*3,0.5)
    for i in np.arange(1,int(DuctilityDemand)*3,0.5):
        MuEstimate.append(function(i))

    #Try Preliminary Ry
    R = [np.interp(DuctilityDemand,MuEstimate,REstimate)]
    # print 'Tryiing %.2f'%R[0]

    # import matplotlib.pylab as plt
    # fig, (ax1) = plt.subplots(nrows=1, ncols=1)
    # ax1.plot(REstimate,MuEstimate,linewidth=2.0)

    from numpy import array
    from scipy.optimize import leastsq

    def residual_function(R):
        return function(R) - DuctilityDemand

    Coeff = leastsq(residual_function,  array(R), args=(), ftol=0.1, xtol=0.1, maxfev=2000)
    Coeff = Coeff[0]
    # print Coeff

    Fy = Sa*m/Coeff[0]
    Dy = Fy/k
    uI, vI, aI, FI, TatMax, Collapse = PlotGMHysteresis_Numerical_ElastoPlastic(GMData*g, Dt, T, Dy, Zeta, Coeff, '%s'%GMid, False)
    DD = uI/Dy
    if round(DD,1)!=DuctilityDemand:
        Coeff[0]=R[0]

    #Show Plot To Verify Results
    if ShowPlot:
        Fy = Sa*m/Coeff[0]
        Dy = Fy/k
        uI, vI, aI, FI, TatMax, Collapse = PlotGMHysteresis_Numerical_ElastoPlastic(GMData*g, Dt, T, Dy, Zeta, Coeff, '%s'%GMid, True)
    else:
        Fy = Sa*m/Coeff[0]
        Dy = Fy/k
        uI, vI, aI, FI, TatMax, Collapse = PlotGMHysteresis_Numerical_ElastoPlastic(GMData*g, Dt, T, Dy, Zeta, Coeff, '%s'%GMid, False)
        DD = uI/Dy

    return Coeff[0], DD

def FindRyForGivenDuctilityDemandThroughInterpolation(GMid, GMData, Dt, T, Zeta, DuctilityDemand, PrintResult=False):
    import numpy as np

    # SDOF Properties
    m = 1
    wn = 2*np.pi/T
    #Stiffness
    k = wn**2.0*m

    #Find Elastic Sa
    g = 386.4
    u, v, Sa = FindSa(GMData*g, Dt, T, Zeta)

    #Define Function
    def function(_R):
        Fy = Sa*m/_R
        Dy = Fy/k

        uI, vI, aI, FI, TatMax, Collapse = PlotGMHysteresis_Numerical_ElastoPlastic(GMData*g, Dt, T, Dy, Zeta, _R, '%s'%GMid)

        DD = uI/Dy

        return DD

    #Try to Estimate A good Starting Point for R
    MuEstimate = []
    REstimate = np.arange(1,DuctilityDemand*10,0.1)
    for i in REstimate:
        DDtemp = function(i)
        MuEstimate.append(DDtemp)
        if PrintResult:
            print(('Ran R: %.2f, Mu: %.2f'%(i, DDtemp)))
        if DDtemp > DuctilityDemand:
            break

    #Try Preliminary Ry
    R = [np.interp(DuctilityDemand,MuEstimate,REstimate[0:len(MuEstimate)])]

    return R[0], DuctilityDemand

def GetGMDurations(Folder_Location):
    Ds = {}

    GMids, GMFiles, Dt, NumPoints, GMData = ExtractGroundMotion(Folder_Location)
    for GMid in GMids:
        Ds[GMid]=FindDs(GMData[GMid], Dt[GMid])

    return Ds

def GetResponseSpectrumSd(Folder_Location):
    Sd = {}
    Tn = {}
    for subdir, dirs, files in os.walk(Folder_Location):
        for file in files:
            if file.startswith('ResponseSpectrum_(')==False:
                continue
            tempSd = []
            tempTn = []
            data = open(Folder_Location+'/'+file,'r')
            lines = data.readlines()
            GMid = lines[0].split()[1]
            for i in range(5,len(lines)):
                tempSd.append(float(lines[i].split()[1]))
                tempTn.append(float(lines[i].split()[0]))

            Sd[GMid]=tempSd
            Tn[GMid]=tempTn
    return Tn, Sd

def GetSacOEta(GMid, GMData, Dt, T, Zeta, mu, PY, PC, kappa, Lambda, c, OpenSeesCommand = os.getcwd()+'//tcl//OpenSees', ShowPlot=False, PlotIDA = False):
    import numpy as np
    g = 386.4

    #Mass
    mass = 1.0
    wn = 2*np.pi/T
    #Stiffness
    k = wn**2.0*mass

    u, v, Sa = FindSa(GMData*g, Dt, T, Zeta)
    SFAdjustment = Sa/g
    Fy = SFAdjustment*mass*g
    Dy = Fy/k
    R = Sa/(Fy/(mass*g))

    #Set Initial Accuracy of Scale Factors, reduce if small datapoints is recieved
    Accuracy = 0.10
    SF = [0.8,0.9]+np.arange(1.0, 30.0, Accuracy).tolist()

    Sas = []
    DuoDy = []
    TimeAtMax = []
    Collapse = []
    NLPeriod = []

    def run(_SF):

        for i in range(0,len(_SF)):
            # uI, vI, aI, FI, TatMax, collapse, ElongatedPeriod = PlotGMHysteresis_OpenSees_Deterioration(GMData*g*_SF[i], Dt, T, Dy, Zeta, PY, Lambda, Lambda*2.0, c, R, mu, kappa, PC, '%s'%GMid, OpenSeesCommand=OpenSeesCommand, Plot=ShowPlot)
            try:
                uI, vI, aI, FI, TatMax, collapse, ElongatedPeriod = PlotGMHysteresis_OpenSees_Deterioration_Peak_Oriented(GMData*g*_SF[i], Dt, T, Dy, Zeta, PY, Lambda, Lambda*2.0, c, R, mu, kappa, PC, '%s'%GMid, OpenSeesCommand=OpenSeesCommand, Plot=ShowPlot)
            except:
                # uI = 1e10
                # vI = 1e10
                # aI = 1e10
                # FI = 1e10
                # TatMax = 1.0
                # collapse = 1
                # ElongatedPeriod = 1e10
                return None, None, None, None, None

            Sas.append(Sa*_SF[i]/g/(Fy/mass/g))
            DuoDy.append(uI/Dy)
            TimeAtMax.append(TatMax)
            Collapse.append(collapse)
            NLPeriod.append(ElongatedPeriod)

            if PlotIDA:
                ind = len(Sas)-1
                print(('Latest Results: %f %f %f %d %f'%(Sas[ind], DuoDy[ind], TimeAtMax[ind], Collapse[ind], NLPeriod[ind])))

            if Collapse[i-1] == 1:
                break

    while len(Sas) < 10:
        Sas = []
        DuoDy = []
        TimeAtMax = []
        Collapse = []
        NLPeriod = []

        run(SF)

        Accuracy = Accuracy/2.0
        SF = [0.8,0.9]+np.arange(1.0, 30.0, Accuracy).tolist()

    Sas = [0] + Sas
    DuoDy = [0] + DuoDy
    TimeAtMax = [0] + TimeAtMax
    Collapse = [0] + Collapse
    NLPeriod = [T] + NLPeriod

    if PlotIDA:
        import matplotlib.pylab as plt

        fig, (ax1) = plt.subplots(nrows=1, ncols=1)
        import matplotlib.pylab as plt
        ax1.plot(DuoDy,Sas,linewidth=2.0)
        ax1.set_xlabel('$\mu$')
        ax1.set_ylabel('$S_a$')
        for i in range(0,len(Sas)):
            if Collapse[i]==1:
                ax1.scatter(DuoDy[i],Sas[i])

        plt.show()

    return Sas, DuoDy, TimeAtMax, Collapse, NLPeriod

def ReadLatLong(GM_Folder):

    GMids, GMFiles, Dt, NumPoints, GMData = ExtractGroundMotion(GM_Folder)

    import OpenSees.OpenSeesPostProcessor as OPP

    LatLonData = OPP.get_string_table_from_data_file(GM_Folder,'LATLON.dat')

    LatDic = {}
    LongDic = {}

    for GMid in GMids:
        #Get 3 or 4 Letter Location Code
        if GMid[0:2]=='00':
            Loc = GMid[2:5]
        else:
            Loc = GMid[1:5]

        #Extract GPS Location
        for i in range(len(LatLonData)):
            if LatLonData[i][0]==Loc:
                LatDic[GMid]=float(LatLonData[i][1])
                LongDic[GMid]=float(LatLonData[i][2])

    return LatDic, LongDic

def GetSAtogram_OLD(GMData, Dt, Zeta, minT=0.1, maxT=10, PeriodDivisions=100, TimeDivisions=10):
    import numpy as np
    # Periods = np.arange(minT,maxT+dT,dT)
    Periods = np.linspace(minT,maxT,PeriodDivisions)
    SaArray = np.zeros((TimeDivisions, len(Periods)))
    GMDataArray = []
    Split = np.linspace(0,len(GMData),TimeDivisions+1)
    for i in range(1,TimeDivisions):
        GMDataArray.append(GMData[0:Split[i]])
    for i in range(len(GMDataArray)):
        for j in range(len(Periods)):
            u, v, a = FindSa(GMDataArray[i], Dt, Periods[j], Zeta)
            SaArray[i,j]=a
    SaDiffArray = np.zeros((TimeDivisions,len(Periods)))

    SaDiffArray[0,:] = np.zeros(len(Periods))
    SaDiffArray[1,:] = SaArray[0,:]
    for i in range(1,len(SaArray)-1):
        SaDiffArray[i+1,:] = SaArray[i,:]-SaArray[i-1,:]

    return [Split*Dt, Periods, SaDiffArray]

def GetSAtogram(GM_Folder, GMid, Zeta, minT=0.1, maxT=10, PeriodDivisions=100, TimeDivisions=10):
    import numpy as np
    import GMHelper
    GMids, GMFiles, Dt, NumPoints, GMData = GMHelper.ExtractGroundMotion(GM_Folder)

    Periods = np.linspace(minT,maxT,PeriodDivisions)
    Duration = np.linspace(0,len(GMData[GMid])*Dt[GMid],TimeDivisions)

    SaArray = []
    for t in Periods:
        a, time = GetSDOFAcc(GM_Folder, GMid, Zeta, t, TimeDivisions)
        SaArray.append(a)
    SaArray = np.array(SaArray)

    return [Duration, Periods, SaArray]

def GetInp(GM_Folder, GMid, Zeta, T1, T2=None):
    import numpy as np
    if T2 == None:
        T2 = 2*T1
    T = np.arange(T1,T2+0.01,0.01)
    Sas = []
    for t in T:
        SA = GetSaThroughSpectrumFile(GM_Folder, GMid, t)
        Sas.append(SA)

    import scipy.stats.mstats as sm
    Saavg = sm.gmean(Sas)

    Inp = Sas[0]*(Saavg/Sas[0])**.4

    return Inp

def GetASAR(GM_Folder, GMid, Zeta, T1, R=0.4):
    import numpy as np
    Xf = 1.-R/100.
    T2 = T1/Xf
    T = np.arange(T1,T2+0.01,0.01)
    Sas = []
    for t in T:
        SA = GetSaThroughSpectrumFile(GM_Folder, GMid, t)
        Sas.append(SA)

    SasOT2 = np.array(Sas)/(np.array(T)**2)
    ASAR = T1/(1-Xf)*np.trapz(SasOT2,np.array(T))

    return ASAR

def GetSa(GM_Folder, GMid, Zeta, T):
    import GMHelper
    GMids, GMFiles, Dt, NumPoints, GMData = GMHelper.ExtractGroundMotion(GM_Folder)

    import numpy as np

    #Using Newmarks: Linear System
    #Using Linear Acceleration Method
    #From Chopra, Page 177

    #Mass
    m = 1
    wn = 2*np.pi/T

    #Stiffness
    k = wn**2.0*m

    #Damping 2%
    c=Zeta*2*m*wn

    gamma = 0.5
    beta = 1./4.

    u = np.zeros(len(GMData[GMid]))
    v = np.zeros(len(GMData[GMid]))
    a = np.zeros(len(GMData[GMid]))
    p = np.array(GMData[GMid])

    a[0]=p[0]/m

    a1 = 1/beta/Dt[GMid]**2*m+gamma/beta/Dt[GMid]*c
    a2 = 1/beta/Dt[GMid]*m+(gamma/beta-1)*c
    a3 = (1./2./beta-1.)*m+Dt[GMid]*(gamma/2./beta-1.)*c

    khat = k + a1

    for i in range(1,len(GMData[GMid])):
        phat = p[i] + a1*u[i-1] + a2*v[i-1] + a3*a[i-1]
        u[i] = phat/khat
        v[i] = gamma/beta/Dt[GMid]*(u[i]-u[i-1]) + (1-gamma/beta)*v[i-1]+Dt[GMid]*(1.-gamma/2./beta)*a[i-1]
        a[i] = 1./beta/Dt[GMid]**2.0*(u[i]-u[i-1])-1/beta/Dt[GMid]*v[i-1]-(1/beta/2-1)*a[i-1]

    a = a-p

    return max(abs(max(u)),abs(min(u))),max(abs(max(v)),abs(min(v))),max(abs(max(a)),abs(min(a)))

def GetSDOFAcc(GM_Folder, GMid, Zeta, T, TimeDivisons=10):
    import GMHelper
    GMids, GMFiles, Dt, NumPoints, GMData = GMHelper.ExtractGroundMotion(GM_Folder)

    import numpy as np

    #Using Newmarks: Linear System
    #Using Linear Acceleration Method
    #From Chopra, Page 177

    #Mass
    m = 1
    wn = 2*np.pi/T

    #Stiffness
    k = wn**2.0*m

    #Damping 2%
    c=Zeta*2*m*wn

    gamma = 0.5
    beta = 1./4.

    u = np.zeros(len(GMData[GMid]))
    v = np.zeros(len(GMData[GMid]))
    a = np.zeros(len(GMData[GMid]))
    p = np.array(GMData[GMid])

    a[0]=p[0]/m

    a1 = 1/beta/Dt[GMid]**2*m+gamma/beta/Dt[GMid]*c
    a2 = 1/beta/Dt[GMid]*m+(gamma/beta-1)*c
    a3 = (1./2./beta-1.)*m+Dt[GMid]*(gamma/2./beta-1.)*c

    khat = k + a1

    for i in range(1,len(GMData[GMid])):
        phat = p[i] + a1*u[i-1] + a2*v[i-1] + a3*a[i-1]
        u[i] = phat/khat
        v[i] = gamma/beta/Dt[GMid]*(u[i]-u[i-1]) + (1-gamma/beta)*v[i-1]+Dt[GMid]*(1.-gamma/2./beta)*a[i-1]
        a[i] = 1./beta/Dt[GMid]**2.0*(u[i]-u[i-1])-1/beta/Dt[GMid]*v[i-1]-(1/beta/2-1)*a[i-1]

    a = a-p

    MaxA = []
    Increment = int(len(a)/TimeDivisons)
    for i in range(1,TimeDivisons+1,1):
        MaxA.append(max(abs(a[(i-1)*Increment:i*Increment])))
    t = Dt[GMid]*np.linspace(0,len(a),TimeDivisons)
    return MaxA, t

def GetAbsAcc(GM_Folder, GMid, Zeta, T):
    import GMHelper
    GMids, GMFiles, Dt, NumPoints, GMData = GMHelper.ExtractGroundMotion(GM_Folder)

    import numpy as np

    #Using Newmarks: Linear System
    #Using Linear Acceleration Method
    #From Chopra, Page 177

    #Mass
    m = 1
    wn = 2*np.pi/T

    #Stiffness
    k = wn**2.0*m

    #Damping 2%
    c=Zeta*2*m*wn

    gamma = 0.5
    beta = 1./4.

    u = np.zeros(len(GMData[GMid]))
    v = np.zeros(len(GMData[GMid]))
    a = np.zeros(len(GMData[GMid]))
    p = np.array(GMData[GMid])

    a[0]=p[0]/m

    a1 = 1/beta/Dt[GMid]**2*m+gamma/beta/Dt[GMid]*c
    a2 = 1/beta/Dt[GMid]*m+(gamma/beta-1)*c
    a3 = (1./2./beta-1.)*m+Dt[GMid]*(gamma/2./beta-1.)*c

    khat = k + a1

    for i in range(1,len(GMData[GMid])):
        phat = p[i] + a1*u[i-1] + a2*v[i-1] + a3*a[i-1]
        u[i] = phat/khat
        v[i] = gamma/beta/Dt[GMid]*(u[i]-u[i-1]) + (1-gamma/beta)*v[i-1]+Dt[GMid]*(1.-gamma/2./beta)*a[i-1]
        a[i] = 1./beta/Dt[GMid]**2.0*(u[i]-u[i-1])-1/beta/Dt[GMid]*v[i-1]-(1/beta/2-1)*a[i-1]

    a = a-p

    return a

def GetAbsAccWithFreeVibration(GM_Folder, GMid, Zeta, T, FreeVibTime=30.):
    import GMHelper
    import numpy as np
    GMids, GMFiles, Dt, NumPoints, GMData = GMHelper.ExtractGroundMotion(GM_Folder)

    GMData = list(GMData[GMid]) + list(np.zeros(int(30/Dt[GMid])))

    #Using Newmarks: Linear System
    #Using Linear Acceleration Method
    #From Chopra, Page 177

    #Mass
    m = 1
    wn = 2*np.pi/T

    #Stiffness
    k = wn**2.0*m

    #Damping 2%
    c=Zeta*2*m*wn

    gamma = 0.5
    beta = 1./4.

    u = np.zeros(len(GMData))
    v = np.zeros(len(GMData))
    a = np.zeros(len(GMData))
    p = np.array(GMData)

    a[0]=p[0]/m

    a1 = 1/beta/Dt[GMid]**2*m+gamma/beta/Dt[GMid]*c
    a2 = 1/beta/Dt[GMid]*m+(gamma/beta-1)*c
    a3 = (1./2./beta-1.)*m+Dt[GMid]*(gamma/2./beta-1.)*c

    khat = k + a1

    for i in range(1,len(GMData)):
        phat = p[i] + a1*u[i-1] + a2*v[i-1] + a3*a[i-1]
        u[i] = phat/khat
        v[i] = gamma/beta/Dt[GMid]*(u[i]-u[i-1]) + (1-gamma/beta)*v[i-1]+Dt[GMid]*(1.-gamma/2./beta)*a[i-1]
        a[i] = 1./beta/Dt[GMid]**2.0*(u[i]-u[i-1])-1/beta/Dt[GMid]*v[i-1]-(1/beta/2-1)*a[i-1]

    a = a-p

    return a

def GetSS2p0(GM_Folder, GMid, Zeta, T1, T2, PeriodDivisions=30, TimeDivisions=30):
    import GMHelper
    import numpy as np
    GMids, GMFiles, Dt, NumPoints, GMData = GMHelper.ExtractGroundMotion(GM_Folder)

    Periods = np.linspace(T1,T2,PeriodDivisions)
    Duration = np.linspace(0,len(GMData[GMid])*Dt[GMid],TimeDivisions)

    dT = (T2-T1)/float(PeriodDivisions)
    dDs = float(len(GMData[GMid])*Dt[GMid])/float(TimeDivisions)

    SaTogram = [] #Period, Duration
    for t in Periods:
        a, time = GetSDOFAcc(GM_Folder, GMid, Zeta, t, TimeDivisons=TimeDivisions)
        SaTogram.append(a)
    SaTogram = np.array(SaTogram)

    maxSaAtT1 = max(SaTogram[0,:]) #t. T
    indexMax = SaTogram[0,:].tolist().index(maxSaAtT1)

    Num = sum(sum(SaTogram[:,indexMax:]*dT*dDs))
    Din = maxSaAtT1*(T2-T1)*(dDs*(TimeDivisions-indexMax))

    SS2p0 = Num/Din

    return SS2p0

def GetSS2p1(GM_Folder, GMid, Zeta, T1, T2, R=1, PeriodDivisions=10, dT=0.05):
    import GMHelper
    import numpy as np
    GMids, GMFiles, Dt, NumPoints, GMData = GMHelper.ExtractGroundMotion(GM_Folder)

    Periods = np.linspace(T1,T2,PeriodDivisions)

    if Periods[1]-Periods[0] > dT:
        PeriodDivisions = int((T2-T1)/dT)+1
        Periods = np.linspace(T1,T2,PeriodDivisions)

    a = abs(GetAbsAcc(GM_Folder, GMid, Zeta, T1))
    SaTarget = max(a)/R

    minDiff = []
    for i in range(len(a)):
        minDiff.append(abs(SaTarget-a[i]))
    indexMin = minDiff.index(min(minDiff))

    # Sas = []
    SasReduced = []
    for t in Periods:
        a = GetAbsAcc(GM_Folder, GMid, Zeta, t)
        # Sas.append(max(abs(a)))
        SasReduced.append(max(abs(a[indexMin:])))
    # Sas = np.array(Sas)
    SasReduced = np.array(SasReduced)

    # Sadiff = Sas-SasReduced

    SS = np.trapz(SasReduced, Periods)/(SasReduced[0]*(T2-T1))

    return SS

def GetSS2p1WithDataArray(AccArray, ArrayPeriods, T1, T2, R=1):
    import numpy as np
    dT = ArrayPeriods[1]-ArrayPeriods[0]

    Periods = np.arange(T1, T2, dT)
    # Find Acc Array
    periodIndex = np.where(ArrayPeriods==(round(T1,2)))[0][0]
    a = AccArray[periodIndex,:]
    SaTarget = max(a)/R

    minDiff = []
    for i in range(len(a)):
        minDiff.append(abs(SaTarget-a[i]))
    indexMin = minDiff.index(min(minDiff))

    SasReduced = []
    for t in Periods:
        periodIndex = np.where(ArrayPeriods==round(t,2))[0][0]
        a = AccArray[periodIndex,:]
        SasReduced.append(max(abs(a[indexMin:])))
    SasReduced = np.array(SasReduced)

    SS = np.trapz(SasReduced, Periods)/(SasReduced[0]*(T2-T1))

    return SS

def GetSS2p5WithDataArray(AccArray, ArrayPeriods, T1, T2, R=1):
    import numpy as np
    dT = ArrayPeriods[1]-ArrayPeriods[0]

    Periods = np.arange(T1, T2, dT)
    # Find Acc Array
    periodIndex = np.where(ArrayPeriods==(round(T1,2)))[0][0]
    a = AccArray[periodIndex,:]
    SaTarget = max(a)/R

    minDiff = []
    for i in range(len(a)):
        minDiff.append(abs(SaTarget-a[i]))
    indexMin = minDiff.index(min(minDiff))

    numDurationPoints = len(a)-indexMin

    SasReduced = []
    for t in Periods:
        periodIndex = np.where(ArrayPeriods==round(t,2))[0][0]
        a = AccArray[periodIndex,:]
        SasReduced.append(sum(abs(a[indexMin:])))
    SasReduced = sum(np.array(SasReduced))*dT

    SS = SasReduced/(SaTarget)

    return SS

def GetTiFromFFT(a, Fs, FilterPeriod, width=2.0):
    import scipy.fftpack as spfft
    import numpy as np
    import math

    NFFT = 2**math.ceil(np.log(len(a))/np.log(2))
    Y = spfft.fft(a,NFFT)
    freq = Fs/2.0*np.linspace(0.00001,1,NFFT/2+1)

    #Apply Filter at centroid
    Center = sum(abs(Y[0:NFFT/2+1])*freq)/sum(abs(Y[0:NFFT/2+1]))

    #Gaussian Filter
    filter = np.exp(-width*(freq-Center)**2.0)

    ffta = abs(Y[0:NFFT/2+1])*filter
    ffta = ffta/max(ffta)

    fftunfiltered = abs(Y[0:NFFT/2+1])
    fftunfiltered = fftunfiltered/max(fftunfiltered)

    peri = 1.0/freq

    fftunfiltered = fftunfiltered[::-1]
    ffta = ffta[::-1]
    peri = peri[::-1]
    freq = freq[::-1]
    filter = filter[::-1]

    Bracket = 0.5
    fftunfilteredbracketed  = np.where(fftunfiltered-Bracket>0,fftunfiltered,0)

    approx_period_uf = 1/(sum(fftunfilteredbracketed*freq)/sum(fftunfilteredbracketed))
    approx_period = 1/(sum(ffta*freq)/sum(ffta))

    class Output:
        def __init__(self):
            self.Tfft = approx_period
            self.TfftUF = approx_period_uf
            self.Coefficients = ffta
            self.CoefficientsUF = fftunfiltered
            self.Periods = peri
            self.Frequency = freq
            self.Filter = filter
    O = Output()

    return O

def GetDsOnBuildingAcc(GM_Folder, GMid, Zeta, T, min=0.05, max=.95):
    import GMHelper

    GMids, GMFiles, Dt, NumPoints, GMData = GMHelper.ExtractGroundMotion(GM_Folder)

    Acc = GetAbsAcc(GM_Folder, GMid, Zeta, T)

    return FindDs(Acc, Dt[GMid], min, max)

def GetGeometricSSThroughSprectrum(Folder_Location, GMid, T1, T2):
    O = GetResponseSpectrumDataPoints(Folder_Location,GMid)
    Period = O.Period
    Acc = O.Sa
    Vel = O.Sv
    Disp = O.Sd

    PeriodTemp = []
    AccTemp = []

    if T1 == T2:
        return 0.0

    for i in range(len(Period)):
        if Period[i] <= T2 and Period[i] >= T1:
            PeriodTemp.append(Period[i])
            AccTemp.append(Acc[i])

    import scipy.stats.mstats as sm

    SS = sm.gmean(AccTemp)/AccTemp[0]

    return SS

def GetSaGeoMean(Folder_Location, GMid, T1, T2):
    O = GetResponseSpectrumDataPoints(Folder_Location,GMid)
    Period = O.Period
    Acc = O.Sa
    Vel = O.Sv
    Disp = O.Sd

    PeriodTemp = []
    AccTemp = []

    if T1 == T2:
        return 0.0

    for i in range(len(Period)):
        if Period[i] <= T2 and Period[i] >= T1:
            PeriodTemp.append(Period[i])
            AccTemp.append(Acc[i])

    import scipy.stats.mstats as sm
    SS = sm.gmean(AccTemp)

    return SS

def GetSaGeoMeanFromSa(Period, Sa, T1, T2):
    Acc = Sa

    PeriodTemp = []
    AccTemp = []

    if T1 == T2:
        return 0.0

    for i in range(len(Period)):
        if Period[i] <= T2 and Period[i] >= T1:
            PeriodTemp.append(Period[i])
            AccTemp.append(Acc[i])

    import scipy.stats.mstats as sm
    SS = sm.gmean(AccTemp)

    return SS

def GetGMLatLong(Folder_Location):
    f = open(Folder_Location+'GMLocations.txt')
    Locations = []
    for line in f.readlines():
        data = line.split()
        temp = {}
        temp['GMid'] = data[0]
        temp['Long'] = float(data[1])
        temp['Lat'] = float(data[2])
        Locations.append(temp)
    return Locations

def GetSaRotD(GMData1, GMData2, Dt, T, Zeta, Theta=15):
    import numpy as np
    Rot = np.arange(0,360.01,Theta)
    GMData = []
    Sa = []
    for i in range(len(Rot)):
        minInd = min(len(GMData1),len(GMData2))
        Data = np.array(GMData1)[:minInd]*np.cos(Rot[i]*2.0*np.pi/360.)+np.array(GMData2)[:minInd]*np.sin(Rot[i]*2.0*np.pi/360.)
        GMData.append(Data)
        u, v, a = FindSa(Data, Dt, T, Zeta)
        Sa.append(a)

    class Output():
        def __init__(self):
            self.SaRotD50 = np.median(Sa)
            self.SaRotD100 = np.max(Sa)
            self.SaRotD0 = np.min(Sa)
            self.Sas = Sa
            self.Dt = Dt
            for i in range(len(Rot)):
                if Sa[i] == self.SaRotD50:
                    self.RotD50 = Rot[i]
                    self.DataRotD50 = GMData[i]
                elif Sa[i] == self.SaRotD100:
                    self.RotD100 = Rot[i]
                    self.DataRotD100 = GMData[i]
                elif Sa[i] == self.SaRotD0:
                    self.RotD0 = Rot[i]
                    self.DataRotD0 = GMData[i]

    return Output()

def GetSaArithmeticMean(GMData1, GMData2, Dt, T, Zeta, Theta=15):
    u1, v1, a1 = FindSa(GMData1, Dt, T, Zeta)
    u2, v2, a2 = FindSa(GMData2, Dt, T, Zeta)
    return (a1*a2)**.5

def GetJapanStationInformation():
    import pandas as pd
    # import numpy as np
    # Stations = pd.read_csv('GroundMotions/StationData/Japan/sitepub_all_en_reformatted.csv')
    # Folder = 'GroundMotions/StationData/Japan/SoilData/'
    # Stations['Vs30'] = 0
    # for root, dirs, files in os.walk(Folder):
    #     for file in files:
    #         print file
    #         f = open(Folder+file,'r')
    #         data = f.readlines()
    #         if data[0].split()[0] == 'N-Value':
    #             data = data[3:]
    #             Vs = []
    #             z = []
    #             for i in range(len(data)):
    #                 try:
    #                     temp = data[i].split()
    #                     if i == 0:
    #                         z.append(0.)
    #                         Vs.append(float(temp[3]))
    #                     Vs.append(float(temp[3]))
    #                     z.append(float(temp[0][:-1]))
    #                 except:
    #                     Vs.append(Vs[-1])
    #                     z.append(30.)
    #                     break
    #         else:
    #             data = data[2:]
    #             Vs = []
    #             z = []
    #             for i in range(len(data)):
    #                 try:
    #                     temp = data[i].split()
    #                     if i == 0:
    #                         z.append(0.)
    #                         Vs.append(float(temp[4][:-1]))
    #                     tempz = float(temp[2][:-1])
    #                     tempVs = float(temp[4][:-1])
    #                 except:
    #                     break
    #                 z.append(tempz)
    #                 Vs.append(tempVs)
    #
    #         sta = file[:-4]
    #         try:
    #             Vs30 = np.mean(np.interp(np.arange(0,31,1),z,Vs))
    #         except:
    #             print 'Exception'
    #             Stations.loc[Stations['Site Code']==sta,'Vs30'] = 0
    #
    #         Stations.loc[Stations['Site Code']==sta,'Vs30'] = Vs30
    #
    # LongLat = np.transpose([Stations['Longitude'].tolist(),Stations['Latitude'].tolist()])
    #
    #
    # import VelocityModel
    # Z1p0 = VelocityModel.GetZValuesJapan(LongLat, 1.0)
    # Z2p5 = VelocityModel.GetZValuesJapan(LongLat, 2.5)
    # Stations['Z1p0'] = 0
    # Stations['Z2p5'] = 0
    # for i in range(len(Stations)):
    #     Stations.iloc[i,11] = Z1p0[i]
    #     Stations.iloc[i,12] = Z2p5[i]
    #
    # Stations.to_csv('GroundMotions/StationData/Japan/StationData.csv')

    # Load CSV File that was Previously Made
    Stations = pd.read_csv('GroundMotions/StationData/Japan/StationData.csv')
    return Stations

def ComputeHousnerSpectralIntensity(Sa, Period, Tstart = 0.1, Tend = 2.5):
    import numpy as np
    Period = np.array(Period)
    Freq = 2.*np.pi/Period
    Sv = np.array(Sa)/Freq*386.4 #Assuming Sa is in G
    ind = np.where((Period>=Tstart) & (Period<=Tend))[0]
    return np.trapz(Sv[ind[0]:ind[-1]], Period[ind[0]:ind[-1]])

def ConvertSacFilesToGMData(SACFolder, GMFolder):
    import obspy
    import numpy as np

    for root, dirs, files in os.walk(SACFolder):
        for file in files:
            if file.endswith('.SAC'):
                st = obspy.read(os.getcwd()+'/GroundMotions/Nisqually/SAC/%s'%file)
                if st[0].stats['channel'].endswith('Z'):
                    continue

                GMData = np.gradient(st[0].data, st[0].stats.delta)/9.81
                Dt = st[0].stats['delta']
                NumPoints = len(GMData)

                ID = file.replace('.SAC','')
                np.savetxt(GMFolder+'SortedEQFile_(%s).dat'%ID, GMData)
                np.savetxt(GMFolder+'DtFile_(%s).dat'%ID, [Dt])
                np.savetxt(GMFolder+'NumPointsFile_(%s).dat'%ID, [NumPoints], fmt='%d')

def GetVs30PNSN(LatLongs):
    import os
    import numpy as np

    Data = np.loadtxt(os.getcwd()+'/GroundMotions/VS30Models/PNSN/waor45.txt', delimiter=',')

    Vs30 = []
    for i in range(len(LatLongs)):
        def round_to(x, acc=0.01):
            return round(x/acc)*acc

        # -124.9375,42.0000
        # -116.4500,49.0000,603
        XGrid = np.arange(-124.9375,-116.44,0.0125)
        YGrid = np.arange(42.0,49.001,0.0125)
        Vs30Grid = Data[:,2].reshape(len(YGrid), len(XGrid))

        lo = round_to(LatLongs[i][1],0.0125)
        la = round_to(LatLongs[i][0],0.0125)

        # print lo
        # print la

        # for j in range(len(Data)):
        #     if Data[j][0] == lo and Data[j][1] == la:
        #         Vs30.append(Data[j][2])
        #         j += 1
        #         break
        #
        # if j == len(Data)-1:
        #     Vs30.append(0)

        # if True:
        try:
            x = abs(XGrid-lo).argmin()
            y = abs(YGrid-la).argmin()
            Vs30.append(Vs30Grid[y][x])
        except:
            Vs30.append(0)

    return Vs30

def GetVs30ForSeattle():
    import os
    import numpy as np

    Data = np.loadtxt(os.getcwd()+'/GroundMotions/VS30Models/PNSN/waor45.txt', delimiter=',')


    def round_to(x, acc=0.01):
        return round(x/acc)*acc

    # -124.9375,42.0000: Bottom Left Coordinate
    # -116.4500,49.0000,603 Top Right Coordinate
    XGrid = np.arange(-124.9375,-116.44,0.0125)
    YGrid = np.arange(42.0,49.001,0.0125)
    Vs30Grid = Data[:,2].reshape(len(YGrid), len(XGrid))

    x,y = np.meshgrid(XGrid,YGrid)

    return x, y, Vs30Grid

def ConvertZ1p0toZ2p5(Z1p0):
    ### According to Cambell and Bozorgnia 2007, PEER Report
    return 0.519 + 3.595*Z1p0

def ConvertZ1p5toZ2p5(Z1p5):
    ### According to Cambell and Bozorgnia 2007, PEER Report
    return 0.636 + 1.549*Z1p5

def ConvertZ2p5toZ1p5(Z2p5):
    ### According to Cambell and Bozorgnia 2007, PEER Report
    return (Z2p5-0.636)/1.549

def BasinAmplficiationCB14(Z2p5=1.):
    import numpy as np
    Period = np.array([0.01,0.02,0.03,0.05,0.075,0.1,0.15,0.2,0.25,0.3,0.4,0.5,0.75,1.0,1.5,2.0,3.0,4.0,5.0,7.5,10.,0.0,-1.0])
    c14 = np.array([-0.007,-0.0167,-0.0422,-0.0663,-0.0794,-0.0294,0.0642,0.0968,0.1441,0.1597,0.141,0.1474,0.1764,0.2593,0.2881,0.3112,0.3478,0.3747,0.3382,0.3754,0.3506,- 0.0064,0.106])
    c15 = np.array([-0.207,-0.199,-0.202,-0.339,-0.404,-0.416,-0.407,-0.311,-0.172,-0.084,0.085,0.233,0.411,0.479,0.566,0.562,0.534,0.522,0.477,0.321,0.174,- 0.202,0.332])
    c16 = np.array([0.39,0.387,0.378,0.295,0.322,0.384,0.417,0.404,0.466,0.528,0.54,0.638,0.776,0.771,0.748,0.763,0.686,0.691,0.67,0.757,0.621,0.393,0.585])
    k3 = np.array([1.839,1.84,1.841,1.843,1.845,1.847,1.852,1.856,1.861,1.865,1.874,1.883,1.906,1.929,1.974,2.019,2.11,2.2,2.291,2.517,2.744,1.839,1.929])

    Sj = 1 #1 for Japan, 0 elsewhere

    import numpy as np
    if Z2p5 <= 1.:
        f_sed = (c14 + c15*Sj)*(Z2p5-1.)
    elif Z2p5 < 3.:
        f_sed = 0*Period
    else:
        f_sed = c16*k3*np.exp(-0.75)*(1.- np.exp(-0.25*(Z2p5-3.)))

    return f_sed[:-2], Period[:-2]

def BasinAmplficiationCB14wVS30(Z2p5=1., VS30 = 760):
    import numpy as np
    Period = np.array([0.01,0.02,0.03,0.05,0.075,0.1,0.15,0.2,0.25,0.3,0.4,0.5,0.75,1.0,1.5,2.0,3.0,4.0,5.0,7.5,10.,0.0,-1.0])
    c14 = np.array([-0.007,-0.0167,-0.0422,-0.0663,-0.0794,-0.0294,0.0642,0.0968,0.1441,0.1597,0.141,0.1474,0.1764,0.2593,0.2881,0.3112,0.3478,0.3747,0.3382,0.3754,0.3506,- 0.0064,0.106])
    c15 = np.array([-0.207,-0.199,-0.202,-0.339,-0.404,-0.416,-0.407,-0.311,-0.172,-0.084,0.085,0.233,0.411,0.479,0.566,0.562,0.534,0.522,0.477,0.321,0.174,- 0.202,0.332])
    c16 = np.array([0.39,0.387,0.378,0.295,0.322,0.384,0.417,0.404,0.466,0.528,0.54,0.638,0.776,0.771,0.748,0.763,0.686,0.691,0.67,0.757,0.621,0.393,0.585])
    k3 = np.array([1.839,1.84,1.841,1.843,1.845,1.847,1.852,1.856,1.861,1.865,1.874,1.883,1.906,1.929,1.974,2.019,2.11,2.2,2.291,2.517,2.744,1.839,1.929])

    Sj = 1 #1 for Japan, 0 elsewhere

    import numpy as np
    if Z2p5 <= 1.:
        f_sed = (c14 + c15*Sj)*(Z2p5-1.)
    elif Z2p5 < 3.:
        f_sed = 0*Period
    else:
        f_sed = c16*k3*np.exp(-0.75)*(1.- np.exp(-0.25*(Z2p5-3.)))

    k1 = [865,865,908,1054,1086,1032,878,748,654,587,503,457,410,400,400,400,400,400,400,400,400,865,400]
    k1 = [- 1.186,- 1.219,- 1.273,- 1.346,- 1.471,- 1.624,- 1.931,- 2.188,- 2.381,- 2.518,- 2.657,- 2.669,- 2.401,- 1.955,- 1.025,- 0.299,0.0,0.0,0.0,0.0,0.0,- 1.186,- 1.955]

    c11 = [1.094,1.149,1.29,1.449,1.535,1.615,1.877,2.069,2.205,2.306,2.398,2.355,1.995,1.447,0.33,- 0.514,- 0.848,- 0.793,- 0.748,- 0.664,- 0.576,1.09,1.713]
    c12 = [2.191,2.189,2.164,2.138,2.446,2.969,3.544,3.707,3.343,3.334,3.544,3.016,2.616,2.47,2.108,1.327,0.601,0.568,0.356,0.075,- 0.027,2.186,2.602]

    f_site_G = c11 * ()

    import pygmm
    A1100 = pygmm.CampbellBozorgnia2014(mag=7, )
    f_site = f_site_g + S_J * f_site_J


    return f_sed[:-2], Period[:-2]

def BasinAmplificationMF14(Z1p4=1000):
    import numpy as np

    Period = np.array([0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.11, 0.12, 0.13, 0.15, 0.17, 0.2, 0.22, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.5, 1.7, 2, 2.2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 6.5, 7, 7.5, 8, 8.5, 9, 9.5, 10
])
    pd = np.array([-0.0043, -0.0205, -0.0335, -0.0396, -0.0383, -0.0315, -0.0236, -0.0176, -0.0088, 0.0072, 0.0235, 0.046, 0.0583, 0.0746, 0.1006, 0.1206, 0.1418, 0.1599, 0.176, 0.2023, 0.2207, 0.237, 0.2532, 0.2744, 0.2917, 0.3062, 0.3175, 0.3391, 0.3552, 0.3759, 0.3846, 0.3916, 0.3996, 0.4085, 0.4108, 0.412, 0.4109, 0.4078, 0.4088, 0.402, 0.391, 0.3783, 0.3671, 0.3553, 0.3438, 0.332, 0.3202
])
    Dlmin = np.array([15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15.62, 17, 17.86, 19.09, 21, 22.75, 24.39, 25.93, 27.4, 30.13, 32.65, 35, 37.22, 39.32, 41.32, 43.23, 45.07, 48.56, 51.84, 56.42, 59.29, 63.37, 69.69, 75.52, 80.96, 86.08, 90.94, 95.57, 100, 100, 100, 100, 100, 100, 100, 100, 100
])
    D0 = np.array([250, 250, 250, 250, 250, 250, 250, 250, 250, 250, 250, 250, 250, 250, 250, 250, 250, 250, 250, 250, 250, 250, 250, 250, 250, 250, 250, 250, 250, 250, 250, 250, 250, 250, 250, 250, 250, 250, 250, 250, 250, 250, 250, 250, 250, 250, 250
])
    Z1p4 = Z1p4*np.ones(len(Period))

    return pd*np.log(np.max([Dlmin, Z1p4],0)/ D0), Period

def FitCollapseFragility(CollapseSa, PropabilityCollapsed):
    import numpy as np
    from scipy import stats
    from scipy.optimize import minimize

    def FragilityBias(MedianSTDLn):
        guassian = stats.norm(loc=MedianSTDLn[0], scale=MedianSTDLn[1])
        FittedFragility = guassian.cdf(np.log(CollapseSa))
        Bias = FittedFragility - np.array(PropabilityCollapsed)
        return np.sum(np.abs(Bias))

    Results = minimize(FragilityBias, np.array([np.exp(np.mean(np.log(CollapseSa))),np.exp(np.std(np.log(CollapseSa)))]), method='nelder-mead', options = {'xtol': 1e-8, 'disp': False})

    return Results.x

def ComputeProbabilityOfExceedance(DataValues, ExceedingValue):
    import numpy as np
    from scipy import stats

    guassian = stats.norm(loc=np.mean(DataValues), scale=np.std(DataValues))
    return 1 - guassian.cdf(ExceedingValue)

def FitFragilityAndComputeProbOfExceedance(DataValues, Probabilities, ExceedingValue):
    import numpy as np
    from scipy import stats
    from scipy.optimize import minimize

    def FragilityBias(MedianSTDLn):
        guassian = stats.norm(loc=MedianSTDLn[0], scale=MedianSTDLn[1])
        FittedFragility = guassian.cdf(np.log(DataValues))
        Bias = FittedFragility - np.array(Probabilities)
        return np.sum(np.abs(Bias))

    Results = minimize(FragilityBias,
                       np.array([np.exp(np.mean(np.log(DataValues))),
                                 np.exp(np.std(np.log(DataValues)))]),
                       method='nelder-mead',
                       options = {'xtol': 1e-8, 'disp': False})

    guassian = stats.norm(loc=Results.x[0], scale=Results.x[1])
    return 1. - guassian.cdf(np.log(ExceedingValue))

