
from matplotlib import animation, rc
rc('animation', html='html5')
import matplotlib.pylab as plt
import os
import numpy as np

def AnimateNormalizedStrains(O, Skip=1):
    #     Ad Normalized Steel Strain / Steel Yielding
    #     import matplotlib.pylab as plt
    #     from matplotlib import animation, rc
    #     from IPython.display import HTML
    #     import matplotlib.pylab as plt
    fig = plt.figure(figsize=(6.5 ,6.5))
    fig.subplots_adjust(hspace=0.3)

    ax1 = plt.subplot(1 ,3 ,1)
    ax1.set_aspect('equal')

    plt.title('$\epsilon/\epsilon_{cu}$', fontsize=12)

    levels = np.linspace(0, 1.0, 11)

    X = O.CoreXLocation
    Y = O.CoreYLocation
    tempZ1 = []
    for i in range(len(O.CoreStrain[0 ,0 ,:])):
        tempZ1.append(O.CoreStrain[: ,: ,i ] /O.CoreCrushingStrain)
    Z = np.max(tempZ1, axis=0)

    pcol1 = ax1.contourf(X, Y, Z, levels, cmap='afmhot_r', alpha=0.5, extend = 'both')
    cbar1 = plt.colorbar(pcol1, shrink=1.0, pad=0.10)
    cbar1.ax.tick_params(labelsize=8)

    # cbar1.ax.set_ylabel('$\epsilon/\epsilon_cu$', fontsize=8)
    # cbar.add_lines(cbar)

    ax1.xaxis.set_tick_params(labelsize=8)
    ax1.yaxis.set_tick_params(labelsize=8)

    ax2 = plt.subplot(1 ,3 ,2)
    ax2.set_aspect('equal')

    plt.title('$\epsilon/\epsilon_{sy}$', fontsize=12)

    # ax = plt.axes()
    # ax.set_aspect('equal')

    levels = np.linspace(0, 1, 11)

    X = O.CoreXLocation
    Y = O.CoreYLocation
    tempZ2 = []
    for i in range(len(O.CoreStrain[0 ,0 ,:])):
        tempZ2.append(O.CoreStrain[: ,: ,i ] /O.CoreYieldingStrain)
    Z = np.max(tempZ2, axis=0)

    pcol2 = ax2.contourf(X, Y, Z, levels, cmap='afmhot_r', alpha=0.5, extend = 'both')
    cbar2 = plt.colorbar(pcol2, shrink=1.0, pad=0.10)
    cbar2.ax.tick_params(labelsize=8)

    # cbar2.ax.set_ylabel('$\epsilon/\epsilon_{s,y}$', fontsize=8)
    # cbar.add_lines(cbar)

    ax2.yaxis.set_visible(False)
    ax2.xaxis.set_tick_params(labelsize=8)

    ax3 = plt.subplot(1 ,3 ,3)
    ax3.set_aspect('equal')

    plt.title('$\epsilon/\epsilon_{su}$', fontsize=12)

    # ax = plt.axes()
    # ax.set_aspect('equal')

    levels = np.linspace(0, 1, 11)

    X = O.CoreXLocation
    Y = O.CoreYLocation
    tempZ3 = []
    for i in range(len(O.CoreStrain[0 ,0 ,:])):
        tempZ3.append(O.CoreStrain[: ,: ,i ] /O.CoreRuptureStrain)
    Z = np.max(tempZ3, axis=0)

    pcol3 = ax3.contourf(X, Y, Z, levels, cmap='afmhot_r', alpha=0.5, extend = 'both')
    cbar3 = plt.colorbar(pcol3, shrink=1.0, pad=0.10)
    cbar3.ax.tick_params(labelsize=8)

    ax3.yaxis.set_visible(False)
    ax3.xaxis.set_tick_params(labelsize=8)

    # cbar3.ax.set_ylabel('$\epsilon/\epsilon_{s,u}$', fontsize=8)

    col = pcol1.collections + pcol2.collections + pcol3.collections
    try:
        for c in col: c.remove()
    except:
        pass
    global DamagedZ1, DamagedZ2, DamagedZ3

    DamagedZ1 = np.array(tempZ1[0] ) *0.0 -1.
    DamagedZ2 = np.array(tempZ1[0]) * 0.0 - 1.
    DamagedZ3 = np.array(tempZ1[0]) * 0.0 - 1.

    def updateDamaged(var, var_step):
        #         if n == 0:
        #             var = np.array(tempZ1[0])*0.0 -1.
        for i in range(len(var_step)):
            for j in range(len(var_step[0])):
                if var_step[i][j] >= 1.0:
                    var[i][j] = 1.

        scat = plt.contourf(X, Y, var, [0, 1.0], cmap='gray_r', alpha=0.25)
        return scat

    def animate(n):
        global DamagedZ1, DamagedZ2, DamagedZ3
        if n == 0:
            DamagedZ1 = np.array(tempZ1[0]) * 0.0 - 1.
            DamagedZ2 = np.array(tempZ1[0]) * 0.0 - 1.
            DamagedZ3 = np.array(tempZ1[0]) * 0.0 - 1.

        n = n * Skip
        global pcol1, pcol2, pcol3, scat1, scat2, scat3
        try:
            col = pcol1.collections + pcol2.collections + pcol3.collections
            for c in col: c.remove()
        except:
            pass

        try:
            col = scat1.collections + scat2.collections + scat3.collections
            for c in col: c.remove()
        except:
            pass

        # ax = plt.subplot(1, 3, 1)
        pcol1 = ax1.contourf(X, Y, tempZ1[n], levels, cmap='afmhot_r', alpha=0.75, extend='both')
        scat1 = updateDamaged(DamagedZ1, tempZ1[n])

        # ax.xaxis.set_visible(False)

        # ax = plt.subplot(1, 3, 2)
        pcol2 = ax2.contourf(X, Y, tempZ2[n], levels, cmap='afmhot_r', alpha=0.75, extend='both')
        scat2 = updateDamaged(DamagedZ2, tempZ2[n])

        ax2.get_xaxis().set_ticks([])

        plt.xlabel('step: %d' % n)

        # ax = plt.subplot(1, 3, 3)
        pcol3 = ax3.contourf(X, Y, tempZ3[n], levels, cmap='afmhot_r', alpha=0.75, extend='both')
        scat3 = updateDamaged(DamagedZ3, tempZ3[n])

        # ax.xaxis.set_visible(False)

        if n >= (NoOfFrames - 1) * Skip - 1.:
            try:
                col = pcol1.collections + pcol2.collections + pcol3.collections
                for c in col: c.remove()
            except:
                pass

    NoOfFrames = int(len(O.CoreStrain[0, 0, :]) / Skip)
    anim = animation.FuncAnimation(fig, animate, frames=NoOfFrames)

    Time = 10.0
    Writer = animation.writers['ffmpeg']
    writer = Writer(fps=NoOfFrames / Time, metadata=dict(artist='Me'), bitrate=12 * 1e6)
    anim.save(os.getcwd() + '/Figures/' + 'ATCWallAnimation-%s.mp4' % O.Archetype.Name, writer=writer)

    return anim


def AnimateWallDeformation(O, Dt, Skip=10, GravitySystemXLoc=None, GammaRacking=1.0, AnimationFileName='Animation'):
    import matplotlib.pylab as plt
    from matplotlib import animation, rc
    from IPython.display import HTML
    from matplotlib.lines import Line2D

    fig = plt.figure(figsize=(6.5, 6.5))

    ax = plt.axes()
    ax.set_aspect('equal')

    X = np.array(O.CoreXLocation)
    Y = np.array(O.CoreYLocation)
    Z = O.CoreStrain[:, :, 0]

    ax = fig.gca()

    plt.xlim([-500, 1000])

    plt.ylim([0, round(1.1 * O.YGrids[-1] / 500) * 500])

    levels = np.linspace(-0.03, 0.03, 10)

    levels = [-0.005, 0, 0.0025, 0.005, 0.01, 0.015, 0.02, 0.03]

    ax.xaxis.set_visible(False)

    ax.set_yticks(O.YGrids)
    ax.set_yticklabels(
        ['%d' % (i - len(O.Archetype.BasementProperties.BasementMass)) for i in range(len(O.YGrids) + 2)])

    import matplotlib.colors
    cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", ["#b0b0b0", "#cb181d"])

    pc = ax.contourf(X, Y, Z, levels, cmap=cmap, alpha=0.5, extend='both')
    cs = ax.contour(X, Y, Z, levels, cmap=cmap, alpha=0.5, extend='both')

    cbar = plt.colorbar(pc, shrink=1.0, pad=0.10)
    cbar.ax.set_ylabel('Strain')
    cbar.add_lines(cs)

    # Compute Location of Core
    AmpFactor = 1.0

    import matplotlib.lines as mlines
    from scipy.stats import norm
    AddKey = True
    if AddKey:
        rot = 20.
        Rotation = int(4 + rot * 0.09 * AmpFactor * 3)
        slab_0 = mlines.Line2D([], [], color='#b0b0b0', marker='o',
                               markersize=15, label='>20%')
        rot = 50.
        Rotation = int(4 + rot * 0.09 * AmpFactor * 3)
        slab_50 = mlines.Line2D([], [], color='#fb6a4a', marker='o',
                                markersize=Rotation, label='>50%')
        rot = 80.
        Rotation = int(4 + rot * 0.09 * AmpFactor * 3)
        slab_100 = mlines.Line2D([], [], color='#cb181d', marker='o',
                                 markersize=Rotation, label='>80%')
        plt.legend(handles=[slab_100, slab_50, slab_0], loc=2, labelspacing=2, frameon=False)

    try:
        for c in pc.collections: c.remove()
    except:
        pass
    try:
        for c in cs.collections: c.remove()
    except:
        pass

    if GravitySystemXLoc == None:
        GravitySystemXLoc = O.XGrids[-1]

    Rotations = np.zeros(len(O.YGrids) - 1)

    ax.axhline(O.YGrids[len(O.Archetype.BasementProperties.BasementMass)], color='#000000', linewidth=1.0,
               linestyle='--', zorder=1)

    global DamagedZ1, scat
    DamagedZ1 = np.array(X) * 0.0 - 1.

    def updateDamaged(var, var_step, X, Y, value=0.002):
        for i in range(len(var_step)):
            for j in range(len(var_step[0])):
                if abs(var_step[i][j]) >= value:
                    var[i][j] = 1.

        scat = ax.contourf(X, Y, var, [0, 1.0], cmap='gray_r', alpha=0.25)
        return scat

    def animate(n):
        n = n * Skip

        global pc, cs, pc2, cs2, pc3, pc4, pc5, Rotations, scat, DamagedZ1

        Z = O.CoreStrain[:, :, n]

        try:
            for c in pc.collections: c.remove()
        except:
            pass

        try:
            for c in cs.collections: c.remove()
        except:
            pass

        try:
            for c in pc2.collections: c.remove()
        except:
            pass

        try:
            for c in cs2.collections: c.remove()
        except:
            pass

        try:
            col = scat.collections
            for c in col: c.remove()
        except:
            pass

        if O.Archetype.CustomSection is not None:
            LWall = []
            for i in range(len(O.Archetype.CustomSection)):
                LWall.append(O.Archetype.CustomSection[i].l_w)
        else:
            if O.CoupledDirection == False:
                LWall = O.Archetype.l_w
            else:
                LWall = (O.Archetype.b_f - O.Archetype.CouplingBeamLength) / 2.

        XLoc = []
        for i in range(len(LWall)):
            XLoc.append(np.linspace(-.5, .5, O.NoOfSamplePoints) * LWall[i])
        XLoc = np.array(XLoc)

        if O.NoOfIntPoints == 3:
            # L_ip = [0.6666667, 0.1666667, 0.6666667]
            L_ip = [0.1666667, 0.6666667, 0.1666667]
        elif O.NoOfIntPoints == 5:
            # L_ip = [0.35555556, 0.27222222, 0.05000000, 0.27222222, 0.35555556]
            L_ip = [0.05000000, 0.27222222, 0.35555556, 0.27222222, 0.05000000]
        else:
            # L_ip = [0.243810, 0.215873, 0.138413, 0.023810, 0.138413, 0.215873, 0.243810]
            L_ip = [0.023810, 0.138413, 0.215873, 0.243810, 0.215873, 0.138413, 0.023810]

        for i in range(1, len(O.YGrids)):
            for k in range(O.NoOfDivisionsPerFloor):
                if O.CoupledDirection == False:
                    colNo = 1 + 3 * (i - 1) * O.NoOfDivisionsPerFloor + 3 * k
                    colNoBefore = colNo - 3

                    if i != 1:
                        dxx = O.CoreDisp[n, colNo] - O.CoreDisp[n, colNoBefore]
                        dyy = O.CoreDisp[n, colNo + 1] - O.CoreDisp[n, colNoBefore + 1]
                    else:
                        dxx = O.CoreDisp[n, colNo]
                        dyy = O.CoreDisp[n, colNo + 1]

                    dx = O.CoreDisp[n, colNo]  # - O.CoreDisp[n, 1+3*(i-1)]
                    dy = O.CoreDisp[n, colNo + 1]  # - O.CoreDisp[n, 2+3*(i-1)]
                    drz = O.CoreDisp[n, colNo + 2]

                    fac = L_ip

                    for j in range(O.NoOfIntPoints):
                        colNo = (i - 1) * O.NoOfIntPoints * O.NoOfDivisionsPerFloor + k * O.NoOfIntPoints + j
                        X[:, colNo] = O.CoreXLocation[:, colNo] + (
                                    dx - dxx * (fac[j] / 2. + np.sum(fac[j + 1:]))) * AmpFactor
                        Y[:, colNo] = O.CoreYLocation[:, colNo] + (
                                    dy - dyy * (fac[j] / 2. + np.sum(fac[j + 1:]))) * AmpFactor

                        # Add Roation
                        Y[:, colNo] = Y[:, colNo] + drz * XLoc[i - 1] * AmpFactor

                else:
                    colNo = 1 + 6 * (i - 1) * O.NoOfDivisionsPerFloor + 6 * k
                    colNoBefore = colNo - 6

                    if i != 1:
                        dxxL = O.CoreDisp[n, colNo] - O.CoreDisp[n, colNoBefore]
                        dyyL = O.CoreDisp[n, colNo + 1] - O.CoreDisp[n, colNoBefore + 1]

                        dxxR = O.CoreDisp[n, colNo + 3] - O.CoreDisp[n, colNoBefore + 3]
                        dyyR = O.CoreDisp[n, colNo + 1 + 3] - O.CoreDisp[n, colNoBefore + 1 + 3]
                    else:
                        dxxL = O.CoreDisp[n, colNo]
                        dyyL = O.CoreDisp[n, colNo + 1]

                        dxxR = O.CoreDisp[n, colNo + 3]
                        dyyR = O.CoreDisp[n, colNo + 1 + 3]

                    dxL = O.CoreDisp[n, colNo]  # - O.CoreDisp[n, 1+3*(i-1)]
                    dyL = O.CoreDisp[n, colNo + 1]  # - O.CoreDisp[n, 2+3*(i-1)]
                    drzL = O.CoreDisp[n, colNo + 2]

                    dxR = O.CoreDisp[n, colNo + 3]  # - O.CoreDisp[n, 1+3*(i-1)]
                    dyR = O.CoreDisp[n, colNo + 1 + 3]  # - O.CoreDisp[n, 2+3*(i-1)]
                    drzR = O.CoreDisp[n, colNo + 2 + 3]

                    fac = L_ip

                    for j in range(O.NoOfIntPoints):
                        if i == 1 and j == 0 and k == 0:
                            continue

                        colNo = (i - 1) * O.NoOfIntPoints * O.NoOfDivisionsPerFloor + j + k * O.NoOfIntPoints
                        X[:O.NoOfSamplePoints, colNo] = O.CoreXLocation[:O.NoOfSamplePoints, colNo] + (
                                    dxL - dxxL * (fac[j] / 2. + np.sum(fac[j + 1:]))) * AmpFactor
                        Y[:O.NoOfSamplePoints, colNo] = O.CoreYLocation[:O.NoOfSamplePoints, colNo] + (
                                    dyL - dyyL * (fac[j] / 2. + np.sum(fac[j + 1:]))) * AmpFactor

                        # Add Rotation
                        Y[:O.NoOfSamplePoints, colNo] = Y[:O.NoOfSamplePoints, colNo] + drzL * XLoc[i - 1] * AmpFactor

                        # Right Columns

                        colNo = (i - 1) * O.NoOfIntPoints * O.NoOfDivisionsPerFloor + j + k * O.NoOfIntPoints
                        X[O.NoOfSamplePoints:O.NoOfSamplePoints * 2., colNo] = O.CoreXLocation[
                                                                               O.NoOfSamplePoints:O.NoOfSamplePoints * 2.,
                                                                               colNo] + (dxR - dxxR * (
                                    fac[j] / 2. + np.sum(fac[j + 1:]))) * AmpFactor
                        Y[O.NoOfSamplePoints:O.NoOfSamplePoints * 2., colNo] = O.CoreYLocation[
                                                                               O.NoOfSamplePoints:O.NoOfSamplePoints * 2.,
                                                                               colNo] + (dyR - dyyR * (
                                    fac[j] / 2. + np.sum(fac[j + 1:]))) * AmpFactor

                        # Add Roation
                        Y[O.NoOfSamplePoints:O.NoOfSamplePoints * 2., colNo] = Y[
                                                                               O.NoOfSamplePoints:O.NoOfSamplePoints * 2.,
                                                                               colNo] + drzR * XLoc[i - 1] * AmpFactor

        # Plotting PDelta Column and Slab-Column Hinges
        if True:
            if n == 0:
                try:
                    DamagedZ1 = np.array(DamagedZ1) * 0.0 - 1.
                    for line in pc3:
                        ax.lines.remove(line)
                    for line in pc4:
                        ax.lines.remove(line)
                    for story in range(1, len(O.Archetype.YGrids)):
                        pc5[story - 1].set_offsets(np.transpose([[],
                                                                 []]))
                    Rotations = np.zeros(len(O.YGrids) - 1)
                except:
                    pass
            if n != 0 and pc3 != []:
                for story in range(1, len(O.Archetype.YGrids)):
                    if story == 1:
                        pc3[story - 1].set_data(
                            [GravitySystemXLoc + 0, GravitySystemXLoc + O.AllDispl[n, story] * AmpFactor],
                            [O.Archetype.YGrids[story - 1], O.Archetype.YGrids[story]])
                    else:
                        pc3[story - 1].set_data([GravitySystemXLoc + O.AllDispl[n, story - 1] * AmpFactor,
                                                 GravitySystemXLoc + O.AllDispl[n, story] * AmpFactor],
                                                [O.Archetype.YGrids[story - 1], O.Archetype.YGrids[story]])

                    pc4[story - 1].set_data([O.Archetype.l_w + O.AllDispl[n, story] * AmpFactor,
                                             GravitySystemXLoc + O.AllDispl[n, story] * AmpFactor],
                                            [O.Archetype.YGrids[story], O.Archetype.YGrids[story]])

                    pc5[story - 1].set_offsets(np.transpose([[O.Archetype.l_w + O.AllDispl[n, story] * AmpFactor + 24,
                                                              GravitySystemXLoc + O.AllDispl[
                                                                  n, story] * AmpFactor - 24],
                                                             [O.Archetype.YGrids[story], O.Archetype.YGrids[story]]]))

                    # Calculate Rotation
                    if story == 1:
                        Rotation = abs((O.AllDispl[n, story]) / (O.YGrids[story])) * 100. * GammaRacking
                        if Rotations[story - 1] < Rotation:
                            Rotations[story - 1] = Rotation
                    else:
                        Rotation = abs((O.AllDispl[n, story] - O.AllDispl[n, story - 1]) / (
                                    O.YGrids[story] - O.YGrids[story - 1])) * 100. * GammaRacking
                        if Rotations[story - 1] < Rotation:
                            Rotations[story - 1] = Rotation

                    rot = norm.cdf(Rotations[story - 1], 5.9, 0.12) * 100.
                    Rotation = int(4 + rot * 0.09 * AmpFactor * 3)

                    pc5[story - 1].set_sizes([Rotation, Rotation])
                    if rot >= 80:
                        pc5[story - 1].set_edgecolors('#cb181d')
                        pc5[story - 1].set_facecolors('#cb181d')
                    elif rot >= 50:
                        pc5[story - 1].set_edgecolors('#fb6a4a')
                        pc5[story - 1].set_facecolors('#fb6a4a')
                    elif rot >= 20:
                        pc5[story - 1].set_edgecolors('#b0b0b0')
                        pc5[story - 1].set_facecolors('#b0b0b0')
                    elif rot >= 20:
                        pc5[story - 1].set_edgecolors('#b0b0b0')
                        pc5[story - 1].set_facecolors('#b0b0b0')
                    else:
                        pc5[story - 1].set_edgecolors('#b0b0b0')
                        pc5[story - 1].set_facecolors('#b0b0b0')

                # Assign Size and Color
            #               print(Rotations)
            else:
                pc3 = []
                pc4 = []
                pc5 = []
                for story in range(1, len(O.Archetype.YGrids)):
                    if story == 1:
                        line, = ax.plot([GravitySystemXLoc + 0, GravitySystemXLoc + O.AllDispl[n, story] * AmpFactor],
                                        [O.Archetype.YGrids[story - 1], O.Archetype.YGrids[story]], color='#b0b0b0',
                                        linewidth=1.5)
                    else:
                        line, = ax.plot([GravitySystemXLoc + O.AllDispl[n, story - 1] * AmpFactor,
                                         GravitySystemXLoc + O.AllDispl[n, story] * AmpFactor],
                                        [O.Archetype.YGrids[story - 1], O.Archetype.YGrids[story]], color='#b0b0b0',
                                        linewidth=1.5)

                    linebeam, = ax.plot([O.Archetype.l_w + O.AllDispl[n, story] * AmpFactor,
                                         GravitySystemXLoc + O.AllDispl[n, story] * AmpFactor],
                                        [O.Archetype.YGrids[story], O.Archetype.YGrids[story]], color='#b0b0b0',
                                        linewidth=1.5)

                    linehinge = ax.scatter([],
                                           [], color='#b0b0b0', linewidth=1.5, zorder=10)

                    pc3.append(line)
                    pc4.append(linebeam)
                    pc5.append(linehinge)

        if O.CoupledDirection == False:
            pc = ax.contourf(X, Y, Z, levels, cmap=cmap, alpha=0.5, extend='both')
            cs = ax.contour(X, Y, Z, levels, cmap=cmap, alpha=0.5, extend='both')

            scat = updateDamaged(DamagedZ1, Z, X, Y, value=0.002)

        else:
            lim = O.NoOfSamplePoints
            pc = ax.contourf(X[:lim], Y[:lim], Z[:lim], levels, cmap=cmap, alpha=0.5, extend='both')
            cs = ax.contour(X[:lim], Y[:lim], Z[:lim], levels, cmap=cmap, alpha=0.5, extend='both')
            pc2 = ax.contourf(X[lim:lim * 2], Y[lim:lim * 2], Z[lim:lim * 2], levels, cmap=cmap, alpha=0.5,
                              extend='both')
            cs2 = ax.contour(X[lim:lim * 2], Y[lim:lim * 2], Z[lim:lim * 2], levels, cmap=cmap, alpha=0.5,
                             extend='both')

        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        #         ax.spines['bottom'].set_visible(False)

        ax.set_title('t = %ds' % (n * Dt))

        return pc, cs

    NoOfFrames = int(len(O.CoreStrain[0, 0, :]) / Skip)
    anim = animation.FuncAnimation(fig, animate, frames=NoOfFrames)

    Time = 60.0
    Writer = animation.writers['ffmpeg']
    writer = Writer(fps=NoOfFrames / Time, metadata=dict(artist='Me'), bitrate=12 * 1e6)
    anim.save(os.getcwd() + '/Figures/' + AnimationFileName + '.mp4', writer=writer)

    return anim

