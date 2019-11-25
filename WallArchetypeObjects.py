# Coded by Marafi
# To Do List:
# Add Concrete Wall Weight
# Add Basement Floors
# Add Using 2.5' Deep by 6' Coupling Beams
# Distribute Forces using MRSA, include 0.85 factor for 2008 designs
####################################################################################
#region Defining Classes
####################################################################################

from __future__ import absolute_import
import numpy as np
import ATCWallArchetypeHelpers as ATCWallHelper
from ATCWallArchetypeObjects import ArchetypeData
from ATCWallArchetypeObjects import CouplingBeam
from ATCWallArchetypeObjects import PlanarWallSection
from ATCWallArchetypeObjects import TWallSection
from ATCWallArchetypeObjects import IWallSection
from six.moves import filter

class Basement:
    def __init__(self, FloorStiffnesses, WallStiffnesses, BasementMass, **kwargs):
        self.FloorStiffnesses = list(FloorStiffnesses)
        self.WallStiffnesses = list(WallStiffnesses)
        self.BasementMass = list(BasementMass)
        self.__dict__.update(kwargs)

####################################################################################
# endregion
####################################################################################

####################################################################################
#region Defining Functions
####################################################################################


####################################################################################
# endregion
####################################################################################

def GetSeattle2008Hazard(Height, R=6, Period = None, IgnoreMinBaseShear = False, Overstrength=1.0):
    Sds = 0.91 * Overstrength
    S1 = 0.529 * Overstrength
    Sd1 = 0.458 * Overstrength
    TL =  6
    I = 1.0
    Cd = 5
    if Period == None:
        CuTa = 1.4 * ASCEHelper.ComputeCuTa((Height / 12.), 0.02, 0.75)
    else:
        CuTa = Period
    Periods = [0.01 ,0.1 ,0.2 ,0.3 ,0.4 ,0.5 ,0.6 ,0.75 ,1 ,2 ,3 ,4 ,5 ,7 ,8 ,9 ,10]
    BasinFactors = [1.235, 1.231, 1.249, 1.340, 1.351, 1.428, 1.477, 1.551, 1.557, 1.583, 1.541, 1.576, 1.581, 1.728, 1.744, 1.703, 1.662]
    # if Height / 12. > 240:
    #     FactorS = np.exp(np.interp(np.log(0.2), np.log(Periods), np.log(BasinFactors))) # Add basin effects
    #     Factor1 = np.exp(np.interp(np.log(1.0), np.log(Periods), np.log(BasinFactors)))  # Add basin effects
    # else:
    #     FactorS = 1.0
    #     Factor1 = 1.0

    FactorS = 1.0
    Factor1 = 1.0
    SaDesign = ASCEHelper.GetDesignSa(CuTa, S1 * Factor1, Sds * FactorS, Sd1 * Factor1, TL, R, I, IgnoreMinBaseShear)
    return SaDesign, Sds, CuTa

def GetSeattle2014Hazard(Height, R=6, Period = None, IgnoreMinBaseShear = False, Overstrength=1.0):
    Sds = 1.12 * Overstrength
    S1 = 0.488 * Overstrength
    Sd1 = 0.488 * Overstrength
    TL =  6
    I = 1.0
    Cd = 5
    if Period == None:
        CuTa = 1.4 * ASCEHelper.ComputeCuTa((Height / 12.), 0.02, 0.75)
    else:
        CuTa = Period
    Periods = [0.01 ,0.1 ,0.2 ,0.3 ,0.4 ,0.5 ,0.6 ,0.75 ,1 ,2 ,3 ,4 ,5 ,7 ,8 ,9 ,10]
    BasinFactors = [1.235, 1.231, 1.249, 1.340, 1.351, 1.428, 1.477, 1.551, 1.557, 1.583, 1.541, 1.576, 1.581, 1.728, 1.744, 1.703, 1.662]
    # if Height / 12. > 240:
    #     FactorS = np.exp(np.interp(np.log(0.2), np.log(Periods), np.log(BasinFactors))) # Add basin effects
    #     Factor1 = np.exp(np.interp(np.log(1.0), np.log(Periods), np.log(BasinFactors)))  # Add basin effects
    # else:
    #     FactorS = 1.0
    #     Factor1 = 1.0
    FactorS = 1.0
    Factor1 = 1.0
    SaDesign = ASCEHelper.GetDesignSa(CuTa, S1 * Factor1, Sds * FactorS, Sd1 * Factor1, TL, R, I, IgnoreMinBaseShear)
    return SaDesign, Sds, CuTa

# Global Variables
# Loading
# DL 150psf Floors  #### All Floors have the same load
# LL 65psf Floors and 20psf Roof
DL_Basements = [155, 155, 155, 230] # Include Basement Wall Loads
LL_Basements = [40, 40, 40, 100] # Check LL
DL = 130 # psf
LL = 50 # psf
DL_Roof = 200 # psf
LL_Roof = 20 # psf
BasementFloorArea = 160. * 160. / 2.
FloorArea = 100. * 100. / 2.
PercentageFloorAreaResistedByWall = 0.5
FirstFloorHeight = 10 * 12.
FloorHeights = 10 * 12.
BasementFloorHeights = 10 * 12.
# Pick Out Prelim. Section Size using Shear
fy = 60.; fu = 105.
fpc_core = 8.
fpc_slabs = 5.
ConcreteDensity = 150.

def CreateArchetype(Basement=None, Use2008Maps = True, Overstrength = 1.0):
    # Defining Story Levels
    YGrids = [0] + np.array(np.arange(0, (NoOfStories) * FloorHeights, FloorHeights) + FloorHeights).tolist()
    # Defining Gravity Loads
    DeadLoads = np.ones(NoOfStories) * DL / 1000.
    DeadLoads[-1] = DeadLoads[-1] * DL_Roof / DL
    LiveLoads = np.ones(NoOfStories) * LL / 1000.
    LiveLoads[-1] = LiveLoads[-1] * LL_Roof / LL
    # Computing Mass of Wall
    WallSelfWeight = []
    i = -1
    for section in Sections:
        i += 1
        if isinstance(section, IWallSection):
            CoreVolume = (section.b_w * section.t_w * 2. + (section.l_w - section.t_w * 2.) * section.t_w) * (
                    YGrids[i + 1] - YGrids[i]) / 12. ** 3.
            EquivalentDL = CoreVolume * ConcreteDensity / 1000.
            WallSelfWeight.append(EquivalentDL)
        elif isinstance(section, PlanarWallSection):
            CoreVolume = section.l_w * section.t_w * (
                    YGrids[i + 1] - YGrids[i]) / 12. ** 3.
            EquivalentDL = CoreVolume * ConcreteDensity / 1000.
            WallSelfWeight.append(EquivalentDL)
    # Defining Mass
    Mass = ( DeadLoads + 0.5 * LiveLoads ) * FloorArea  # Compute
    Mass = Mass + np.array(WallSelfWeight)  # Adding Wall Self Weight
    WallTribArea = FloorArea * PercentageFloorAreaResistedByWall
    WallGravityLoad = WallTribArea * DeadLoads + np.array(WallSelfWeight)
    WallDeadLoads = DeadLoads * FloorArea + np.array(WallSelfWeight)
    WallLiveLoads = LiveLoads * FloorArea
    PDeltaGravityLoad = Mass - WallGravityLoad
    if Basement is not None:
        Height = YGrids[-1] - YGrids[len(Basement.FloorStiffnesses)]
    else:
        Height = YGrids[-1]
    # Seismic Hazard
    R = 6; Cd = 5
    if Use2008Maps:
        SaDesign, Sds, CuTa = GetSeattle2008Hazard(Height, R=R, Overstrength = Overstrength)
    else:
        SaDesign, Sds, CuTa = GetSeattle2014Hazard(Height, R=R, Overstrength = Overstrength)

    if Basement is not None:
        archetypename = ArchetypeData(Name, YGrids, R, CuTa, Length, Thickness, None, None, None,
                                      fpc_core, fy, fu, PDeltaGravityLoad, Mass, WallGravityLoad,
                                      None, None, None, Sections, CuTa=CuTa, SaDesign=SaDesign,
                                      Cd=Cd, BasementProperties=Basement,
                                      WallDeadLoads = list(WallDeadLoads), WallLiveLoads = list(WallLiveLoads), Sds = Sds)
    else:
        archetypename = ArchetypeData(Name, YGrids, R, CuTa, Length, Thickness, None, None, None,
                                    fpc_core, fy, fu, PDeltaGravityLoad, Mass, WallGravityLoad,
                                    None, None, None, Sections, CuTa=CuTa, SaDesign=SaDesign,
                                    Cd=Cd, WallDeadLoads = list(WallDeadLoads), WallLiveLoads = list(WallLiveLoads), Sds = Sds)
    return archetypename

BasementFloorStiffnesses = np.array([8200, 8200, 8200, 10100]) * 0.5
BasementWallStiffnesses = np.array([0.0496e9, 0.0496e9, 0.0496e9, 0.0496e9, ]) * 0.5
BasementMass = (np.array(DL_Basements) + 0.5 * np.array(LL_Basements)) * ( BasementFloorArea - FloorArea ) / 1000.
Basements = Basement(BasementFloorStiffnesses, BasementWallStiffnesses, BasementMass)
BasementFloorStiffnesses = np.array([8200, 8200, 10100]) * 0.5
BasementWallStiffnesses = np.array([0.0496e9, 0.0496e9, 0.0496e9, ]) * 0.5
BasementMass = (np.array(DL_Basements[1:]) + 0.5 * np.array(LL_Basements[1:])) * ( BasementFloorArea - FloorArea ) / 1000.
Basements3Levels = Basement(BasementFloorStiffnesses, BasementWallStiffnesses, BasementMass)
BasementFloorStiffnesses = np.array([8200, 10100]) * 0.5
BasementWallStiffnesses = np.array([0.0496e9, 0.0496e9 ]) * 0.5
BasementMass = (np.array(DL_Basements[2:]) + 0.5 * np.array(LL_Basements[2:])) * ( BasementFloorArea - FloorArea ) / 1000.
Basements2Levels = Basement(BasementFloorStiffnesses, BasementWallStiffnesses, BasementMass)

####################################################################################
#region Defining Archetype
####################################################################################

Archetypes = []
import ASCEHelper

############################### Performance Group #1 ###############################
# 2008 Maps

#region Archetype S4H08SEA and S4H08SEAWB
Name = 'S4H08SEA'
# print 'Importing Archetype: ' + Name
# Compute Seismic Weight
NoOfStories = 4
YGrids = [0] + np.array(np.arange(0,(NoOfStories)*13*12, 13*12)+15*12).tolist()
DeadLoads = np.ones(NoOfStories) * DL / 1000.
DeadLoads[-1] = DeadLoads[-1] * DL_Roof / DL
LiveLoads = np.ones(NoOfStories) * LL / 1000.
LiveLoads[-1] = LiveLoads[-1] * LL_Roof / LL
MassPerSqFt = DL / 1000.
Mass = np.ones(NoOfStories) * MassPerSqFt * FloorArea
Mass[-1] = FloorArea * DL_Roof / 1000. # Adjust for Roof Weight
WallTribArea = FloorArea * 0.5
WeightPerSqFt = DL
BuildingWeight = np.ones(NoOfStories) * WeightPerSqFt * FloorArea
BuildingWeight[-1] = 152. / 1000. * FloorArea # Adjust for Roof Weight
# Seismic Hazard
R = 6; Cd = 5
SaDesign, Sds, CuTa = GetSeattle2008Hazard(YGrids[-1], R=R)
Thickness = 14.
Length = 14. * 12.
Long_Spacing = 4
NoOfCols = 10
BarSize = 8.
Ag = ( (NoOfCols - 1) * Long_Spacing + 6 ) * Thickness
Rho = ( NoOfCols  * 2 + 2 ) * np.pi * ( BarSize / 2. / 8.) ** 2. / Ag
# print Rho
Section1 = PlanarWallSection(Length, Thickness,
                             (NoOfCols - 1) * Long_Spacing + 6,
                             (NoOfCols - 1) * Long_Spacing + 6, BarSize,
                             [3] + (np.ones(NoOfCols - 2) * 2.).tolist() + [3],
                             [3] + (np.ones(NoOfCols - 2) * 2.).tolist() + [3],
                             0.255, 4.037, fpc_core, fy, fu, 3, 4., NoOfCols, 3)
NoOfCols = 6
BarSize = 8.0
Ag = ( (NoOfCols - 1) * Long_Spacing + 6 ) * Thickness
Rho = ( NoOfCols * 2 + 2 ) * np.pi * ( BarSize / 2. / 8.) ** 2. / Ag
# print Rho
Section2 = PlanarWallSection(Length, Thickness,
                             (NoOfCols - 1) * Long_Spacing + 6,
                             (NoOfCols - 1) * Long_Spacing + 6, BarSize,
                             [3] + (np.ones(NoOfCols - 2) * 2.).tolist() + [3],
                             [3] + (np.ones(NoOfCols - 2) * 2.).tolist() + [3],
                             0.255, 4.037, fpc_core, fy, fu, 3, 4., 8, 3)
Section3 = PlanarWallSection(Length, Thickness, 0, 0, 10.173,
                             [],
                             [],
                             0.255, 4.037, fpc_core, fy, fu, None, None, None, None)
Sections = [Section1, Section1, Section2, Section2]
S4H08SEA = CreateArchetype()
Archetypes.append(S4H08SEA)
Name = 'S4H08SEAWB'
NoOfStories = 6
Sections = [
            Section1, Section1,
            Section1, Section1, Section2, Section2
            ]
S4H08SEAWB = CreateArchetype(Basements2Levels)
Archetypes.append(S4H08SEAWB)
#endregion

#region Archetype S8H08SEA and S8H08SEAWB
Name = 'S8H08SEA'
# print 'Importing Archetype: ' + Name
#### Input Variables
NoOfStories = 8
Thickness = 14.
Length = 16. * 12.
Flange_Thickness = 8*12. # Assume 6' Long Core
Long_Spacing = 4
BarSize = 7.0
Rho = 0.9 #In Fraction
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section1 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, 3., 4., Spacing, Thickness - 3.5)
BarSize = 5.0
Rho = 0.55 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section2 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, 3., 4., Spacing*2, Thickness - 3.5)
BarSize = 4.0
Rho = 0.25 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section3 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, None, None, None, None)
Sections = [Section1, Section1, Section1,
            Section2, Section2, Section2,
            Section3, Section3,
            ]
S8H08SEA = CreateArchetype()
Archetypes.append(S8H08SEA)
Name = 'S8H08SEAWB'
NoOfStories = 11
Sections = [
            Section1, Section1, Section1,
            Section1, Section1, Section1,
            Section2, Section2, Section2,
            Section3, Section3,
            ]
S8H08SEAWB = CreateArchetype(Basements3Levels)
Archetypes.append(S8H08SEAWB)
#endregion

#region Archetype S12H08SEA and S12H08SEAWB
Name = 'S12H08SEA'
# print 'Importing Archetype: ' + Name
#### Input Variables
NoOfStories = 12
Thickness = 14.
Length = 20. * 12.
Flange_Thickness = 10.0*12. # Assume 6' Long Core
Long_Spacing = 4
BarSize = 5.0
Rho = 0.50 #In Fraction
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section1 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, 3., 4., Spacing, Thickness - 3.5)
BarSize = 5.0
Rho = 0.50 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section2 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, 3., 4., Spacing*2, Thickness - 3.5)
ThicknessBelow = float(Thickness)
Thickness = 14.
Length = Length - (ThicknessBelow - Thickness) * 2.
BarSize = 4.0
Rho = 0.35 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section3 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu,  None, None, None, None)
BarSize = 4.0
Rho = 0.25 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section4 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, None, None, None, None)
Sections = [Section1, Section1, Section1,
            Section2, Section2, Section2,
            Section3, Section3, Section3,
            Section4, Section4, Section4,
            ]
S12H08SEA = CreateArchetype()
Archetypes.append(S12H08SEA)
Name = 'S12H08SEAWB'
NoOfStories = 16
Sections = [
            Section1, Section1, Section1, Section1,
            Section1, Section1, Section1,
            Section2, Section2, Section2,
            Section3, Section3, Section3,
            Section4, Section4, Section4,
            ]
S12H08SEAWB = CreateArchetype(Basements)
Archetypes.append(S12H08SEAWB)
#endregion

#region Archetype S16H08SEA and S16H08SEAWB
Name = 'S16H08SEA'
#### Input Variables
NoOfStories = 16
Thickness = 14.
Length = 22. * 12.
Flange_Thickness = 11.*12. # Assume 6' Long Core
Long_Spacing = 4
BarSize = 6.0
Rho = 0.5 #In Fraction
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section1 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, 3., 4., Spacing, Thickness - 3.5)
BarSize = 5.0
Rho = 0.5 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section2 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, 3., 4., Spacing*2, Thickness - 3.5)
ThicknessBelow = float(Thickness)
Thickness = 14.
Length = Length - (ThicknessBelow - Thickness) * 2.
BarSize = 4.0
Rho = 0.25 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section3 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu,  None, None, None, None)
BarSize = 4.0
Rho = 0.25 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section4 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, None, None, None, None)
Sections = [Section1, Section1, Section1, Section1,
            Section2, Section2, Section2, Section2,
            Section3, Section3, Section3, Section3,
            Section4, Section4, Section4, Section4,
            ]
S16H08SEA = CreateArchetype()
Archetypes.append(S16H08SEA)
Name = 'S16H08SEAWB'
NoOfStories = 20 # Include Basement Floors Here
Sections = [Section1, Section1, Section1, Section1,
            Section1, Section1, Section1, Section1,
            Section2, Section2, Section2, Section2,
            Section3, Section3, Section3, Section3,
            Section4, Section4, Section4, Section4,
            ]
S16H08SEAWB = CreateArchetype(Basements)
Archetypes.append(S16H08SEAWB)
#endregion

#region Archetype S20H08SEA and S20H08SEAWB
Name = 'S20H08SEA'
# print 'Importing Archetype: ' + Name
#### Input Variables
NoOfStories = 20
Thickness = 14.
Length = 24. * 12.
Flange_Thickness = 12*12. # Assume 6' Long Core
Long_Spacing = 4
BarSize = 6.0
Rho = 0.5 #In Fraction
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section1 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, 3., 4., Spacing, Thickness - 3.5)
BarSize = 5.0
Rho = 0.5 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section2 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, 3., 4., Spacing*2, Thickness - 3.5)
ThicknessBelow = float(Thickness)
Thickness = 14.
Length = Length - (ThicknessBelow - Thickness) * 2.
BarSize = 4.0
Rho =0.35 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section3 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, None, None, None, None)
BarSize = 4.0
Rho = 0.25 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section4 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, None, None, None, None)
ThicknessBelow = float(Thickness)
Thickness = 14.
Length = Length - (ThicknessBelow - Thickness) * 2.
BarSize = 4.0
Rho = 0.25 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section5 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, None, None, None, None)
Sections = [Section1, Section1, Section1, Section1,
            Section2, Section2, Section2, Section2,
            Section3, Section3, Section3, Section3,
            Section4, Section4, Section4, Section4,
            Section5, Section5, Section5, Section5,
            ]
S20H08SEA = CreateArchetype()
Archetypes.append(S20H08SEA)
Name = 'S20H08SEAWB'
NoOfStories = 24 # Include Basement Floors Here
Sections = [Section1, Section1, Section1, Section1,
            Section1, Section1, Section1, Section1,
            Section2, Section2, Section2, Section2,
            Section3, Section3, Section3, Section3,
            Section4, Section4, Section4, Section4,
            Section5, Section5, Section5, Section5,
            ]
S20H08SEAWB = CreateArchetype(Basements)
Archetypes.append(S20H08SEAWB)
#endregion

#region Archetype S24H08SEA and S24H08SEAWB
Name = 'S24H08SEA'
# print 'Importing Archetype: ' + Name
#### Input Variables
NoOfStories = 24
Thickness = 18
Length = 26. * 12.
Flange_Thickness = 13.*12. # Assume 6' Long Core
Long_Spacing = 4
BarSize = 7.0
Rho = 1.0 #In Fraction
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section1 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, 3., 4., Spacing, Thickness - 3.5)
BarSize = 6.0
Rho = 0.75 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section2 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, 3., 4., Spacing*2., Thickness - 3.5)
ThicknessBelow = float(Thickness)
Thickness = 18.
Length = Length - (ThicknessBelow - Thickness) * 2.
BarSize = 6.0
Rho = 0.60 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section3 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, 3., 4., Spacing*2., Thickness - 3.5)
BarSize = 5.0
Rho = 0.50 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section4 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, 3., 4., Spacing*2., Thickness - 3.5)
ThicknessBelow = float(Thickness)
Thickness = 14.
Length = Length - (ThicknessBelow - Thickness) * 2.
BarSize = 5.0
Rho = 0.50 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section5 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, 3., 4., Spacing*2., Thickness - 3.5)
BarSize = 5.0
Rho = 0.50 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section6 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, None, None, None, None)
Sections = [Section1, Section1, Section1, Section1,
            Section2, Section2, Section2, Section2,
            Section3, Section3, Section3, Section3,
            Section4, Section4, Section4, Section4,
            Section5, Section5, Section5, Section5,
            Section6, Section6, Section6, Section6,
            ]
S24H08SEA = CreateArchetype()
Archetypes.append(S24H08SEA)
Name = 'S24H08SEAWB'
NoOfStories = 28
Sections = [Section1, Section1, Section1, Section1,
            Section1, Section1, Section1, Section1,
            Section2, Section2, Section2, Section2,
            Section3, Section3, Section3, Section3,
            Section4, Section4, Section4, Section4,
            Section5, Section5, Section5, Section5,
            Section6, Section6, Section6, Section6,
            ]
S24H08SEAWB = CreateArchetype(Basements)
Archetypes.append(S24H08SEAWB)
#endregion

#region Archetype S28H08SEA and S28H08SEAWB
Name = 'S28H08SEA'
# print 'Importing Archetype: ' + Name
#### Input Variables
NoOfStories = 28
Thickness = 18.
Length = 28. * 12.
Flange_Thickness = 14*12. # Assume 6' Long Core
Long_Spacing = 4
BarSize = 7.0
Rho = 0.85 #In Fraction
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section1 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, 3., 4., Spacing, Thickness - 3.5)
BarSize = 6.0
Rho = 0.6 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section2 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu,  3., 4., Spacing*2, Thickness - 3.5)

ThicknessBelow = float(Thickness)
Thickness = 18.
Length = Length - (ThicknessBelow - Thickness) * 2.
BarSize = 6.0
Rho = 0.5 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section3 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu,  3., 4., Spacing*2, Thickness - 3.5)
BarSize = 6.0
Rho = 0.5 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section4 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu,  3., 4., Spacing*2, Thickness - 3.5)

ThicknessBelow = float(Thickness)
Thickness = 16.
Length = Length - (ThicknessBelow - Thickness) * 2.
BarSize = 6.0
Rho = 0.5 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section5 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu,  3., 4., Spacing*2, Thickness - 3.5)
BarSize = 6.0
Rho = 0.5 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section6 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu,  3., 4., Spacing*2, Thickness - 3.5)

ThicknessBelow = float(Thickness)
Thickness = 16.
Length = Length - (ThicknessBelow - Thickness) * 2.
BarSize = 6.0
Rho = 0.5 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section7 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu,  3., 4., Spacing*2, Thickness - 3.5)

Sections = [
    Section1, Section1, Section1, Section1,
    Section2, Section2, Section2, Section2,
    Section3, Section3, Section3, Section3,
    Section4, Section4, Section4, Section4,
    Section5, Section5, Section5, Section5,
    Section6, Section6, Section6, Section6,
    Section7, Section7, Section7, Section7,
            ]
S28H08SEA = CreateArchetype()
Archetypes.append(S28H08SEA)

Name = 'S28H08SEAWB'
NoOfStories = 32
Sections = [
    Section1, Section1, Section1, Section1,
    Section1, Section1, Section1, Section1,
    Section2, Section2, Section2, Section2,
    Section3, Section3, Section3, Section3,
    Section4, Section4, Section4, Section4,
    Section5, Section5, Section5, Section5,
    Section6, Section6, Section6, Section6,
    Section7, Section7, Section7, Section7,
            ]
S28H08SEAWB = CreateArchetype(Basements)
Archetypes.append(S28H08SEAWB)
#endregion

#region Archetype S32H08SEA and S32H08SEAWB
Name = 'S32H08SEA'
# print 'Importing Archetype: ' + Name
#### Input Variables
NoOfStories = 32
Thickness = 20.
Length = 30. * 12.
Flange_Thickness = 15*12. # Assume 6' Long Core
Long_Spacing = 4
BarSize = 6.0
Rho = 0.75 #In Fraction
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section1 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, 3., 4., Spacing, Thickness - 3.5)
BarSize = 5.0
Rho = 0.5 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section2 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, 3., 4., Spacing*2., Thickness - 3.5)
ThicknessBelow = float(Thickness)
Thickness = 20.
Length = Length - (ThicknessBelow - Thickness) * 2.
BarSize = 5.0
Rho = 0.5 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section3 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu,  3., 4., Spacing*2, Thickness - 3.5)

BarSize = 5.0
Rho = 0.5 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section4 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu,  3., 4., Spacing*2, Thickness - 3.5)

ThicknessBelow = float(Thickness)
Thickness = 18.
Length = Length - (ThicknessBelow - Thickness) * 2.
BarSize = 5.0
Rho = 0.5 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section5 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu,  3., 4., Spacing*2, Thickness - 3.5)

BarSize = 5.0
Rho = 0.5 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section6 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu,  3., 4., Spacing*2, Thickness - 3.5)

ThicknessBelow = float(Thickness)
Thickness = 18.
Length = Length - (ThicknessBelow - Thickness) * 2.
BarSize = 5.0
Rho = 0.5 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section7 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu,  3., 4., Spacing*2, Thickness - 3.5)

BarSize = 5.0
Rho = 0.5 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section8 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu,  3., 4., Spacing*2, Thickness - 3.5)

Sections = [
    Section1, Section1, Section1, Section1,
    Section2, Section2, Section2, Section2,
    Section3, Section3, Section3, Section3,
    Section4, Section4, Section4, Section4,
    Section5, Section5, Section5, Section5,
    Section6, Section6, Section6, Section6,
    Section7, Section7, Section7, Section7,
    Section8, Section8, Section8, Section8,
            ]
S32H08SEA = CreateArchetype()
Archetypes.append(S32H08SEA)

Name = 'S32H08SEAWB'
NoOfStories = 36
Sections = [
    Section1, Section1, Section1, Section1,
    Section1, Section1, Section1, Section1,
    Section2, Section2, Section2, Section2,
    Section3, Section3, Section3, Section3,
    Section4, Section4, Section4, Section4,
    Section5, Section5, Section5, Section5,
    Section6, Section6, Section6, Section6,
    Section7, Section7, Section7, Section7,
    Section8, Section8, Section8, Section8,
            ]
S32H08SEAWB = CreateArchetype(Basements)
Archetypes.append(S32H08SEAWB)
#endregion

#region Archetype S36H08SEA and S36H08SEAWB
Name = 'S36H08SEA'
# print 'Importing Archetype: ' + Name
#### Input Variables
NoOfStories = 36
Thickness = 22.
Length = 32. * 12.
Flange_Thickness = 16*12. # Assume 6' Long Core
Long_Spacing = 4
BarSize = 6.0
Rho = 0.6 #In Fraction
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section1 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, 3., 4., Spacing, Thickness - 3.5)
BarSize = 6.0
Rho = 0.5 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section2 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, 3., 4., Spacing*2, Thickness - 3.5)

ThicknessBelow = float(Thickness)
Thickness = 22.
Length = Length - (ThicknessBelow - Thickness) * 2.
BarSize = 6.0
Rho = 0.5 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section3 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, 3., 4., Spacing*2, Thickness - 3.5)

BarSize = 6.0
Rho = 0.5 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section4 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, 3., 4., Spacing*2, Thickness - 3.5)

ThicknessBelow = float(Thickness)
Thickness = 16.
Length = Length - (ThicknessBelow - Thickness) * 2.
BarSize = 6.0
Rho = 0.5 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section5 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu,  3., 4., Spacing*2, Thickness - 3.5)

BarSize = 6.0
Rho = 0.5 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section6 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu,  3., 4., Spacing*2, Thickness - 3.5)

ThicknessBelow = float(Thickness)
Thickness = 16.
Length = Length - (ThicknessBelow - Thickness) * 2.
BarSize = 6.0
Rho = 0.5 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section7 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu,  3., 4., Spacing*2, Thickness - 3.5)

BarSize = 6.0
Rho = 0.5 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section8 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu,  3., 4., Spacing*2, Thickness - 3.5)

ThicknessBelow = float(Thickness)
Thickness = 16.
Length = Length - (ThicknessBelow - Thickness) * 2.
BarSize = 6.0
Rho = 0.5 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section9 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu,  3., 4., Spacing*2, Thickness - 3.5)

Sections = [
    Section1, Section1, Section1, Section1,
    Section2, Section2, Section2, Section2,
    Section3, Section3, Section3, Section3,
    Section4, Section4, Section4, Section4,
    Section5, Section5, Section5, Section5,
    Section6, Section6, Section6, Section6,
    Section7, Section7, Section7, Section7,
    Section8, Section8, Section8, Section8,
    Section9, Section9, Section9, Section9,
            ]
S36H08SEA = CreateArchetype()
Archetypes.append(S36H08SEA)
Name = 'S36H08SEAWB'
NoOfStories = 40
Sections = [
    Section1, Section1, Section1, Section1,
    Section1, Section1, Section1, Section1,
    Section2, Section2, Section2, Section2,
    Section3, Section3, Section3, Section3,
    Section4, Section4, Section4, Section4,
    Section5, Section5, Section5, Section5,
    Section6, Section6, Section6, Section6,
    Section7, Section7, Section7, Section7,
    Section8, Section8, Section8, Section8,
    Section9, Section9, Section9, Section9,
    ]
S36H08SEAWB = CreateArchetype(Basements)
Archetypes.append(S36H08SEAWB)
#endregion

#region Archetype S40H08SEA and S40H08SEAWB
Name = 'S40H08SEA'
# print 'Importing Archetype: ' + Name
#### Input Variables
NoOfStories = 40
Thickness = 24.
Length = 34. * 12.
Flange_Thickness = 17.*12. # Assume 6' Long Core
Long_Spacing = 4
BarSize = 6.0
Rho = 0.6 #In Fraction
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section1 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, 3., 4., Spacing, Thickness - 3.5)
BarSize = 6.0
Rho = 0.6 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section2 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, 3., 4., Spacing*2, Thickness - 3.5)
ThicknessBelow = float(Thickness)
Thickness = 24.
Length = Length - (ThicknessBelow - Thickness) * 2.
BarSize = 6.0
Rho = 0.5 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section3 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu,  3., 4., Spacing*2, Thickness - 3.5)

BarSize = 6.0
Rho = 0.5 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section4 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu,  3., 4., Spacing*2, Thickness - 3.5)

ThicknessBelow = float(Thickness)
Thickness = 18.
Length = Length - (ThicknessBelow - Thickness) * 2.
BarSize = 6.0
Rho = 0.5 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section5 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu,  3., 4., Spacing*2, Thickness - 3.5)

BarSize = 6.0
Rho = 0.5 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section6 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu,  3., 4., Spacing*2, Thickness - 3.5)

ThicknessBelow = float(Thickness)
Thickness = 18.
Length = Length - (ThicknessBelow - Thickness) * 2.
BarSize = 6.0
Rho = 0.5 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section7 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu,  3., 4., Spacing*2, Thickness - 3.5)

BarSize = 6.0
Rho = 0.5 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section8 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu,  3., 4., Spacing*2, Thickness - 3.5)

ThicknessBelow = float(Thickness)
Thickness = 16.
Length = Length - (ThicknessBelow - Thickness) * 2.
BarSize = 6.0
Rho = 0.5 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section9 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu,  3., 4., Spacing*2, Thickness - 3.5)

BarSize = 6.0
Rho = 0.5 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section10 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu,  3., 4., Spacing*2, Thickness - 3.5)

Sections = [
    Section1, Section1, Section1, Section1,
    Section2, Section2, Section2, Section2,
    Section3, Section3, Section3, Section3,
    Section4, Section4, Section4, Section4,
    Section5, Section5, Section5, Section5,
    Section6, Section6, Section6, Section6,
    Section7, Section7, Section7, Section7,
    Section8, Section8, Section8, Section8,
    Section9, Section9, Section9, Section9,
    Section10, Section10, Section10, Section10,
            ]
S40H08SEA = CreateArchetype()
Archetypes.append(S40H08SEA)
Name = 'S40H08SEAWB'
NoOfStories = 44
Sections = [
    Section1, Section1, Section1, Section1,
    Section1, Section1, Section1, Section1,
    Section2, Section2, Section2, Section2,
    Section3, Section3, Section3, Section3,
    Section4, Section4, Section4, Section4,
    Section5, Section5, Section5, Section5,
    Section6, Section6, Section6, Section6,
    Section7, Section7, Section7, Section7,
    Section8, Section8, Section8, Section8,
    Section9, Section9, Section9, Section9,
    Section10, Section10, Section10, Section10,
            ]
S40H08SEAWB = CreateArchetype(Basements)
Archetypes.append(S40H08SEAWB)
#endregion

# 2014 Maps

#region Archetype S4H14SEA and S4H14SEAWB
Name = 'S4H14SEA'
# print 'Importing Archetype: ' + Name
# Compute Seismic Weight
NoOfStories = 4
YGrids = [0] + np.array(np.arange(0,(NoOfStories)*13*12, 13*12)+15*12).tolist()
DeadLoads = np.ones(NoOfStories) * DL / 1000.
DeadLoads[-1] = DeadLoads[-1] * DL_Roof / DL
LiveLoads = np.ones(NoOfStories) * LL / 1000.
LiveLoads[-1] = LiveLoads[-1] * LL_Roof / LL
MassPerSqFt = DL / 1000.
Mass = np.ones(NoOfStories) * MassPerSqFt * FloorArea
Mass[-1] = FloorArea * DL_Roof / 1000. # Adjust for Roof Weight
WallTribArea = FloorArea * 0.5
WeightPerSqFt = DL
BuildingWeight = np.ones(NoOfStories) * WeightPerSqFt * FloorArea
BuildingWeight[-1] = 152. / 1000. * FloorArea # Adjust for Roof Weight
# Seismic Hazard
R = 6; Cd = 5
SaDesign, Sds, CuTa = GetSeattle2008Hazard(YGrids[-1], R=R)
Thickness = 18.
Length = 16. * 12.
Long_Spacing = 4
NoOfCols = 13
BarSize = 8.
Ag = ( (NoOfCols - 1) * Long_Spacing + 6 ) * Thickness
Rho = ( NoOfCols  * 2 + 2 ) * np.pi * ( BarSize / 2. / 8.) ** 2. / Ag
# print Rho
Section1 = PlanarWallSection(Length, Thickness,
                             (NoOfCols - 1) * Long_Spacing + 6,
                             (NoOfCols - 1) * Long_Spacing + 6, BarSize,
                             [3] + (np.ones(NoOfCols - 2) * 2.).tolist() + [3],
                             [3] + (np.ones(NoOfCols - 2) * 2.).tolist() + [3],
                             0.255, 4.037, fpc_core, fy, fu, 3, 4., NoOfCols, 3)
NoOfCols = 8
BarSize = 8.0
Ag = ( (NoOfCols - 1) * Long_Spacing + 6 ) * Thickness
Rho = ( NoOfCols * 2 + 2 ) * np.pi * ( BarSize / 2. / 8.) ** 2. / Ag
# print Rho
Section2 = PlanarWallSection(Length, Thickness,
                             (NoOfCols - 1) * Long_Spacing + 6,
                             (NoOfCols - 1) * Long_Spacing + 6, BarSize,
                             [3] + (np.ones(NoOfCols - 2) * 2.).tolist() + [3],
                             [3] + (np.ones(NoOfCols - 2) * 2.).tolist() + [3],
                             0.255, 4.037, fpc_core, fy, fu, 3, 4., 8, 3)
Section3 = PlanarWallSection(Length, Thickness, 0, 0, 10.173,
                             [],
                             [],
                             0.255, 4.037, fpc_core, fy, fu, None, None, None, None)
Sections = [Section1, Section1, Section2, Section2]
S4H14SEA = CreateArchetype(Use2008Maps = False)
Archetypes.append(S4H14SEA)
Name = 'S4H14SEAWB'
NoOfStories = 6
Sections = [
            Section1, Section1,
            Section1, Section1, Section2, Section2
            ]
S4H14SEAWB = CreateArchetype(Basements2Levels, Use2008Maps = False)
Archetypes.append(S4H14SEAWB)
#endregion

#region Archetype S8H14SEA and S8H14SEAWB
Name = 'S8H14SEA'
# print 'Importing Archetype: ' + Name
#### Input Variables
NoOfStories = 8
Thickness = 16.
Length = 18. * 12.
Flange_Thickness = 9.*12. # Assume 6' Long Core
Long_Spacing = 4
BarSize = 7.0
Rho = 0.95 #In Fraction
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section1 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, 3., 4., Spacing, Thickness - 3.5)
BarSize = 5.0
Rho = 0.70 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section2 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, 3., 4., Spacing*2, Thickness - 3.5)
BarSize = 4.0
Rho = 0.25 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section3 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, None, None, None, None)
Sections = [Section1, Section1, Section1,
            Section2, Section2, Section2,
            Section3, Section3,
            ]
S8H14SEA = CreateArchetype(Use2008Maps = False)
Archetypes.append(S8H14SEA)
Name = 'S8H14SEAWB'
NoOfStories = 11
Sections = [
            Section1, Section1, Section1,
            Section1, Section1, Section1,
            Section2, Section2, Section2,
            Section3, Section3,
            ]
S8H14SEAWB = CreateArchetype(Basements3Levels, Use2008Maps = False)
Archetypes.append(S8H14SEAWB)
#endregion

#region Archetype S12H14SEA and S12H14SEAWB
Name = 'S12H14SEA'
# print 'Importing Archetype: ' + Name
#### Input Variables
NoOfStories = 12
Thickness = 18.
Length = 20. * 12.
Flange_Thickness = 10.*12. # Assume 6' Long Core
Long_Spacing = 4
BarSize = 6.0
Rho = 0.85 #In Fraction
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section1 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, 3., 4., Spacing, Thickness - 3.5)
BarSize = 6.0
Rho = 0.6 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section2 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, 3., 4., Spacing*2, Thickness - 3.5)
ThicknessBelow = float(Thickness)
Thickness = 18.
Length = Length - (ThicknessBelow - Thickness) * 2.
BarSize = 4.0
Rho = 0.40 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section3 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu,  None, None, None, None)
BarSize = 4.0
Rho = 0.25 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section4 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, None, None, None, None)
Sections = [Section1, Section1, Section1,
            Section2, Section2, Section2,
            Section3, Section3, Section3,
            Section4, Section4, Section4,
            ]
S12H14SEA = CreateArchetype(Use2008Maps = False)
Archetypes.append(S12H14SEA)
Name = 'S12H14SEAWB'
NoOfStories = 16
Sections = [
            Section1, Section1, Section1, Section1,
            Section1, Section1, Section1,
            Section2, Section2, Section2,
            Section3, Section3, Section3,
            Section4, Section4, Section4,
            ]
S12H14SEAWB = CreateArchetype(Basements, Use2008Maps = False)
Archetypes.append(S12H14SEAWB)
#endregion

#region Archetype S16H14SEA and S16H14SEAWB
Name = 'S16H14SEA'
#### Input Variables
NoOfStories = 16
Thickness = 22.
Length = 24. * 12.
Flange_Thickness = 12.*12. # Assume 6' Long Core
Long_Spacing = 4
BarSize = 7.0
Rho = 0.6 #In Fraction
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section1 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, 3., 4., Spacing, Thickness - 3.5)
BarSize = 6.0
Rho = 0.5 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section2 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, 3., 4., Spacing*2, Thickness - 3.5)
ThicknessBelow = float(Thickness)
Thickness = 22.
Length = Length - (ThicknessBelow - Thickness) * 2.
BarSize = 5.0
Rho = 0.40 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section3 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu,  None, None, None, None)
BarSize = 4.0
Rho = 0.25 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section4 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, None, None, None, None)
Sections = [Section1, Section1, Section1, Section1,
            Section2, Section2, Section2, Section2,
            Section3, Section3, Section3, Section3,
            Section4, Section4, Section4, Section4,
            ]
S16H14SEA = CreateArchetype(Use2008Maps = False)
Archetypes.append(S16H14SEA)
Name = 'S16H14SEAWB'
NoOfStories = 20 # Include Basement Floors Here
Sections = [Section1, Section1, Section1, Section1,
            Section1, Section1, Section1, Section1,
            Section2, Section2, Section2, Section2,
            Section3, Section3, Section3, Section3,
            Section4, Section4, Section4, Section4,
            ]
S16H14SEAWB = CreateArchetype(Basements, False)
Archetypes.append(S16H14SEAWB)
#endregion

#region Archetype S20H14SEA and S20H14SEAWB
Name = 'S20H14SEA'
# print 'Importing Archetype: ' + Name
#### Input Variables
NoOfStories = 20
Thickness = 24.
Length = 26. * 12.
Flange_Thickness = 13*12. # Assume 6' Long Core
Long_Spacing = 4
BarSize = 7.0
Rho = 0.55 #In Fraction
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section1 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, 3., 4., Spacing, Thickness - 3.5)
BarSize = 5.0
Rho = 0.5 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section2 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, 3., 4., Spacing*2, Thickness - 3.5)
ThicknessBelow = float(Thickness)
Thickness = 24.
Length = Length - (ThicknessBelow - Thickness) * 2.
BarSize = 5.0
Rho =0.45 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section3 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, None, None, None, None)

BarSize = 4.0
Rho = 0.25 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section4 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, None, None, None, None)
ThicknessBelow = float(Thickness)
Thickness = 20.
Length = Length - (ThicknessBelow - Thickness) * 2.
BarSize = 4.0
Rho = 0.25 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section5 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, None, None, None, None)
Sections = [Section1, Section1, Section1, Section1,
            Section2, Section2, Section2, Section2,
            Section3, Section3, Section3, Section3,
            Section4, Section4, Section4, Section4,
            Section5, Section5, Section5, Section5,
            ]
S20H14SEA = CreateArchetype(Use2008Maps = False)
Archetypes.append(S20H14SEA)
Name = 'S20H14SEAWB'
NoOfStories = 24 # Include Basement Floors Here
Sections = [Section1, Section1, Section1, Section1,
            Section1, Section1, Section1, Section1,
            Section2, Section2, Section2, Section2,
            Section3, Section3, Section3, Section3,
            Section4, Section4, Section4, Section4,
            Section5, Section5, Section5, Section5,
            ]
S20H14SEAWB = CreateArchetype(Basements, False)
Archetypes.append(S20H14SEAWB)
#endregion

#region Archetype S24H14SEA and S24H14SEAWB
Name = 'S24H14SEA'
# print 'Importing Archetype: ' + Name
#### Input Variables
NoOfStories = 24
Thickness = 26.
Length = 28. * 12.
Flange_Thickness = 14.*12. # Assume 6' Long Core
Long_Spacing = 4
BarSize = 8.0
Rho = 1.1 #In Fraction
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section1 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, 3., 4., Spacing, Thickness - 3.5)
BarSize = 7.0
Rho = 0.75 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section2 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu,  3., 4., Spacing*2, Thickness - 3.5)

ThicknessBelow = float(Thickness)
Thickness = 26.
Length = Length - (ThicknessBelow - Thickness) * 2.
BarSize = 6.0
Rho = 0.6 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section3 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu,  3., 4., Spacing*2, Thickness - 3.5)

BarSize = 6.0
Rho = 0.5 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section4 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu,  3., 4., Spacing*2, Thickness - 3.5)

ThicknessBelow = float(Thickness)
Thickness = 22.
Length = Length - (ThicknessBelow - Thickness) * 2.
BarSize = 6.0
Rho = 0.5 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section5 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu,  3., 4., Spacing*2, Thickness - 3.5)

BarSize = 6.0
Rho = 0.5 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section6 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu,  3., 4., Spacing*2, Thickness - 3.5)

Sections = [Section1, Section1, Section1, Section1,
            Section2, Section2, Section2, Section2,
            Section3, Section3, Section3, Section3,
            Section4, Section4, Section4, Section4,
            Section5, Section5, Section5, Section5,
            Section6, Section6, Section6, Section6,
            ]
S24H14SEA = CreateArchetype(Use2008Maps = False)
Archetypes.append(S24H14SEA)
Name = 'S24H14SEAWB'
NoOfStories = 28
Sections = [Section1, Section1, Section1, Section1,
            Section1, Section1, Section1, Section1,
            Section2, Section2, Section2, Section2,
            Section3, Section3, Section3, Section3,
            Section4, Section4, Section4, Section4,
            Section5, Section5, Section5, Section5,
            Section6, Section6, Section6, Section6,
            ]
S24H14SEAWB = CreateArchetype(Basements, False)
Archetypes.append(S24H14SEAWB)
#endregion

#region Archetype S28H14SEA and S28H14SEAWB
Name = 'S28H14SEA'
# print 'Importing Archetype: ' + Name
#### Input Variables
NoOfStories = 28
Thickness = 28.
Length = 30. * 12.
Flange_Thickness = 15*12. # Assume 6' Long Core
Long_Spacing = 4
BarSize = 7.0
Rho = 0.95 #In Fraction
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section1 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, 3., 4., Spacing, Thickness - 3.5)
BarSize = 7.0
Rho = 0.7 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section2 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, 3., 4., Spacing*2, Thickness - 3.5)
ThicknessBelow = float(Thickness)
Thickness = 28.
Length = Length - (ThicknessBelow - Thickness) * 2.
BarSize = 6.0
Rho = 0.6 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section3 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, 3., 4., Spacing*2, Thickness - 3.5)

BarSize = 6.0
Rho = 0.5 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section4 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu,  3., 4., Spacing*2, Thickness - 3.5)

ThicknessBelow = float(Thickness)
Thickness = 24.
Length = Length - (ThicknessBelow - Thickness) * 2.
BarSize = 6.0
Rho = 0.5 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section5 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu,  3., 4., Spacing*2, Thickness - 3.5)

BarSize = 6.0
Rho = 0.5 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section6 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu,  3., 4., Spacing*2, Thickness - 3.5)

BarSize = 6.0
Rho = 0.5 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section7 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu,  3., 4., Spacing*2, Thickness - 3.5)

Sections = [
    Section1, Section1, Section1, Section1,
    Section2, Section2, Section2, Section2,
    Section3, Section3, Section3, Section3,
    Section4, Section4, Section4, Section4,
    Section5, Section5, Section5, Section5,
    Section6, Section6, Section6, Section6,
    Section7, Section7, Section7, Section7,
            ]
S28H14SEA = CreateArchetype(Use2008Maps = False)
Archetypes.append(S28H14SEA)

Name = 'S28H14SEAWB'
NoOfStories = 32
Sections = [
    Section1, Section1, Section1, Section1,
    Section1, Section1, Section1, Section1,
    Section2, Section2, Section2, Section2,
    Section3, Section3, Section3, Section3,
    Section4, Section4, Section4, Section4,
    Section5, Section5, Section5, Section5,
    Section6, Section6, Section6, Section6,
    Section7, Section7, Section7, Section7,
            ]
S28H14SEAWB = CreateArchetype(Basements, False)
Archetypes.append(S28H14SEAWB)
#endregion

#region Archetype S32H14SEA and S32H14SEAWB
Name = 'S32H14SEA'
# print 'Importing Archetype: ' + Name
#### Input Variables
NoOfStories = 32
Thickness = 30.
Length = 32. * 12.
Flange_Thickness = 16.*12. # Assume 6' Long Core
Long_Spacing = 4
BarSize = 7.0
Rho = 0.95 #In Fraction
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section1 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, 3., 4., Spacing, Thickness - 3.5)
BarSize = 7.0
Rho = 0.8 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section2 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, 3., 4., Spacing*2, Thickness - 3.5)
ThicknessBelow = float(Thickness)
Thickness = 30.
Length = Length - (ThicknessBelow - Thickness) * 2.
BarSize = 7.0
Rho = 0.7 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section3 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, 3., 4., Spacing*2, Thickness - 3.5)
BarSize = 5.0
Rho = 0.5 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section4 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, 3., 4., Spacing*2, Thickness - 3.5)
ThicknessBelow = float(Thickness)
Thickness = 26.
Length = Length - (ThicknessBelow - Thickness) * 2.
BarSize = 6.0
Rho = 0.5 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section5 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu,  3., 4., Spacing*2, Thickness - 3.5)

BarSize = 6.0
Rho = 0.5 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section6 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu,  3., 4., Spacing*2, Thickness - 3.5)

ThicknessBelow = float(Thickness)
Thickness = 26.
Length = Length - (ThicknessBelow - Thickness) * 2.
BarSize = 6.0
Rho = 0.5 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section7 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu,  3., 4., Spacing*2, Thickness - 3.5)

BarSize = 6.0
Rho = 0.5 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section8 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu,  3., 4., Spacing*2, Thickness - 3.5)
Sections = [
    Section1, Section1, Section1, Section1,
    Section2, Section2, Section2, Section2,
    Section3, Section3, Section3, Section3,
    Section4, Section4, Section4, Section4,
    Section5, Section5, Section5, Section5,
    Section6, Section6, Section6, Section6,
    Section7, Section7, Section7, Section7,
    Section8, Section8, Section8, Section8,
            ]
S32H14SEA = CreateArchetype(Use2008Maps = False)
Archetypes.append(S32H14SEA)

Name = 'S32H14SEAWB'
NoOfStories = 36
Sections = [
    Section1, Section1, Section1, Section1,
    Section1, Section1, Section1, Section1,
    Section2, Section2, Section2, Section2,
    Section3, Section3, Section3, Section3,
    Section4, Section4, Section4, Section4,
    Section5, Section5, Section5, Section5,
    Section6, Section6, Section6, Section6,
    Section7, Section7, Section7, Section7,
    Section8, Section8, Section8, Section8,
            ]
S32H14SEAWB = CreateArchetype(Basements, False)
Archetypes.append(S32H14SEAWB)
#endregion

#region Archetype S36H14SEA and S36H14SEAWB
Name = 'S36H14SEA'
# print 'Importing Archetype: ' + Name
#### Input Variables
NoOfStories = 36
Thickness = 32.
Length = 34. * 12.
Flange_Thickness = 17.*12. # Assume 6' Long Core
Long_Spacing = 4
BarSize = 8.0
Rho = 1.1 #In Fraction
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section1 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, 3., 4., Spacing, Thickness - 3.5)
BarSize = 8.0
Rho = 0.8 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section2 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, 3., 4., Spacing*2, Thickness - 3.5)

ThicknessBelow = float(Thickness)
Thickness = 32.
Length = Length - (ThicknessBelow - Thickness) * 2.
BarSize = 7.0
Rho = 0.7 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section3 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, 3., 4., Spacing*2, Thickness - 3.5)

BarSize = 7.0
Rho = 0.6 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section4 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, 3., 4., Spacing*2, Thickness - 3.5)

ThicknessBelow = float(Thickness)
Thickness = 28.
Length = Length - (ThicknessBelow - Thickness) * 2.
BarSize = 6.0
Rho = 0.5 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section5 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu,  3., 4., Spacing*2, Thickness - 3.5)
BarSize = 6.0
Rho = 0.5 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section6 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu,  3., 4., Spacing*2, Thickness - 3.5)
ThicknessBelow = float(Thickness)
Thickness = 28.
Length = Length - (ThicknessBelow - Thickness) * 2.
BarSize = 6.0
Rho = 0.5 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section7 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu,  3., 4., Spacing*2, Thickness - 3.5)
BarSize = 6.0
Rho = 0.5 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section8 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu,  3., 4., Spacing*2, Thickness - 3.5)
BarSize = 6.0
Rho = 0.5 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section9 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu,  3., 4., Spacing*2, Thickness - 3.5)
Sections = [
    Section1, Section1, Section1, Section1,
    Section2, Section2, Section2, Section2,
    Section3, Section3, Section3, Section3,
    Section4, Section4, Section4, Section4,
    Section5, Section5, Section5, Section5,
    Section6, Section6, Section6, Section6,
    Section7, Section7, Section7, Section7,
    Section8, Section8, Section8, Section8,
    Section9, Section9, Section9, Section9,
            ]
S36H14SEA = CreateArchetype(Use2008Maps = False)
Archetypes.append(S36H14SEA)
Name = 'S36H14SEAWB'
NoOfStories = 40
Sections = [
    Section1, Section1, Section1, Section1,
    Section1, Section1, Section1, Section1,
    Section2, Section2, Section2, Section2,
    Section3, Section3, Section3, Section3,
    Section4, Section4, Section4, Section4,
    Section5, Section5, Section5, Section5,
    Section6, Section6, Section6, Section6,
    Section7, Section7, Section7, Section7,
    Section8, Section8, Section8, Section8,
    Section9, Section9, Section9, Section9,
    ]
S36H14SEAWB = CreateArchetype(Basements, False)
Archetypes.append(S36H14SEAWB)
#endregion

#region Archetype S40H14SEA and S40H14SEAWB
Name = 'S40H14SEA'
# print 'Importing Archetype: ' + Name
#### Input Variables
NoOfStories = 40
Thickness = 34.
Length = 36. * 12.
Flange_Thickness = 18.*12. # Assume 6' Long Core
Long_Spacing = 4
BarSize = 8.0
Rho = 1.2 #In Fraction
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section1 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, 3., 4., Spacing, Thickness - 3.5)
BarSize = 8.0
Rho = 1.0 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section2 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, 3., 4., Spacing*2, Thickness - 3.5)
ThicknessBelow = float(Thickness)
Thickness = 34.
Length = Length - (ThicknessBelow - Thickness) * 2.
BarSize = 7.0
Rho = 0.8 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section3 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, 3., 4., Spacing*2, Thickness - 3.5)
BarSize = 6.0
Rho = 0.8 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section4 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, 3., 4., Spacing*2, Thickness - 3.5)
ThicknessBelow = float(Thickness)
Thickness = 28.
Length = Length - (ThicknessBelow - Thickness) * 2.
BarSize = 5.0
Rho = 0.7 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section5 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, 3., 4., Spacing*2, Thickness - 3.5)
BarSize = 6.0
Rho = 0.5 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section6 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu,  3., 4., Spacing*2, Thickness - 3.5)
ThicknessBelow = float(Thickness)
Thickness = 28.
Length = Length - (ThicknessBelow - Thickness) * 2.
BarSize = 6.0
Rho = 0.5 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section7 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu,  3., 4., Spacing*2, Thickness - 3.5)
BarSize = 6.0
Rho = 0.5 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section8 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu,  3., 4., Spacing*2, Thickness - 3.5)
ThicknessBelow = float(Thickness)
Thickness = 24.
Length = Length - (ThicknessBelow - Thickness) * 2.
BarSize = 6.0
Rho = 0.5 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section9 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu,  3., 4., Spacing*2, Thickness - 3.5)
BarSize = 6.0
Rho = 0.5 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section10 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu,  3., 4., Spacing*2, Thickness - 3.5)
Sections = [
    Section1, Section1, Section1, Section1,
    Section2, Section2, Section2, Section2,
    Section3, Section3, Section3, Section3,
    Section4, Section4, Section4, Section4,
    Section5, Section5, Section5, Section5,
    Section6, Section6, Section6, Section6,
    Section7, Section7, Section7, Section7,
    Section8, Section8, Section8, Section8,
    Section9, Section9, Section9, Section9,
    Section10, Section10, Section10, Section10,
            ]
S40H14SEA = CreateArchetype(Use2008Maps = False)
Archetypes.append(S40H14SEA)
Name = 'S40H14SEAWB'
NoOfStories = 44
Sections = [
    Section1, Section1, Section1, Section1,
    Section1, Section1, Section1, Section1,
    Section2, Section2, Section2, Section2,
    Section3, Section3, Section3, Section3,
    Section4, Section4, Section4, Section4,
    Section5, Section5, Section5, Section5,
    Section6, Section6, Section6, Section6,
    Section7, Section7, Section7, Section7,
    Section8, Section8, Section8, Section8,
    Section9, Section9, Section9, Section9,
    Section10, Section10, Section10, Section10,
            ]
S40H14SEAWB = CreateArchetype(Basements, False)
Archetypes.append(S40H14SEAWB)
#endregion

############################### Performance Group #2 ###############################
##### 2008 Maps ######

#region Archetype S4H08SEAPG2 and S4H08SEAWBPG2
Name = 'S4H08SEAPG2'
# print 'Importing Archetype: ' + Name
# Compute Seismic Weight
NoOfStories = 4
YGrids = [0] + np.array(np.arange(0,(NoOfStories)*13*12, 13*12)+15*12).tolist()
DeadLoads = np.ones(NoOfStories) * DL / 1000.
DeadLoads[-1] = DeadLoads[-1] * DL_Roof / DL
LiveLoads = np.ones(NoOfStories) * LL / 1000.
LiveLoads[-1] = LiveLoads[-1] * LL_Roof / LL
MassPerSqFt = DL / 1000.
Mass = np.ones(NoOfStories) * MassPerSqFt * FloorArea
Mass[-1] = FloorArea * DL_Roof / 1000. # Adjust for Roof Weight
WallTribArea = FloorArea * 0.5
WeightPerSqFt = DL
BuildingWeight = np.ones(NoOfStories) * WeightPerSqFt * FloorArea
BuildingWeight[-1] = 152. / 1000. * FloorArea # Adjust for Roof Weight
# Seismic Hazard
R = 6; Cd = 5
SaDesign, Sds, CuTa = GetSeattle2008Hazard(YGrids[-1], R=R)
Thickness = 14.
Length = 10. * 12.
Long_Spacing = 4
NoOfCols = 14
BarSize = 8.
Ag = ( (NoOfCols - 1) * Long_Spacing + 6 ) * Thickness
Rho = ( NoOfCols  * 2 + 2 ) * np.pi * ( BarSize / 2. / 8.) ** 2. / Ag
# print Rho
Section1 = PlanarWallSection(Length, Thickness,
                             (NoOfCols - 1) * Long_Spacing + 6,
                             (NoOfCols - 1) * Long_Spacing + 6, BarSize,
                             [3] + (np.ones(NoOfCols - 2) * 2.).tolist() + [3],
                             [3] + (np.ones(NoOfCols - 2) * 2.).tolist() + [3],
                             0.255, 4.037, fpc_core, fy, fu, 3, 4., NoOfCols, 3)
NoOfCols = 6
BarSize = 8.0
Ag = ( (NoOfCols - 1) * Long_Spacing + 6 ) * Thickness
Rho = ( NoOfCols * 2 + 2 ) * np.pi * ( BarSize / 2. / 8.) ** 2. / Ag
# print Rho
Section2 = PlanarWallSection(Length, Thickness,
                             (NoOfCols - 1) * Long_Spacing + 6,
                             (NoOfCols - 1) * Long_Spacing + 6, BarSize,
                             [3] + (np.ones(NoOfCols - 2) * 2.).tolist() + [3],
                             [3] + (np.ones(NoOfCols - 2) * 2.).tolist() + [3],
                             0.255, 4.037, fpc_core, fy, fu, 3, 4., 8, 3)
Section3 = PlanarWallSection(Length, Thickness, 0, 0, 10.173,
                             [],
                             [],
                             0.255, 4.037, fpc_core, fy, fu, None, None, None, None)
Sections = [Section1, Section1, Section2, Section2]
S4H08SEAPG2 = CreateArchetype()
Archetypes.append(S4H08SEAPG2)
Name = 'S4H08SEAWBPG2'
NoOfStories = 6
Sections = [
            Section1, Section1,
            Section1, Section1, Section2, Section2
            ]
S4H08SEAWBPG2 = CreateArchetype(Basements2Levels)
Archetypes.append(S4H08SEAWBPG2)
#endregion

#region Archetype S8H08SEAPG2 and S8H08SEAWBPG2
Name = 'S8H08SEAPG2'
# print 'Importing Archetype: ' + Name
#### Input Variables
NoOfStories = 8
Thickness = 20.
Length = 11. * 12.
Flange_Thickness = 5.5*12. # Assume 6' Long Core
Long_Spacing = 4
BarSize = 8.0
Rho = 2.0 #In Fraction
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section1 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, 3., 4., Spacing, Thickness - 3.5)
BarSize = 8.0
Rho = 1.1 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section2 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, 3., 4., Spacing*2, Thickness - 3.5)
BarSize = 4.0
Rho = 0.25 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section3 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, None, None, None, None)
Sections = [Section1, Section1, Section1,
            Section2, Section2, Section2,
            Section3, Section3,
            ]
S8H08SEAPG2 = CreateArchetype()
Archetypes.append(S8H08SEAPG2)
Name = 'S8H08SEAWBPG2'
NoOfStories = 11
Sections = [
            Section1, Section1, Section1,
            Section1, Section1, Section1,
            Section2, Section2, Section2,
            Section3, Section3,
            ]
S8H08SEAWBPG2 = CreateArchetype(Basements3Levels)
Archetypes.append(S8H08SEAWBPG2)
#endregion

#region Archetype S12H08SEAPG2 and S12H08SEAWBPG2
Name = 'S12H08SEAPG2'
# print 'Importing Archetype: ' + Name
#### Input Variables
NoOfStories = 12
Thickness = 20.
Length = 14. * 12.
Flange_Thickness = 7*12. # Assume 6' Long Core
Long_Spacing = 4
BarSize = 8.0
Rho = 1.6 #In Fraction
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section1 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, 3., 4., Spacing, Thickness - 3.5)
BarSize = 8.0
Rho = 1.0 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section2 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, 3., 4., Spacing*2, Thickness - 3.5)
ThicknessBelow = float(Thickness)
Thickness = 16.
Length = Length - (ThicknessBelow - Thickness) * 2.
BarSize = 4.0
Rho = 0.45 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section3 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu,  None, None, None, None)
BarSize = 4.0
Rho = 0.25 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section4 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, None, None, None, None)
Sections = [Section1, Section1, Section1,
            Section2, Section2, Section2,
            Section3, Section3, Section3,
            Section4, Section4, Section4,
            ]
S12H08SEAPG2 = CreateArchetype()
Archetypes.append(S12H08SEAPG2)
Name = 'S12H08SEAWBPG2'
NoOfStories = 16
Sections = [
            Section1, Section1, Section1, Section1,
            Section1, Section1, Section1,
            Section2, Section2, Section2,
            Section3, Section3, Section3,
            Section4, Section4, Section4,
            ]
S12H08SEAWBPG2 = CreateArchetype(Basements)
Archetypes.append(S12H08SEAWBPG2)
#endregion

#region Archetype S16H08SEAPG2 and S16H08SEAWBPG2
Name = 'S16H08SEAPG2'
#### Input Variables
NoOfStories = 16
Thickness = 22.
Length = 16. * 12.
Flange_Thickness = 8.*12. # Assume 6' Long Core
Long_Spacing = 4
BarSize = 8.0
Rho = 1.4 #In Fraction
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section1 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, 3., 4., Spacing, Thickness - 3.5)
BarSize = 8.0
Rho = 1.0 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section2 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, 3., 4., Spacing*2, Thickness - 3.5)
ThicknessBelow = float(Thickness)
Thickness = 16.
Length = Length - (ThicknessBelow - Thickness) * 2.
BarSize = 4.0
Rho = 0.35 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section3 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu,  None, None, None, None)
BarSize = 4.0
Rho = 0.25 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section4 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, None, None, None, None)
Sections = [Section1, Section1, Section1, Section1,
            Section2, Section2, Section2, Section2,
            Section3, Section3, Section3, Section3,
            Section4, Section4, Section4, Section4,
            ]
S16H08SEAPG2 = CreateArchetype()
Archetypes.append(S16H08SEAPG2)
Name = 'S16H08SEAWBPG2'
NoOfStories = 20 # Include Basement Floors Here
Sections = [Section1, Section1, Section1, Section1,
            Section1, Section1, Section1, Section1,
            Section2, Section2, Section2, Section2,
            Section3, Section3, Section3, Section3,
            Section4, Section4, Section4, Section4,
            ]
S16H08SEAWBPG2 = CreateArchetype(Basements)
Archetypes.append(S16H08SEAWBPG2)
#endregion

#region Archetype S20H08SEAPG2 and S20H08SEAWBPG2
Name = 'S20H08SEAPG2'
# print 'Importing Archetype: ' + Name
#### Input Variables
NoOfStories = 20
Thickness = 24.
Length = 18. * 12.
Flange_Thickness = 9*12. # Assume 6' Long Core
Long_Spacing = 4
BarSize = 8.0
Rho = 1.2 #In Fraction
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section1 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, 3., 4., Spacing, Thickness - 3.5)
BarSize = 8.0
Rho = 0.9 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section2 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, 3., 4., Spacing*2, Thickness - 3.5)
ThicknessBelow = float(Thickness)
Thickness = 18.
Length = Length - (ThicknessBelow - Thickness) * 2.
BarSize = 4.0
Rho = 0.5 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section3 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, 3., 4., Spacing*2, Thickness - 3.5)
BarSize = 4.0
Rho = 0.25 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section4 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, None, None, None, None)
ThicknessBelow = float(Thickness)
Thickness = 18.
Length = Length - (ThicknessBelow - Thickness) * 2.
BarSize = 4.0
Rho = 0.25 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section5 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, None, None, None, None)
Sections = [Section1, Section1, Section1, Section1,
            Section2, Section2, Section2, Section2,
            Section3, Section3, Section3, Section3,
            Section4, Section4, Section4, Section4,
            Section5, Section5, Section5, Section5,
            ]
S20H08SEA = CreateArchetype()
Archetypes.append(S20H08SEA)
Name = 'S20H08SEAWBPG2'
NoOfStories = 24 # Include Basement Floors Here
Sections = [Section1, Section1, Section1, Section1,
            Section1, Section1, Section1, Section1,
            Section2, Section2, Section2, Section2,
            Section3, Section3, Section3, Section3,
            Section4, Section4, Section4, Section4,
            Section5, Section5, Section5, Section5,
            ]
S20H08SEAWBPG2 = CreateArchetype(Basements)
Archetypes.append(S20H08SEAWBPG2)
#endregion

#region Archetype S24H08SEAPG2 and S24H08SEAWBPG2
Name = 'S24H08SEAPG2'
# print 'Importing Archetype: ' + Name
#### Input Variables
NoOfStories = 24
Thickness = 28
Length = 21. * 12.
Flange_Thickness = 10.5*12. # Assume 6' Long Core
Long_Spacing = 4
BarSize = 5.0
Rho = 0.7 #In Fraction
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section1 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, 3., 4., Spacing, Thickness - 3.5)
BarSize = 5.0
Rho = 0.5 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section2 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, 3., 4., Spacing*2., Thickness - 3.5)
ThicknessBelow = float(Thickness)
Thickness = 18.
Length = Length - (ThicknessBelow - Thickness) * 2.
BarSize = 4.0
Rho = 0.5 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section3 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, 3., 4., Spacing*2., Thickness - 3.5)
BarSize = 4.0
Rho = 0.25 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section4 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu,   None, None, None, None)
ThicknessBelow = float(Thickness)
Thickness = 18.
Length = Length - (ThicknessBelow - Thickness) * 2.
BarSize = 4.0
Rho = 0.25 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section5 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, None, None, None, None)
BarSize = 4.0
Rho = 0.25 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section6 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, None, None, None, None)
Sections = [Section1, Section1, Section1, Section1,
            Section2, Section2, Section2, Section2,
            Section3, Section3, Section3, Section3,
            Section4, Section4, Section4, Section4,
            Section5, Section5, Section5, Section5,
            Section6, Section6, Section6, Section6,
            ]
S24H08SEAPG2 = CreateArchetype()
Archetypes.append(S24H08SEAPG2)
Name = 'S24H08SEAWBPG2'
NoOfStories = 28
Sections = [Section1, Section1, Section1, Section1,
            Section1, Section1, Section1, Section1,
            Section2, Section2, Section2, Section2,
            Section3, Section3, Section3, Section3,
            Section4, Section4, Section4, Section4,
            Section5, Section5, Section5, Section5,
            Section6, Section6, Section6, Section6,
            ]
S24H08SEAWBPG2 = CreateArchetype(Basements)
Archetypes.append(S24H08SEAWBPG2)
#endregion

##### 2014 Maps ######

#region Archetype S4H14SEAPG2 and S4H14SEAWBPG2
Name = 'S4H14SEAPG2'
# print 'Importing Archetype: ' + Name
# Compute Seismic Weight
NoOfStories = 4
YGrids = [0] + np.array(np.arange(0,(NoOfStories)*13*12, 13*12)+15*12).tolist()
DeadLoads = np.ones(NoOfStories) * DL / 1000.
DeadLoads[-1] = DeadLoads[-1] * DL_Roof / DL
LiveLoads = np.ones(NoOfStories) * LL / 1000.
LiveLoads[-1] = LiveLoads[-1] * LL_Roof / LL
MassPerSqFt = DL / 1000.
Mass = np.ones(NoOfStories) * MassPerSqFt * FloorArea
Mass[-1] = FloorArea * DL_Roof / 1000. # Adjust for Roof Weight
WallTribArea = FloorArea * 0.5
WeightPerSqFt = DL
BuildingWeight = np.ones(NoOfStories) * WeightPerSqFt * FloorArea
BuildingWeight[-1] = 152. / 1000. * FloorArea # Adjust for Roof Weight
# Seismic Hazard
R = 6; Cd = 5
SaDesign, Sds, CuTa = GetSeattle2008Hazard(YGrids[-1], R=R)
Thickness = 18.
Length = 12. * 12.
Long_Spacing = 4
NoOfCols = 12
BarSize = 10.
Ag = ( (NoOfCols - 1) * Long_Spacing + 6 ) * Thickness
Rho = ( NoOfCols  * 2 + 2 ) * np.pi * ( BarSize / 2. / 8.) ** 2. / Ag
# print Rho
Section1 = PlanarWallSection(Length, Thickness,
                             (NoOfCols - 1) * Long_Spacing + 6,
                             (NoOfCols - 1) * Long_Spacing + 6, BarSize,
                             [3] + (np.ones(NoOfCols - 2) * 2.).tolist() + [3],
                             [3] + (np.ones(NoOfCols - 2) * 2.).tolist() + [3],
                             0.255, 4.037, fpc_core, fy, fu, 3, 4., NoOfCols, 3)
NoOfCols = 10
BarSize = 8.0
Ag = ( (NoOfCols - 1) * Long_Spacing + 6 ) * Thickness
Rho = ( NoOfCols * 2 + 2 ) * np.pi * ( BarSize / 2. / 8.) ** 2. / Ag
# print Rho
Section2 = PlanarWallSection(Length, Thickness,
                             (NoOfCols - 1) * Long_Spacing + 6,
                             (NoOfCols - 1) * Long_Spacing + 6, BarSize,
                             [3] + (np.ones(NoOfCols - 2) * 2.).tolist() + [3],
                             [3] + (np.ones(NoOfCols - 2) * 2.).tolist() + [3],
                             0.255, 4.037, fpc_core, fy, fu, 3, 4., 8, 3)
Section3 = PlanarWallSection(Length, Thickness, 0, 0, 10.173,
                             [],
                             [],
                             0.255, 4.037, fpc_core, fy, fu, None, None, None, None)
Sections = [Section1, Section1, Section2, Section2]
S4H14SEAPG2 = CreateArchetype(Use2008Maps = False)
Archetypes.append(S4H14SEAPG2)
Name = 'S4H14SEAWBPG2'
NoOfStories = 6
Sections = [
            Section1, Section1,
            Section1, Section1, Section2, Section2
            ]
S4H14SEAWBPG2 = CreateArchetype(Basements2Levels, Use2008Maps = False)
Archetypes.append(S4H14SEAWBPG2)
#endregion

#region Archetype S8H14SEAPG2 and S8H14SEAWBPG2
Name = 'S8H14SEAPG2'
# print 'Importing Archetype: ' + Name
#### Input Variables
NoOfStories = 8
Thickness = 24.
Length = 12. * 12.
Flange_Thickness = 6*12. # Assume 6' Long Core
Long_Spacing = 4
BarSize = 8.0
Rho = 2.0 #In Fraction
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section1 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, 3., 4., Spacing, Thickness - 3.5)
BarSize = 8.0
Rho = 1.0 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section2 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, 3., 4., Spacing*2, Thickness - 3.5)
BarSize = 4.0
Rho = 0.25 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section3 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, None, None, None, None)
Sections = [Section1, Section1, Section1,
            Section2, Section2, Section2,
            Section3, Section3,
            ]
S8H14SEAPG2 = CreateArchetype(Use2008Maps = False)
Archetypes.append(S8H14SEAPG2)
Name = 'S8H14SEAWBPG2'
NoOfStories = 11
Sections = [
            Section1, Section1, Section1,
            Section1, Section1, Section1,
            Section2, Section2, Section2,
            Section3, Section3,
            ]
S8H14SEAWBPG2 = CreateArchetype(Basements3Levels, Use2008Maps = False)
Archetypes.append(S8H14SEAWBPG2)
#endregion

#region Archetype S12H14SEAPG2 and S12H14SEAWBPG2
Name = 'S12H14SEAPG2'
# print 'Importing Archetype: ' + Name
#### Input Variables
NoOfStories = 12
Thickness = 24.
Length = 15. * 12.
Flange_Thickness = 7.5*12. # Assume 6' Long Core
Long_Spacing = 4
BarSize = 8.0
Rho = 1.6 #In Fraction
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section1 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, 3., 4., Spacing, Thickness - 3.5)
BarSize = 8.0
Rho = 1.2 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section2 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, 3., 4., Spacing*2, Thickness - 3.5)
ThicknessBelow = float(Thickness)
Thickness = 18.
Length = Length - (ThicknessBelow - Thickness) * 2.
BarSize = 4.0
Rho = 0.7 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section3 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu,  3., 4., Spacing*2, Thickness - 3.5)
BarSize = 4.0
Rho = 0.25 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section4 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, None, None, None, None)
Sections = [Section1, Section1, Section1,
            Section2, Section2, Section2,
            Section3, Section3, Section3,
            Section4, Section4, Section4,
            ]
S12H14SEAPG2 = CreateArchetype(Use2008Maps = False)
Archetypes.append(S12H14SEAPG2)
Name = 'S12H14SEAWBPG2'
NoOfStories = 16
Sections = [
            Section1, Section1, Section1, Section1,
            Section1, Section1, Section1,
            Section2, Section2, Section2,
            Section3, Section3, Section3,
            Section4, Section4, Section4,
            ]
S12H14SEAWBPG2 = CreateArchetype(Basements, Use2008Maps = False)
Archetypes.append(S12H14SEAWBPG2)
#endregion

#region Archetype S16H14SEAPG2 and S16H14SEAWBPG2
Name = 'S16H14SEAPG2'
#### Input Variables
NoOfStories = 16
Thickness = 28.
Length = 17. * 12.
Flange_Thickness = 8.5*12. # Assume 6' Long Core
Long_Spacing = 4
BarSize = 8.0
Rho = 1.5 #In Fraction
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section1 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, 3., 4., Spacing, Thickness - 3.5)
BarSize = 8.0
Rho = 1.0 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section2 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, 3., 4., Spacing*2, Thickness - 3.5)
ThicknessBelow = float(Thickness)
Thickness = 20.
Length = Length - (ThicknessBelow - Thickness) * 2.
BarSize = 5.0
Rho = 0.60 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section3 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, 3., 4., Spacing*2, Thickness - 3.5)
BarSize = 4.0
Rho = 0.25 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section4 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, None, None, None, None)
Sections = [Section1, Section1, Section1, Section1,
            Section2, Section2, Section2, Section2,
            Section3, Section3, Section3, Section3,
            Section4, Section4, Section4, Section4,
            ]
S16H14SEAPG2 = CreateArchetype(Use2008Maps = False)
Archetypes.append(S16H14SEAPG2)
Name = 'S16H14SEAWBPG2'
NoOfStories = 20 # Include Basement Floors Here
Sections = [Section1, Section1, Section1, Section1,
            Section1, Section1, Section1, Section1,
            Section2, Section2, Section2, Section2,
            Section3, Section3, Section3, Section3,
            Section4, Section4, Section4, Section4,
            ]
S16H14SEAWBPG2 = CreateArchetype(Basements, False)
Archetypes.append(S16H14SEAWBPG2)
#endregion

#region Archetype S20H14SEAPG2 and S20H14SEAWBPG2
Name = 'S20H14SEAPG2'
# print 'Importing Archetype: ' + Name
#### Input Variables
NoOfStories = 20
Thickness = 30.
Length = 19. * 12.
Flange_Thickness = 9.5*12. # Assume 6' Long Core
Long_Spacing = 4
BarSize = 9.0
Rho = 1.4 #In Fraction
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section1 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, 3., 4., Spacing, Thickness - 3.5)
BarSize = 8.0
Rho = 0.95 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section2 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, 3., 4., Spacing*2, Thickness - 3.5)
ThicknessBelow = float(Thickness)
Thickness = 22.
Length = Length - (ThicknessBelow - Thickness) * 2.
BarSize = 6.0
Rho =0.7 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section3 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, 3., 4., Spacing*2, Thickness - 3.5)

BarSize = 4.0
Rho = 0.25 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section4 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, None, None, None, None)
ThicknessBelow = float(Thickness)
Thickness = 22.
Length = Length - (ThicknessBelow - Thickness) * 2.
BarSize = 4.0
Rho = 0.25 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section5 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, None, None, None, None)
Sections = [Section1, Section1, Section1, Section1,
            Section2, Section2, Section2, Section2,
            Section3, Section3, Section3, Section3,
            Section4, Section4, Section4, Section4,
            Section5, Section5, Section5, Section5,
            ]
S20H14SEAPG2 = CreateArchetype(Use2008Maps = False)
Archetypes.append(S20H14SEAPG2)
Name = 'S20H14SEAWBPG2'
NoOfStories = 24 # Include Basement Floors Here
Sections = [Section1, Section1, Section1, Section1,
            Section1, Section1, Section1, Section1,
            Section2, Section2, Section2, Section2,
            Section3, Section3, Section3, Section3,
            Section4, Section4, Section4, Section4,
            Section5, Section5, Section5, Section5,
            ]
S20H14SEAWBPG2 = CreateArchetype(Basements, False)
Archetypes.append(S20H14SEAWBPG2)
#endregion

#region Archetype S24H14SEAPG2 and S24H14SEAWBPG2
Name = 'S24H14SEAPG2'
# print 'Importing Archetype: ' + Name
#### Input Variables
NoOfStories = 24
Thickness = 32.
Length = 21. * 12.
Flange_Thickness = 10.5*12. # Assume 6' Long Core
Long_Spacing = 4
BarSize = 9.0
Rho = 1.3 #In Fraction
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section1 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, 3., 4., Spacing, Thickness - 3.5)
BarSize = 8.0
Rho = 1.1 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section2 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, 3., 4., Spacing*2., Thickness - 3.5)
ThicknessBelow = float(Thickness)
Thickness = 26.
Length = Length - (ThicknessBelow - Thickness) * 2.
BarSize = 7.0
Rho = 0.8 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section3 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, 3., 4., Spacing*2., Thickness - 3.5)
BarSize = 4.0
Rho = 0.35 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section4 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu,   None, None, None, None)
ThicknessBelow = float(Thickness)
Thickness = 26.
Length = Length - (ThicknessBelow - Thickness) * 2.
BarSize = 4.0
Rho = 0.25 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section5 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, None, None, None, None)
BarSize = 4.0
Rho = 0.25 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section6 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, None, None, None, None)
Sections = [Section1, Section1, Section1, Section1,
            Section2, Section2, Section2, Section2,
            Section3, Section3, Section3, Section3,
            Section4, Section4, Section4, Section4,
            Section5, Section5, Section5, Section5,
            Section6, Section6, Section6, Section6,
            ]
S24H14SEAPG2 = CreateArchetype(Use2008Maps = False)
Archetypes.append(S24H14SEAPG2)
Name = 'S24H14SEAWBPG2'
NoOfStories = 28
Sections = [Section1, Section1, Section1, Section1,
            Section1, Section1, Section1, Section1,
            Section2, Section2, Section2, Section2,
            Section3, Section3, Section3, Section3,
            Section4, Section4, Section4, Section4,
            Section5, Section5, Section5, Section5,
            Section6, Section6, Section6, Section6,
            ]
S24H14SEAWBPG2 = CreateArchetype(Basements, False)
Archetypes.append(S24H14SEAWBPG2)
#endregion

#region Archetype S24H14SEAPG2 and S24H14SEAWBPG2
Name = 'S24H14SEAPG2TEST'
# print 'Importing Archetype: ' + Name
#### Input Variables
NoOfStories = 24
Thickness = 32.
Length = 21. * 12.
Flange_Thickness = 10.5*12. # Assume 6' Long Core
Long_Spacing = 4
BarSize = 9.0
Rho = 1.3 #In Fraction
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section1 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, 3., 4., Spacing, Thickness - 3.5)
BarSize = 8.0
Rho = 1.1 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section2 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, 3., 4., Spacing*2., Thickness - 3.5)
ThicknessBelow = float(Thickness)
Thickness = 26.
Length = Length - (ThicknessBelow - Thickness) * 2.
BarSize = 7.0
Rho = 0.8 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section3 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, 3., 4., Spacing*2., Thickness - 3.5)
BarSize = 4.0
Rho = 0.35 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section4 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu,   None, None, None, None)
ThicknessBelow = float(Thickness)
Thickness = 26.
Length = Length - (ThicknessBelow - Thickness) * 2.
BarSize = 4.0
Rho = 0.25 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section5 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, None, None, None, None)
BarSize = 4.0
Rho = 0.25 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section6 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, None, None, None, None)
Sections = [Section1, Section1, Section1, Section1,
            Section1, Section1, Section1, Section1,
            Section1, Section1, Section1, Section1,
            Section1, Section1, Section1, Section1,
            Section1, Section1, Section1, Section1,
            Section1, Section1, Section1, Section1,
            ]
S24H14SEAPG2TEST = CreateArchetype(Use2008Maps = False)
Archetypes.append(S24H14SEAPG2TEST)
Name = 'S24H14SEAWBPG2TEST'
NoOfStories = 28
Sections = [Section1, Section1, Section1, Section1,
            Section1, Section1, Section1, Section1,
            Section1, Section1, Section1, Section1,
            Section1, Section1, Section1, Section1,
            Section1, Section1, Section1, Section1,
            Section1, Section1, Section1, Section1,
            Section1, Section1, Section1, Section1,
            ]
S24H14SEAWBPG2TEST = CreateArchetype(Basements, False)
Archetypes.append(S24H14SEAWBPG2TEST)
#endregion

##### 2014 Maps ######
# PG3 : 25% Over-strength on ASCE 7 loads

#region Archetype S4H14SEAPG3 and S4H14SEAWBPG3
Name = 'S4H14SEAPG3'
# print 'Importing Archetype: ' + Name
# Compute Seismic Weight
NoOfStories = 4
YGrids = [0] + np.array(np.arange(0,(NoOfStories)*13*12, 13*12)+15*12).tolist()
DeadLoads = np.ones(NoOfStories) * DL / 1000.
DeadLoads[-1] = DeadLoads[-1] * DL_Roof / DL
LiveLoads = np.ones(NoOfStories) * LL / 1000.
LiveLoads[-1] = LiveLoads[-1] * LL_Roof / LL
MassPerSqFt = DL / 1000.
Mass = np.ones(NoOfStories) * MassPerSqFt * FloorArea
Mass[-1] = FloorArea * DL_Roof / 1000. # Adjust for Roof Weight
WallTribArea = FloorArea * 0.5
WeightPerSqFt = DL
BuildingWeight = np.ones(NoOfStories) * WeightPerSqFt * FloorArea
BuildingWeight[-1] = 152. / 1000. * FloorArea # Adjust for Roof Weight

# Seismic Hazard
R = 6; Cd = 5
SaDesign, Sds, CuTa = GetSeattle2014Hazard(YGrids[-1], R=R, Overstrength = 1.25)


Thickness = 20.
Length = 13. * 12.
Long_Spacing = 4
NoOfCols = 16
BarSize = 10.
Ag = ( (NoOfCols - 1) * Long_Spacing + 6 ) * Thickness
Rho = ( NoOfCols  * 2 + 2 ) * np.pi * ( BarSize / 2. / 8.) ** 2. / Ag
# print Rho
Section1 = PlanarWallSection(Length, Thickness,
                             (NoOfCols - 1) * Long_Spacing + 6,
                             (NoOfCols - 1) * Long_Spacing + 6, BarSize,
                             [3] + (np.ones(NoOfCols - 2) * 2.).tolist() + [3],
                             [3] + (np.ones(NoOfCols - 2) * 2.).tolist() + [3],
                             0.255, 4.037, fpc_core, fy, fu, 3, 4., NoOfCols, 3)
NoOfCols = 16
BarSize = 8.0
Ag = ( (NoOfCols - 1) * Long_Spacing + 6 ) * Thickness
Rho = ( NoOfCols * 2 + 2 ) * np.pi * ( BarSize / 2. / 8.) ** 2. / Ag
# print Rho
Section2 = PlanarWallSection(Length, Thickness,
                             (NoOfCols - 1) * Long_Spacing + 6,
                             (NoOfCols - 1) * Long_Spacing + 6, BarSize,
                             [3] + (np.ones(NoOfCols - 2) * 2.).tolist() + [3],
                             [3] + (np.ones(NoOfCols - 2) * 2.).tolist() + [3],
                             0.255, 4.037, fpc_core, fy, fu, 3, 4., 8, 3)
Section3 = PlanarWallSection(Length, Thickness, 0, 0, 10.173,
                             [],
                             [],
                             0.255, 4.037, fpc_core, fy, fu, None, None, None, None)
Sections = [Section1, Section1, Section2, Section2]
S4H14SEAPG3 = CreateArchetype(Use2008Maps = False, Overstrength = 1.25)
Archetypes.append(S4H14SEAPG3)
Name = 'S4H14SEAWBPG3'
NoOfStories = 6
Sections = [
            Section1, Section1,
            Section1, Section1, Section2, Section2
            ]
S4H14SEAWBPG3 = CreateArchetype(Basements2Levels, Use2008Maps = False, Overstrength = 1.25)
Archetypes.append(S4H14SEAWBPG3)
#endregion

#region Archetype S8H14SEAPG3 and S8H14SEAWBPG3
Name = 'S8H14SEAPG3'
# print 'Importing Archetype: ' + Name
#### Input Variables
NoOfStories = 8
Thickness = 24.
Length = 14. * 12.
Flange_Thickness = 7*12. # Assume 6' Long Core
Long_Spacing = 4
BarSize = 8.0
Rho = 2.0 #In Fraction
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section1 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, 3., 4., Spacing, Thickness - 3.5)
BarSize = 8.0
Rho = 1.1 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section2 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, 3., 4., Spacing*2, Thickness - 3.5)
BarSize = 4.0
Rho = 0.30 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section3 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, None, None, None, None)
Sections = [Section1, Section1, Section1,
            Section2, Section2, Section2,
            Section3, Section3,
            ]
S8H14SEAPG3 = CreateArchetype(Use2008Maps = False, Overstrength = 1.25)
Archetypes.append(S8H14SEAPG3)
Name = 'S8H14SEAWBPG3'
NoOfStories = 11
Sections = [
            Section1, Section1, Section1,
            Section1, Section1, Section1,
            Section2, Section2, Section2,
            Section3, Section3,
            ]
S8H14SEAWBPG3 = CreateArchetype(Basements3Levels, Use2008Maps = False, Overstrength = 1.25)
Archetypes.append(S8H14SEAWBPG3)
#endregion

#region Archetype S12H14SEAPG3 and S12H14SEAWBPG3
Name = 'S12H14SEAPG3'
# print 'Importing Archetype: ' + Name
#### Input Variables
NoOfStories = 12
Thickness = 24.
Length = 18. * 12.
Flange_Thickness = 9*12. # Assume 6' Long Core
Long_Spacing = 4
BarSize = 8.0
Rho = 1.55 #In Fraction
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section1 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, 3., 4., Spacing, Thickness - 3.5)
BarSize = 8.0
Rho = 1.25 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section2 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, 3., 4., Spacing*2, Thickness - 3.5)
ThicknessBelow = float(Thickness)
Thickness = 18.
Length = Length - (ThicknessBelow - Thickness) * 2.
BarSize = 4.0
Rho = 0.75 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section3 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu,  3., 4., Spacing*2, Thickness - 3.5)
BarSize = 4.0
Rho = 0.25 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section4 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, None, None, None, None)
Sections = [Section1, Section1, Section1,
            Section2, Section2, Section2,
            Section3, Section3, Section3,
            Section4, Section4, Section4,
            ]
S12H14SEAPG3 = CreateArchetype(Use2008Maps = False, Overstrength = 1.25)
Archetypes.append(S12H14SEAPG3)
Name = 'S12H14SEAWBPG3'
NoOfStories = 16
Sections = [
            Section1, Section1, Section1, Section1,
            Section1, Section1, Section1,
            Section2, Section2, Section2,
            Section3, Section3, Section3,
            Section4, Section4, Section4,
            ]
S12H14SEAWBPG3 = CreateArchetype(Basements, Use2008Maps = False, Overstrength = 1.25)
Archetypes.append(S12H14SEAWBPG3)
#endregion

#region Archetype S16H14SEAPG3 and S16H14SEAWBPG3
Name = 'S16H14SEAPG3'
#### Input Variables
NoOfStories = 16
Thickness = 34.
Length = 22. * 12.
Flange_Thickness = 11*12. # Assume 6' Long Core
Long_Spacing = 4
BarSize = 8.0
Rho = 0.9 #In Fraction
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section1 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, 3., 4., Spacing, Thickness - 3.5)
BarSize = 7.0
Rho = 0.8 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section2 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, 3., 4., Spacing*2, Thickness - 3.5)
ThicknessBelow = float(Thickness)
Thickness = 24.
Length = Length - (ThicknessBelow - Thickness) * 2.
BarSize = 5.0
Rho = 0.60 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section3 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, 3., 4., Spacing*2, Thickness - 3.5)
BarSize = 4.0
Rho = 0.25 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section4 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, None, None, None, None)
Sections = [Section1, Section1, Section1, Section1,
            Section2, Section2, Section2, Section2,
            Section3, Section3, Section3, Section3,
            Section4, Section4, Section4, Section4,
            ]
S16H14SEAPG3 = CreateArchetype(Use2008Maps = False, Overstrength = 1.25)
Archetypes.append(S16H14SEAPG3)
Name = 'S16H14SEAWBPG3'
NoOfStories = 20 # Include Basement Floors Here
Sections = [Section1, Section1, Section1, Section1,
            Section1, Section1, Section1, Section1,
            Section2, Section2, Section2, Section2,
            Section3, Section3, Section3, Section3,
            Section4, Section4, Section4, Section4,
            ]
S16H14SEAWBPG3 = CreateArchetype(Basements, False, Overstrength = 1.25)
Archetypes.append(S16H14SEAWBPG3)
#endregion

#region Archetype S20H14SEAPG3 and S20H14SEAWBPG3
Name = 'S20H14SEAPG3'
# print 'Importing Archetype: ' + Name
#### Input Variables
NoOfStories = 20
Thickness = 40.
Length = 26. * 12.
Flange_Thickness = 13.*12. # Assume 6' Long Core

Long_Spacing = 4
BarSize = 8.0
Rho = 0.675 #In Fraction
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section1 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, 3., 4., Spacing, Thickness - 3.5)
BarSize = 7.0
Rho = 0.6 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section2 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, 3., 4., Spacing*2, Thickness - 3.5)
ThicknessBelow = float(Thickness)
Thickness = 26.
Length = Length - (ThicknessBelow - Thickness) * 2.
BarSize = 6.0
Rho =0.6 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section3 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, 3., 4., Spacing*2, Thickness - 3.5)

BarSize = 4.0
Rho = 0.25 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section4 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, None, None, None, None)
ThicknessBelow = float(Thickness)
Thickness = 22.
Length = Length - (ThicknessBelow - Thickness) * 2.
BarSize = 4.0
Rho = 0.25 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section5 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, None, None, None, None)
Sections = [Section1, Section1, Section1, Section1,
            Section2, Section2, Section2, Section2,
            Section3, Section3, Section3, Section3,
            Section4, Section4, Section4, Section4,
            Section5, Section5, Section5, Section5,
            ]
S20H14SEAPG3 = CreateArchetype(Use2008Maps = False, Overstrength = 1.25)
Archetypes.append(S20H14SEAPG3)
Name = 'S20H14SEAWBPG3'
NoOfStories = 24 # Include Basement Floors Here
Sections = [Section1, Section1, Section1, Section1,
            Section1, Section1, Section1, Section1,
            Section2, Section2, Section2, Section2,
            Section3, Section3, Section3, Section3,
            Section4, Section4, Section4, Section4,
            Section5, Section5, Section5, Section5,
            ]
S20H14SEAWBPG3 = CreateArchetype(Basements, False, Overstrength = 1.25)
Archetypes.append(S20H14SEAWBPG3)
#endregion

#region Archetype S24H14SEAPG3 and S24H14SEAWBPG3
Name = 'S24H14SEAPG3'
# print 'Importing Archetype: ' + Name
#### Input Variables
NoOfStories = 24
Thickness = 44.
Length = 30. * 12.
Flange_Thickness = 15*12. # Assume 6' Long Core
Long_Spacing = 4
BarSize = 7.0
Rho = 0.525
#In Fraction
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section1 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, 3., 4., Spacing, Thickness - 3.5)
BarSize = 7.0
Rho = 0.525 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section2 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, 3., 4., Spacing*2., Thickness - 3.5)
ThicknessBelow = float(Thickness)
Thickness = 30.
Length = Length - (ThicknessBelow - Thickness) * 2.
BarSize = 6.0
Rho = 0.55 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section3 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, 3., 4., Spacing*2., Thickness - 3.5)
BarSize = 4.0
Rho = 0.30 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section4 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu,   None, None, None, None)
ThicknessBelow = float(Thickness)
Thickness = 24.
Length = Length - (ThicknessBelow - Thickness) * 2.
BarSize = 4.0
Rho = 0.25 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section5 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, None, None, None, None)
BarSize = 4.0
Rho = 0.25 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section6 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, None, None, None, None)
Sections = [Section1, Section1, Section1, Section1,
            Section2, Section2, Section2, Section2,
            Section3, Section3, Section3, Section3,
            Section4, Section4, Section4, Section4,
            Section5, Section5, Section5, Section5,
            Section6, Section6, Section6, Section6,
            ]
S24H14SEAPG3 = CreateArchetype(Use2008Maps = False, Overstrength = 1.25)
Archetypes.append(S24H14SEAPG3)
Name = 'S24H14SEAWBPG3'
NoOfStories = 28
Sections = [Section1, Section1, Section1, Section1,
            Section1, Section1, Section1, Section1,
            Section2, Section2, Section2, Section2,
            Section3, Section3, Section3, Section3,
            Section4, Section4, Section4, Section4,
            Section5, Section5, Section5, Section5,
            Section6, Section6, Section6, Section6,
            ]
S24H14SEAWBPG3 = CreateArchetype(Basements, False, Overstrength = 1.25)
Archetypes.append(S24H14SEAWBPG3)
#endregion

# PG4 : 50% Over-strength on ASCE 7 loads

#region Archetype S4H14SEAPG4 and S4H14SEAWBPG4
Name = 'S4H14SEAPG4'
# print 'Importing Archetype: ' + Name
# Compute Seismic Weight
NoOfStories = 4
YGrids = [0] + np.array(np.arange(0,(NoOfStories)*13*12, 13*12)+15*12).tolist()
DeadLoads = np.ones(NoOfStories) * DL / 1000.
DeadLoads[-1] = DeadLoads[-1] * DL_Roof / DL
LiveLoads = np.ones(NoOfStories) * LL / 1000.
LiveLoads[-1] = LiveLoads[-1] * LL_Roof / LL
MassPerSqFt = DL / 1000.
Mass = np.ones(NoOfStories) * MassPerSqFt * FloorArea
Mass[-1] = FloorArea * DL_Roof / 1000. # Adjust for Roof Weight
WallTribArea = FloorArea * 0.5
WeightPerSqFt = DL
BuildingWeight = np.ones(NoOfStories) * WeightPerSqFt * FloorArea
BuildingWeight[-1] = 152. / 1000. * FloorArea # Adjust for Roof Weight

# Seismic Hazard
R = 6; Cd = 5
SaDesign, Sds, CuTa = GetSeattle2014Hazard(YGrids[-1], R=R, Overstrength = 1.50)


Thickness = 22.
Length = 15. * 12.
Long_Spacing = 4
NoOfCols = 18
BarSize = 10.
Ag = ( (NoOfCols - 1) * Long_Spacing + 6 ) * Thickness
Rho = ( NoOfCols  * 2 + 2 ) * np.pi * ( BarSize / 2. / 8.) ** 2. / Ag
# print Rho
Section1 = PlanarWallSection(Length, Thickness,
                             (NoOfCols - 1) * Long_Spacing + 6,
                             (NoOfCols - 1) * Long_Spacing + 6, BarSize,
                             [3] + (np.ones(NoOfCols - 2) * 2.).tolist() + [3],
                             [3] + (np.ones(NoOfCols - 2) * 2.).tolist() + [3],
                             0.255, 4.037, fpc_core, fy, fu, 3, 4., NoOfCols, 3)
NoOfCols = 18
BarSize = 8.0
Ag = ( (NoOfCols - 1) * Long_Spacing + 6 ) * Thickness
Rho = ( NoOfCols * 2 + 2 ) * np.pi * ( BarSize / 2. / 8.) ** 2. / Ag
# print Rho
Section2 = PlanarWallSection(Length, Thickness,
                             (NoOfCols - 1) * Long_Spacing + 6,
                             (NoOfCols - 1) * Long_Spacing + 6, BarSize,
                             [3] + (np.ones(NoOfCols - 2) * 2.).tolist() + [3],
                             [3] + (np.ones(NoOfCols - 2) * 2.).tolist() + [3],
                             0.255, 4.037, fpc_core, fy, fu, 3, 4., 8, 3)
Section3 = PlanarWallSection(Length, Thickness, 0, 0, 10.173,
                             [],
                             [],
                             0.255, 4.037, fpc_core, fy, fu, None, None, None, None)
Sections = [Section1, Section1, Section2, Section2]
S4H14SEAPG4 = CreateArchetype(Use2008Maps = False, Overstrength = 1.50)
Archetypes.append(S4H14SEAPG4)
Name = 'S4H14SEAWBPG4'
NoOfStories = 6
Sections = [
            Section1, Section1,
            Section1, Section1, Section2, Section2
            ]
S4H14SEAWBPG4 = CreateArchetype(Basements2Levels, Use2008Maps = False, Overstrength = 1.50)
Archetypes.append(S4H14SEAWBPG4)
#endregion

#region Archetype S8H14SEAPG4 and S8H14SEAWBPG4
Name = 'S8H14SEAPG4'
# print 'Importing Archetype: ' + Name
#### Input Variables
NoOfStories = 8
Thickness = 26.
Length = 15. * 12.
Flange_Thickness = 7.5*12. # Assume 6' Long Core
Long_Spacing = 4
BarSize = 10.0
Rho = 2.0 #In Fraction
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section1 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, 3., 4., Spacing, Thickness - 3.5)
BarSize = 9.0
Rho = 1.3 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section2 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, 3., 4., Spacing*2, Thickness - 3.5)
BarSize = 5.0
Rho = 0.40 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section3 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, None, None, None, None)
Sections = [Section1, Section1, Section1,
            Section2, Section2, Section2,
            Section3, Section3,
            ]
S8H14SEAPG4 = CreateArchetype(Use2008Maps = False, Overstrength = 1.50)
Archetypes.append(S8H14SEAPG4)
Name = 'S8H14SEAWBPG4'
NoOfStories = 11
Sections = [
            Section1, Section1, Section1,
            Section1, Section1, Section1,
            Section2, Section2, Section2,
            Section3, Section3,
            ]
S8H14SEAWBPG4 = CreateArchetype(Basements3Levels, Use2008Maps = False, Overstrength = 1.50)
Archetypes.append(S8H14SEAWBPG4)
#endregion

#region Archetype S12H14SEAPG4 and S12H14SEAWBPG4
Name = 'S12H14SEAPG4'
# print 'Importing Archetype: ' + Name
#### Input Variables
NoOfStories = 12
Thickness = 30.
Length = 18. * 12.
Flange_Thickness = 9.0*12. # Assume 6' Long Core
Long_Spacing = 4
BarSize = 9.0
Rho = 1.70 #In Fraction
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section1 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, 3., 4., Spacing, Thickness - 3.5)
BarSize = 8.0
Rho = 1.35 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section2 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, 3., 4., Spacing*2, Thickness - 3.5)
ThicknessBelow = float(Thickness)
Thickness = 22.
Length = Length - (ThicknessBelow - Thickness) * 2.
BarSize = 6.0
Rho = 0.9 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section3 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu,  3., 4., Spacing*2, Thickness - 3.5)
BarSize = 4.0
Rho = 0.35 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section4 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, None, None, None, None)
Sections = [Section1, Section1, Section1,
            Section2, Section2, Section2,
            Section3, Section3, Section3,
            Section4, Section4, Section4,
            ]
S12H14SEAPG4 = CreateArchetype(Use2008Maps = False, Overstrength = 1.50)
Archetypes.append(S12H14SEAPG4)
Name = 'S12H14SEAWBPG4'
NoOfStories = 16
Sections = [
            Section1, Section1, Section1, Section1,
            Section1, Section1, Section1,
            Section2, Section2, Section2,
            Section3, Section3, Section3,
            Section4, Section4, Section4,
            ]
S12H14SEAWBPG4 = CreateArchetype(Basements, Use2008Maps = False, Overstrength = 1.50)
Archetypes.append(S12H14SEAWBPG4)
#endregion

#region Archetype S16H14SEAPG4 and S16H14SEAWBPG4
Name = 'S16H14SEAPG4'
#### Input Variables
NoOfStories = 16
Thickness = 34.
Length = 23. * 12.
Flange_Thickness = 11.5*12. # Assume 6' Long Core
Long_Spacing = 4
BarSize = 9.0
Rho = 1.2 #In Fraction
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section1 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, 3., 4., Spacing, Thickness - 3.5)
BarSize = 8.0
Rho = 0.9 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section2 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, 3., 4., Spacing*2, Thickness - 3.5)
ThicknessBelow = float(Thickness)
Thickness = 26.
Length = Length - (ThicknessBelow - Thickness) * 2.
BarSize = 5.0
Rho = 0.65 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section3 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, 3., 4., Spacing*2, Thickness - 3.5)
BarSize = 4.0
Rho = 0.25 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section4 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, None, None, None, None)
Sections = [Section1, Section1, Section1, Section1,
            Section2, Section2, Section2, Section2,
            Section3, Section3, Section3, Section3,
            Section4, Section4, Section4, Section4,
            ]
S16H14SEAPG4 = CreateArchetype(Use2008Maps = False, Overstrength = 1.50)
Archetypes.append(S16H14SEAPG4)
Name = 'S16H14SEAWBPG4'
NoOfStories = 20 # Include Basement Floors Here
Sections = [Section1, Section1, Section1, Section1,
            Section1, Section1, Section1, Section1,
            Section2, Section2, Section2, Section2,
            Section3, Section3, Section3, Section3,
            Section4, Section4, Section4, Section4,
            ]
S16H14SEAWBPG4 = CreateArchetype(Basements, False, Overstrength = 1.50)
Archetypes.append(S16H14SEAWBPG4)
#endregion

#region Archetype S20H14SEAPG4 and S20H14SEAWBPG4
Name = 'S20H14SEAPG4'
# print 'Importing Archetype: ' + Name
#### Input Variables
NoOfStories = 20
Thickness = 44.
Length = 27. * 12.
Flange_Thickness = 13.5*12. # Assume 6' Long Core

Long_Spacing = 4
BarSize = 8.0
Rho = 0.825 #In Fraction
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section1 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, 3., 4., Spacing, Thickness - 3.5)
BarSize = 8.0
Rho = 0.7 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section2 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, 3., 4., Spacing*2, Thickness - 3.5)
ThicknessBelow = float(Thickness)
Thickness = 30.
Length = Length - (ThicknessBelow - Thickness) * 2.
BarSize = 8.0
Rho =0.70 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section3 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, 3., 4., Spacing*2, Thickness - 3.5)

BarSize = 5.0
Rho = 0.4 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section4 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, None, None, None, None)
ThicknessBelow = float(Thickness)
Thickness = 22.
Length = Length - (ThicknessBelow - Thickness) * 2.
BarSize = 4.0
Rho = 0.25 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section5 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, None, None, None, None)
Sections = [Section1, Section1, Section1, Section1,
            Section2, Section2, Section2, Section2,
            Section3, Section3, Section3, Section3,
            Section4, Section4, Section4, Section4,
            Section5, Section5, Section5, Section5,
            ]
S20H14SEAPG4 = CreateArchetype(Use2008Maps = False, Overstrength = 1.50)
Archetypes.append(S20H14SEAPG4)
Name = 'S20H14SEAWBPG4'
NoOfStories = 24 # Include Basement Floors Here
Sections = [Section1, Section1, Section1, Section1,
            Section1, Section1, Section1, Section1,
            Section2, Section2, Section2, Section2,
            Section3, Section3, Section3, Section3,
            Section4, Section4, Section4, Section4,
            Section5, Section5, Section5, Section5,
            ]
S20H14SEAWBPG4 = CreateArchetype(Basements, False, Overstrength = 1.50)
Archetypes.append(S20H14SEAWBPG4)
#endregion

#region Archetype S24H14SEAPG4 and S24H14SEAWBPG4
Name = 'S24H14SEAPG4'
# print 'Importing Archetype: ' + Name
#### Input Variables
NoOfStories = 24
Thickness = 50.
Length = 31. * 12.
Flange_Thickness = 15.5*12. # Assume 6' Long Core
Long_Spacing = 4
BarSize = 8.0
Rho = 0.7
#In Fraction
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section1 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, 3., 4., Spacing, Thickness - 3.5)
BarSize = 8.0
Rho = 0.6 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section2 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, 3., 4., Spacing*2., Thickness - 3.5)
ThicknessBelow = float(Thickness)
Thickness = 36.
Length = Length - (ThicknessBelow - Thickness) * 2.
BarSize = 7.0
Rho = 0.65 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section3 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, 3., 4., Spacing*2., Thickness - 3.5)
BarSize = 5.0
Rho = 0.40 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section4 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu,   None, None, None, None)
ThicknessBelow = float(Thickness)
Thickness = 26.
Length = Length - (ThicknessBelow - Thickness) * 2.
BarSize = 4.0
Rho = 0.25 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section5 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, None, None, None, None)
BarSize = 4.0
Rho = 0.25 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section6 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, None, None, None, None)
Sections = [Section1, Section1, Section1, Section1,
            Section2, Section2, Section2, Section2,
            Section3, Section3, Section3, Section3,
            Section4, Section4, Section4, Section4,
            Section5, Section5, Section5, Section5,
            Section6, Section6, Section6, Section6,
            ]
S24H14SEAPG4 = CreateArchetype(Use2008Maps = False, Overstrength = 1.50)
Archetypes.append(S24H14SEAPG4)
Name = 'S24H14SEAWBPG4'
NoOfStories = 28
Sections = [Section1, Section1, Section1, Section1,
            Section1, Section1, Section1, Section1,
            Section2, Section2, Section2, Section2,
            Section3, Section3, Section3, Section3,
            Section4, Section4, Section4, Section4,
            Section5, Section5, Section5, Section5,
            Section6, Section6, Section6, Section6,
            ]
S24H14SEAWBPG4 = CreateArchetype(Basements, False, Overstrength = 1.50)
Archetypes.append(S24H14SEAWBPG4)
#endregion

# PG5: 1.5% Drift Limit

#region Archetype S4H14SEAPG5 and S4H14SEAWBPG5
Name = 'S4H14SEAPG5'
# print 'Importing Archetype: ' + Name
# Compute Seismic Weight
NoOfStories = 4
YGrids = [0] + np.array(np.arange(0,(NoOfStories)*13*12, 13*12)+15*12).tolist()
DeadLoads = np.ones(NoOfStories) * DL / 1000.
DeadLoads[-1] = DeadLoads[-1] * DL_Roof / DL
LiveLoads = np.ones(NoOfStories) * LL / 1000.
LiveLoads[-1] = LiveLoads[-1] * LL_Roof / LL
MassPerSqFt = DL / 1000.
Mass = np.ones(NoOfStories) * MassPerSqFt * FloorArea
Mass[-1] = FloorArea * DL_Roof / 1000. # Adjust for Roof Weight
WallTribArea = FloorArea * 0.5
WeightPerSqFt = DL
BuildingWeight = np.ones(NoOfStories) * WeightPerSqFt * FloorArea
BuildingWeight[-1] = 152. / 1000. * FloorArea # Adjust for Roof Weight
# Seismic Hazard
R = 6; Cd = 5
SaDesign, Sds, CuTa = GetSeattle2008Hazard(YGrids[-1], R=R)

Thickness = 24.
Length = 13. * 12.
Long_Spacing = 4
NoOfCols = 9
BarSize = 10.
Ag = ( (NoOfCols - 1) * Long_Spacing + 6 ) * Thickness
Rho = ( NoOfCols  * 2 + 2 ) * np.pi * ( BarSize / 2. / 8.) ** 2. / Ag
# print Rho
Section1 = PlanarWallSection(Length, Thickness,
                             (NoOfCols - 1) * Long_Spacing + 6,
                             (NoOfCols - 1) * Long_Spacing + 6, BarSize,
                             [3] + (np.ones(NoOfCols - 2) * 2.).tolist() + [3],
                             [3] + (np.ones(NoOfCols - 2) * 2.).tolist() + [3],
                             0.255, 4.037, fpc_core, fy, fu, 3, 4., NoOfCols, 3)
NoOfCols = 9
BarSize = 8.0
Ag = ( (NoOfCols - 1) * Long_Spacing + 6 ) * Thickness
Rho = ( NoOfCols * 2 + 2 ) * np.pi * ( BarSize / 2. / 8.) ** 2. / Ag
# print Rho
Section2 = PlanarWallSection(Length, Thickness,
                             (NoOfCols - 1) * Long_Spacing + 6,
                             (NoOfCols - 1) * Long_Spacing + 6, BarSize,
                             [3] + (np.ones(NoOfCols - 2) * 2.).tolist() + [3],
                             [3] + (np.ones(NoOfCols - 2) * 2.).tolist() + [3],
                             0.255, 4.037, fpc_core, fy, fu, 3, 4., 8, 3)
Section3 = PlanarWallSection(Length, Thickness, 0, 0, 10.173,
                             [],
                             [],
                             0.255, 4.037, fpc_core, fy, fu, None, None, None, None)
Sections = [Section1, Section1, Section2, Section2]
S4H14SEAPG5 = CreateArchetype(Use2008Maps = False)
Archetypes.append(S4H14SEAPG5)
Name = 'S4H14SEAWBPG5'
NoOfStories = 6
Sections = [
            Section1, Section1,
            Section1, Section1, Section2, Section2
            ]
S4H14SEAWBPG5 = CreateArchetype(Basements2Levels, Use2008Maps = False)
Archetypes.append(S4H14SEAWBPG5)
#endregion

#region Archetype S8H14SEAPG5 and S8H14SEAWBPG5
Name = 'S8H14SEAPG5'
# print 'Importing Archetype: ' + Name
#### Input Variables
NoOfStories = 8
Thickness = 24.
Length = 14. * 12.
Flange_Thickness = 7*12. # Assume 6' Long Core
Long_Spacing = 4
BarSize = 8.0
Rho = 1.25 #In Fraction
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section1 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, 3., 4., Spacing, Thickness - 3.5)
BarSize = 8.0
Rho = 0.8 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section2 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, 3., 4., Spacing*2, Thickness - 3.5)
BarSize = 4.0
Rho = 0.25 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section3 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, None, None, None, None)
Sections = [Section1, Section1, Section1,
            Section2, Section2, Section2,
            Section3, Section3,
            ]
S8H14SEAPG5 = CreateArchetype(Use2008Maps = False)
Archetypes.append(S8H14SEAPG5)
Name = 'S8H14SEAWBPG5'
NoOfStories = 11
Sections = [
            Section1, Section1, Section1,
            Section1, Section1, Section1,
            Section2, Section2, Section2,
            Section3, Section3,
            ]
S8H14SEAWBPG5 = CreateArchetype(Basements3Levels, Use2008Maps = False)
Archetypes.append(S8H14SEAWBPG5)
#endregion

#region Archetype S12H14SEAPG5 and S12H14SEAWBPG5
Name = 'S12H14SEAPG5'
# print 'Importing Archetype: ' + Name
#### Input Variables
NoOfStories = 12
Thickness = 26.
Length = 17. * 12.
Flange_Thickness = 8.5*12. # Assume 6' Long Core
Long_Spacing = 4
BarSize = 8.0
Rho = 1.025 #In Fraction
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section1 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, 3., 4., Spacing, Thickness - 3.5)
BarSize = 7.0
Rho = 0.80 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section2 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, 3., 4., Spacing*2, Thickness - 3.5)
ThicknessBelow = float(Thickness)
Thickness = 18.
Length = Length - (ThicknessBelow - Thickness) * 2.
BarSize = 5.0
Rho = 0.60 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section3 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu,  3., 4., Spacing*2, Thickness - 3.5)
BarSize = 4.0
Rho = 0.25 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section4 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, None, None, None, None)
Sections = [Section1, Section1, Section1,
            Section2, Section2, Section2,
            Section3, Section3, Section3,
            Section4, Section4, Section4,
            ]
S12H14SEAPG5 = CreateArchetype(Use2008Maps = False)
Archetypes.append(S12H14SEAPG5)
Name = 'S12H14SEAWBPG5'
NoOfStories = 16
Sections = [
            Section1, Section1, Section1, Section1,
            Section1, Section1, Section1,
            Section2, Section2, Section2,
            Section3, Section3, Section3,
            Section4, Section4, Section4,
            ]
S12H14SEAWBPG5 = CreateArchetype(Basements, Use2008Maps = False)
Archetypes.append(S12H14SEAWBPG5)
#endregion

#region Archetype S16H14SEAPG5 and S16H14SEAWBPG5
Name = 'S16H14SEAPG5'
#### Input Variables
NoOfStories = 16
Thickness = 32.
Length = 20. * 12.
Flange_Thickness = 10*12. # Assume 6' Long Core
Long_Spacing = 4
BarSize = 8.0
Rho = 0.725 #In Fraction
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section1 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, 3., 4., Spacing, Thickness - 3.5)
BarSize = 8.0
Rho = 0.6 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section2 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, 3., 4., Spacing*2, Thickness - 3.5)
ThicknessBelow = float(Thickness)
Thickness = 18.
Length = Length - (ThicknessBelow - Thickness) * 2.
BarSize = 5.0
Rho = 0.50 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section3 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, 3., 4., Spacing*2, Thickness - 3.5)
BarSize = 4.0
Rho = 0.25 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section4 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, None, None, None, None)
Sections = [Section1, Section1, Section1, Section1,
            Section2, Section2, Section2, Section2,
            Section3, Section3, Section3, Section3,
            Section4, Section4, Section4, Section4,
            ]
S16H14SEAPG5 = CreateArchetype(Use2008Maps = False)
Archetypes.append(S16H14SEAPG5)
Name = 'S16H14SEAWBPG5'
NoOfStories = 20 # Include Basement Floors Here
Sections = [Section1, Section1, Section1, Section1,
            Section1, Section1, Section1, Section1,
            Section2, Section2, Section2, Section2,
            Section3, Section3, Section3, Section3,
            Section4, Section4, Section4, Section4,
            ]
S16H14SEAWBPG5 = CreateArchetype(Basements, False)
Archetypes.append(S16H14SEAWBPG5)
#endregion

#region Archetype S20H14SEAPG5 and S20H14SEAWBPG5
Name = 'S20H14SEAPG5'
# print 'Importing Archetype: ' + Name
#### Input Variables
NoOfStories = 20
Thickness = 36.
Length = 23. * 12.
Flange_Thickness = 11.5*12. # Assume 6' Long Core
Long_Spacing = 4
BarSize = 7.0
Rho = 0.525 #In Fraction
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section1 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, 3., 4., Spacing, Thickness - 3.5)
BarSize = 7.0
Rho = 0.525 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section2 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, 3., 4., Spacing*2, Thickness - 3.5)
ThicknessBelow = float(Thickness)
Thickness = 20.
Length = Length - (ThicknessBelow - Thickness) * 2.
BarSize = 6.0
Rho =0.5 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section3 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, 3., 4., Spacing*2, Thickness - 3.5)

BarSize = 4.0
Rho = 0.25 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section4 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, None, None, None, None)
ThicknessBelow = float(Thickness)
Thickness = 22.
Length = Length - (ThicknessBelow - Thickness) * 2.
BarSize = 4.0
Rho = 0.25 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section5 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, None, None, None, None)
Sections = [Section1, Section1, Section1, Section1,
            Section2, Section2, Section2, Section2,
            Section3, Section3, Section3, Section3,
            Section4, Section4, Section4, Section4,
            Section5, Section5, Section5, Section5,
            ]
S20H14SEAPG5 = CreateArchetype(Use2008Maps = False)
Archetypes.append(S20H14SEAPG5)
Name = 'S20H14SEAWBPG5'
NoOfStories = 24 # Include Basement Floors Here
Sections = [Section1, Section1, Section1, Section1,
            Section1, Section1, Section1, Section1,
            Section2, Section2, Section2, Section2,
            Section3, Section3, Section3, Section3,
            Section4, Section4, Section4, Section4,
            Section5, Section5, Section5, Section5,
            ]
S20H14SEAWBPG5 = CreateArchetype(Basements, False)
Archetypes.append(S20H14SEAWBPG5)
#endregion

#region Archetype S24H14SEAPG5 and S24H14SEAWBPG5
Name = 'S24H14SEAPG5'
# print 'Importing Archetype: ' + Name
#### Input Variables
NoOfStories = 24
Thickness = 40.
Length = 25. * 12.
Flange_Thickness = 12.5*12. # Assume 6' Long Core
Long_Spacing = 4
BarSize = 7.0
Rho = 0.50 #In Fraction
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section1 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, 3., 4., Spacing, Thickness - 3.5)
BarSize = 7.0
Rho = 0.501 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section2 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, 3., 4., Spacing*2., Thickness - 3.5)
ThicknessBelow = float(Thickness)
Thickness = 26.
Length = Length - (ThicknessBelow - Thickness) * 2.
BarSize = 7.0
Rho = 0.55 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section3 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, 3., 4., Spacing*2., Thickness - 3.5)
BarSize = 4.0
Rho = 0.35 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section4 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu,   None, None, None, None)
ThicknessBelow = float(Thickness)
Thickness = 20.
Length = Length - (ThicknessBelow - Thickness) * 2.
BarSize = 4.0
Rho = 0.25 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section5 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, None, None, None, None)
BarSize = 4.0
Rho = 0.25 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section6 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, None, None, None, None)
Sections = [Section1, Section1, Section1, Section1,
            Section2, Section2, Section2, Section2,
            Section3, Section3, Section3, Section3,
            Section4, Section4, Section4, Section4,
            Section5, Section5, Section5, Section5,
            Section6, Section6, Section6, Section6,
            ]
S24H14SEAPG5 = CreateArchetype(Use2008Maps = False)
Archetypes.append(S24H14SEAPG5)
Name = 'S24H14SEAWBPG5'
NoOfStories = 28
Sections = [Section1, Section1, Section1, Section1,
            Section1, Section1, Section1, Section1,
            Section2, Section2, Section2, Section2,
            Section3, Section3, Section3, Section3,
            Section4, Section4, Section4, Section4,
            Section5, Section5, Section5, Section5,
            Section6, Section6, Section6, Section6,
            ]
S24H14SEAWBPG5 = CreateArchetype(Basements, False)
Archetypes.append(S24H14SEAWBPG5)
#endregion

# PG6: 1.25% Drift Limit

#region Archetype S4H14SEAPG6 and S4H14SEAWBPG6
Name = 'S4H14SEAPG6'
# print 'Importing Archetype: ' + Name
# Compute Seismic Weight
NoOfStories = 4
YGrids = [0] + np.array(np.arange(0,(NoOfStories)*13*12, 13*12)+15*12).tolist()
DeadLoads = np.ones(NoOfStories) * DL / 1000.
DeadLoads[-1] = DeadLoads[-1] * DL_Roof / DL
LiveLoads = np.ones(NoOfStories) * LL / 1000.
LiveLoads[-1] = LiveLoads[-1] * LL_Roof / LL
MassPerSqFt = DL / 1000.
Mass = np.ones(NoOfStories) * MassPerSqFt * FloorArea
Mass[-1] = FloorArea * DL_Roof / 1000. # Adjust for Roof Weight
WallTribArea = FloorArea * 0.5
WeightPerSqFt = DL
BuildingWeight = np.ones(NoOfStories) * WeightPerSqFt * FloorArea
BuildingWeight[-1] = 152. / 1000. * FloorArea # Adjust for Roof Weight
# Seismic Hazard
R = 6; Cd = 5
SaDesign, Sds, CuTa = GetSeattle2008Hazard(YGrids[-1], R=R)

Thickness = 28.
Length = 14. * 12.
Long_Spacing = 4
NoOfCols = 10
BarSize = 9.
Ag = ( (NoOfCols - 1) * Long_Spacing + 6 ) * Thickness
Rho = ( NoOfCols  * 2 + 2 ) * np.pi * ( BarSize / 2. / 8.) ** 2. / Ag
# print Rho
Section1 = PlanarWallSection(Length, Thickness,
                             (NoOfCols - 1) * Long_Spacing + 6,
                             (NoOfCols - 1) * Long_Spacing + 6, BarSize,
                             [3] + (np.ones(NoOfCols - 2) * 2.).tolist() + [3],
                             [3] + (np.ones(NoOfCols - 2) * 2.).tolist() + [3],
                             0.255, 4.037, fpc_core, fy, fu, 3, 4., NoOfCols, 3)
NoOfCols = 8
BarSize = 8.0
Ag = ( (NoOfCols - 1) * Long_Spacing + 6 ) * Thickness
Rho = ( NoOfCols * 2 + 2 ) * np.pi * ( BarSize / 2. / 8.) ** 2. / Ag
# print Rho
Section2 = PlanarWallSection(Length, Thickness,
                             (NoOfCols - 1) * Long_Spacing + 6,
                             (NoOfCols - 1) * Long_Spacing + 6, BarSize,
                             [3] + (np.ones(NoOfCols - 2) * 2.).tolist() + [3],
                             [3] + (np.ones(NoOfCols - 2) * 2.).tolist() + [3],
                             0.255, 4.037, fpc_core, fy, fu, 3, 4., 8, 3)
Section3 = PlanarWallSection(Length, Thickness, 0, 0, 10.173,
                             [],
                             [],
                             0.255, 4.037, fpc_core, fy, fu, None, None, None, None)
Sections = [Section1, Section1, Section2, Section2]
S4H14SEAPG6 = CreateArchetype(Use2008Maps = False)
Archetypes.append(S4H14SEAPG6)
Name = 'S4H14SEAWBPG6'
NoOfStories = 6
Sections = [
            Section1, Section1,
            Section1, Section1, Section2, Section2
            ]
S4H14SEAWBPG6 = CreateArchetype(Basements2Levels, Use2008Maps = False)
Archetypes.append(S4H14SEAWBPG6)
#endregion

#region Archetype S8H14SEAPG6 and S8H14SEAWBPG6
Name = 'S8H14SEAPG6'
# print 'Importing Archetype: ' + Name
#### Input Variables
NoOfStories = 8
Thickness = 24.
Length = 15. * 12.
Flange_Thickness = 7.5*12. # Assume 6' Long Core
Long_Spacing = 4
BarSize = 7.0
Rho = 0.975 #In Fraction
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section1 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, 3., 4., Spacing, Thickness - 3.5)
BarSize = 7.0
Rho = 0.65 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section2 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, 3., 4., Spacing*2, Thickness - 3.5)
BarSize = 4.0
Rho = 0.25 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section3 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, None, None, None, None)
Sections = [Section1, Section1, Section1,
            Section2, Section2, Section2,
            Section3, Section3,
            ]
S8H14SEAPG6 = CreateArchetype(Use2008Maps = False)
Archetypes.append(S8H14SEAPG6)
Name = 'S8H14SEAWBPG6'
NoOfStories = 11
Sections = [
            Section1, Section1, Section1,
            Section1, Section1, Section1,
            Section2, Section2, Section2,
            Section3, Section3,
            ]
S8H14SEAWBPG6 = CreateArchetype(Basements3Levels, Use2008Maps = False)
Archetypes.append(S8H14SEAWBPG6)
#endregion

#region Archetype S12H14SEAPG6 and S12H14SEAWBPG6
Name = 'S12H14SEAPG6'
# print 'Importing Archetype: ' + Name
#### Input Variables
NoOfStories = 12
Thickness = 28.
Length = 19. * 12.
Flange_Thickness = 9.5*12. # Assume 6' Long Core
Long_Spacing = 4
BarSize = 7.0
Rho = 0.55 #In Fraction
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section1 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, 3., 4., Spacing, Thickness - 3.5)
BarSize = 6.0
Rho = 0.50 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section2 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, 3., 4., Spacing*2, Thickness - 3.5)
ThicknessBelow = float(Thickness)
Thickness = 16.
Length = Length - (ThicknessBelow - Thickness) * 2.
BarSize = 5.0
Rho = 0.50 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section3 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu,  3., 4., Spacing*2, Thickness - 3.5)
BarSize = 4.0
Rho = 0.25 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section4 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, None, None, None, None)
Sections = [Section1, Section1, Section1,
            Section2, Section2, Section2,
            Section3, Section3, Section3,
            Section4, Section4, Section4,
            ]
S12H14SEAPG6 = CreateArchetype(Use2008Maps = False)
Archetypes.append(S12H14SEAPG6)
Name = 'S12H14SEAWBPG6'
NoOfStories = 16
Sections = [
            Section1, Section1, Section1, Section1,
            Section1, Section1, Section1,
            Section2, Section2, Section2,
            Section3, Section3, Section3,
            Section4, Section4, Section4,
            ]
S12H14SEAWBPG6 = CreateArchetype(Basements, Use2008Maps = False)
Archetypes.append(S12H14SEAWBPG6)
#endregion

#region Archetype S16H14SEAPG6 and S16H14SEAWBPG6
Name = 'S16H14SEAPG6'
#### Input Variables
NoOfStories = 16
Thickness = 32.
Length = 22. * 12.
Flange_Thickness = 11.*12. # Assume 6' Long Core
Long_Spacing = 4
BarSize = 6.0
Rho = 0.50 #In Fraction
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section1 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, 3., 4., Spacing, Thickness - 3.5)
BarSize = 6.0
Rho = 0.5 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section2 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, 3., 4., Spacing*2, Thickness - 3.5)
ThicknessBelow = float(Thickness)
Thickness = 16.
Length = Length - (ThicknessBelow - Thickness) * 2.
BarSize = 5.0
Rho = 0.50 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section3 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, 3., 4., Spacing*2, Thickness - 3.5)
BarSize = 4.0
Rho = 0.25 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section4 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, None, None, None, None)
Sections = [Section1, Section1, Section1, Section1,
            Section2, Section2, Section2, Section2,
            Section3, Section3, Section3, Section3,
            Section4, Section4, Section4, Section4,
            ]
S16H14SEAPG6 = CreateArchetype(Use2008Maps = False)
Archetypes.append(S16H14SEAPG6)
Name = 'S16H14SEAWBPG6'
NoOfStories = 20 # Include Basement Floors Here
Sections = [Section1, Section1, Section1, Section1,
            Section1, Section1, Section1, Section1,
            Section2, Section2, Section2, Section2,
            Section3, Section3, Section3, Section3,
            Section4, Section4, Section4, Section4,
            ]
S16H14SEAWBPG6 = CreateArchetype(Basements, False)
Archetypes.append(S16H14SEAWBPG6)
#endregion

#region Archetype S20H14SEAPG6 and S20H14SEAWBPG6
Name = 'S20H14SEAPG6'
# print 'Importing Archetype: ' + Name
#### Input Variables
NoOfStories = 20
Thickness = 24.
Length = 25. * 12.
Flange_Thickness = 12.5*12. # Assume 6' Long Core
Long_Spacing = 4
BarSize = 7.0
Rho = 0.501 #In Fraction
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section1 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, 3., 4., Spacing, Thickness - 3.5)
BarSize = 7.0
Rho = 0.501 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section2 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, 3., 4., Spacing*2, Thickness - 3.5)
ThicknessBelow = float(Thickness)
Thickness = 20.
Length = Length - (ThicknessBelow - Thickness) * 2.
BarSize = 6.0
Rho =0.501 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section3 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, 3., 4., Spacing*2, Thickness - 3.5)

BarSize = 4.0
Rho = 0.30 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section4 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, None, None, None, None)
ThicknessBelow = float(Thickness)
Thickness = 16.
Length = Length - (ThicknessBelow - Thickness) * 2.
BarSize = 4.0
Rho = 0.25 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section5 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, None, None, None, None)
Sections = [Section1, Section1, Section1, Section1,
            Section2, Section2, Section2, Section2,
            Section3, Section3, Section3, Section3,
            Section4, Section4, Section4, Section4,
            Section5, Section5, Section5, Section5,
            ]
S20H14SEAPG6 = CreateArchetype(Use2008Maps = False)
Archetypes.append(S20H14SEAPG6)
Name = 'S20H14SEAWBPG6'
NoOfStories = 24 # Include Basement Floors Here
Sections = [Section1, Section1, Section1, Section1,
            Section1, Section1, Section1, Section1,
            Section2, Section2, Section2, Section2,
            Section3, Section3, Section3, Section3,
            Section4, Section4, Section4, Section4,
            Section5, Section5, Section5, Section5,
            ]
S20H14SEAWBPG6 = CreateArchetype(Basements, False)
Archetypes.append(S20H14SEAWBPG6)
#endregion

#region Archetype S24H14SEAPG6 and S24H14SEAWBPG6
Name = 'S24H14SEAPG6'
# print 'Importing Archetype: ' + Name
#### Input Variables
NoOfStories = 24
Thickness = 20.
Length = 28. * 12.
Flange_Thickness = 14.*12. # Assume 6' Long Core
Long_Spacing = 4
BarSize = 5.0
Rho = 0.50 #In Fraction
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section1 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, 3., 4., Spacing, Thickness - 3.5)
BarSize = 5.0
Rho = 0.501 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section2 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, 3., 4., Spacing*2., Thickness - 3.5)
ThicknessBelow = float(Thickness)
Thickness = 20.
Length = Length - (ThicknessBelow - Thickness) * 2.
BarSize = 5.0
Rho = 0.501 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section3 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, 3., 4., Spacing*2., Thickness - 3.5)
BarSize = 4.0
Rho = 0.25 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section4 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu,   None, None, None, None)
ThicknessBelow = float(Thickness)
Thickness = 16.
Length = Length - (ThicknessBelow - Thickness) * 2.
BarSize = 4.0
Rho = 0.25 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section5 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, None, None, None, None)
BarSize = 4.0
Rho = 0.25 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section6 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, None, None, None, None)
Sections = [Section1, Section1, Section1, Section1,
            Section2, Section2, Section2, Section2,
            Section3, Section3, Section3, Section3,
            Section4, Section4, Section4, Section4,
            Section5, Section5, Section5, Section5,
            Section6, Section6, Section6, Section6,
            ]
S24H14SEAPG6 = CreateArchetype(Use2008Maps = False)
Archetypes.append(S24H14SEAPG6)
Name = 'S24H14SEAWBPG6'
NoOfStories = 28
Sections = [Section1, Section1, Section1, Section1,
            Section1, Section1, Section1, Section1,
            Section2, Section2, Section2, Section2,
            Section3, Section3, Section3, Section3,
            Section4, Section4, Section4, Section4,
            Section5, Section5, Section5, Section5,
            Section6, Section6, Section6, Section6,
            ]
S24H14SEAWBPG6 = CreateArchetype(Basements, False)
Archetypes.append(S24H14SEAWBPG6)
#endregion


# PG7 : 25% Over-strength on ASCE 7 loads

#region Archetype S4H14SEAPG7 and S4H14SEAWBPG7
Name = 'S4H14SEAPG7'
# print 'Importing Archetype: ' + Name
# Compute Seismic Weight
NoOfStories = 4
YGrids = [0] + np.array(np.arange(0,(NoOfStories)*13*12, 13*12)+15*12).tolist()
DeadLoads = np.ones(NoOfStories) * DL / 1000.
DeadLoads[-1] = DeadLoads[-1] * DL_Roof / DL
LiveLoads = np.ones(NoOfStories) * LL / 1000.
LiveLoads[-1] = LiveLoads[-1] * LL_Roof / LL
MassPerSqFt = DL / 1000.
Mass = np.ones(NoOfStories) * MassPerSqFt * FloorArea
Mass[-1] = FloorArea * DL_Roof / 1000. # Adjust for Roof Weight
WallTribArea = FloorArea * 0.5
WeightPerSqFt = DL
BuildingWeight = np.ones(NoOfStories) * WeightPerSqFt * FloorArea
BuildingWeight[-1] = 152. / 1000. * FloorArea # Adjust for Roof Weight

# Seismic Hazard
R = 6; Cd = 5
SaDesign, Sds, CuTa = GetSeattle2014Hazard(YGrids[-1], R=R, Overstrength = 1.25)


Thickness = 24.
Length = 16. * 12.
Long_Spacing = 4
NoOfCols = 10
BarSize = 10.
Ag = ( (NoOfCols - 1) * Long_Spacing + 6 ) * Thickness
Rho = ( NoOfCols  * 2 + 2 ) * np.pi * ( BarSize / 2. / 8.) ** 2. / Ag

# print Rho
Section1 = PlanarWallSection(Length, Thickness,
                             (NoOfCols - 1) * Long_Spacing + 6,
                             (NoOfCols - 1) * Long_Spacing + 6, BarSize,
                             [3] + (np.ones(NoOfCols - 2) * 2.).tolist() + [3],
                             [3] + (np.ones(NoOfCols - 2) * 2.).tolist() + [3],
                             0.255, 4.037, fpc_core, fy, fu, 3, 4., NoOfCols, 3)
NoOfCols = 10
BarSize = 9.0
Ag = ( (NoOfCols - 1) * Long_Spacing + 6 ) * Thickness
Rho = ( NoOfCols * 2 + 2 ) * np.pi * ( BarSize / 2. / 8.) ** 2. / Ag

# print Rho
Section2 = PlanarWallSection(Length, Thickness,
                             (NoOfCols - 1) * Long_Spacing + 6,
                             (NoOfCols - 1) * Long_Spacing + 6, BarSize,
                             [3] + (np.ones(NoOfCols - 2) * 2.).tolist() + [3],
                             [3] + (np.ones(NoOfCols - 2) * 2.).tolist() + [3],
                             0.255, 4.037, fpc_core, fy, fu, 3, 4., 8, 3)
Section3 = PlanarWallSection(Length, Thickness, 0, 0, 10.173,
                             [],
                             [],
                             0.255, 4.037, fpc_core, fy, fu, None, None, None, None)
Sections = [Section1, Section1, Section2, Section2]
S4H14SEAPG7 = CreateArchetype(Use2008Maps = False, Overstrength = 1.25)
Archetypes.append(S4H14SEAPG7)
Name = 'S4H14SEAWBPG7'
NoOfStories = 6
Sections = [
            Section1, Section1,
            Section1, Section1, Section2, Section2
            ]
S4H14SEAWBPG7 = CreateArchetype(Basements2Levels, Use2008Maps = False, Overstrength = 1.25)
Archetypes.append(S4H14SEAWBPG7)
#endregion

#region Archetype S8H14SEAPG7 and S8H14SEAWBPG7
Name = 'S8H14SEAPG7'
# print 'Importing Archetype: ' + Name
#### Input Variables
NoOfStories = 8
Thickness = 24.
Length = 16. * 12.
Flange_Thickness = 8*12. # Assume 6' Long Core
Long_Spacing = 4
BarSize = 7.0
Rho = 1.25 #In Fraction
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section1 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, 3., 4., Spacing, Thickness - 3.5)
BarSize = 7.0
Rho = 0.85 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section2 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, 3., 4., Spacing*2, Thickness - 3.5)
BarSize = 4.0
Rho = 0.30 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section3 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, None, None, None, None)
Sections = [Section1, Section1, Section1,
            Section2, Section2, Section2,
            Section3, Section3,
            ]
S8H14SEAPG7 = CreateArchetype(Use2008Maps = False, Overstrength = 1.25)
Archetypes.append(S8H14SEAPG7)
Name = 'S8H14SEAWBPG7'
NoOfStories = 11
Sections = [
            Section1, Section1, Section1,
            Section1, Section1, Section1,
            Section2, Section2, Section2,
            Section3, Section3,
            ]
S8H14SEAWBPG7 = CreateArchetype(Basements3Levels, Use2008Maps = False, Overstrength = 1.25)
Archetypes.append(S8H14SEAWBPG7)
#endregion

#region Archetype S12H14SEAPG7 and S12H14SEAWBPG7
Name = 'S12H14SEAPG7'
# print 'Importing Archetype: ' + Name
#### Input Variables
NoOfStories = 12
Thickness = 26.
Length = 20. * 12.
Flange_Thickness = 10*12. # Assume 6' Long Core
Long_Spacing = 4
BarSize = 8.0
Rho = 1.0 #In Fraction
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section1 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, 3., 4., Spacing, Thickness - 3.5)
BarSize = 8.0
Rho = 0.8 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section2 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, 3., 4., Spacing*2, Thickness - 3.5)
ThicknessBelow = float(Thickness)
Thickness = 18.
Length = Length - (ThicknessBelow - Thickness) * 2.
BarSize = 4.0
Rho = 0.75 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section3 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu,  3., 4., Spacing*2, Thickness - 3.5)
BarSize = 4.0
Rho = 0.25 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section4 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, None, None, None, None)
Sections = [Section1, Section1, Section1,
            Section2, Section2, Section2,
            Section3, Section3, Section3,
            Section4, Section4, Section4,
            ]
S12H14SEAPG7 = CreateArchetype(Use2008Maps = False, Overstrength = 1.25)
Archetypes.append(S12H14SEAPG7)
Name = 'S12H14SEAWBPG7'
NoOfStories = 16
Sections = [
            Section1, Section1, Section1, Section1,
            Section1, Section1, Section1,
            Section2, Section2, Section2,
            Section3, Section3, Section3,
            Section4, Section4, Section4,
            ]
S12H14SEAWBPG7 = CreateArchetype(Basements, Use2008Maps = False, Overstrength = 1.25)
Archetypes.append(S12H14SEAWBPG7)
#endregion

#region Archetype S16H14SEAPG7 and S16H14SEAWBPG7
Name = 'S16H14SEAPG7'
#### Input Variables
NoOfStories = 16
Thickness = 32.
Length = 24. * 12.
Flange_Thickness = 12*12. # Assume 6' Long Core
Long_Spacing = 4
BarSize = 8.0
Rho = 0.75 #In Fraction
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section1 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, 3., 4., Spacing, Thickness - 3.5)
BarSize = 7.0
Rho = 0.55 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section2 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, 3., 4., Spacing*2, Thickness - 3.5)
ThicknessBelow = float(Thickness)
Thickness = 22.
Length = Length - (ThicknessBelow - Thickness) * 2.
BarSize = 5.0
Rho = 0.50 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section3 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, 3., 4., Spacing*2, Thickness - 3.5)
BarSize = 4.0
Rho = 0.25 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section4 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, None, None, None, None)
Sections = [Section1, Section1, Section1, Section1,
            Section2, Section2, Section2, Section2,
            Section3, Section3, Section3, Section3,
            Section4, Section4, Section4, Section4,
            ]
S16H14SEAPG7 = CreateArchetype(Use2008Maps = False, Overstrength = 1.25)
Archetypes.append(S16H14SEAPG7)
Name = 'S16H14SEAWBPG7'
NoOfStories = 20 # Include Basement Floors Here
Sections = [Section1, Section1, Section1, Section1,
            Section1, Section1, Section1, Section1,
            Section2, Section2, Section2, Section2,
            Section3, Section3, Section3, Section3,
            Section4, Section4, Section4, Section4,
            ]
S16H14SEAWBPG7 = CreateArchetype(Basements, False, Overstrength = 1.25)
Archetypes.append(S16H14SEAWBPG7)
#endregion

#region Archetype S20H14SEAPG7 and S20H14SEAWBPG7
Name = 'S20H14SEAPG7'
# print 'Importing Archetype: ' + Name
#### Input Variables
NoOfStories = 20
Thickness = 40.
Length = 28. * 12.
Flange_Thickness = 14.*12. # Assume 6' Long Core

Long_Spacing = 4
BarSize = 8.0
Rho = 0.50 #In Fraction
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section1 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, 3., 4., Spacing, Thickness - 3.5)
BarSize = 7.0
Rho = 0.5 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section2 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, 3., 4., Spacing*2, Thickness - 3.5)
ThicknessBelow = float(Thickness)
Thickness = 24.
Length = Length - (ThicknessBelow - Thickness) * 2.
BarSize = 6.0
Rho =0.5 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section3 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, 3., 4., Spacing*2, Thickness - 3.5)

BarSize = 4.0
Rho = 0.25 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section4 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, None, None, None, None)
ThicknessBelow = float(Thickness)
Thickness = 22.
Length = Length - (ThicknessBelow - Thickness) * 2.
BarSize = 4.0
Rho = 0.25 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section5 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, None, None, None, None)
Sections = [Section1, Section1, Section1, Section1,
            Section2, Section2, Section2, Section2,
            Section3, Section3, Section3, Section3,
            Section4, Section4, Section4, Section4,
            Section5, Section5, Section5, Section5,
            ]
S20H14SEAPG7 = CreateArchetype(Use2008Maps = False, Overstrength = 1.25)
Archetypes.append(S20H14SEAPG7)
Name = 'S20H14SEAWBPG7'
NoOfStories = 24 # Include Basement Floors Here
Sections = [Section1, Section1, Section1, Section1,
            Section1, Section1, Section1, Section1,
            Section2, Section2, Section2, Section2,
            Section3, Section3, Section3, Section3,
            Section4, Section4, Section4, Section4,
            Section5, Section5, Section5, Section5,
            ]
S20H14SEAWBPG7 = CreateArchetype(Basements, False, Overstrength = 1.25)
Archetypes.append(S20H14SEAWBPG7)
#endregion

#region Archetype S24H14SEAPG7 and S24H14SEAWBPG7
Name = 'S24H14SEAPG7'
# print 'Importing Archetype: ' + Name
#### Input Variables
NoOfStories = 24
Thickness = 24.
Length = 35. * 12.
Flange_Thickness = 17.5*12. # Assume 6' Long Core
Long_Spacing = 4
BarSize = 6.0
Rho = 0.501
#In Fraction
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section1 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, 3., 4., Spacing, Thickness - 3.5)
BarSize = 6.0
Rho = 0.501 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section2 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, 3., 4., Spacing*2., Thickness - 3.5)
ThicknessBelow = float(Thickness)
Thickness = 24.
Length = Length - (ThicknessBelow - Thickness) * 2.
BarSize = 6.0
Rho = 0.501 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section3 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, 3., 4., Spacing*2., Thickness - 3.5)
BarSize = 4.0
Rho = 0.30 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section4 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu,   None, None, None, None)
ThicknessBelow = float(Thickness)
Thickness = 24.
Length = Length - (ThicknessBelow - Thickness) * 2.
BarSize = 4.0
Rho = 0.25 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section5 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, None, None, None, None)
BarSize = 4.0
Rho = 0.25 #In percentages
Abar = np.pi*(BarSize / 8. / 2.)**2.
Spacing = Abar * 2. / Thickness / Rho * 100
Section6 = IWallSection(Length, Flange_Thickness, Thickness, Rho, BarSize,
                        fpc_core, fy, fu, None, None, None, None)
Sections = [Section1, Section1, Section1, Section1,
            Section2, Section2, Section2, Section2,
            Section3, Section3, Section3, Section3,
            Section4, Section4, Section4, Section4,
            Section5, Section5, Section5, Section5,
            Section6, Section6, Section6, Section6,
            ]
S24H14SEAPG7 = CreateArchetype(Use2008Maps = False, Overstrength = 1.25)
Archetypes.append(S24H14SEAPG7)
Name = 'S24H14SEAWBPG7'
NoOfStories = 28
Sections = [Section1, Section1, Section1, Section1,
            Section1, Section1, Section1, Section1,
            Section2, Section2, Section2, Section2,
            Section3, Section3, Section3, Section3,
            Section4, Section4, Section4, Section4,
            Section5, Section5, Section5, Section5,
            Section6, Section6, Section6, Section6,
            ]
S24H14SEAWBPG7 = CreateArchetype(Basements, False, Overstrength = 1.25)
Archetypes.append(S24H14SEAWBPG7)
#endregion


####################################################################################
# endregion
####################################################################################

# import ATCWallArchetypeHelpers as WallHelper
# for arch in Archetypes:
#     print WallHelper.GetAxialLoadRatio(arch)

def GetArchetypeByName(Name):
    return [x for x in Archetypes if x.Name == Name][0]
