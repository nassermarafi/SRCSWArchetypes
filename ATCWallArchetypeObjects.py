
####################################################################################
#region Defining Classes
####################################################################################

from __future__ import absolute_import
import numpy as np
import ATCWallArchetypeHelpers as ATCWallHelper
import ASCEHelper
from six.moves import range

class ArchetypeData:
    def __init__(self, Name, YGrids, R, T1, l_w, t, b_f, rho, rho_t, fpc, fy, fu, GravityLoad, Mass, WallGravityLoad,
                 BoundaryElement=None, CouplingBeams=None, CouplingBeamLength=None, CustomWallSection=None, **kwargs):
        self.Name = Name
        self.YGrids = YGrids
        self.R = R
        self.T1 = T1
        self.l_w = l_w
        self.t = t
        self.b_f = b_f
        self.rho = rho
        self.rho_t = rho_t
        self.fpc = fpc
        self.fy = fy
        self.fu = fu
        self.GravityLoad = list(GravityLoad)
        self.Mass = list(Mass)
        self.WallGravityLoad = list(WallGravityLoad)
        self.BoundaryElement = BoundaryElement

        self.CustomSection = CustomWallSection

        self.fce = 1.3 * fpc  # TBI Table 7.1
        self.fye = 1.17 * fy  # TBI Table 7.1

        self.CouplingBeams = CouplingBeams
        self.CouplingBeamLength = CouplingBeamLength

        self.__dict__.update(kwargs)

class CouplingBeam:
    def __init__(self, b, h, fpc, fy, NoOfBarsX, NoOfBarsY, BarDia, DiagonalReinf=True, FaceBarSize=None,
                 FaceBarSpacing=None, TieSpacing = 6, TieDia = 5, no_of_ties_x=None, no_of_ties_y=None, **kwargs):
        """
        :param b: 
        :param h: 
        :param NoOfBarsX: 
        :param NoOfBarsY: 
        :param fpc: fpc of conc in ksi
        :param fy: fy of steel in ksi
        :param BarDia:
        :param DiagonalReinf:
        """

        self.b = b
        self.h = h
        self.fpc = fpc
        self.fy = fy
        self.NoOfBarsX = NoOfBarsX
        self.NoOfBarsY = NoOfBarsY
        self.BarDia = BarDia
        self.TieSpacing = TieSpacing
        self.TieDia = TieDia
        self.DiagonalReinf = DiagonalReinf

        self.FaceBarSize = FaceBarSize
        self.FaceBarSpacing = FaceBarSpacing

        self.no_of_ties_x = no_of_ties_x
        self.no_of_ties_y = no_of_ties_y


        self.fce = 1.3 * fpc  # TBI Table 7.1
        self.fye = 1.17 * fy  # TBI Table 7.1

        self.__dict__.update(kwargs)

class PlanarWallSection:
    def __init__(self, l_w, t_w, left_boundary, right_boundary, boundary_bar_size, left_reinf_layout,
                 right_reinf_layout,
                 web_rho, bar_size_web, fpc, fy, fu, boundary_tie_spacing, boundary_tie_bar,
                 boundary_tie_x_no, boundary_tie_y_no, **kwargs):
        self.l_w = l_w
        self.t_w = t_w
        self.right_boundary = right_boundary
        self.left_boundary = left_boundary
        self.boundary_bar_size = boundary_bar_size
        self.right_reinf_layout = right_reinf_layout
        self.left_reinf_layout = left_reinf_layout
        self.web_rho = web_rho
        self.bar_size_web = bar_size_web

        self.fu = fu
        self.fpc = fpc
        self.fy = fy
        self.boundary_tie_spacing = boundary_tie_spacing
        self.boundary_tie_bar = boundary_tie_bar
        self.boundary_tie_x_no = boundary_tie_x_no
        self.boundary_tie_y_no = boundary_tie_y_no

        self.fce = 1.3 * fpc  # TBI Table 7.1
        self.fye = 1.17 * fy  # TBI Table 7.1

        self.__dict__.update(kwargs)

class TWallSection:
    def __init__(self, l_w, b_w, t_w, rho, bar_size, fpc, fy, fu, tie_vertical_spacing, tie_bar_size,
                 tie_x_spacing, tie_y_spacing, flange_is_left = True,
                 boundary_length = None, boundary_reinf_layout = None, **kwargs):
        self.l_w = l_w
        self.t_w = t_w
        self.b_w = b_w

        self.flange_is_left = flange_is_left
        self.boundary_length = boundary_length
        self.boundary_reinf_layout = boundary_reinf_layout

        self.rho = rho
        self.bar_size = bar_size

        self.fu = fu
        self.fpc = fpc
        self.fy = fy

        self.tie_vertical_spacing = tie_vertical_spacing
        self.tie_bar_size = tie_bar_size
        self.tie_x_spacing = tie_x_spacing # along thickness
        self.tie_y_spacing = tie_y_spacing # along length

        self.fce = 1.3 * fpc  # TBI Table 7.1
        self.fye = 1.17 * fy  # TBI Table 7.1

        self.__dict__.update(kwargs)

class IWallSection:
    def __init__(self, l_w, b_w, t_w, rho, bar_size, fpc, fy, fu, tie_vertical_spacing,
                 tie_bar_size, tie_x_spacing, tie_y_spacing, **kwargs):
        self.l_w = l_w
        self.t_w = t_w
        self.b_w = b_w

        self.rho = rho
        self.bar_size = bar_size

        self.fu = fu
        self.fpc = fpc
        self.fy = fy

        self.tie_vertical_spacing = tie_vertical_spacing
        self.tie_bar_size = tie_bar_size
        self.tie_x_spacing = tie_x_spacing # along thickness
        self.tie_y_spacing = tie_y_spacing # along length

        self.fce = 1.3 * fpc  # TBI Table 7.1
        self.fye = 1.17 * fy  # TBI Table 7.1

        self.__dict__.update(kwargs)

####################################################################################
# endregion
####################################################################################

####################################################################################
#region Defining Archetypes
####################################################################################

### D Max Designs

# 8 Story Wall

#region ATCWALLS4O0S100ELF

# Geometry
NoOfStories = 4
YGrids = [0] + np.array(np.arange(0,(NoOfStories)*13*12, 13*12)+15*12).tolist()

# Loading
# DL 175psf Floors and 140psf Roof
# LL 65psf Floors and 20psf Roof

MassPerSqFt = (199.5)/1000. # use this 200 psf
FloorArea = 60.*120.
Mass = np.ones(NoOfStories) * MassPerSqFt * FloorArea
Mass[-1] = FloorArea*154.9/1000. # Adjust for Roof Weight

WallTribArea = 30.*38.
WeightPerSqFt = 212./1000.

BuildingWeight = np.ones(NoOfStories) * WeightPerSqFt * FloorArea
BuildingWeight[-1] = 152. / 1000. * FloorArea # Adjust for Roof Weight

PDeltaColumnGravityLoad = -1 * BuildingWeight * (FloorArea - WallTribArea) / FloorArea

WallGravityLoad = BuildingWeight * WallTribArea / FloorArea * -1
# WallGravityLoad = Mass * WallTribArea / FloorArea * -1 # This is what Kamal used in Atena

# Material Properties
fy = 60
fu = 105
fpc = 5.0

Section1 = PlanarWallSection(24.*12., 24., 46., 46., 9.0,
                             [3, 2, 2, 2, 2,
                              2, 2, 2, 2, 2,
                              3],
                             [3, 2, 2, 2, 2,
                              2, 2, 2, 2, 2,
                              3],
                             0.255, 4.037, fpc, fy, fu, 3, 4., 11, 3)

Section2 = PlanarWallSection(24.*12., 24., 38., 38., 8.0,
                             [3, 2, 2, 2, 2,
                              2, 2, 3],
                             [3, 2, 2, 2, 2,
                              2, 2, 3],
                             0.255, 4.037, fpc, fy, fu, 3, 4., 8, 3)


ATCWALLS4O0S100ELF = ArchetypeData('ATCWALLS4O0S100ELF', YGrids, None, None, None, None, None, None, None, fpc, fy, fu, PDeltaColumnGravityLoad,
                              Mass, WallGravityLoad,
                                   CustomWallSection=[Section1, Section1,
                                                     Section2, Section2,],
                                   heff = 41.*12)

#endregion

#region ATCWALLS4O4S50ELF

# Geometry
NoOfStories = 4
YGrids = [0] + np.array(np.arange(0,(NoOfStories)*13*12, 13*12)+15*12).tolist()

# Loading
# DL 175psf Floors and 140psf Roof
# LL 65psf Floors and 20psf Roof

MassPerSqFt = (199.5)/1000. # use this 200 psf
FloorArea = 60.*120.
Mass = np.ones(NoOfStories) * MassPerSqFt * FloorArea
Mass[-1] = FloorArea*154.9/1000. # Adjust for Roof Weight

WallTribArea = 30.*38.
WeightPerSqFt = 212./1000.

BuildingWeight = np.ones(NoOfStories) * WeightPerSqFt * FloorArea
BuildingWeight[-1] = 152. / 1000. * FloorArea # Adjust for Roof Weight

PDeltaColumnGravityLoad = -1 * BuildingWeight * (FloorArea - WallTribArea) / FloorArea

WallGravityLoad = BuildingWeight * WallTribArea / FloorArea * -1
# WallGravityLoad = Mass * WallTribArea / FloorArea * -1 # This is what Kamal used in Atena

# Material Properties
fy = 60.
fu = 105.
fpc = 5.

def FindConfinementRebarSpacing(NoOfBars):
    RebarLayout = []
    NoOfCols = int(int(NoOfBars-2)/2)
    for i in range(NoOfCols):
        if i == 0 or i == NoOfCols - 1:
            RebarLayout.append(3)
        else:
            RebarLayout.append(2)

    return RebarLayout

LLBE = 54
NoOfBarsLBE = 22
BarSizeLBE = 11.28
LBE_S_strirrup = 3.5

LSBE = 30
NoOfBarsSBE = 12
BarSizeSBE = 8.0
SBE_S_strirrup = 3.0

Section1A = PlanarWallSection(9.0*12., 24., LLBE, LSBE, BarSizeLBE,
                              FindConfinementRebarSpacing(NoOfBarsLBE),
                              FindConfinementRebarSpacing(NoOfBarsSBE),
                              0.255, 4., fpc, fy, fu, LBE_S_strirrup, 4., len(FindConfinementRebarSpacing(NoOfBarsLBE)), 3)

Section1B = PlanarWallSection(9.0*12., 24., LSBE, LLBE, BarSizeLBE,
                              FindConfinementRebarSpacing(NoOfBarsSBE),
                              FindConfinementRebarSpacing(NoOfBarsLBE),
                              0.255, 4., fpc, fy, fu, LBE_S_strirrup, 4., len(FindConfinementRebarSpacing(NoOfBarsSBE)), 3)

CB1 = CouplingBeam(24., 5.*12., fpc, fy, 4, 2, 11.274, True, 4., 6., 2.5, 4, 3, 10) #Last Row needs to be fixed
CB2 = CouplingBeam(24., 3.*12., fpc, fy, 4, 2, 11.274, True, 4., 6., 2.5, 4, 3, 10) #Last Row needs to be fixed
CB3 = CouplingBeam(24., 3.*12., fpc, fy, 4, 2, 8.0, True, 4., 6., 2.5, 4, 3, 10) #Last Row needs to be fixed

CouplingBeams = [CB1, CB2, CB2, CB3]
CouplingBeamLengths = [6.0*12., 6.0*12, 6.0*12, 6.0*12]

ATCWALLS4O4S50ELF = ArchetypeData('ATCWALLS4O4S50ELF', YGrids, None, None, None, None, None, None, None, fpc, fy, fu, PDeltaColumnGravityLoad,
                                Mass, WallGravityLoad,
                                CustomWallSection=[[Section1A, Section1B],
                                                   [Section1A, Section1B],
                                                   [Section1A, Section1B],
                                                   [Section1A, Section1B],
                                                   ],
                                CouplingBeamLength=CouplingBeamLengths,
                                CouplingBeams=CouplingBeams, heff = 41.*12.)

#endregion

# 8 Story Wall

#region ATCWALLS8O0S100ELF

# Geometry
NoOfStories = 8
YGrids = [0] + np.array(np.arange(0,(NoOfStories)*13*12, 13*12)+15*12).tolist()

# Loading
# DL 175psf Floors and 140psf Roof
# LL 65psf Floors and 20psf Roof

MassPerSqFt = (199.5)/1000. # use this 200 psf
FloorArea = 60.*120.
Mass = np.ones(NoOfStories) * MassPerSqFt * FloorArea
Mass[-1] = FloorArea*154.9/1000. # Adjust for Roof Weight

WallTribArea = 30.*44.
WeightPerSqFt = 212./1000.

BuildingWeight = np.ones(NoOfStories) * WeightPerSqFt * FloorArea
BuildingWeight[-1] = 152. / 1000. * FloorArea # Adjust for Roof Weight

PDeltaColumnGravityLoad = -1 * BuildingWeight * (FloorArea - WallTribArea) / FloorArea

WallGravityLoad = BuildingWeight * WallTribArea / FloorArea * -1
# WallGravityLoad = Mass * WallTribArea / FloorArea * -1 # This is what Kamal used in Atena

# Material Properties
fy = 60
fu = 105
fpc = 5.0

Section1 = PlanarWallSection(30.*12., 24., 72, 72, 10.173,
                             [3, 2, 2, 2, 2,
                              2, 2, 2, 2, 2,
                              2, 2, 2, 2, 2, 3],
                             [3, 2, 2, 2, 2,
                              2, 2, 2, 2, 2,
                              2, 2, 2, 2, 2, 3],
                             0.255, 4.037, fpc, fy, fu, 3, 4., 16, 3)

Section2 = PlanarWallSection(30.*12., 24., 36, 36, 10.173,
                             [3, 2, 2, 2, 2,
                              2, 2, 3],
                             [3, 2, 2, 2, 2,
                              2, 2, 3],
                             0.255, 4.037, fpc, fy, fu, 3, 4., 8, 3)

Section3 = PlanarWallSection(30.*12., 24., 0, 0, 10.173,
                             [],
                             [],
                             0.255, 4.037, fpc, fy, fu, None, None, None, None)

ATCWALLS8O0S100ELF = ArchetypeData('ATCWALLS8O0S100ELF', YGrids, None, None, None, None, None, None, None, fpc, fy, fu, PDeltaColumnGravityLoad,
                              Mass, WallGravityLoad,
                                   CustomWallSection=[Section1, Section1, Section1,
                                                     Section2, Section2, Section2,
                                                     Section3, Section3],
                                   heff = 73.57*12.)

#endregion

#region ATCWALLS8O0S100MRSA

# Geometry
NoOfStories = 8
YGrids = [0] + np.array(np.arange(0,(NoOfStories)*13*12, 13*12)+15*12).tolist()

# Loading
# DL 175psf Floors and 140psf Roof
# LL 65psf Floors and 20psf Roof

MassPerSqFt = (199.5)/1000. # use this 200 psf
FloorArea = 60.*120.
Mass = np.ones(NoOfStories) * MassPerSqFt * FloorArea
Mass[-1] = FloorArea*154.9/1000. # Adjust for Roof Weight

WallTribArea = 30.*44.
WeightPerSqFt = 212./1000.

BuildingWeight = np.ones(NoOfStories) * WeightPerSqFt * FloorArea
BuildingWeight[-1] = 152. / 1000. * FloorArea # Adjust for Roof Weight

PDeltaColumnGravityLoad = -1 * BuildingWeight * (FloorArea - WallTribArea) / FloorArea

WallGravityLoad = BuildingWeight * WallTribArea / FloorArea * -1
# WallGravityLoad = Mass * WallTribArea / FloorArea * -1 # This is what Kamal used in Atena

# Material Properties
fy = 60
fu = 105
fpc = 5.0

Section1 = PlanarWallSection(30.*12., 24., 72, 72, 9.027,
                             [3, 2, 2, 2, 2,
                              2, 2, 2, 2, 2,
                              2, 2, 2, 2, 3],
                             [3, 2, 2, 2, 2,
                              2, 2, 2, 2, 2,
                              2, 2, 2, 2, 3],
                             0.255, 4.037, fpc, fy, fu, 3, 4., 15, 3)

Section2 = PlanarWallSection(30.*12., 24., 36, 36, 9.027,
                             [3, 2, 2, 2, 2,
                              3],
                             [3, 2, 2, 2, 2,
                              3],
                             0.255, 4.037, fpc, fy, fu, 3, 4., 6, 3)

Section3 = PlanarWallSection(30.*12., 24., 0, 0, 9.027,
                             [],
                             [],
                             0.255, 4.037, fpc, fy, fu, None, None, None, None)

ATCWALLS8O0S100MRSA = ArchetypeData('ATCWALLS8O0S100MRSA', YGrids, None, None, None, None, None, None, None, fpc, fy, fu, PDeltaColumnGravityLoad,
                              Mass, WallGravityLoad,
                                   CustomWallSection=[Section1, Section1, Section1,
                                                     Section2, Section2, Section2,
                                                     Section3, Section3],
                                   heff = 73.57*12.)

#endregion

#region ATCWALLS8O1S50ELF

# Geometry
NoOfStories = 8
YGrids = [0] + np.array(np.arange(0,(NoOfStories)*13*12, 13*12)+15*12).tolist()

# Loading
# DL 175psf Floors and 140psf Roof
# LL 65psf Floors and 20psf Roof

MassPerSqFt = (199.5)/1000. # use this 200 psf
FloorArea = 60.*120.
Mass = np.ones(NoOfStories) * MassPerSqFt * FloorArea
Mass[-1] = FloorArea*154.9/1000. # Adjust for Roof Weight

WallTribArea = 30.*44.
WeightPerSqFt = 212./1000.

BuildingWeight = np.ones(NoOfStories) * WeightPerSqFt * FloorArea
BuildingWeight[-1] = 152. / 1000. * FloorArea # Adjust for Roof Weight

PDeltaColumnGravityLoad = -1 * BuildingWeight * (FloorArea - WallTribArea) / FloorArea

WallGravityLoad = BuildingWeight * WallTribArea / FloorArea * -1
# WallGravityLoad = Mass * WallTribArea / FloorArea * -1 # This is what Kamal used in Atena

# Material Properties
fy = 60.
fu = 105.
fpc = 5.

Section1A = PlanarWallSection(11.25*12., 24., 60, 30, 11.274,
                             [3, 2, 2, 2, 2,
                              2, 2, 2, 2, 2,
                              2, 2, 2, 2, 2, 3],
                             [3, 2, 2, 2, 2, 3],
                             0.255, 4., fpc, fy, fu, 3., 4., 16, 3)

Section1B = PlanarWallSection(11.25*12., 24., 30, 60, 11.274,
                             [3, 2, 2, 2, 2, 3],
                             [3, 2, 2, 2, 2,
                              2, 2, 2, 2, 2,
                              2, 2, 2, 2, 2, 3],
                             0.255, 4., fpc, fy, fu, 3., 4., 6, 3)

Section1 = PlanarWallSection(30.*12., 24., 60, 60, 10.173,
                             [3, 2, 2, 2, 2,
                              2, 2, 2, 2, 2,
                              2, 2, 2, 2, 2, 3],
                             [3, 2, 2, 2, 2,
                              2, 2, 2, 2, 2,
                              2, 2, 2, 2, 2, 3],
                             0.255, 4., fpc, fy, fu, 3, 4., 16, 3)

Section2 = PlanarWallSection(30.*12., 24., 36, 36, 10.173,
                             [3, 2, 2, 2, 2,
                              2, 2, 3],
                             [3, 2, 2, 2, 2,
                              2, 2, 3],
                             0.255, 4., fpc, fy, fu, 3, 4., 8, 3)

Section3 = PlanarWallSection(30.*12., 24., 0, 0, 10.173,
                             [],
                             [],
                             0.255, 4., fpc, fy, fu, None, None, None, None)

CouplingBeams = [None, None, None, None, None, None, None, None]
CouplingBeamLengths = [7.5*12., None, None, None, None, None, None, None]

ATCWALLS8O1S50ELF = ArchetypeData('ATCWALLS8O1S50ELF', YGrids, None, None, None, None, None, None, None, fpc, fy, fu, PDeltaColumnGravityLoad,
                                Mass, WallGravityLoad,
                                CustomWallSection=[[Section1A, Section1B], Section1, Section1,
                                                     Section2, Section2, Section2,
                                                     Section3, Section3],
                                CouplingBeamLength=CouplingBeamLengths,
                                CouplingBeams=CouplingBeams, heff = 73.57*12.)

#endregion

#region ATCWALLS8O1S50MRSA

# Geometry
NoOfStories = 8
YGrids = [0] + np.array(np.arange(0,(NoOfStories)*13*12, 13*12)+15*12).tolist()

# Loading
# DL 175psf Floors and 140psf Roof
# LL 65psf Floors and 20psf Roof

MassPerSqFt = (199.5)/1000. # use this 200 psf
FloorArea = 60.*120.
Mass = np.ones(NoOfStories) * MassPerSqFt * FloorArea
Mass[-1] = FloorArea*154.9/1000. # Adjust for Roof Weight

WallTribArea = 30.*44.
WeightPerSqFt = 212./1000.

BuildingWeight = np.ones(NoOfStories) * WeightPerSqFt * FloorArea
BuildingWeight[-1] = 152. / 1000. * FloorArea # Adjust for Roof Weight

PDeltaColumnGravityLoad = -1 * BuildingWeight * (FloorArea - WallTribArea) / FloorArea

WallGravityLoad = BuildingWeight * WallTribArea / FloorArea * -1
# WallGravityLoad = Mass * WallTribArea / FloorArea * -1 # This is what Kamal used in Atena


# Material Properties
fy = 60.
fu = 105.
fpc = 5.

Section1A = PlanarWallSection(11.25*12., 24., 60, 30, 10.173,
                             [3, 2, 2, 2, 2,
                              2, 2, 2, 2, 2,
                              2, 2, 2, 2, 2, 3],
                             [3, 2, 2, 2, 2, 3],
                             0.255, 4., fpc, fy, fu, 3., 4., 16, 3)

Section1B = PlanarWallSection(11.25*12., 24., 30, 60, 10.173,
                             [3, 2, 2, 2, 2, 3],
                             [3, 2, 2, 2, 2,
                              2, 2, 2, 2, 2,
                              2, 2, 2, 2, 2, 3],
                             0.255, 4., fpc, fy, fu, 3., 4., 6, 3)

Section1 = PlanarWallSection(30.*12., 24., 60, 60, 10.173,
                             [3, 2, 2, 2, 2,
                              2, 2, 2, 2, 2,
                              2, 2, 2, 2, 2, 3],
                             [3, 2, 2, 2, 2,
                              2, 2, 2, 2, 2,
                              2, 2, 2, 2, 2, 3],
                             0.255, 4., fpc, fy, fu, 3, 4., 16, 3)

Section2 = PlanarWallSection(30.*12., 24., 36, 36, 9.027,
                             [3, 2, 2, 2, 2,
                              3],
                             [3, 2, 2, 2, 2,
                              3],
                             0.255, 4., fpc, fy, fu, 3, 4., 6, 3)

Section3 = PlanarWallSection(30.*12., 24., 0, 0, 9.027,
                             [],
                             [],
                             0.255, 4., fpc, fy, fu, None, None, None, None)

CouplingBeams = [None, None, None, None, None, None, None, None]
CouplingBeamLengths = [7.5*12., None, None, None, None, None, None, None]

ATCWALLS8O1S50MRSA = ArchetypeData('ATCWALLS8O1S50MRSA', YGrids, None, None, None, None, None, None, None, fpc, fy, fu, PDeltaColumnGravityLoad,
                                Mass, WallGravityLoad,
                                CustomWallSection=[[Section1A, Section1B], Section1, Section1,
                                                     Section2, Section2, Section2,
                                                     Section3, Section3],
                                CouplingBeamLength=CouplingBeamLengths,
                                CouplingBeams=CouplingBeams, heff = 73.57*12.)

#endregion

#region ATCWALLS8O1S75ELF

# Geometry
NoOfStories = 8
YGrids = [0] + np.array(np.arange(0,(NoOfStories)*13*12, 13*12)+15*12).tolist()

# Loading
# DL 175psf Floors and 140psf Roof
# LL 65psf Floors and 20psf Roof

MassPerSqFt = (199.5)/1000. # use this 200 psf
FloorArea = 60.*120.
Mass = np.ones(NoOfStories) * MassPerSqFt * FloorArea
Mass[-1] = FloorArea*154.9/1000. # Adjust for Roof Weight

WallTribArea = 30.*44.
WeightPerSqFt = 212./1000.

BuildingWeight = np.ones(NoOfStories) * WeightPerSqFt * FloorArea
BuildingWeight[-1] = 152. / 1000. * FloorArea # Adjust for Roof Weight

PDeltaColumnGravityLoad = -1 * BuildingWeight * (FloorArea - WallTribArea) / FloorArea

WallGravityLoad = BuildingWeight * WallTribArea / FloorArea * -1
# WallGravityLoad = Mass * WallTribArea / FloorArea * -1 # This is what Kamal used in Atena


# Material Properties
fy = 60.
fu = 105.
fpc = 5.

Section1A = PlanarWallSection(13.75*12., 24., 60, 30, 11.274,
                             [3, 2, 2, 2, 2,
                              2, 2, 2, 2, 2,
                              2, 2, 2, 2, 3],
                             [3, 2, 2, 2, 2, 3],
                             0.255, 4., fpc, fy, fu, 3., 4., 15, 3)

Section1B = PlanarWallSection(13.75*12., 24., 30, 60, 11.274,
                             [3, 2, 2, 2, 2, 3],
                             [3, 2, 2, 2, 2,
                              2, 2, 2, 2, 2,
                              2, 2, 2, 2, 3],
                             0.255, 4., fpc, fy, fu, 3., 4., 6, 3)

Section1 = PlanarWallSection(30.*12., 24., 60, 60, 10.173,
                             [3, 2, 2, 2, 2,
                              2, 2, 2, 2, 2,
                              2, 2, 2, 2, 2, 3],
                             [3, 2, 2, 2, 2,
                              2, 2, 2, 2, 2,
                              2, 2, 2, 2, 2, 3],
                             0.255, 4., fpc, fy, fu, 3, 4., 16, 3)

Section2 = PlanarWallSection(30.*12., 24., 36, 36, 10.173,
                             [3, 2, 2, 2, 2,
                              2, 2, 3],
                             [3, 2, 2, 2, 2,
                              2, 2, 3],
                             0.255, 4., fpc, fy, fu, 3, 4., 8, 3)

Section3 = PlanarWallSection(30.*12., 24., 0, 0, 10.173,
                             [],
                             [],
                             0.255, 4., fpc, fy, fu, None, None, None, None)

CouplingBeams = [None, None, None, None, None, None, None, None]
CouplingBeamLengths = [2.5*12., None, None, None, None, None, None, None]

ATCWALLS8O1S75ELF = ArchetypeData('ATCWALLS8O1S75ELF', YGrids, None, None, None, None, None, None, None, fpc, fy, fu, PDeltaColumnGravityLoad,
                                Mass, WallGravityLoad,
                                CustomWallSection=[[Section1A, Section1B], Section1, Section1,
                                                     Section2, Section2, Section2,
                                                     Section3, Section3],
                                CouplingBeamLength=CouplingBeamLengths,
                                CouplingBeams=CouplingBeams, heff = 73.57*12.)

#endregion

#region ATCWALLS8O1S75MRSA

# Geometry
NoOfStories = 8
YGrids = [0] + np.array(np.arange(0,(NoOfStories)*13*12, 13*12)+15*12).tolist()

# Loading
# DL 175psf Floors and 140psf Roof
# LL 65psf Floors and 20psf Roof

MassPerSqFt = (199.5)/1000. # use this 200 psf
FloorArea = 60.*120.
Mass = np.ones(NoOfStories) * MassPerSqFt * FloorArea
Mass[-1] = FloorArea*154.9/1000. # Adjust for Roof Weight

WallTribArea = 30.*44.
WeightPerSqFt = 212./1000.

BuildingWeight = np.ones(NoOfStories) * WeightPerSqFt * FloorArea
BuildingWeight[-1] = 152. / 1000. * FloorArea # Adjust for Roof Weight

PDeltaColumnGravityLoad = -1 * BuildingWeight * (FloorArea - WallTribArea) / FloorArea

WallGravityLoad = BuildingWeight * WallTribArea / FloorArea * -1
# WallGravityLoad = Mass * WallTribArea / FloorArea * -1 # This is what Kamal used in Atena


# Material Properties
fy = 60.
fu = 105.
fpc = 5.

Section1A = PlanarWallSection(13.75*12., 24., 60, 30, 10.173,
                             [3, 2, 2, 2, 2,
                              2, 2, 2, 2, 2,
                              2, 2, 2, 2, 3],
                             [3, 2, 2, 2, 2, 3],
                             0.255, 4., fpc, fy, fu, 3., 4., 15, 3)

Section1B = PlanarWallSection(13.75*12., 24., 30, 60, 10.173,
                             [3, 2, 2, 2, 2, 3],
                             [3, 2, 2, 2, 2,
                              2, 2, 2, 2, 2,
                              2, 2, 2, 2, 3],
                             0.255, 4., fpc, fy, fu, 3., 4., 6, 3)

Section1 = PlanarWallSection(30.*12., 24., 60, 60, 10.173,
                             [3, 2, 2, 2, 2,
                              2, 2, 2, 2, 2,
                              2, 2, 2, 2, 3],
                             [3, 2, 2, 2, 2,
                              2, 2, 2, 2, 2,
                              2, 2, 2, 2, 3],
                             0.255, 4., fpc, fy, fu, 3, 4., 15, 3)

Section2 = PlanarWallSection(30.*12., 24., 36, 36, 9.027,
                             [3, 2, 2, 2, 2,
                              3],
                             [3, 2, 2, 2, 2,
                              3],
                             0.255, 4., fpc, fy, fu, 3, 4., 6, 3)

Section3 = PlanarWallSection(30.*12., 24., 0, 0, 10.173,
                             [],
                             [],
                             0.255, 4., fpc, fy, fu, None, None, None, None)

CouplingBeams = [None, None, None, None, None, None, None, None]
CouplingBeamLengths = [2.5*12., None, None, None, None, None, None, None]

ATCWALLS8O1S75MRSA = ArchetypeData('ATCWALLS8O1S75MRSA', YGrids, None, None, None, None, None, None, None, fpc, fy, fu, PDeltaColumnGravityLoad,
                                Mass, WallGravityLoad,
                                CustomWallSection=[[Section1A, Section1B], Section1, Section1,
                                                     Section2, Section2, Section2,
                                                     Section3, Section3],
                                CouplingBeamLength=CouplingBeamLengths,
                                CouplingBeams=CouplingBeams, heff = 73.57*12.)

#endregion

#region ATCWALLS8O2S50ELF

# Geometry
NoOfStories = 8
YGrids = [0] + np.array(np.arange(0,(NoOfStories)*13*12, 13*12)+15*12).tolist()

# Loading
# DL 175psf Floors and 140psf Roof
# LL 65psf Floors and 20psf Roof

MassPerSqFt = (199.5)/1000. # use this 200 psf
FloorArea = 60.*120.
Mass = np.ones(NoOfStories) * MassPerSqFt * FloorArea
Mass[-1] = FloorArea*154.9/1000. # Adjust for Roof Weight

WallTribArea = 30.*44.
WeightPerSqFt = 212./1000.

BuildingWeight = np.ones(NoOfStories) * WeightPerSqFt * FloorArea
BuildingWeight[-1] = 152. / 1000. * FloorArea # Adjust for Roof Weight

PDeltaColumnGravityLoad = -1 * BuildingWeight * (FloorArea - WallTribArea) / FloorArea

WallGravityLoad = BuildingWeight * WallTribArea / FloorArea * -1
# WallGravityLoad = Mass * WallTribArea / FloorArea * -1 # This is what Kamal used in Atena

# Material Properties
fy = 60.
fu = 105.
fpc = 5.

Section1A = PlanarWallSection(11.25*12., 24., 60, 30, 11.274,
                             [3, 2, 2, 2, 2,
                              2, 2, 2, 2, 2,
                              2, 2, 2, 2, 2,
                              2, 3],
                             [3, 2, 2, 2, 2,
                              3],
                             0.255, 4., fpc, fy, fu, 3., 4., 17, 3)

Section1B = PlanarWallSection(11.25*12., 24., 30, 60, 11.274,
                             [3, 2, 2, 2, 2, 3],
                             [3, 2, 2, 2, 2,
                              2, 2, 2, 2, 2,
                              2, 2, 2, 2, 2,
                              2, 3],
                              0.255, 4., fpc, fy, fu, 3., 4., 6, 3)

Section1 = PlanarWallSection(30.*12., 24., 60, 60, 11.274,
                             [3, 2, 2, 2, 2,
                              2, 2, 2, 2, 2,
                              2, 2, 2, 2, 2,
                              2, 3],
                             [3, 2, 2, 2, 2,
                              2, 2, 2, 2, 2,
                              2, 2, 2, 2, 2,
                              2, 3],
                             0.255, 4., fpc, fy, fu, 3, 4., 17, 3)

Section2 = PlanarWallSection(30.*12., 24., 30, 30, 9.027,
                             [3, 2, 2, 2, 2,
                              3],
                             [3, 2, 2, 2, 2,
                              3],
                              0.255, 4., fpc, fy, fu, 3, 4., 6, 3)

Section3 = PlanarWallSection(30.*12., 24., 0, 0, 9.027,
                             [],
                             [],
                             0.255, 4., fpc, fy, fu, None, None, None, None)

CB1 = CouplingBeam(24., 60., fpc, fy, 4, 2, 11.274, True, 4., 6., 2.5, 4, 3, 10)

CouplingBeams = [CB1, None, None, None, None, None, None, None]
CouplingBeamLengths = [7.5*12., 7.5*12, None, None, None, None, None, None]

ATCWALLS8O2S50ELF = ArchetypeData('ATCWALLS8O2S50ELF', YGrids, None, None, None, None, None, None, None, fpc, fy, fu, PDeltaColumnGravityLoad,
                                Mass, WallGravityLoad,
                                CustomWallSection=[[Section1A, Section1B], [Section1A, Section1B], Section1,
                                                     Section2, Section2, Section2,
                                                     Section3, Section3],
                                CouplingBeamLength=CouplingBeamLengths,
                                CouplingBeams=CouplingBeams, heff = 73.57*12.)

#endregion

#region ATCWALLS8O2S50MRSA

# Geometry
NoOfStories = 8
YGrids = [0] + np.array(np.arange(0,(NoOfStories)*13*12, 13*12)+15*12).tolist()

# Loading
# DL 175psf Floors and 140psf Roof
# LL 65psf Floors and 20psf Roof

MassPerSqFt = (199.5)/1000. # use this 200 psf
FloorArea = 60.*120.
Mass = np.ones(NoOfStories) * MassPerSqFt * FloorArea
Mass[-1] = FloorArea*154.9/1000. # Adjust for Roof Weight

WallTribArea = 30.*44.
WeightPerSqFt = 212./1000.

BuildingWeight = np.ones(NoOfStories) * WeightPerSqFt * FloorArea
BuildingWeight[-1] = 152. / 1000. * FloorArea # Adjust for Roof Weight

PDeltaColumnGravityLoad = -1 * BuildingWeight * (FloorArea - WallTribArea) / FloorArea

WallGravityLoad = BuildingWeight * WallTribArea / FloorArea * -1
# WallGravityLoad = Mass * WallTribArea / FloorArea * -1 # This is what Kamal used in Atena

# Material Properties
fy = 60.
fu = 105.
fpc = 5.

Section1A = PlanarWallSection(11.25*12., 24., 60, 30, 11.274,
                             [3, 2, 2, 2, 2,
                              2, 2, 2, 2, 2,
                              2, 2, 2, 2, 3],
                             [3, 2, 2, 2, 2,
                              3],
                             0.255, 4., fpc, fy, fu, 3., 4., 15, 3)

Section1B = PlanarWallSection(11.25*12., 24., 30, 60, 11.274,
                             [3, 2, 2, 2, 2,
                              3],
                             [3, 2, 2, 2, 2,
                              2, 2, 2, 2, 2,
                              2, 2, 2, 2, 3],
                              0.255, 4., fpc, fy, fu, 3., 4., 6, 3)

Section1 = PlanarWallSection(30.*12., 24., 60, 60, 11.274,
                             [3, 2, 2, 2, 2,
                              2, 2, 2, 2, 2,
                              2, 2, 2, 2, 3],
                             [3, 2, 2, 2, 2,
                              2, 2, 2, 2, 2,
                              2, 2, 2, 2, 3],
                             0.255, 4., fpc, fy, fu, 3, 4., 15, 3)

Section2 = PlanarWallSection(30.*12., 24., 30, 30, 8.0,
                             [3, 2, 2, 2, 3],
                             [3, 2, 2, 2, 3],
                              0.255, 4., fpc, fy, fu, 3, 4., 6, 3)

Section3 = PlanarWallSection(30.*12., 24., 0, 0, 9.027,
                             [],
                             [],
                             0.255, 4., fpc, fy, fu, None, None, None, None)

CB1 = CouplingBeam(24., 60., fpc, fy, 4, 2, 11.274, True, 4., 6., 2.5, 4, 3, 10)

CouplingBeams = [CB1, None, None, None, None, None, None, None]
CouplingBeamLengths = [7.5*12., 7.5*12, None, None, None, None, None, None]

ATCWALLS8O2S50MRSA = ArchetypeData('ATCWALLS8O2S50MRSA', YGrids, None, None, None, None, None, None, None, fpc, fy, fu, PDeltaColumnGravityLoad,
                                Mass, WallGravityLoad,
                                CustomWallSection=[[Section1A, Section1B], [Section1A, Section1B], Section1,
                                                     Section2, Section2, Section2,
                                                     Section3, Section3],
                                CouplingBeamLength=CouplingBeamLengths,
                                CouplingBeams=CouplingBeams, heff = 73.57*12.)

#endregion

#region ATCWALLS8O2S75ELF

# Geometry
NoOfStories = 8
YGrids = [0] + np.array(np.arange(0,(NoOfStories)*13*12, 13*12)+15*12).tolist()

# Loading
# DL 175psf Floors and 140psf Roof
# LL 65psf Floors and 20psf Roof

MassPerSqFt = (199.5)/1000. # use this 200 psf
FloorArea = 60.*120.
Mass = np.ones(NoOfStories) * MassPerSqFt * FloorArea
Mass[-1] = FloorArea*154.9/1000. # Adjust for Roof Weight

WallTribArea = 30.*44.
WeightPerSqFt = 212./1000.

BuildingWeight = np.ones(NoOfStories) * WeightPerSqFt * FloorArea
BuildingWeight[-1] = 152. / 1000. * FloorArea # Adjust for Roof Weight

PDeltaColumnGravityLoad = -1 * BuildingWeight * (FloorArea - WallTribArea) / FloorArea

WallGravityLoad = BuildingWeight * WallTribArea / FloorArea * -1
# WallGravityLoad = Mass * WallTribArea / FloorArea * -1 # This is what Kamal used in Atena

# Material Properties
fy = 60.
fu = 105.
fpc = 5.

Section1A = PlanarWallSection(13.75*12., 24., 60, 30, 11.274,
                             [3, 2, 2, 2, 2,
                              2, 2, 2, 2, 2,
                              2, 2, 2, 2, 2, 3],
                             [3, 2, 2, 2, 2, 3],
                             0.255, 4., fpc, fy, fu, 3., 4., 16, 3)

Section1B = PlanarWallSection(13.75*12., 24., 30, 60, 11.274,
                             [3, 2, 2, 2, 2, 3],
                             [3, 2, 2, 2, 2,
                              2, 2, 2, 2, 2,
                              2, 2, 2, 2, 2, 3],
                             0.255, 4., fpc, fy, fu, 3., 4., 6, 3)

Section1 = PlanarWallSection(30.*12., 24., 60, 60, 11.274,
                             [3, 2, 2, 2, 2,
                              2, 2, 2, 2, 2,
                              2, 2, 2, 2, 2, 3],
                             [3, 2, 2, 2, 2,
                              2, 2, 2, 2, 2,
                              2, 2, 2, 2, 2, 3],
                             0.255, 4., fpc, fy, fu, 3, 4., 16, 3)

Section2 = PlanarWallSection(30.*12., 24., 30, 30, 9.027,
                             [3, 2, 2, 2, 2, 2, 2, 3],
                             [3, 2, 2, 2, 2, 2, 2, 3],
                             0.255, 4., fpc, fy, fu, 3, 4., 8, 3)

Section3 = PlanarWallSection(30.*12., 24., 0, 0, 9.027,
                             [],
                             [],
                             0.255, 4., fpc, fy, fu, None, None, None, None)

CB1 = CouplingBeam(24., 60., fpc, fy, 4, 2, 10.173, True, 4., 6., 3., 5, 3, 10)

CouplingBeams = [CB1, None, None, None, None, None, None, None]
CouplingBeamLengths = [2.5*12., 2.5*12, None, None, None, None, None, None]

ATCWALLS8O2S75ELF = ArchetypeData('ATCWALLS8O2S75ELF', YGrids, None, None, None, None, None, None, None, fpc, fy, fu, PDeltaColumnGravityLoad,
                                Mass, WallGravityLoad,
                                CustomWallSection=[[Section1A, Section1B], [Section1A, Section1B], Section1,
                                                     Section2, Section2, Section2,
                                                     Section3, Section3],
                                CouplingBeamLength=CouplingBeamLengths,
                                CouplingBeams=CouplingBeams, heff = 73.57*12.)

#endregion

#region ATCWALLS8O2S75MRSA

# Geometry
NoOfStories = 8
YGrids = [0] + np.array(np.arange(0,(NoOfStories)*13*12, 13*12)+15*12).tolist()

# Loading
# DL 175psf Floors and 140psf Roof
# LL 65psf Floors and 20psf Roof

MassPerSqFt = (199.5)/1000. # use this 200 psf
FloorArea = 60.*120.
Mass = np.ones(NoOfStories) * MassPerSqFt * FloorArea
Mass[-1] = FloorArea*154.9/1000. # Adjust for Roof Weight

WallTribArea = 30.*44.
WeightPerSqFt = 212./1000.

BuildingWeight = np.ones(NoOfStories) * WeightPerSqFt * FloorArea
BuildingWeight[-1] = 152. / 1000. * FloorArea # Adjust for Roof Weight

PDeltaColumnGravityLoad = -1 * BuildingWeight * (FloorArea - WallTribArea) / FloorArea

WallGravityLoad = BuildingWeight * WallTribArea / FloorArea * -1
# WallGravityLoad = Mass * WallTribArea / FloorArea * -1 # This is what Kamal used in Atena

# Material Properties
fy = 60.
fu = 105.
fpc = 5.

Section1A = PlanarWallSection(13.75*12., 24., 60, 30, 10.173,
                             [3, 2, 2, 2, 2,
                              2, 2, 2, 2, 2,
                              2, 2, 2, 3],
                             [3, 2, 2, 2, 2,
                              3],
                             0.255, 4., fpc, fy, fu, 3., 4., 14, 3)

Section1B = PlanarWallSection(13.75*12., 24., 30, 60, 10.173,
                             [3, 2, 2, 2, 2,
                              3],
                             [3, 2, 2, 2, 2,
                              2, 2, 2, 2, 2,
                              2, 2, 2, 3],
                             0.255, 4., fpc, fy, fu, 3., 4., 6, 3)

Section1 = PlanarWallSection(30.*12., 24., 60, 60, 10.173,
                             [3, 2, 2, 2, 2,
                              2, 2, 2, 2, 2,
                              2, 2, 2, 3],
                             [3, 2, 2, 2, 2,
                              2, 2, 2, 2, 2,
                              2, 2, 2, 3],
                             0.255, 4., fpc, fy, fu, 3, 4., 14, 3)

Section2 = PlanarWallSection(30.*12., 24., 36, 36, 8.0,
                             [3, 2, 2, 2, 3],
                             [3, 2, 2, 2, 3],
                             0.255, 4., fpc, fy, fu, 3, 4., 5, 3)

Section3 = PlanarWallSection(30.*12., 24., 0, 0, 8.0,
                             [],
                             [],
                             0.255, 4., fpc, fy, fu, None, None, None, None)

CB1 = CouplingBeam(24., 60., fpc, fy, 4, 2, 10.173, True, 4., 6., 3., 5, 3, 10) #drawings say 11 but i think its a mistake.

CouplingBeams = [CB1, None, None, None, None, None, None, None]
CouplingBeamLengths = [2.5*12., 2.5*12, None, None, None, None, None, None]

ATCWALLS8O2S75MRSA = ArchetypeData('ATCWALLS8O2S75MRSA', YGrids, None, None, None, None, None, None, None, fpc, fy, fu, PDeltaColumnGravityLoad,
                                Mass, WallGravityLoad,
                                CustomWallSection=[[Section1A, Section1B], [Section1A, Section1B], Section1,
                                                     Section2, Section2, Section2,
                                                     Section3, Section3],
                                CouplingBeamLength=CouplingBeamLengths,
                                CouplingBeams=CouplingBeams, heff = 73.57*12.)

#endregion

#region ATCWALLS8O8S50ELF

# Geometry
NoOfStories = 8
YGrids = [0] + np.array(np.arange(0,(NoOfStories)*13*12, 13*12)+15*12).tolist()

# Loading
# DL 175psf Floors and 140psf Roof
# LL 65psf Floors and 20psf Roof

MassPerSqFt = (199.5)/1000. # use this 200 psf
FloorArea = 60.*120.
Mass = np.ones(NoOfStories) * MassPerSqFt * FloorArea
Mass[-1] = FloorArea*154.9/1000. # Adjust for Roof Weight

WallTribArea = 30.*44.
WeightPerSqFt = 212./1000.

BuildingWeight = np.ones(NoOfStories) * WeightPerSqFt * FloorArea
BuildingWeight[-1] = 152. / 1000. * FloorArea # Adjust for Roof Weight

PDeltaColumnGravityLoad = -1 * BuildingWeight * (FloorArea - WallTribArea) / FloorArea

WallGravityLoad = BuildingWeight * WallTribArea / FloorArea * -1
# WallGravityLoad = Mass * WallTribArea / FloorArea * -1 # This is what Kamal used in Atena

# Material Properties
fy = 60.
fu = 105.
fpc = 5.

def FindConfinementRebarSpacing(NoOfBars):
    RebarLayout = []
    NoOfCols = int(int(NoOfBars-2)/2)
    for i in range(NoOfCols):
        if i == 0 or i == NoOfCols - 1:
            RebarLayout.append(3)
        else:
            RebarLayout.append(2)

    return RebarLayout

LLBE = 72
NoOfBarsLBE = 42.
BarSizeLBE = 11.28
LBE_S_strirrup = 3.5

LSBE = 30
NoOfBarsSBE = 12
BarSizeSBE = 11.28
SBE_S_strirrup = 3.5

Section1A = PlanarWallSection(11.25*12., 24., LLBE, LSBE, BarSizeLBE,
                              np.ones(14)*3,
                              FindConfinementRebarSpacing(NoOfBarsSBE),
                              0.255, 4., fpc, fy, fu, LBE_S_strirrup, 4., 13, 3)

Section1B = PlanarWallSection(11.25*12., 24., LSBE, LLBE, BarSizeLBE,
                              FindConfinementRebarSpacing(NoOfBarsSBE),
                              np.ones(14) * 3,
                              0.255, 4., fpc, fy, fu, LBE_S_strirrup, 4., len(FindConfinementRebarSpacing(NoOfBarsSBE)), 3)

LLBE = 42
NoOfBarsLBE = 22
BarSizeLBE = 11.28
LBE_S_strirrup = 3.5

LSBE = 30
NoOfBarsSBE = 12.
BarSizeSBE = 5.0
SBE_S_strirrup = 3.0

Section2A = PlanarWallSection(11.25*12., 24., LLBE, LSBE, BarSizeLBE,
                              FindConfinementRebarSpacing(NoOfBarsLBE),
                              FindConfinementRebarSpacing(NoOfBarsSBE),
                              0.255, 4., fpc, fy, fu, LBE_S_strirrup, 4., len(FindConfinementRebarSpacing(NoOfBarsLBE)), 3)

Section2B = PlanarWallSection(11.25*12., 24., LSBE, LLBE, BarSizeLBE,
                              FindConfinementRebarSpacing(NoOfBarsSBE),
                              FindConfinementRebarSpacing(NoOfBarsLBE),
                              0.255, 4., fpc, fy, fu, LBE_S_strirrup, 4., len(FindConfinementRebarSpacing(NoOfBarsSBE)), 3)

LLBE = 30
NoOfBarsLBE = 14
BarSizeLBE = 10.173
LBE_S_strirrup = 2.5

LSBE = 24
NoOfBarsSBE = 12.
BarSizeSBE = 5.0
SBE_S_strirrup = 2.5

Section3A = PlanarWallSection(11.25*12., 24., LLBE, LSBE, BarSizeLBE,
                              FindConfinementRebarSpacing(NoOfBarsLBE),
                              FindConfinementRebarSpacing(NoOfBarsSBE),
                              0.255, 4., fpc, fy, fu, LBE_S_strirrup, 4., len(FindConfinementRebarSpacing(NoOfBarsLBE)), 3)

Section3B = PlanarWallSection(11.25*12., 24., LSBE, LLBE, BarSizeLBE,
                              FindConfinementRebarSpacing(NoOfBarsSBE),
                              FindConfinementRebarSpacing(NoOfBarsLBE),
                              0.255, 4., fpc, fy, fu, LBE_S_strirrup, 4., len(FindConfinementRebarSpacing(NoOfBarsSBE)), 3)

CB1 = CouplingBeam(24., 5*12., fpc, fy, 4, 3, 10.16, True, 4., 6., 2.5, 4, 3, 10)
CB2 = CouplingBeam(24., 3*12., fpc, fy, 4, 3, 11.28, True, 4., 6., 2.5, 4, 3, 10)
CB3 = CouplingBeam(24., 3*12., fpc, fy, 4, 3, 10.16, True, 4., 6., 2.5, 4, 3, 10)
CB4 = CouplingBeam(24., 3*12., fpc, fy, 4, 2, 10.16, True, 4., 6., 2.5, 4, 3, 10)
CB5 = CouplingBeam(24., 3*12., fpc, fy, 4, 2, 10.16, True, 4., 6., 2.5, 4, 3, 10) #Last Row needs to be fixed

CouplingBeams = [CB1, CB2, CB2, CB2, CB3, CB3, CB4, CB5]
CouplingBeamLengths = [7.5*12., 7.5*12, 7.5*12, 7.5*12, 7.5*12, 7.5*12, 7.5*12, 7.5*12]

ATCWALLS8O8S50ELF = ArchetypeData('ATCWALLS8O8S50ELF', YGrids, None, None, None, None, None, None, None, fpc, fy, fu, PDeltaColumnGravityLoad,
                                Mass, WallGravityLoad,
                                CustomWallSection=[[Section1A, Section1B],
                                                   [Section1A, Section1B],
                                                   [Section1A, Section1B],
                                                   [Section2A, Section2B],
                                                   [Section2A, Section2B],
                                                   [Section2A, Section2B],

                                                   [Section3A, Section3B],
                                                   [Section3A, Section3B],
                                                   ],
                                CouplingBeamLength=CouplingBeamLengths,
                                CouplingBeams=CouplingBeams, heff = 73.57*12.)

#endregion

#region ATCWALLS8O8S50ELF48

# Geometry
NoOfStories = 8
YGrids = [0] + np.array(np.arange(0,(NoOfStories)*13*12, 13*12)+15*12).tolist()

# Loading
# DL 175psf Floors and 140psf Roof
# LL 65psf Floors and 20psf Roof

MassPerSqFt = (199.5)/1000. # use this 200 psf
FloorArea = 60.*120.
Mass = np.ones(NoOfStories) * MassPerSqFt * FloorArea
Mass[-1] = FloorArea*154.9/1000. # Adjust for Roof Weight

WallTribArea = 30.*44.
WeightPerSqFt = 212./1000.

BuildingWeight = np.ones(NoOfStories) * WeightPerSqFt * FloorArea
BuildingWeight[-1] = 152. / 1000. * FloorArea # Adjust for Roof Weight

PDeltaColumnGravityLoad = -1 * BuildingWeight * (FloorArea - WallTribArea) / FloorArea

WallGravityLoad = BuildingWeight * WallTribArea / FloorArea * -1
# WallGravityLoad = Mass * WallTribArea / FloorArea * -1 # This is what Kamal used in Atena

# Material Properties
fy = 60.
fu = 105.
fpc = 5.

def FindConfinementRebarSpacing(NoOfBars):
    RebarLayout = []
    NoOfCols = int(int(NoOfBars-2)/2)
    for i in range(NoOfCols):
        if i == 0 or i == NoOfCols - 1:
            RebarLayout.append(3)
        else:
            RebarLayout.append(2)

    return RebarLayout

LLBE = 72
NoOfBarsLBE = 42.
BarSizeLBE = 11.28
LBE_S_strirrup = 3.5

LSBE = 30
NoOfBarsSBE = 12
BarSizeSBE = 11.28
SBE_S_strirrup = 3.5

Section1A = PlanarWallSection(11.25*12., 48., LLBE, LSBE, BarSizeLBE,
                              np.ones(14)*3,
                              FindConfinementRebarSpacing(NoOfBarsSBE),
                              0.255, 4., fpc, fy, fu, LBE_S_strirrup, 4., 13, 3)

Section1B = PlanarWallSection(11.25*12., 48., LSBE, LLBE, BarSizeLBE,
                              FindConfinementRebarSpacing(NoOfBarsSBE),
                              np.ones(14) * 3,
                              0.255, 4., fpc, fy, fu, LBE_S_strirrup, 4., len(FindConfinementRebarSpacing(NoOfBarsSBE)), 3)

LLBE = 42
NoOfBarsLBE = 22
BarSizeLBE = 11.28
LBE_S_strirrup = 3.5

LSBE = 30
NoOfBarsSBE = 12.
BarSizeSBE = 5.0
SBE_S_strirrup = 3.0

Section2A = PlanarWallSection(11.25*12., 48., LLBE, LSBE, BarSizeLBE,
                              FindConfinementRebarSpacing(NoOfBarsLBE),
                              FindConfinementRebarSpacing(NoOfBarsSBE),
                              0.255, 4., fpc, fy, fu, LBE_S_strirrup, 4., len(FindConfinementRebarSpacing(NoOfBarsLBE)), 3)

Section2B = PlanarWallSection(11.25*12., 48., LSBE, LLBE, BarSizeLBE,
                              FindConfinementRebarSpacing(NoOfBarsSBE),
                              FindConfinementRebarSpacing(NoOfBarsLBE),
                              0.255, 4., fpc, fy, fu, LBE_S_strirrup, 4., len(FindConfinementRebarSpacing(NoOfBarsSBE)), 3)

LLBE = 30
NoOfBarsLBE = 14
BarSizeLBE = 10.173
LBE_S_strirrup = 2.5

LSBE = 24
NoOfBarsSBE = 12.
BarSizeSBE = 5.0
SBE_S_strirrup = 2.5

Section3A = PlanarWallSection(11.25*12., 48., LLBE, LSBE, BarSizeLBE,
                              FindConfinementRebarSpacing(NoOfBarsLBE),
                              FindConfinementRebarSpacing(NoOfBarsSBE),
                              0.255, 4., fpc, fy, fu, LBE_S_strirrup, 4., len(FindConfinementRebarSpacing(NoOfBarsLBE)), 3)

Section3B = PlanarWallSection(11.25*12., 48., LSBE, LLBE, BarSizeLBE,
                              FindConfinementRebarSpacing(NoOfBarsSBE),
                              FindConfinementRebarSpacing(NoOfBarsLBE),
                              0.255, 4., fpc, fy, fu, LBE_S_strirrup, 4., len(FindConfinementRebarSpacing(NoOfBarsSBE)), 3)

CB1 = CouplingBeam(24., 5*12., fpc, fy, 4, 3, 10.16, True, 4., 6., 2.5, 4, 3, 10)
CB2 = CouplingBeam(24., 3*12., fpc, fy, 4, 3, 11.28, True, 4., 6., 2.5, 4, 3, 10)
CB3 = CouplingBeam(24., 3*12., fpc, fy, 4, 3, 10.16, True, 4., 6., 2.5, 4, 3, 10)
CB4 = CouplingBeam(24., 3*12., fpc, fy, 4, 2, 10.16, True, 4., 6., 2.5, 4, 3, 10)
CB5 = CouplingBeam(24., 3*12., fpc, fy, 4, 2, 10.16, True, 4., 6., 2.5, 4, 3, 10) #Last Row needs to be fixed

CouplingBeams = [CB1, CB2, CB2, CB2, CB3, CB3, CB4, CB5]
CouplingBeamLengths = [7.5*12., 7.5*12, 7.5*12, 7.5*12, 7.5*12, 7.5*12, 7.5*12, 7.5*12]

ATCWALLS8O8S50ELF48 = ArchetypeData('ATCWALLS8O8S50ELF48', YGrids, None, None, None, None, None, None, None, fpc, fy, fu, PDeltaColumnGravityLoad,
                                Mass, WallGravityLoad,
                                CustomWallSection=[[Section1A, Section1B],
                                                   [Section1A, Section1B],
                                                   [Section1A, Section1B],
                                                   [Section2A, Section2B],
                                                   [Section2A, Section2B],
                                                   [Section2A, Section2B],

                                                   [Section3A, Section3B],
                                                   [Section3A, Section3B],
                                                   ],
                                CouplingBeamLength=CouplingBeamLengths,
                                CouplingBeams=CouplingBeams, heff = 73.57*12.)

#endregion

# 5th Story Opening

#region ATCWALLS8O1AT5S50ELF

# Geometry
NoOfStories = 8
YGrids = [0] + np.array(np.arange(0,(NoOfStories)*13*12, 13*12)+15*12).tolist()

# Loading
# DL 175psf Floors and 140psf Roof
# LL 65psf Floors and 20psf Roof

MassPerSqFt = (199.5)/1000. # use this 200 psf
FloorArea = 60.*120.
Mass = np.ones(NoOfStories) * MassPerSqFt * FloorArea
Mass[-1] = FloorArea*154.9/1000. # Adjust for Roof Weight

WallTribArea = 30.*44.
WeightPerSqFt = 212./1000.

BuildingWeight = np.ones(NoOfStories) * WeightPerSqFt * FloorArea
BuildingWeight[-1] = 152. / 1000. * FloorArea # Adjust for Roof Weight

PDeltaColumnGravityLoad = -1 * BuildingWeight * (FloorArea - WallTribArea) / FloorArea

WallGravityLoad = BuildingWeight * WallTribArea / FloorArea * -1
# WallGravityLoad = Mass * WallTribArea / FloorArea * -1 # This is what Kamal used in Atena


# Material Properties
fy = 60.
fu = 105.
fpc = 5.

Section1 = PlanarWallSection(30.*12., 24., 72, 72, 11.274,
                             [3, 2, 2, 2, 2,
                              2, 2, 2, 2, 2,
                              2, 2, 2, 3],
                             [3, 2, 2, 2, 2,
                              2, 2, 2, 2, 2,
                              2, 2, 2, 3],
                             0.255, 4., fpc, fy, fu, 3., 4., 14, 3)

Section2 = PlanarWallSection(30.*12., 24., 36, 36, 10.173,
                             [3, 2, 2, 2, 2,
                              2, 2, 3],
                             [3, 2, 2, 2, 2,
                              2, 2, 3],
                             0.255, 4., fpc, fy, fu, 3, 4., 8, 3)

Section2A = PlanarWallSection(11.25*12., 24., 36, 12, 9.027,
                             [3, 2, 2, 2, 2,
                              2, 2, 3],
                             [3, 2, 3],
                             0.255, 4., fpc, fy, fu, 3., 4., 8, 3)

Section2B = PlanarWallSection(11.25*12., 24., 12, 36, 9.027,
                              [3, 2, 3],
                              [3, 2, 2, 2, 2,
                               2, 2, 3],
                             0.255, 4., fpc, fy, fu, 3., 4., 3, 3)

Section3 = PlanarWallSection(30.*12., 24., 0, 0, 10.173,
                             [],
                             [],
                             0.255, 4., fpc, fy, fu, None, None, None, None)

CouplingBeams = [None, None, None, None, None, None, None, None]
CouplingBeamLengths = [None, None, None, None, 7.5*12., None, None, None]

ATCWALLS8O1AT5S50ELF = ArchetypeData('ATCWALLS8O1AT5S50ELF', YGrids, None, None, None, None, None, None, None, fpc, fy, fu, PDeltaColumnGravityLoad,
                                Mass, WallGravityLoad,
                                CustomWallSection=[ Section1, Section1, Section1,
                                                    Section2, [Section2A, Section2B], Section2,
                                                    Section3, Section3],
                                CouplingBeamLength=CouplingBeamLengths,
                                CouplingBeams=CouplingBeams, heff = 73.57*12.)

#endregion

# 12 Story Wall

#region ATCWALLS12O0S100ELF
NoOfStories = 12
YGrids = [0] + np.array(np.arange(0,(NoOfStories)*13*12, 13*12)+15*12).tolist()

# Loading
# DL 175psf Floors and 140psf Roof
# LL 65psf Floors and 20psf Roof

MassPerSqFt = (199.5)/1000. # use this 200 psf
FloorArea = 120.*120./4.
Mass = np.ones(NoOfStories) * MassPerSqFt * FloorArea
Mass[-1] = FloorArea*154.9/1000. # Adjust for Roof Weight

WallTribArea = 30.*44.
WeightPerSqFt = 212./1000.

BuildingWeight = np.ones(NoOfStories) * WeightPerSqFt * FloorArea
BuildingWeight[-1] = 152. / 1000. * FloorArea # Adjust for Roof Weight

PDeltaColumnGravityLoad = -1 * BuildingWeight * (FloorArea - WallTribArea) / FloorArea

WallGravityLoad = BuildingWeight * WallTribArea / FloorArea * -1
# WallGravityLoad = Mass * WallTribArea / FloorArea * -1 # This is what Kamal used in Atena

# Material Properties
fy = 60
fu = 105
fpc = 5

Section1 = PlanarWallSection(30.*12., 24., 72, 72, 10.127,
                             [3, 2, 2, 2, 2,
                              2, 2, 2, 2, 2,
                              2, 2, 3],
                             [3, 2, 2, 2, 2,
                              2, 2, 2, 2, 2,
                              2, 2, 3],
                             0.255, 4., fpc, fy, fu, 3.5, 4., 13, 3)

Section2 = PlanarWallSection(30.*12., 24., 42, 42, 10.127,
                             [3, 2, 2, 2,
                              2, 2, 2, 3],
                             [3, 2, 2, 2,
                              2, 2, 2, 3],
                             0.255, 4., fpc, fy, fu, 3.5, 4., 8, 3)

Section3 = PlanarWallSection(30.*12., 24., 30, 30, 8.,
                             [3, 2, 2, 2, 3],
                             [3, 2, 2, 2, 3],
                             0.255, 4., fpc, fy, fu, 3, 4., 5, 3)

Section4 = PlanarWallSection(30.*12., 24., 0, 0, 10.273,
                             [],
                             [],
                             0.255, 4., fpc, fy, fu, None, None, None, None)

ATCWALLS12O0S100ELF = ArchetypeData('ATCWALLS12O0S100ELF', YGrids, None, None, None, None, None, None, None,
                                    fpc, fy, fu, PDeltaColumnGravityLoad,
                              Mass, WallGravityLoad,
                                    CustomWallSection=[Section1, Section1, Section1,
                                                     Section2, Section2, Section2,
                                                     Section3, Section3, Section4,
                                                     Section4, Section4, Section4],
                                    heff = 108.*12.)

#endregion

#region ATCWALLS12O0S100MRSA
NoOfStories = 12
YGrids = [0] + np.array(np.arange(0,(NoOfStories)*13*12, 13*12)+15*12).tolist()

# Loading
# DL 175psf Floors and 140psf Roof
# LL 65psf Floors and 20psf Roof

MassPerSqFt = (199.5)/1000. # use this 200 psf
FloorArea = 120.*120./4.
Mass = np.ones(NoOfStories) * MassPerSqFt * FloorArea
Mass[-1] = FloorArea*154.9/1000. # Adjust for Roof Weight

WallTribArea = 30.*44.
WeightPerSqFt = 212./1000.

BuildingWeight = np.ones(NoOfStories) * WeightPerSqFt * FloorArea
BuildingWeight[-1] = 152. / 1000. * FloorArea # Adjust for Roof Weight

PDeltaColumnGravityLoad = -1 * BuildingWeight * (FloorArea - WallTribArea) / FloorArea

WallGravityLoad = BuildingWeight * WallTribArea / FloorArea * -1
# WallGravityLoad = Mass * WallTribArea / FloorArea * -1 # This is what Kamal used in Atena

# Material Properties
fy = 60
fu = 105
fpc = 5

Section1 = PlanarWallSection(30.*12., 24., 54, 54, 6.992,
                             [3, 2, 2, 2, 2,
                              2, 2, 2, 2, 3],
                             [3, 2, 2, 2, 2,
                              2, 2, 2, 2, 3],
                             0.255, 4., fpc, fy, fu, 3.5, 4., 10, 3)

Section2 = PlanarWallSection(30.*12., 24., 36, 36, 8.00,
                             [3, 2, 2, 2, 2,
                              3],
                             [3, 2, 2, 2, 2,
                              3],
                             0.255, 4., fpc, fy, fu, 3.5, 4., 6, 3)

Section3 = PlanarWallSection(30.*12., 24., 0, 0, 10.273,
                             [],
                             [],
                             0.255, 4., fpc, fy, fu, None, None, None, None)

ATCWALLS12O0S100MRSA = ArchetypeData('ATCWALLS12O0S100MRSA', YGrids, None, None, None, None, None, None, None,
                                    fpc, fy, fu, PDeltaColumnGravityLoad,
                              Mass, WallGravityLoad,
                                    CustomWallSection=[Section1, Section1, Section1,
                                                     Section2, Section2, Section2,
                                                     Section3, Section3, Section3,
                                                     Section3, Section3, Section3],
                                    heff = 108.*12.)

#endregion

#region ATCWALLS12O1S50ELF

# Geometry
NoOfStories = 12
YGrids = [0] + np.array(np.arange(0,(NoOfStories)*13*12, 13*12)+15*12).tolist()

# Loading
# DL 175psf Floors and 140psf Roof
# LL 65psf Floors and 20psf Roof

MassPerSqFt = (199.5)/1000. # use this 200 psf
FloorArea = 120.*120./4.
Mass = np.ones(NoOfStories) * MassPerSqFt * FloorArea
Mass[-1] = FloorArea*154.9/1000. # Adjust for Roof Weight

WallTribArea = 30.*44.
WeightPerSqFt = 212./1000.

BuildingWeight = np.ones(NoOfStories) * WeightPerSqFt * FloorArea
BuildingWeight[-1] = 152. / 1000. * FloorArea # Adjust for Roof Weight

PDeltaColumnGravityLoad = -1 * BuildingWeight * (FloorArea - WallTribArea) / FloorArea

WallGravityLoad = BuildingWeight * WallTribArea / FloorArea * -1
# WallGravityLoad = Mass * WallTribArea / FloorArea * -1 # This is what Kamal used in Atena

# Material Properties
fy = 60.
fu = 105.
fpc = 5.

Section1A = PlanarWallSection(11.25*12., 24., 72., 30., 10.173,
                             [3, 2, 2, 2, 2,
                              2, 2, 2, 2, 2,
                              2, 2, 2, 2, 2,
                              3],
                             [3, 2, 2, 2, 2,
                              3],
                             0.255, 4., fpc, fy, fu, 3.5, 4., 16, 3)

Section1B = PlanarWallSection(11.25*12., 24., 30, 72., 10.173,
                             [3, 2, 2, 2, 2,
                              3],
                             [3, 2, 2, 2, 2,
                              2, 2, 2, 2, 2,
                              2, 2, 2, 2, 2,
                              3],
                             0.255, 4., fpc, fy, fu, 3.0, 4., 6, 3)

Section1 = PlanarWallSection(30.*12., 24., 72, 72, 10.173,
                             [3, 2, 2, 2, 2,
                              2, 2, 2, 2, 2,
                              2, 2, 2, 2, 2,
                              3],
                             [3, 2, 2, 2, 2,
                              2, 2, 2, 2, 2,
                              2, 2, 2, 2, 2,
                              3],
                             0.255, 4., fpc, fy, fu, 3, 4., 16, 3)

Section2 = PlanarWallSection(30.*12., 24., 42, 42, 10.173,
                             [3, 2, 2, 2,
                              2, 2, 2, 3],
                             [3, 2, 2, 2,
                              2, 2, 2, 3],
                             0.255, 4., fpc, fy, fu, 3.5, 4., 8, 3)

Section3 = PlanarWallSection(30.*12., 24., 30, 30, 8.023,
                             [3, 2, 2, 2, 3],
                             [3, 2, 2, 2, 3],
                             0.255, 4., fpc, fy, fu, 3.0, 4., 5, 3)

Section4 = PlanarWallSection(30.*12., 24., 0, 0, 10.173,
                             [],
                             [],
                             0.255, 4., fpc, fy, fu, None, None, None, None)

CouplingBeams = [None, None, None, None, None, None, None, None]
CouplingBeamLengths = [7.5*12., None, None, None, None, None, None, None]

ATCWALLS12O1S50ELF = ArchetypeData('ATCWALLS12O1S50ELF', YGrids, None, None, None, None, None, None, None, fpc, fy, fu, PDeltaColumnGravityLoad,
                                Mass, WallGravityLoad,
                                CustomWallSection=[[Section1A, Section1B], Section1, Section1,
                                                     Section2, Section2, Section2,
                                                     Section3, Section3, Section4,
                                                     Section4, Section4, Section4],
                                CouplingBeamLength=CouplingBeamLengths,
                                CouplingBeams=CouplingBeams, heff = 108.*12.)

#endregion

#region ATCWALLS12O1S50MRSA

# Geometry
NoOfStories = 12
YGrids = [0] + np.array(np.arange(0,(NoOfStories)*13*12, 13*12)+15*12).tolist()

# Loading
# DL 175psf Floors and 140psf Roof
# LL 65psf Floors and 20psf Roof

MassPerSqFt = (199.5)/1000. # use this 200 psf
FloorArea = 120.*120./4.
Mass = np.ones(NoOfStories) * MassPerSqFt * FloorArea
Mass[-1] = FloorArea*154.9/1000. # Adjust for Roof Weight

WallTribArea = 30.*44.
WeightPerSqFt = 212./1000.

BuildingWeight = np.ones(NoOfStories) * WeightPerSqFt * FloorArea
BuildingWeight[-1] = 152. / 1000. * FloorArea # Adjust for Roof Weight

PDeltaColumnGravityLoad = -1 * BuildingWeight * (FloorArea - WallTribArea) / FloorArea

WallGravityLoad = BuildingWeight * WallTribArea / FloorArea * -1
# WallGravityLoad = Mass * WallTribArea / FloorArea * -1 # This is what Kamal used in Atena

# Material Properties
fy = 60.
fu = 105.
fpc = 5.

Section1A = PlanarWallSection(11.25*12., 24., 54., 30., 10.173,
                             [3, 2, 2, 2, 2,
                              2, 2, 2, 2, 2,
                              2, 3],
                             [3, 2, 2, 2, 3],
                             0.255, 4., fpc, fy, fu, 3.5, 4., 12, 3)

Section1B = PlanarWallSection(11.25*12., 24., 30, 54., 10.173,
                             [3, 2, 2, 2, 3],
                             [3, 2, 2, 2, 2,
                              2, 2, 2, 2, 2,
                              2, 3],
                             0.255, 4., fpc, fy, fu, 3.0, 4., 5, 3)

Section1 = PlanarWallSection(30.*12., 24., 54., 54., 10.173,
                             [3, 2, 2, 2, 2,
                              2, 2, 2, 2, 2,
                              2, 3],
                             [3, 2, 2, 2, 2,
                              2, 2, 2, 2, 2,
                              2, 3],
                             0.255, 4., fpc, fy, fu, 3.5, 4., 12, 3)

Section2 = PlanarWallSection(30.*12., 24., 36., 36., 8.0,
                             [3, 2, 2, 2, 2,
                              3],
                             [3, 2, 2, 2, 2,
                              3],
                             0.255, 4., fpc, fy, fu, 3.0, 4., 6, 3)


Section3 = PlanarWallSection(30.*12., 24., 0, 0, 10.173,
                             [],
                             [],
                             0.255, 4., fpc, fy, fu, None, None, None, None)

CouplingBeams = [None, None, None, None, None, None, None, None]
CouplingBeamLengths = [7.5*12., None, None, None, None, None, None, None]

ATCWALLS12O1S50MRSA = ArchetypeData('ATCWALLS12O1S50MRSA', YGrids, None, None, None, None, None, None, None, fpc, fy, fu, PDeltaColumnGravityLoad,
                                Mass, WallGravityLoad,
                                CustomWallSection=[[Section1A, Section1B], Section1, Section1,
                                                     Section2, Section2, Section2,
                                                     Section3, Section3, Section3,
                                                   Section3, Section3, Section3],
                                CouplingBeamLength=CouplingBeamLengths,
                                CouplingBeams=CouplingBeams, heff = 108.*12.)

#endregion

#region ATCWALLS12O1S75ELF

# Geometry
NoOfStories = 12
YGrids = [0] + np.array(np.arange(0,(NoOfStories)*13*12, 13*12)+15*12).tolist()

# Loading
# DL 175psf Floors and 140psf Roof
# LL 65psf Floors and 20psf Roof

MassPerSqFt = (199.5)/1000. # use this 200 psf
FloorArea = 120.*120./4.
Mass = np.ones(NoOfStories) * MassPerSqFt * FloorArea
Mass[-1] = FloorArea*154.9/1000. # Adjust for Roof Weight

WallTribArea = 30.*44.
WeightPerSqFt = 212./1000.

BuildingWeight = np.ones(NoOfStories) * WeightPerSqFt * FloorArea
BuildingWeight[-1] = 152. / 1000. * FloorArea # Adjust for Roof Weight

PDeltaColumnGravityLoad = -1 * BuildingWeight * (FloorArea - WallTribArea) / FloorArea

WallGravityLoad = BuildingWeight * WallTribArea / FloorArea * -1
# WallGravityLoad = Mass * WallTribArea / FloorArea * -1 # This is what Kamal used in Atena

# Material Properties
fy = 60.
fu = 105.
fpc = 5.

Section1A = PlanarWallSection(13.75*12., 24., 72., 30., 10.173,
                             [3, 2, 2, 2, 2,
                              2, 2, 2, 2, 2,
                              2, 2, 2, 2, 3],
                             [3, 2, 2, 2, 2,
                              3],
                             0.255, 4., fpc, fy, fu, 3.5, 4., 15, 3)

Section1B = PlanarWallSection(13.75*12., 24., 30, 72., 10.173,
                             [3, 2, 2, 2, 2,
                              3],
                             [3, 2, 2, 2, 2,
                              2, 2, 2, 2, 2,
                              2, 2, 2, 2, 3],
                             0.255, 4., fpc, fy, fu, 3.0, 4., 6, 3)

Section1 = PlanarWallSection(30.*12., 24., 72, 72, 10.173,
                             [3, 2, 2, 2, 2,
                              2, 2, 2, 2, 2,
                              2, 2, 2, 2, 3],
                             [3, 2, 2, 2, 2,
                              2, 2, 2, 2, 2,
                              2, 2, 2, 2, 3],
                             0.255, 4., fpc, fy, fu, 3.5, 4., 15, 3)

Section2 = PlanarWallSection(30.*12., 24., 42, 42, 10.173,
                             [3, 2, 2, 2,
                              2, 2, 2, 3],
                             [3, 2, 2, 2,
                              2, 2, 2, 3],
                             0.255, 4., fpc, fy, fu, 3.5, 4., 8, 3)

Section3 = PlanarWallSection(30.*12., 24., 30, 30, 8.023,
                             [3, 2, 2, 2, 3],
                             [3, 2, 2, 2, 3],
                             0.255, 4., fpc, fy, fu, 3.0, 4., 5, 3)

Section4 = PlanarWallSection(30.*12., 24., 0, 0, 10.173,
                             [],
                             [],
                             0.255, 4., fpc, fy, fu, None, None, None, None)

CouplingBeams = [None, None, None, None, None, None, None, None]
CouplingBeamLengths = [2.5*12., None, None, None, None, None, None, None]

ATCWALLS12O1S75ELF = ArchetypeData('ATCWALLS12O1S75ELF', YGrids, None, None, None, None, None, None, None, fpc, fy, fu, PDeltaColumnGravityLoad,
                                Mass, WallGravityLoad,
                                CustomWallSection=[[Section1A, Section1B], Section1, Section1,
                                                     Section2, Section2, Section2,
                                                     Section3, Section3, Section4,
                                                     Section4, Section4, Section4],
                                CouplingBeamLength=CouplingBeamLengths,
                                CouplingBeams=CouplingBeams, heff = 108.*12.)

#endregion

#region ATCWALLS12O1S75MRSA

# Geometry
NoOfStories = 12
YGrids = [0] + np.array(np.arange(0,(NoOfStories)*13*12, 13*12)+15*12).tolist()

# Loading
# DL 175psf Floors and 140psf Roof
# LL 65psf Floors and 20psf Roof

MassPerSqFt = (199.5)/1000. # use this 200 psf
FloorArea = 120.*120./4.
Mass = np.ones(NoOfStories) * MassPerSqFt * FloorArea
Mass[-1] = FloorArea*154.9/1000. # Adjust for Roof Weight

WallTribArea = 30.*44.
WeightPerSqFt = 212./1000.

BuildingWeight = np.ones(NoOfStories) * WeightPerSqFt * FloorArea
BuildingWeight[-1] = 152. / 1000. * FloorArea # Adjust for Roof Weight

PDeltaColumnGravityLoad = -1 * BuildingWeight * (FloorArea - WallTribArea) / FloorArea

WallGravityLoad = BuildingWeight * WallTribArea / FloorArea * -1
# WallGravityLoad = Mass * WallTribArea / FloorArea * -1 # This is what Kamal used in Atena

# Material Properties
fy = 60.
fu = 105.
fpc = 5.

Section1A = PlanarWallSection(13.75*12., 24., 54., 30., 10.173,
                             [3, 2, 2, 2, 2,
                              2, 2, 2, 2, 2,
                              3],
                             [3, 2, 2, 2, 3],
                             0.255, 4., fpc, fy, fu, 3.5, 4., 15, 3)

Section1B = PlanarWallSection(13.75*12., 24., 30, 54., 10.173,
                             [3, 2, 2, 2, 3],
                             [3, 2, 2, 2, 2,
                              2, 2, 2, 2, 2,
                              3],
                             0.255, 4., fpc, fy, fu, 3.0, 4., 5, 3)

Section1 = PlanarWallSection(30.*12., 24., 54., 54., 10.173,
                             [3, 2, 2, 2, 2,
                              2, 2, 2, 2, 2,
                              3],
                             [3, 2, 2, 2, 2,
                              2, 2, 2, 2, 2,
                              3],
                             0.255, 4., fpc, fy, fu, 3.5, 4., 15, 3)

Section2 = PlanarWallSection(30.*12., 24., 36., 36., 8.00,
                             [3, 2, 2, 2,
                              2, 2, 2, 3],
                             [3, 2, 2, 2,
                              2, 2, 2, 3],
                             0.255, 4., fpc, fy, fu, 3.5, 4., 8, 3)

Section3 = PlanarWallSection(30.*12., 24., 0, 0, 10.173,
                             [],
                             [],
                             0.255, 4., fpc, fy, fu, None, None, None, None)

CouplingBeams = [None, None, None, None, None, None, None, None]
CouplingBeamLengths = [2.5*12., None, None, None, None, None, None, None]

ATCWALLS12O1S75MRSA = ArchetypeData('ATCWALLS12O1S75MRSA', YGrids, None, None, None, None, None, None, None, fpc, fy, fu, PDeltaColumnGravityLoad,
                                Mass, WallGravityLoad,
                                CustomWallSection=[[Section1A, Section1B], Section1, Section1,
                                                   Section2, Section2, Section2,
                                                   Section3, Section3, Section3,
                                                   Section3, Section3, Section3],
                                CouplingBeamLength=CouplingBeamLengths,
                                CouplingBeams=CouplingBeams, heff = 108.*12.)

#endregion

#region ATCWALLS12O2S50ELF

# Geometry
NoOfStories = 12
YGrids = [0] + np.array(np.arange(0,(NoOfStories)*13*12, 13*12)+15*12).tolist()

# Loading
# DL 175psf Floors and 140psf Roof
# LL 65psf Floors and 20psf Roof

MassPerSqFt = (199.5)/1000. # use this 200 psf
FloorArea = 120.*120./4.
Mass = np.ones(NoOfStories) * MassPerSqFt * FloorArea
Mass[-1] = FloorArea*154.9/1000. # Adjust for Roof Weight

WallTribArea = 30.*44.
WeightPerSqFt = 212./1000.

BuildingWeight = np.ones(NoOfStories) * WeightPerSqFt * FloorArea
BuildingWeight[-1] = 152. / 1000. * FloorArea # Adjust for Roof Weight

PDeltaColumnGravityLoad = -1 * BuildingWeight * (FloorArea - WallTribArea) / FloorArea

WallGravityLoad = BuildingWeight * WallTribArea / FloorArea * -1
# WallGravityLoad = Mass * WallTribArea / FloorArea * -1 # This is what Kamal used in Atena

# Material Properties
fy = 60.
fu = 105.
fpc = 5.

Section1A = PlanarWallSection(11.25*12., 24., 72, 30, 10.173,
                             [3, 2, 2, 2, 2,
                              2, 2, 2, 2, 2,
                              2, 2, 2, 2, 2,
                              3],
                             [3, 2, 2, 2, 2,
                              3],
                             0.255, 4., fpc, fy, fu, 3.5, 4., 16, 3)

Section1B = PlanarWallSection(11.25*12., 24., 30, 72., 10.173,
                             [3, 2, 2, 2, 2,
                              3],
                             [3, 2, 2, 2, 2,
                              2, 2, 2, 2, 2,
                              2, 2, 2, 2, 2,
                              3],
                             0.255, 4., fpc, fy, fu, 3., 4., 6, 3)

Section1 = PlanarWallSection(30.*12., 24., 72., 72., 10.173,
                             [3, 2, 2, 2, 2,
                              2, 2, 2, 2, 2,
                              2, 2, 2, 2, 2,
                              3],
                             [3, 2, 2, 2, 2,
                              2, 2, 2, 2, 2,
                              2, 2, 2, 2, 2,
                              3],
                             0.255, 4., fpc, fy, fu, 3.0, 4., 16, 3)

Section2 = PlanarWallSection(30.*12., 24., 42, 42, 10.173,
                             [3, 2, 2, 2, 2,
                              2, 2, 3],
                             [3, 2, 2, 2, 2,
                              2, 2, 3],
                             0.255, 4., fpc, fy, fu, 3.5, 4., 8, 3)

Section3 = PlanarWallSection(30.*12., 24., 30., 30, 8.023,
                             [3, 2, 2, 2, 3],
                             [3, 2, 2, 2, 3],
                             0.255, 4., fpc, fy, fu, 3., 4., 5, 3)

Section4 = PlanarWallSection(30.*12., 24., 0, 0, 10.173,
                             [],
                             [],
                             0.255, 4., fpc, fy, fu, None, None, None, None)

CB1 = CouplingBeam(24., 60., fpc, fy, 4, 2, 9.027, True, 4., 6., 2.5, 4, 3, 10)

CouplingBeams = [CB1, None, None, None, None, None, None, None]
CouplingBeamLengths = [7.5*12., 7.5*12, None, None, None, None, None, None]

ATCWALLS12O2S50ELF = ArchetypeData('ATCWALLS12O2S50ELF', YGrids, None, None, None, None, None, None, None, fpc, fy, fu, PDeltaColumnGravityLoad,
                                Mass, WallGravityLoad,
                                CustomWallSection=[[Section1A, Section1B], [Section1A, Section1B], Section1,
                                                   Section1, Section2, Section2,
                                                   Section3, Section3, Section4,
                                                   Section4, Section4, Section4],
                                CouplingBeamLength=CouplingBeamLengths,
                                CouplingBeams=CouplingBeams, heff = 108.*12.)

#endregion

#region ATCWALLS12O2S50MRSA

# Geometry
NoOfStories = 12
YGrids = [0] + np.array(np.arange(0,(NoOfStories)*13*12, 13*12)+15*12).tolist()

# Loading
# DL 175psf Floors and 140psf Roof
# LL 65psf Floors and 20psf Roof

MassPerSqFt = (199.5)/1000. # use this 200 psf
FloorArea = 120.*120./4.
Mass = np.ones(NoOfStories) * MassPerSqFt * FloorArea
Mass[-1] = FloorArea*154.9/1000. # Adjust for Roof Weight

WallTribArea = 30.*44.
WeightPerSqFt = 212./1000.

BuildingWeight = np.ones(NoOfStories) * WeightPerSqFt * FloorArea
BuildingWeight[-1] = 152. / 1000. * FloorArea # Adjust for Roof Weight

PDeltaColumnGravityLoad = -1 * BuildingWeight * (FloorArea - WallTribArea) / FloorArea

WallGravityLoad = BuildingWeight * WallTribArea / FloorArea * -1
# WallGravityLoad = Mass * WallTribArea / FloorArea * -1 # This is what Kamal used in Atena

# Material Properties
fy = 60.
fu = 105.
fpc = 5.

Section1A = PlanarWallSection(11.25*12., 24., 54., 30., 9.027,
                             [3, 2, 2, 2, 2,
                              2, 2, 2, 2, 3],
                             [3, 2, 2, 2, 3],
                             0.255, 4., fpc, fy, fu, 3.5, 4., 10, 3)

Section1B = PlanarWallSection(11.25*12., 24., 30, 54., 9.027,
                             [3, 2, 2, 2, 3],
                             [3, 2, 2, 2, 2,
                              2, 2, 2, 2, 3],
                             0.255, 4., fpc, fy, fu, 3., 4., 5, 3)

Section1 = PlanarWallSection(30.*12., 24., 54., 54., 9.027,
                             [3, 2, 2, 2, 2,
                              2, 2, 2, 2, 3],
                             [3, 2, 2, 2, 2,
                              2, 2, 2, 2, 3],
                             0.255, 4., fpc, fy, fu, 3.0, 4., 10, 3)

Section2 = PlanarWallSection(30.*12., 24., 36, 36, 8.0,
                             [3, 2, 2, 2, 2,
                              3],
                             [3, 2, 2, 2, 2,
                              3],
                             0.255, 4., fpc, fy, fu, 3.0, 4., 8, 3)

Section3 = PlanarWallSection(30.*12., 24., 0, 0, 10.173,
                             [],
                             [],
                             0.255, 4., fpc, fy, fu, None, None, None, None)

CB1 = CouplingBeam(24., 60., fpc, fy, 4, 2, 9.027, True, 4., 6., 2.5, 4, 3, 10)

CouplingBeams = [CB1, None, None, None, None, None, None, None]
CouplingBeamLengths = [7.5*12., 7.5*12, None, None, None, None, None, None]

ATCWALLS12O2S50MRSA = ArchetypeData('ATCWALLS12O2S50MRSA', YGrids, None, None, None, None, None, None, None, fpc, fy, fu, PDeltaColumnGravityLoad,
                                Mass, WallGravityLoad,
                                CustomWallSection=[[Section1A, Section1B], [Section1A, Section1B], Section1,
                                                   Section1, Section2, Section2,
                                                   Section3, Section3, Section3,
                                                   Section3, Section3, Section3],
                                CouplingBeamLength=CouplingBeamLengths,
                                CouplingBeams=CouplingBeams, heff = 108.*12.)

#endregion

#region ATCWALLS12O2S75ELF

# Geometry
NoOfStories = 12
YGrids = [0] + np.array(np.arange(0,(NoOfStories)*13*12, 13*12)+15*12).tolist()

# Loading
# DL 175psf Floors and 140psf Roof
# LL 65psf Floors and 20psf Roof

MassPerSqFt = (199.5)/1000. # use this 200 psf
FloorArea = 120.*120./4.
Mass = np.ones(NoOfStories) * MassPerSqFt * FloorArea
Mass[-1] = FloorArea*154.9/1000. # Adjust for Roof Weight

WallTribArea = 30.*44.
WeightPerSqFt = 212./1000.

BuildingWeight = np.ones(NoOfStories) * WeightPerSqFt * FloorArea
BuildingWeight[-1] = 152. / 1000. * FloorArea # Adjust for Roof Weight

PDeltaColumnGravityLoad = -1 * BuildingWeight * (FloorArea - WallTribArea) / FloorArea

WallGravityLoad = BuildingWeight * WallTribArea / FloorArea * -1
# WallGravityLoad = Mass * WallTribArea / FloorArea * -1 # This is what Kamal used in Atena

# Material Properties
fy = 60.
fu = 105.
fpc = 5.

Section1A = PlanarWallSection(13.75*12., 24., 72, 30, 11.,
                             [3, 2, 2, 2, 2,
                              2, 2, 2, 2, 2,
                              2, 2, 2, 2, 3],
                             [3, 2, 2, 2, 2,
                              3],
                             0.255, 4., fpc, fy, fu, 3.5, 4., 15, 3)

Section1B = PlanarWallSection(13.75*12., 24., 30, 72., 11.,
                             [3, 2, 2, 2, 2,
                              3],
                              [3, 2, 2, 2, 2,
                               2, 2, 2, 2, 2,
                               2, 2, 2, 2, 3],
                             0.255, 4., fpc, fy, fu, 3., 4., 6, 3)

Section1 = PlanarWallSection(30.*12., 24., 72., 72., 11.,
                             [3, 2, 2, 2, 2,
                              2, 2, 2, 2, 2,
                              2, 2, 2, 2, 3],
                             [3, 2, 2, 2, 2,
                              2, 2, 2, 2, 2,
                              2, 2, 2, 2, 3],
                             0.255, 4., fpc, fy, fu, 3.5, 4., 15, 3)

Section2 = PlanarWallSection(30.*12., 24., 42, 42, 10.173,
                             [3, 2, 2, 2, 2,
                              2, 2, 3],
                             [3, 2, 2, 2, 2,
                              2, 2, 3],
                             0.255, 4., fpc, fy, fu, 3.0, 4., 8, 3)

Section3 = PlanarWallSection(30.*12., 24., 30, 30, 8.023,
                             [3, 2, 2, 2, 3],
                             [3, 2, 2, 2, 3],
                             0.255, 4., fpc, fy, fu, 3.0, 4., 6, 3)

Section4 = PlanarWallSection(30.*12., 24., 0, 0, 10.173,
                             [],
                             [],
                             0.255, 4., fpc, fy, fu, None, None, None, None)

CB1 = CouplingBeam(24., 60., fpc, fy, 4, 2, 8.023, True, 4., 6., 2.5, 4, 3, 10)

CouplingBeams = [CB1, None, None, None, None, None, None, None]
CouplingBeamLengths = [2.5*12., 2.5*12, None, None, None, None, None, None]

ATCWALLS12O2S75ELF = ArchetypeData('ATCWALLS12O2S75ELF', YGrids, None, None, None, None, None, None, None, fpc, fy, fu, PDeltaColumnGravityLoad,
                                Mass, WallGravityLoad,
                                CustomWallSection=[[Section1A, Section1B], [Section1A, Section1B], Section1,
                                                   Section1, Section2, Section2,
                                                   Section3, Section3, Section4,
                                                   Section4, Section4, Section4],
                                CouplingBeamLength=CouplingBeamLengths,
                                CouplingBeams=CouplingBeams, heff = 108.*12.)

#endregion

#region ATCWALLS12O2S75MRSA

# Geometry
NoOfStories = 12
YGrids = [0] + np.array(np.arange(0,(NoOfStories)*13*12, 13*12)+15*12).tolist()

# Loading
# DL 175psf Floors and 140psf Roof
# LL 65psf Floors and 20psf Roof

MassPerSqFt = (199.5)/1000. # use this 200 psf
FloorArea = 120.*120./4.
Mass = np.ones(NoOfStories) * MassPerSqFt * FloorArea
Mass[-1] = FloorArea*154.9/1000. # Adjust for Roof Weight

WallTribArea = 30.*44.
WeightPerSqFt = 212./1000.

BuildingWeight = np.ones(NoOfStories) * WeightPerSqFt * FloorArea
BuildingWeight[-1] = 152. / 1000. * FloorArea # Adjust for Roof Weight

PDeltaColumnGravityLoad = -1 * BuildingWeight * (FloorArea - WallTribArea) / FloorArea

WallGravityLoad = BuildingWeight * WallTribArea / FloorArea * -1
# WallGravityLoad = Mass * WallTribArea / FloorArea * -1 # This is what Kamal used in Atena

# Material Properties
fy = 60.
fu = 105.
fpc = 5.

Section1A = PlanarWallSection(13.75*12., 24., 54, 30, 9.027,
                             [3, 2, 2, 2, 2,
                              2, 2, 2, 2, 3],
                             [3, 2, 2, 2, 3],
                              0.255, 4., fpc, fy, fu, 3.5, 4., 10, 3)

Section1B = PlanarWallSection(13.75*12., 24., 30, 54., 9.027,
                              [3, 2, 2, 2, 3],
                              [3, 2, 2, 2, 2,
                               2, 2, 2, 2, 3],
                               0.255, 4., fpc, fy, fu, 3., 4., 6, 3)

Section1 = PlanarWallSection(30.*12., 24., 54., 54., 9.027,
                             [3, 2, 2, 2, 2,
                              2, 2, 2, 2, 3],
                             [3, 2, 2, 2, 2,
                              2, 2, 2, 2, 3],
                             0.255, 4., fpc, fy, fu, 3.5, 4., 15, 3)

Section2 = PlanarWallSection(30.*12., 24., 36., 36., 8.0,
                             [3, 2, 2, 2, 2,
                              3],
                             [3, 2, 2, 2, 2,
                              3],
                             0.255, 4., fpc, fy, fu, 3.0, 4., 8, 3)

Section3 = PlanarWallSection(30.*12., 24., 0, 0, 10.173,
                             [],
                             [],
                             0.255, 4., fpc, fy, fu, None, None, None, None)

CB1 = CouplingBeam(24., 60., fpc, fy, 4, 2, 8.023, True, 4., 6., 3.0, 4, 3, 10)

CouplingBeams = [CB1, None, None, None, None, None, None, None]
CouplingBeamLengths = [2.5*12., 2.5*12, None, None, None, None, None, None]

ATCWALLS12O2S75MRSA = ArchetypeData('ATCWALLS12O2S75MRSA', YGrids, None, None, None, None, None, None, None, fpc, fy, fu, PDeltaColumnGravityLoad,
                                Mass, WallGravityLoad,
                                CustomWallSection=[[Section1A, Section1B], [Section1A, Section1B], Section1,
                                                   Section1, Section2, Section2,
                                                   Section3, Section3, Section3,
                                                   Section3, Section3, Section3],
                                CouplingBeamLength=CouplingBeamLengths,
                                CouplingBeams=CouplingBeams, heff = 108.*12.)

#endregion

#region ATCWALLS12O12S50ELF

# Geometry
NoOfStories = 12
YGrids = [0] + np.array(np.arange(0,(NoOfStories)*13*12, 13*12)+15*12).tolist()

# Loading
# DL 175psf Floors and 140psf Roof
# LL 65psf Floors and 20psf Roof

MassPerSqFt = (199.5)/1000. # use this 200 psf
FloorArea = 120.*120./4.
Mass = np.ones(NoOfStories) * MassPerSqFt * FloorArea
Mass[-1] = FloorArea*154.9/1000. # Adjust for Roof Weight

WallTribArea = 30.*44.
WeightPerSqFt = 212./1000.

BuildingWeight = np.ones(NoOfStories) * WeightPerSqFt * FloorArea
BuildingWeight[-1] = 152. / 1000. * FloorArea # Adjust for Roof Weight

PDeltaColumnGravityLoad = -1 * BuildingWeight * (FloorArea - WallTribArea) / FloorArea

WallGravityLoad = BuildingWeight * WallTribArea / FloorArea * -1
# WallGravityLoad = Mass * WallTribArea / FloorArea * -1 # This is what Kamal used in Atena

# Material Properties
fy = 60.
fu = 105.
fpc = 5.


def FindConfinementRebarSpacing(NoOfBars):
    RebarLayout = []
    NoOfCols = int(int(NoOfBars-2)/2)
    for i in range(NoOfCols):
        if i == 0 or i == NoOfCols - 1:
            RebarLayout.append(3)
        else:
            RebarLayout.append(2)

    return RebarLayout

LLBE = 72.
NoOfBarsLBE = 42.
BarSizeLBE = 10.173
LBE_S_strirrup = 3.5

LSBE = 30
NoOfBarsSBE = 12.
BarSizeSBE = 7.0
SBE_S_strirrup = 3.5

Section1A = PlanarWallSection(11.25*12., 24., LLBE, LSBE, BarSizeLBE,
                              np.ones(14)*3.,#FindConfinementRebarSpacing(NoOfBarsLBE),
                              FindConfinementRebarSpacing(NoOfBarsSBE),
                              0.255, 4., fpc, fy, fu, LBE_S_strirrup, 4., len(FindConfinementRebarSpacing(NoOfBarsLBE)), 3)

Section1B = PlanarWallSection(11.25*12., 24., LSBE, LLBE, BarSizeLBE,
                              FindConfinementRebarSpacing(NoOfBarsSBE),
                              np.ones(14) * 3.,#FindConfinementRebarSpacing(NoOfBarsLBE),
                              0.255, 4., fpc, fy, fu, LBE_S_strirrup, 4., len(FindConfinementRebarSpacing(NoOfBarsSBE)), 3)

LLBE = 54
NoOfBarsLBE = 39.
BarSizeLBE = 9.0
LBE_S_strirrup = 3.5

LSBE = 30
NoOfBarsSBE = 12.
BarSizeSBE = 7.0 #Should be #9 bar but i converted it to a #10
SBE_S_strirrup = 3.0

Section2A = PlanarWallSection(11.25*12., 24., LLBE, LSBE, BarSizeLBE,
                              np.ones(13) * 3.,#FindConfinementRebarSpacing(NoOfBarsLBE),
                              FindConfinementRebarSpacing(NoOfBarsSBE),
                              0.255, 4., fpc, fy, fu, LBE_S_strirrup, 4., len(FindConfinementRebarSpacing(NoOfBarsLBE)), 3)

Section2B = PlanarWallSection(11.25*12., 24., LSBE, LLBE, BarSizeLBE,
                              FindConfinementRebarSpacing(NoOfBarsSBE),
                              np.ones(13) * 3.,  #FindConfinementRebarSpacing(NoOfBarsLBE),
                              0.255, 4., fpc, fy, fu, LBE_S_strirrup, 4., len(FindConfinementRebarSpacing(NoOfBarsSBE)), 3)

LLBE = 36
NoOfBarsLBE = 24.
BarSizeLBE = 9.
LBE_S_strirrup = 3.0

LSBE = 30.
NoOfBarsSBE = 12.
BarSizeSBE = 7.0
SBE_S_strirrup = 2.5

Section3A = PlanarWallSection(11.25*12., 24., LLBE, LSBE, BarSizeLBE,
                              FindConfinementRebarSpacing(NoOfBarsLBE),
                              FindConfinementRebarSpacing(NoOfBarsSBE),
                              0.255, 4., fpc, fy, fu, LBE_S_strirrup, 4., len(FindConfinementRebarSpacing(NoOfBarsLBE)), 3)

Section3B = PlanarWallSection(11.25*12., 24., LSBE, LLBE, BarSizeLBE,
                              FindConfinementRebarSpacing(NoOfBarsSBE),
                              FindConfinementRebarSpacing(NoOfBarsLBE),
                              0.255, 4., fpc, fy, fu, LBE_S_strirrup, 4., len(FindConfinementRebarSpacing(NoOfBarsSBE)), 3)

LLBE = 16.
NoOfBarsLBE = 9.
BarSizeLBE = 4.0
LBE_S_strirrup = 2.0

LSBE = 30.
NoOfBarsSBE = 12.
BarSizeSBE = 7.0
SBE_S_strirrup = 2.5

Section4A = PlanarWallSection(11.25*12., 24., LLBE, LSBE, BarSizeLBE,
                              FindConfinementRebarSpacing(NoOfBarsLBE),
                              FindConfinementRebarSpacing(NoOfBarsSBE),
                              0.255, 4., fpc, fy, fu, LBE_S_strirrup, 4., len(FindConfinementRebarSpacing(NoOfBarsLBE)), 3)

Section4B = PlanarWallSection(11.25*12., 24., LSBE, LLBE, BarSizeLBE,
                              FindConfinementRebarSpacing(NoOfBarsSBE),
                              FindConfinementRebarSpacing(NoOfBarsLBE),
                              0.255, 4., fpc, fy, fu, LBE_S_strirrup, 4., len(FindConfinementRebarSpacing(NoOfBarsSBE)), 3)

CB1 = CouplingBeam(24., 5*12., fpc, fy, 4, 2, 10.16, True, 4., 6., 2.5, 4, 3, 10)
CB2 = CouplingBeam(24., 3*12., fpc, fy, 4, 2, 11.27, True, 4., 6., 2.5, 4, 3, 10)
CB3 = CouplingBeam(24., 3*12., fpc, fy, 4, 2, 10.17, True, 4., 6., 2.5, 4, 3, 10)
CB4 = CouplingBeam(24., 3*12., fpc, fy, 3, 2, 10.17, True, 4., 6., 2.5, 4, 3, 10)

CouplingBeams = [CB1, CB2, CB2, CB2, CB2, CB3, CB3, CB3, CB4, CB4, CB4, CB4]
CouplingBeamLengths = [7.5*12., 7.5*12, 7.5*12, 7.5*12, 7.5*12, 7.5*12, 7.5*12, 7.5*12, 7.5*12, 7.5*12, 7.5*12, 7.5*12]

ATCWALLS12O12S50ELF = ArchetypeData('ATCWALLS12O12S50ELF', YGrids, None, None, None, None, None, None, None, fpc, fy, fu, PDeltaColumnGravityLoad,
                                Mass, WallGravityLoad,
                                CustomWallSection=[[Section1A, Section1B],
                                                   [Section1A, Section1B],
                                                   [Section1A, Section1B],
                                                   [Section2A, Section2B],
                                                   [Section2A, Section2B],
                                                   [Section2A, Section2B],

                                                   [Section3A, Section3B],
                                                   [Section3A, Section3B],
                                                   [Section3A, Section3B],
                                                   [Section4A, Section4B],
                                                   [Section4A, Section4B],
                                                   [Section4A, Section4B],
                                                   ],
                                CouplingBeamLength=CouplingBeamLengths,
                                CouplingBeams=CouplingBeams, heff = 108.*12.)

#endregion

# 20 Story Wall

#region ATCWALLS20O0S100ELF
NoOfStories = 20
YGrids = [0] + np.array(np.arange(0,(NoOfStories)*13*12, 13*12)+15*12).tolist()

# Loading
# DL 175psf Floors and 140psf Roof
# LL 65psf Floors and 20psf Roof

MassPerSqFt = (210.)/1000. # use this 200 psf
FloorArea = 120.*120./2.
Mass = np.ones(NoOfStories) * MassPerSqFt * FloorArea
Mass[-1] = FloorArea*105./1000. # Adjust for Roof Weight

WallTribArea = 30.*15.
WeightPerSqFt = 210./1000.

BuildingWeight = np.ones(NoOfStories) * WeightPerSqFt * FloorArea
BuildingWeight[-1] = 105. / 1000. * FloorArea # Adjust for Roof Weight

PDeltaColumnGravityLoad = -1 * BuildingWeight * (FloorArea - WallTribArea) / FloorArea

WallGravityLoad = BuildingWeight * WallTribArea / FloorArea * -1
# WallGravityLoad = Mass * WallTribArea / FloorArea * -1 # This is what Kamal used in Atena

# Material Properties
fy = 60
fu = 105
fpc = 5

Section1 = IWallSection(402., 96., 18., 1.10, 10.127,
                        fpc, fy, fu,
                        None, None, None, None)

Section2 = IWallSection(402., 72., 18., 1.05, 10.127,
                        fpc, fy, fu,
                        None, None, None, None)

l_w = 402.
t_w = 18.
Rho = 1.22
bar_size = 10.173
NoOfBarRows = np.ceil(l_w/2.*t_w*Rho/100./(np.pi/4.*(bar_size/8.)**2.*2.))
tie_spacing = ATCWallHelper.GetTieSpacing(t_w, bar_size/8., l_w/NoOfBarRows/2.)
Section3 = PlanarWallSection(l_w, t_w, l_w/2., l_w/2., bar_size,
                            np.ones(int(NoOfBarRows))*2., np.ones(int(NoOfBarRows))*2., 0.0, 4.,
                            fpc, fy, fu,
                            tie_spacing, 4.0, 18, 2)

l_w = 360.
t_w = 16.
Rho = 1.04
bar_size = 9.027
NoOfBarRows = np.ceil(l_w/2.*t_w*Rho/100./(np.pi/4.*(bar_size/8.)**2.*2.))
tie_spacing = ATCWallHelper.GetTieSpacing(t_w, bar_size/8., l_w/NoOfBarRows/2.)
Section4 = PlanarWallSection(l_w, t_w, l_w/2., l_w/2., bar_size,
                            np.ones(int(NoOfBarRows))*2., np.ones(int(NoOfBarRows))*2., 0.0, 4.,
                            fpc, fy, fu,
                            tie_spacing, 4.0, 18, 2)

l_w = 360.
t_w = 14.
Rho = 0.95
bar_size = 9.027
NoOfBarRows = np.ceil(l_w/2.*t_w*Rho/100./(np.pi/4.*(bar_size/8.)**2.*2.))
tie_spacing = ATCWallHelper.GetTieSpacing(t_w, bar_size/8., l_w/NoOfBarRows/2.)
Section5 = PlanarWallSection(l_w, t_w, l_w/2., l_w/2., bar_size,
                            np.ones(int(NoOfBarRows))*2., np.ones(int(NoOfBarRows))*2., 0.0, 4.,
                            fpc, fy, fu,
                            tie_spacing, 4.0, 18, 2)

ATCWALLS20O0S100ELF = ArchetypeData('ATCWALLS20O0S100ELF', YGrids, None, None, None, None, None, None, None,
                                    fpc, fy, fu, PDeltaColumnGravityLoad,
                              Mass, WallGravityLoad,
                                    CustomWallSection=[Section1, Section1, Section1, Section1, Section1,
                                                       Section1, Section2, Section2, Section2, Section2,
                                                       Section3, Section3, Section3, Section3, Section4,
                                                       Section4, Section4, Section4, Section5, Section5,],
                                    heff = 108.*12.)

#endregion

#region ATCWALLS20O1S50ELF
NoOfStories = 20
YGrids = [0] + np.array(np.arange(0,(NoOfStories)*13*12, 13*12)+15*12).tolist()

# Loading
# DL 175psf Floors and 140psf Roof
# LL 65psf Floors and 20psf Roof

MassPerSqFt = (210.)/1000. # use this 200 psf
FloorArea = 120.*120./4.
Mass = np.ones(NoOfStories) * MassPerSqFt * FloorArea
Mass[-1] = FloorArea*105./1000. # Adjust for Roof Weight

WallTribArea = 30.*30.
WeightPerSqFt = 210./1000.

BuildingWeight = np.ones(NoOfStories) * WeightPerSqFt * FloorArea
BuildingWeight[-1] = 105. / 1000. * FloorArea # Adjust for Roof Weight

PDeltaColumnGravityLoad = -1 * BuildingWeight * (FloorArea - WallTribArea) / FloorArea

WallGravityLoad = BuildingWeight * WallTribArea / FloorArea * -1
# WallGravityLoad = Mass * WallTribArea / FloorArea * -1 # This is what Kamal used in Atena

# Material Properties
fy = 60
fu = 105
fpc = 5

OpeningWdith = 16.*12.
Section1A = TWallSection((402.-OpeningWdith)/2., 96., 18., 1.10, 10.127,
                        fpc, fy, fu,
                        None, None, None, None, True)

Section1B = TWallSection((402.-OpeningWdith)/2., 96., 18., 1.10, 10.127,
                        fpc, fy, fu,
                        None, None, None, None, False)

Section1 = IWallSection(402., 96., 18., 1.10, 10.127,
                        fpc, fy, fu,
                        None, None, None, None)

Section2 = IWallSection(402., 72., 18., 1.05, 10.127,
                        fpc, fy, fu,
                        None, None, None, None)

l_w = 402.
t_w = 18.
Rho = 1.22
bar_size = 10.173
NoOfBarRows = np.ceil(l_w/2.*t_w*Rho/100./(np.pi/4.*(bar_size/8.)**2.*2.))
tie_spacing = ATCWallHelper.GetTieSpacing(t_w, bar_size/8., l_w/NoOfBarRows/2.)
Section3 = PlanarWallSection(l_w, t_w, l_w/2., l_w/2., bar_size,
                            np.ones(int(NoOfBarRows))*2., np.ones(int(NoOfBarRows))*2., 0.0, 4.,
                            fpc, fy, fu,
                            tie_spacing, 4.0, 18, 2)

l_w = 360.
t_w = 16.
Rho = 1.04
bar_size = 9.027
NoOfBarRows = np.ceil(l_w/2.*t_w*Rho/100./(np.pi/4.*(bar_size/8.)**2.*2.))
tie_spacing = ATCWallHelper.GetTieSpacing(t_w, bar_size/8., l_w/NoOfBarRows/2.)
Section4 = PlanarWallSection(l_w, t_w, l_w/2., l_w/2., bar_size,
                            np.ones(int(NoOfBarRows))*2., np.ones(int(NoOfBarRows))*2., 0.0, 4.,
                            fpc, fy, fu,
                            tie_spacing, 4.0, 18, 2)

l_w = 360.
t_w = 14.
Rho = 0.95
bar_size = 9.027
NoOfBarRows = np.ceil(l_w/2.*t_w*Rho/100./(np.pi/4.*(bar_size/8.)**2.*2.))
tie_spacing = ATCWallHelper.GetTieSpacing(t_w, bar_size/8., l_w/NoOfBarRows/2.)
Section5 = PlanarWallSection(l_w, t_w, l_w/2., l_w/2., bar_size,
                            np.ones(int(NoOfBarRows))*2., np.ones(int(NoOfBarRows))*2., 0.0, 4.,
                            fpc, fy, fu,
                            tie_spacing, 4.0, 18, 2)

CouplingBeams = [None, None, None, None, None, None, None, None]
CouplingBeamLengths = [OpeningWdith, None, None, None, None, None, None, None]

ATCWALLS20O1S50ELF = ArchetypeData('ATCWALLS20O1S50ELF', YGrids, None, None, None, None, None, None, None,
                                    fpc, fy, fu, PDeltaColumnGravityLoad,
                                    Mass, WallGravityLoad,
                                    CustomWallSection=[[Section1A, Section1B], Section1, Section1, Section1, Section1,
                                                       Section1, Section2, Section2, Section2, Section2,
                                                       Section3, Section3, Section3, Section3, Section4,
                                                       Section4, Section4, Section4, Section5, Section5,],
                                   CouplingBeamLength=CouplingBeamLengths,
                                   CouplingBeams=CouplingBeams,
                                   heff=108. * 12.)

#endregion

#region ATCWALLS20O1S75ELF
NoOfStories = 20
YGrids = [0] + np.array(np.arange(0,(NoOfStories)*13*12, 13*12)+15*12).tolist()

# Loading
# DL 175psf Floors and 140psf Roof
# LL 65psf Floors and 20psf Roof

MassPerSqFt = (210.)/1000. # use this 200 psf
FloorArea = 120.*120./4.
Mass = np.ones(NoOfStories) * MassPerSqFt * FloorArea
Mass[-1] = FloorArea*105./1000. # Adjust for Roof Weight

WallTribArea = 30.*30.
WeightPerSqFt = 210./1000.

BuildingWeight = np.ones(NoOfStories) * WeightPerSqFt * FloorArea
BuildingWeight[-1] = 105. / 1000. * FloorArea # Adjust for Roof Weight

PDeltaColumnGravityLoad = -1 * BuildingWeight * (FloorArea - WallTribArea) / FloorArea

WallGravityLoad = BuildingWeight * WallTribArea / FloorArea * -1
# WallGravityLoad = Mass * WallTribArea / FloorArea * -1 # This is what Kamal used in Atena

# Material Properties
fy = 60
fu = 105
fpc = 5

OpeningWdith = 8.*12.
Section1A = TWallSection((402.-OpeningWdith)/2., 96., 18., 1.10, 10.127,
                        fpc, fy, fu,
                        None, None, None, None, True)

Section1B = TWallSection((402.-OpeningWdith)/2., 96., 18., 1.10, 10.127,
                        fpc, fy, fu,
                        None, None, None, None, False)

Section1 = IWallSection(402., 96., 18., 1.10, 10.127,
                        fpc, fy, fu,
                        None, None, None, None)

Section2 = IWallSection(402., 72., 18., 1.05, 10.127,
                        fpc, fy, fu,
                        None, None, None, None)

l_w = 402.
t_w = 18.
Rho = 1.22
bar_size = 10.173
NoOfBarRows = np.ceil(l_w/2.*t_w*Rho/100./(np.pi/4.*(bar_size/8.)**2.*2.))
tie_spacing = ATCWallHelper.GetTieSpacing(t_w, bar_size/8., l_w/NoOfBarRows/2.)
Section3 = PlanarWallSection(l_w, t_w, l_w/2., l_w/2., bar_size,
                            np.ones(int(NoOfBarRows))*2., np.ones(int(NoOfBarRows))*2., 0.0, 4.,
                            fpc, fy, fu,
                            tie_spacing, 4.0, 18, 2)

l_w = 360.
t_w = 16.
Rho = 1.04
bar_size = 9.027
NoOfBarRows = np.ceil(l_w/2.*t_w*Rho/100./(np.pi/4.*(bar_size/8.)**2.*2.))
tie_spacing = ATCWallHelper.GetTieSpacing(t_w, bar_size/8., l_w/NoOfBarRows/2.)
Section4 = PlanarWallSection(l_w, t_w, l_w/2., l_w/2., bar_size,
                            np.ones(int(NoOfBarRows))*2., np.ones(int(NoOfBarRows))*2., 0.0, 4.,
                            fpc, fy, fu,
                            tie_spacing, 4.0, 18, 2)

l_w = 360.
t_w = 14.
Rho = 0.95
bar_size = 9.027
NoOfBarRows = np.ceil(l_w/2.*t_w*Rho/100./(np.pi/4.*(bar_size/8.)**2.*2.))
tie_spacing = ATCWallHelper.GetTieSpacing(t_w, bar_size/8., l_w/NoOfBarRows/2.)
Section5 = PlanarWallSection(l_w, t_w, l_w/2., l_w/2., bar_size,
                            np.ones(int(NoOfBarRows))*2., np.ones(int(NoOfBarRows))*2., 0.0, 4.,
                            fpc, fy, fu,
                            tie_spacing, 4.0, 18, 2)

CouplingBeams = [None, None, None, None, None, None, None, None]
CouplingBeamLengths = [OpeningWdith, None, None, None, None, None, None, None]

ATCWALLS20O1S75ELF = ArchetypeData('ATCWALLS20O1S75ELF', YGrids, None, None, None, None, None, None, None,
                                    fpc, fy, fu, PDeltaColumnGravityLoad,
                                    Mass, WallGravityLoad,
                                    CustomWallSection=[[Section1A, Section1B], Section1, Section1, Section1, Section1,
                                                       Section1, Section2, Section2, Section2, Section2,
                                                       Section3, Section3, Section3, Section3, Section4,
                                                       Section4, Section4, Section4, Section5, Section5,],
                                   CouplingBeamLength=CouplingBeamLengths,
                                   CouplingBeams=CouplingBeams,
                                   heff=108. * 12.)

#endregion

### B Max Designs

# 8 Story B Max Designs

#region ATCWALLS8O0S100ELFBMAX

# Geometry
NoOfStories = 8
YGrids = [0] + np.array(np.arange(0,(NoOfStories)*13*12, 13*12)+15*12).tolist()

# Loading
# DL 175psf Floors and 140psf Roof
# LL 65psf Floors and 20psf Roof

MassPerSqFt = (199.5)/1000. # use this 200 psf
FloorArea = 60.*120.
Mass = np.ones(NoOfStories) * MassPerSqFt * FloorArea
Mass[-1] = FloorArea*154.9/1000. # Adjust for Roof Weight

WallTribArea = 30.*34.
WeightPerSqFt = 212./1000.

BuildingWeight = np.ones(NoOfStories) * WeightPerSqFt * FloorArea
BuildingWeight[-1] = 152. / 1000. * FloorArea # Adjust for Roof Weight

PDeltaColumnGravityLoad = -1 * BuildingWeight * (FloorArea - WallTribArea) / FloorArea

WallGravityLoad = BuildingWeight * WallTribArea / FloorArea * -1
# WallGravityLoad = Mass * WallTribArea / FloorArea * -1 # This is what Kamal used in Atena

# Material Properties
fy = 60
fu = 105
fpc = 5.0

Section1 = PlanarWallSection(20.*12., 24., 0, 0, 10.173,
                             [],
                             [],
                             0.25, 4.037, fpc, fy, fu, None, None, None, None)

ATCWALLS8O0S100ELFBMAX = ArchetypeData('ATCWALLS8O0S100ELFBMAX', YGrids, None, None, None, None, None, None, None, fpc, fy, fu, PDeltaColumnGravityLoad,
                                        Mass, WallGravityLoad,
                                        CustomWallSection=[Section1, Section1, Section1,
                                                          Section1, Section1, Section1,
                                                          Section1, Section1],
                                        heff = 73.57*12.)

#endregion

#region ATCWALLS8O1S50ELFBMAX

# Geometry
NoOfStories = 8
YGrids = [0] + np.array(np.arange(0,(NoOfStories)*13*12, 13*12)+15*12).tolist()

# Loading
# DL 175psf Floors and 140psf Roof
# LL 65psf Floors and 20psf Roof

MassPerSqFt = (199.5)/1000. # use this 200 psf
FloorArea = 60.*120.
Mass = np.ones(NoOfStories) * MassPerSqFt * FloorArea
Mass[-1] = FloorArea*154.9/1000. # Adjust for Roof Weight

WallTribArea = 30.*34.
WeightPerSqFt = 212./1000.

BuildingWeight = np.ones(NoOfStories) * WeightPerSqFt * FloorArea
BuildingWeight[-1] = 152. / 1000. * FloorArea # Adjust for Roof Weight

PDeltaColumnGravityLoad = -1 * BuildingWeight * (FloorArea - WallTribArea) / FloorArea

WallGravityLoad = BuildingWeight * WallTribArea / FloorArea * -1
# WallGravityLoad = Mass * WallTribArea / FloorArea * -1 # This is what Kamal used in Atena


# Material Properties
fy = 60.
fu = 105.
fpc = 5.

Section1A = PlanarWallSection(7.5*12., 24., 0, 0, 0, [], [],
                              0.54, 5., fpc, fy, fu,  None, None, None, None)

Section1B = PlanarWallSection(7.5*12., 24., 0, 0, 0, [], [],
                              0.54, 5., fpc, fy, fu,  None, None, None, None)

Section1 = PlanarWallSection(20.*12., 24., 0, 0, 0,
                             [],
                             [],
                             0.25, 4.037, fpc, fy, fu, None, None, None, None)

CouplingBeams = [None, None, None, None, None, None, None, None]
CouplingBeamLengths = [5.0*12., None, None, None, None, None, None, None]

ATCWALLS8O1S50ELFBMAX = ArchetypeData('ATCWALLS8O1S50ELFBMAX', YGrids, None, None, None, None, None, None, None, fpc, fy, fu, PDeltaColumnGravityLoad,
                                Mass, WallGravityLoad,
                                CustomWallSection=[[Section1A, Section1B], Section1, Section1,
                                                   Section1, Section1, Section1,
                                                   Section1, Section1],
                                CouplingBeamLength=CouplingBeamLengths,
                                CouplingBeams=CouplingBeams, heff = 73.57*12.)

#endregion

#region ATCWALLS8O1S75ELFBMAX

# Geometry
NoOfStories = 8
YGrids = [0] + np.array(np.arange(0,(NoOfStories)*13*12, 13*12)+15*12).tolist()

# Loading
# DL 175psf Floors and 140psf Roof
# LL 65psf Floors and 20psf Roof

MassPerSqFt = (199.5)/1000. # use this 200 psf
FloorArea = 60.*120.
Mass = np.ones(NoOfStories) * MassPerSqFt * FloorArea
Mass[-1] = FloorArea*154.9/1000. # Adjust for Roof Weight

WallTribArea = 30.*34.
WeightPerSqFt = 212./1000.

BuildingWeight = np.ones(NoOfStories) * WeightPerSqFt * FloorArea
BuildingWeight[-1] = 152. / 1000. * FloorArea # Adjust for Roof Weight

PDeltaColumnGravityLoad = -1 * BuildingWeight * (FloorArea - WallTribArea) / FloorArea

WallGravityLoad = BuildingWeight * WallTribArea / FloorArea * -1
# WallGravityLoad = Mass * WallTribArea / FloorArea * -1 # This is what Kamal used in Atena


# Material Properties
fy = 60.
fu = 105.
fpc = 5.

Section1A = PlanarWallSection(9.25*12., 24., 0, 0, 0, [], [],
                              0.54, 4., fpc, fy, fu,  None, None, None, None)

Section1B = PlanarWallSection(9.25*12., 24., 0, 0, 0, [], [],
                              0.54, 4., fpc, fy, fu,  None, None, None, None)

Section1 = PlanarWallSection(20.*12., 24., 0, 0, 0,
                             [],
                             [],
                             0.25, 4.037, fpc, fy, fu, None, None, None, None)

CouplingBeams = [None, None, None, None, None, None, None, None]
CouplingBeamLengths = [1.5*12., None, None, None, None, None, None, None]

ATCWALLS8O1S75ELFBMAX = ArchetypeData('ATCWALLS8O1S75ELFBMAX', YGrids, None, None, None, None, None, None, None, fpc, fy, fu, PDeltaColumnGravityLoad,
                                Mass, WallGravityLoad,
                                CustomWallSection=[[Section1A, Section1B], Section1, Section1,
                                                   Section1, Section1, Section1,
                                                   Section1, Section1],
                                CouplingBeamLength=CouplingBeamLengths,
                                CouplingBeams=CouplingBeams, heff = 73.57*12.)

#endregion

#region ATCWALLS8O2S50ELFBMAX

# Geometry
NoOfStories = 8
YGrids = [0] + np.array(np.arange(0,(NoOfStories)*13*12, 13*12)+15*12).tolist()

# Loading
# DL 175psf Floors and 140psf Roof
# LL 65psf Floors and 20psf Roof

MassPerSqFt = (199.5)/1000. # use this 200 psf
FloorArea = 60.*120.
Mass = np.ones(NoOfStories) * MassPerSqFt * FloorArea
Mass[-1] = FloorArea*154.9/1000. # Adjust for Roof Weight

WallTribArea = 30.*34.
WeightPerSqFt = 212./1000.

BuildingWeight = np.ones(NoOfStories) * WeightPerSqFt * FloorArea
BuildingWeight[-1] = 152. / 1000. * FloorArea # Adjust for Roof Weight

PDeltaColumnGravityLoad = -1 * BuildingWeight * (FloorArea - WallTribArea) / FloorArea

WallGravityLoad = BuildingWeight * WallTribArea / FloorArea * -1
# WallGravityLoad = Mass * WallTribArea / FloorArea * -1 # This is what Kamal used in Atena


# Material Properties
fy = 60.
fu = 105.
fpc = 5.

Section1A = PlanarWallSection(7.5*12., 24., 0, 0, 0, [], [],
                              0.54, 5., fpc, fy, fu,  None, None, None, None)

Section1B = PlanarWallSection(7.5*12., 24., 0, 0, 0, [], [],
                              0.54, 5., fpc, fy, fu,  None, None, None, None)

Section1 = PlanarWallSection(20.*12., 24., 0, 0, 0,
                             [],
                             [],
                             0.25, 4.037, fpc, fy, fu, None, None, None, None)

CB1 = CouplingBeam(24., 60., fpc, fy, 2, 2, 0.1, True, 4., 6., 2.5, 4, 3, 10)

CouplingBeams = [CB1, None, None, None, None, None, None, None]
CouplingBeamLengths = [5.0*12., 5.0*12., None, None, None, None, None, None]

ATCWALLS8O2S50ELFBMAX = ArchetypeData('ATCWALLS8O2S50ELFBMAX', YGrids, None, None, None, None, None, None, None, fpc, fy, fu, PDeltaColumnGravityLoad,
                                Mass, WallGravityLoad,
                                CustomWallSection=[[Section1A, Section1B], [Section1A, Section1B], Section1,
                                                   Section1, Section1, Section1,
                                                   Section1, Section1],
                                CouplingBeamLength=CouplingBeamLengths,
                                CouplingBeams=CouplingBeams, heff = 73.57*12.)

#endregion

#region ATCWALLS8O2S75ELFBMAX

# Geometry
NoOfStories = 8
YGrids = [0] + np.array(np.arange(0,(NoOfStories)*13*12, 13*12)+15*12).tolist()

# Loading
# DL 175psf Floors and 140psf Roof
# LL 65psf Floors and 20psf Roof

MassPerSqFt = (199.5)/1000. # use this 200 psf
FloorArea = 60.*120.
Mass = np.ones(NoOfStories) * MassPerSqFt * FloorArea
Mass[-1] = FloorArea*154.9/1000. # Adjust for Roof Weight

WallTribArea = 30.*34.
WeightPerSqFt = 212./1000.

BuildingWeight = np.ones(NoOfStories) * WeightPerSqFt * FloorArea
BuildingWeight[-1] = 152. / 1000. * FloorArea # Adjust for Roof Weight

PDeltaColumnGravityLoad = -1 * BuildingWeight * (FloorArea - WallTribArea) / FloorArea

WallGravityLoad = BuildingWeight * WallTribArea / FloorArea * -1
# WallGravityLoad = Mass * WallTribArea / FloorArea * -1 # This is what Kamal used in Atena


# Material Properties
fy = 60.
fu = 105.
fpc = 5.

Section1A = PlanarWallSection(9.25*12., 24., 0, 0, 0, [], [],
                              0.54, 4., fpc, fy, fu,  None, None, None, None)

Section1B = PlanarWallSection(9.25*12., 24., 0, 0, 0, [], [],
                              0.54, 4., fpc, fy, fu,  None, None, None, None)

Section1 = PlanarWallSection(20.*12., 24., 0, 0, 0,
                             [],
                             [],
                             0.25, 4.037, fpc, fy, fu, None, None, None, None)

CB1 = CouplingBeam(24., 60., fpc, fy, 2, 2, 0.1, True, 4., 6., 2.5, 4, 3, 10)

CouplingBeams = [CB1, None, None, None, None, None, None, None]
CouplingBeamLengths = [1.5*12., 1.5*12., None, None, None, None, None, None]

ATCWALLS8O2S75ELFBMAX = ArchetypeData('ATCWALLS8O2S75ELFBMAX', YGrids, None, None, None, None, None, None, None, fpc, fy, fu, PDeltaColumnGravityLoad,
                                Mass, WallGravityLoad,
                                CustomWallSection=[[Section1A, Section1B], [Section1A, Section1B], Section1,
                                                   Section1, Section1, Section1,
                                                   Section1, Section1],
                                CouplingBeamLength=CouplingBeamLengths,
                                CouplingBeams=CouplingBeams, heff = 73.57*12.)

#endregion

## 12 Story BMax Designs

#region ATCWALLS12O0S100ELFBMAX

# Geometry
NoOfStories = 12
YGrids = [0] + np.array(np.arange(0,(NoOfStories)*13*12, 13*12)+15*12).tolist()

# Loading
# DL 175psf Floors and 140psf Roof
# LL 65psf Floors and 20psf Roof

MassPerSqFt = (199.5)/1000. # use this 200 psf
FloorArea = 60.*60.
Mass = np.ones(NoOfStories) * MassPerSqFt * FloorArea
Mass[-1] = FloorArea*154.9/1000. # Adjust for Roof Weight

WallTribArea = 30.*34.
WeightPerSqFt = 212./1000.

BuildingWeight = np.ones(NoOfStories) * WeightPerSqFt * FloorArea
BuildingWeight[-1] = 152. / 1000. * FloorArea # Adjust for Roof Weight

PDeltaColumnGravityLoad = -1 * BuildingWeight * (FloorArea - WallTribArea) / FloorArea

WallGravityLoad = BuildingWeight * WallTribArea / FloorArea * -1
# WallGravityLoad = Mass * WallTribArea / FloorArea * -1 # This is what Kamal used in Atena

# Material Properties
fy = 60
fu = 105
fpc = 5.0

Section1 = PlanarWallSection(20.*12., 24., 0, 0, 10.173,
                             [],
                             [],
                             0.25, 4.037, fpc, fy, fu, None, None, None, None)

ATCWALLS12O0S100ELFBMAX = ArchetypeData('ATCWALLS12O0S100ELFBMAX', YGrids, None, None, None, None, None, None, None, fpc, fy, fu, PDeltaColumnGravityLoad,
                                        Mass, WallGravityLoad,
                                        CustomWallSection=[Section1, Section1, Section1,
                                                          Section1, Section1, Section1,
                                                          Section1, Section1, Section1,
                                                          Section1, Section1, Section1,
                                                           ],
                                        heff = 73.57*12.)

#endregion

#region ATCWALLS12O1S50ELFBMAX

# Geometry
NoOfStories = 12
YGrids = [0] + np.array(np.arange(0,(NoOfStories)*13*12, 13*12)+15*12).tolist()

# Loading
# DL 175psf Floors and 140psf Roof
# LL 65psf Floors and 20psf Roof

MassPerSqFt = (199.5)/1000. # use this 200 psf
FloorArea = 60.*60.
Mass = np.ones(NoOfStories) * MassPerSqFt * FloorArea
Mass[-1] = FloorArea*154.9/1000. # Adjust for Roof Weight

WallTribArea = 30.*34.
WeightPerSqFt = 212./1000.

BuildingWeight = np.ones(NoOfStories) * WeightPerSqFt * FloorArea
BuildingWeight[-1] = 152. / 1000. * FloorArea # Adjust for Roof Weight

PDeltaColumnGravityLoad = -1 * BuildingWeight * (FloorArea - WallTribArea) / FloorArea

WallGravityLoad = BuildingWeight * WallTribArea / FloorArea * -1
# WallGravityLoad = Mass * WallTribArea / FloorArea * -1 # This is what Kamal used in Atena

# Material Properties
fy = 60
fu = 105
fpc = 5.0

Section1A = PlanarWallSection(7.5*12., 24., 0, 0, 0, [], [],
                              0.55, 5., fpc, fy, fu,  None, None, None, None)

Section1B = PlanarWallSection(7.5*12., 24., 0, 0, 0, [], [],
                              0.55, 5., fpc, fy, fu,  None, None, None, None)

Section1 = PlanarWallSection(20.*12., 24., 0, 0, 0,
                             [],
                             [],
                             0.25, 4.037, fpc, fy, fu, None, None, None, None)

CouplingBeams = [None, None, None, None, None, None, None, None]
CouplingBeamLengths = [5.0*12., None, None, None, None, None, None, None]

ATCWALLS12O1S50ELFBMAX = ArchetypeData('ATCWALLS12O1S50ELFBMAX', YGrids, None, None, None, None, None, None, None, fpc, fy, fu, PDeltaColumnGravityLoad,
                                Mass, WallGravityLoad,
                                CustomWallSection=[[Section1A, Section1B], Section1, Section1,
                                                   Section1, Section1, Section1,
                                                   Section1, Section1, Section1,
                                                   Section1, Section1, Section1,],
                                CouplingBeamLength=CouplingBeamLengths,
                                CouplingBeams=CouplingBeams, heff = 73.57*12.)

#endregion

#region ATCWALLS12O1S75ELFBMAX


# Geometry
NoOfStories = 12
YGrids = [0] + np.array(np.arange(0,(NoOfStories)*13*12, 13*12)+15*12).tolist()

# Loading
# DL 175psf Floors and 140psf Roof
# LL 65psf Floors and 20psf Roof

MassPerSqFt = (199.5)/1000. # use this 200 psf
FloorArea = 60.*60.
Mass = np.ones(NoOfStories) * MassPerSqFt * FloorArea
Mass[-1] = FloorArea*154.9/1000. # Adjust for Roof Weight

WallTribArea = 30.*34.
WeightPerSqFt = 212./1000.

BuildingWeight = np.ones(NoOfStories) * WeightPerSqFt * FloorArea
BuildingWeight[-1] = 152. / 1000. * FloorArea # Adjust for Roof Weight

PDeltaColumnGravityLoad = -1 * BuildingWeight * (FloorArea - WallTribArea) / FloorArea

WallGravityLoad = BuildingWeight * WallTribArea / FloorArea * -1
# WallGravityLoad = Mass * WallTribArea / FloorArea * -1 # This is what Kamal used in Atena

# Material Properties
fy = 60
fu = 105
fpc = 5.0

Section1A = PlanarWallSection(9.25*12., 24., 0, 0, 0, [], [],
                              0.55, 4., fpc, fy, fu,  None, None, None, None)

Section1B = PlanarWallSection(9.25*12., 24., 0, 0, 0, [], [],
                              0.55, 4., fpc, fy, fu,  None, None, None, None)

Section1 = PlanarWallSection(20.*12., 24., 0, 0, 0,
                             [],
                             [],
                             0.25, 4.037, fpc, fy, fu, None, None, None, None)

CouplingBeams = [None, None, None, None, None, None, None, None]
CouplingBeamLengths = [1.5*12., None, None, None, None, None, None, None]

ATCWALLS12O1S75ELFBMAX = ArchetypeData('ATCWALLS12O1S75ELFBMAX', YGrids, None, None, None, None, None, None, None, fpc, fy, fu, PDeltaColumnGravityLoad,
                                Mass, WallGravityLoad,
                                CustomWallSection=[[Section1A, Section1B], Section1, Section1,
                                                   Section1, Section1, Section1,
                                                   Section1, Section1, Section1,
                                                   Section1, Section1, Section1,],
                                CouplingBeamLength=CouplingBeamLengths,
                                CouplingBeams=CouplingBeams, heff = 73.57*12.)

#endregion

#region ATCWALLS12O2S50ELFBMAX

# Geometry
NoOfStories = 12
YGrids = [0] + np.array(np.arange(0,(NoOfStories)*13*12, 13*12)+15*12).tolist()

# Loading
# DL 175psf Floors and 140psf Roof
# LL 65psf Floors and 20psf Roof

MassPerSqFt = (199.5)/1000. # use this 200 psf
FloorArea = 60.*60.
Mass = np.ones(NoOfStories) * MassPerSqFt * FloorArea
Mass[-1] = FloorArea*154.9/1000. # Adjust for Roof Weight

WallTribArea = 30.*44.
WeightPerSqFt = 212./1000.

BuildingWeight = np.ones(NoOfStories) * WeightPerSqFt * FloorArea
BuildingWeight[-1] = 152. / 1000. * FloorArea # Adjust for Roof Weight

PDeltaColumnGravityLoad = -1 * BuildingWeight * (FloorArea - WallTribArea) / FloorArea

WallGravityLoad = BuildingWeight * WallTribArea / FloorArea * -1
# WallGravityLoad = Mass * WallTribArea / FloorArea * -1 # This is what Kamal used in Atena


# Material Properties
fy = 60.
fu = 105.
fpc = 5.

Section1A = PlanarWallSection(7.5*12., 24., 0, 0, 0, [], [],
                              0.55, 5., fpc, fy, fu,  None, None, None, None)

Section1B = PlanarWallSection(7.5*12., 24., 0, 0, 0, [], [],
                              0.55, 5., fpc, fy, fu,  None, None, None, None)

Section1 = PlanarWallSection(20.*12., 24., 0, 0, 0,
                             [],
                             [],
                             0.25, 4.037, fpc, fy, fu, None, None, None, None)

CB1 = CouplingBeam(24., 60., fpc, fy, 2, 2, 0.1, True, 4., 6., 2.5, 4, 3, 10)

CouplingBeams = [CB1, None, None, None, None, None, None, None]
CouplingBeamLengths = [5.0*12., 5.0*12., None, None, None, None, None, None]

ATCWALLS12O2S50ELFBMAX = ArchetypeData('ATCWALLS12O2S50ELFBMAX', YGrids, None, None, None, None, None, None, None, fpc, fy, fu, PDeltaColumnGravityLoad,
                                Mass, WallGravityLoad,
                                CustomWallSection=[[Section1A, Section1B], [Section1A, Section1B], Section1,
                                                   Section1, Section1, Section1,
                                                   Section1, Section1, Section1,
                                                   Section1, Section1, Section1,],
                                CouplingBeamLength=CouplingBeamLengths,
                                CouplingBeams=CouplingBeams, heff = 73.57*12.)

#endregion

#region ATCWALLS12O2S75ELFBMAX

# Geometry
NoOfStories = 12
YGrids = [0] + np.array(np.arange(0,(NoOfStories)*13*12, 13*12)+15*12).tolist()

# Loading
# DL 175psf Floors and 140psf Roof
# LL 65psf Floors and 20psf Roof

MassPerSqFt = (199.5)/1000. # use this 200 psf
FloorArea = 60.*60.
Mass = np.ones(NoOfStories) * MassPerSqFt * FloorArea
Mass[-1] = FloorArea*154.9/1000. # Adjust for Roof Weight

WallTribArea = 30.*44.
WeightPerSqFt = 212./1000.

BuildingWeight = np.ones(NoOfStories) * WeightPerSqFt * FloorArea
BuildingWeight[-1] = 152. / 1000. * FloorArea # Adjust for Roof Weight

PDeltaColumnGravityLoad = -1 * BuildingWeight * (FloorArea - WallTribArea) / FloorArea

WallGravityLoad = BuildingWeight * WallTribArea / FloorArea * -1
# WallGravityLoad = Mass * WallTribArea / FloorArea * -1 # This is what Kamal used in Atena


# Material Properties
fy = 60.
fu = 105.
fpc = 5.

Section1A = PlanarWallSection(9.25*12., 24., 0, 0, 0, [], [],
                              0.55, 4., fpc, fy, fu,  None, None, None, None)

Section1B = PlanarWallSection(9.25*12., 24., 0, 0, 0, [], [],
                              0.55, 4., fpc, fy, fu,  None, None, None, None)

Section1 = PlanarWallSection(20.*12., 24., 0, 0, 0,
                             [],
                             [],
                             0.25, 4.037, fpc, fy, fu, None, None, None, None)

CB1 = CouplingBeam(24., 60., fpc, fy, 2, 2, 0.1, True, 4., 6., 2.5, 4, 3, 10)

CouplingBeams = [CB1, None, None, None, None, None, None, None]
CouplingBeamLengths = [1.5*12., 1.5*12., None, None, None, None, None, None]

ATCWALLS12O2S75ELFBMAX = ArchetypeData('ATCWALLS12O2S75ELFBMAX', YGrids, None, None, None, None, None, None, None, fpc, fy, fu, PDeltaColumnGravityLoad,
                                Mass, WallGravityLoad,
                                CustomWallSection=[[Section1A, Section1B], [Section1A, Section1B], Section1,
                                                   Section1, Section1, Section1,
                                                   Section1, Section1, Section1,
                                                   Section1, Section1, Section1,],
                                CouplingBeamLength=CouplingBeamLengths,
                                CouplingBeams=CouplingBeams, heff = 73.57*12.)

#endregion

# Seattle Archetypes

####################################################################################
#endregion
####################################################################################

Archetypes = [  ATCWALLS4O0S100ELF,
                ATCWALLS4O4S50ELF,

                ATCWALLS8O0S100ELF,
                ATCWALLS8O1S50ELF,
                ATCWALLS8O1S75ELF,
                ATCWALLS8O2S50ELF,
                ATCWALLS8O2S75ELF,
                ATCWALLS8O8S50ELF,

                ATCWALLS12O0S100ELF,
                ATCWALLS12O1S50ELF,
                ATCWALLS12O1S75ELF,
                ATCWALLS12O2S50ELF,
                ATCWALLS12O2S75ELF,
                ATCWALLS12O12S50ELF,

                # ATCWALLS20O0S100ELF,
                # ATCWALLS20O1S50ELF,
                # ATCWALLS20O1S75ELF,

                ATCWALLS8O0S100MRSA,
                ATCWALLS8O1S50MRSA,
                ATCWALLS8O1S75MRSA,
                ATCWALLS8O2S50MRSA,
                ATCWALLS8O2S75MRSA,

                ATCWALLS12O0S100MRSA,
                ATCWALLS12O1S50MRSA,
                ATCWALLS12O1S75MRSA,
                ATCWALLS12O2S50MRSA,
                ATCWALLS12O2S75MRSA,

                ATCWALLS8O0S100ELFBMAX,
                ATCWALLS8O1S50ELFBMAX,
                ATCWALLS8O1S75ELFBMAX,
                ATCWALLS8O2S50ELFBMAX,
                ATCWALLS8O2S75ELFBMAX,

                ATCWALLS12O0S100ELFBMAX,
                ATCWALLS12O1S50ELFBMAX,
                ATCWALLS12O1S75ELFBMAX,
                ATCWALLS12O2S50ELFBMAX,
                ATCWALLS12O2S75ELFBMAX,

                ATCWALLS8O1AT5S50ELF,
                ATCWALLS8O8S50ELF48,

               ]
