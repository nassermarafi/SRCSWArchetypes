####################################################################################
#region Libraries
####################################################################################

from __future__ import absolute_import
import numpy as np
import os
import OpenSeesAPI
import ATCWallArchetypeHelpers
import ATCWallArchetypeObjects as ATCWallObjects
from six.moves import map
from six.moves import range

####################################################################################
#endregion
####################################################################################

####################################################################################
#region Defining Classes
####################################################################################

# from WallArchetypes.ArchetypeBuilder import ArchetypeData
# from WallArchetypes.ArchetypeBuilder import CouplingBeam
# from WallArchetypes.ArchetypeBuilder import PlanarWallSection

class ArchetypeData:
    def __init__(self, Name, YGrids, R, T1, l_w, t, b_f, rho, rho_t, fpc, fy, fu, GravityLoad, Mass, WallGravityLoad,
                 BoundaryElement = None, CouplingBeams = None, CouplingBeamLength = None, CustomWallSection=None, **kwargs):
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
        self.GravityLoad = GravityLoad
        self.Mass = Mass
        self.WallGravityLoad = WallGravityLoad
        self.BoundaryElement = BoundaryElement

        self.CustomSection = CustomWallSection

        self.fce = 1.3 * fpc #TBI Table 7.1
        self.fye = 1.17 * fy #TBI Table 7.1

        self.CouplingBeams = CouplingBeams
        self.CouplingBeamLength = CouplingBeamLength

        self.__dict__.update(kwargs)

class CouplingBeam:
    def __init__(self, b, h, fpc, fy, NoOfBarsX, NoOfBarsY, BarDia, DiagonalReinf=True, **kwargs):
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
        self.TieSpacing = 6
        self.TieDia = 5
        self.DiagonalReinf = DiagonalReinf

        self.fce = 1.3*fpc #TBI Table 7.1
        self.fye = 1.17*fy #TBI Table 7.1

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

####################################################################################

# endregion
####################################################################################

####################################################################################
#region Defining Methods
####################################################################################

# from WallArchetypes.ArchetypeBuilder import AnalyzeArchetypeSingleLineElement

def AnalyzeArchetypeSingleLineElement(Archetype, GMData, Dt, SupressOutput=True, Viewer=False,
                                          Animation=False, TimeHistory=True, PushOver=False, OpenSeesCommand='OpenSeesSP',
                                          T1=None, T2=None, EnhancedOutput=False, PDeltaColumn=True,
                                          CyclicStatic=False, DriftHistory=None, POModalForces=True,
                                          heff=None, ConstantWallAxialForce=False, ApplyPDeltaLoad=True,
                                          DebugMode=False, RegularizeSteel=False, MaxPORoofDrift=0.05,
                                          PlotSections=False, RegularizeFracture=False, NoOfIterations=10000,
                                          Tolerance=1e-6, CrushingStrengthRatio=0.2, POELFForces=False, CuTa=None,
                                          GfccOGfc=1.75, ConcreteMaterialModel='Concrete02', ConfinementModel='SaatRazvi',
                                          SteelMaterialModel='Steel02', SteelUltimateStrainTension=0.2,
                                          SteelPostYieldStiffness=0.006, WallAxialLoadMultiplier=1.0,
                                          UseForceBasedElements=False, DivisionsPerStory=6,
                                          GfcOfpc=2.0, NoOfIntPoints=5, UnconfinedBeta = 0.01, Regularized = True,
                                          WallThicknessMultipler = 1.0, FBE_Tolerance = 1e-6,
                                          ModalDamping = False, Zeta = 0.02, UseTBIZeta = False,
                                          IncludeSupplementalRayleigh = True,
                                          Options=None, TrackPeriod = False,
                                          HHTTransientIntegrator=False,
                                          OutputTag = '', TrackDrifts = False, **kwargs
                                          ):
    ### Adding Options to Global Variables
    OtherOptions = {}
    if Options is not None:
        for key in Options:
            if key not in locals():
                if type(Options[key]) == str:
                    exec(("%s = \'%s\'"%(key, Options[key])), globals(), locals())
                elif type(Options[key]) == float:
                    exec(("%s = %f" % (key, Options[key])), globals(), locals())
                elif type(Options[key]) == bool:
                    exec(("%s = %s"% (key, Options[key])), globals(), locals())
                OtherOptions[key] = Options[key]
            else:
                if key == 'SupressOutput':
                    SupressOutput = Options[key]
                elif key == 'Viewer':
                    Viewer = Options[key]
                elif key == 'Animation':
                    Animation = Options[key]
                elif key == 'TimeHistory':
                    TimeHistory = Options[key]
                elif key == 'PushOver':
                    PushOver = Options[key]
                elif key == 'OpenSeesCommand':
                    OpenSeesCommand = Options[key]
                elif key == 'T1':
                    T1 = Options[key]
                elif key == 'T2':
                    T2 = Options[key]
                elif key == 'EnhancedOutput':
                    EnhancedOutput = Options[key]
                elif key == 'PDeltaColumn':
                    PDeltaColumn = Options[key]
                elif key == 'CyclicStatic':
                    CyclicStatic = Options[key]
                elif key == 'DriftHistory':
                    DriftHistory = Options[key]
                elif key == 'POModalForces':
                    POModalForces = Options[key]
                elif key == 'heff':
                    heff = Options[key]
                elif key == 'ConstantWallAxialForce':
                    ConstantWallAxialForce = Options[key]
                elif key == 'ApplyPDeltaLoad':
                    ApplyPDeltaLoad = Options[key]
                elif key == 'DebugMode':
                    DebugMode = Options[key]
                elif key == 'RegularizeSteel':
                    RegularizeSteel = Options[key]
                elif key == 'MaxPORoofDrift':
                    MaxPORoofDrift = Options[key]
                elif key == 'PlotSections':
                    PlotSections = Options[key]
                elif key == 'RegularizeFracture':
                    RegularizeFracture = Options[key]
                elif key == 'NoOfIterations':
                    NoOfIterations = Options[key]
                elif key == 'Tolerance':
                    Tolerance = Options[key]
                elif key == 'CrushingStrengthRatio':
                    CrushingStrengthRatio = Options[key]
                elif key == 'POELFForces':
                    POELFForces = Options[key]
                elif key == 'CuTa':
                    CuTa = Options[key]
                elif key == 'GfccOGfc':
                    GfccOGfc = Options[key]
                elif key == 'ConcreteMaterialModel':
                    ConcreteMaterialModel = Options[key]
                elif key == 'SteelMaterialModel':
                    SteelMaterialModel = Options[key]
                elif key == 'SteelUltimateStrainTension':
                    SteelUltimateStrainTension = Options[key]
                elif key == 'SteelPostYieldStiffness':
                    SteelPostYieldStiffness = Options[key]
                elif key == 'WallAxialLoadMultiplier':
                    WallAxialLoadMultiplier = Options[key]
                elif key == 'UseForceBasedElements':
                    UseForceBasedElements = Options[key]
                elif key == 'DivisionsPerStory':
                    DivisionsPerStory = Options[key]
                elif key == 'GfcOfpc':
                    GfcOfpc = Options[key]
                elif key == 'NoOfIntPoints':
                    NoOfIntPoints = Options[key]
                elif key == 'UnconfinedBeta':
                    UnconfinedBeta = Options[key]
                elif key == 'Regularized':
                    Regularized = Options[key]
                elif key == 'WallThicknessMultipler':
                    WallThicknessMultipler = Options[key]
                elif key == 'FBE_Tolerance':
                    FBE_Tolerance = Options[key]
                elif key == 'ModalDamping':
                    ModalDamping = Options[key]
                elif key == 'Zeta':
                    Zeta = Options[key]
                elif key == 'ConfinementModel':
                    ConfinementModel = Options[key]
                elif key == 'TrackPeriod':
                    TrackPeriod = Options[key]
                elif key == 'HHTTransientIntegrator':
                    HHTTransientIntegrator = Options[key]
                elif key == 'OutputTag':
                    OutputTag = Options[key]
                elif key == 'TrackDrifts':
                    TrackDrifts = Options[key]

    for key, value in kwargs.items():
        OtherOptions[key] = value

    #region ########################## Initializing ##########################
    import time
    import uuid
    randomnumber = str(uuid.uuid4()).replace('-', '').upper()
    timestamp = ''
    if not(DebugMode):
        timestamp = time.strftime("%y%m%d-%H%M%S-") + randomnumber
    ModelName = 'PWA'
    FileName = '%s.tcl' % (ModelName)

    Subfolder = '/PWA' + timestamp

    import platform
    if platform.system() == 'Windows':
        TCLFileDirectory = os.getcwd() + '/tcl%s/' % Subfolder
        ResultDirectory = 'Results/'
    else:
        TCLFileDirectory = '/tmp/%s/' % Subfolder# '/dev/shm%s/' % Subfolder
        ResultDirectory = 'Results/'

    if not os.path.exists(TCLFileDirectory):  # Make Directory is unavailable
        os.makedirs(TCLFileDirectory)
    if not os.path.exists(TCLFileDirectory + ResultDirectory):  # Make Directory is unavailable
        os.makedirs(TCLFileDirectory + ResultDirectory)

    OData = OpenSeesAPI.Database.Collector(OpenSeesCommand, TCLFileDirectory, FileName)

    # Save OtherOption in OData so that I can be referenced if needed
    OData._OtherOptions = OtherOptions

    #endregion

    #region ########################## Setup and Source Definition ##########################
    OData.AddObject(OpenSeesAPI.TCL.CodeTitle('Initialization'))
    OData.AddObject(OpenSeesAPI.Model.BasicBuilder(2, 3))  # Start 3d Model
    OData.AddObject(OpenSeesAPI.Output.LogFile(OData.Executable.LogFileName))  # Start Log File

    # endregion

    #region ########################## Define Building Geometry, Nodes and Constraints ##########################
    OData.AddObject(OpenSeesAPI.TCL.CodeTitle('Geometry Setup'))

    # Create Grid Nodes
    XGrids = [0,400]
    YGrids = Archetype.YGrids

    # Define Nodes
    CoreNodes = []
    SplitNodes = []
    PDeltaColumnNode = []
    NoOfDivisionsPerFloor = int(DivisionsPerStory)
    for i in range(0, len(YGrids)):
        if i !=0:
            if Archetype.CustomSection is not None:
                if type(Archetype.CustomSection[i-1]) is list: #check if split column
                    l_wall = Archetype.CustomSection[i-1][0].l_w + Archetype.CustomSection[i-1][1].l_w + Archetype.CouplingBeamLength[i-1]
                else:
                    l_wall = Archetype.CustomSection[i-1].l_w
            else:
                l_wall = Archetype.l_w[i-1]
        else:
            if Archetype.CustomSection is not None:
                if type(Archetype.CustomSection[0]) is list: #check if split column
                    l_wall = Archetype.CustomSection[0][0].l_w + Archetype.CustomSection[0][1].l_w + Archetype.CouplingBeamLength[0]
                else:
                    l_wall = Archetype.CustomSection[0].l_w
            else:
                l_wall = Archetype.l_w[0]

        if i != 0 and NoOfDivisionsPerFloor > 1:
            div = np.linspace(0,1,NoOfDivisionsPerFloor+1) * (YGrids[i] - YGrids[i-1])
            for y in div[1:-1]:
                CoreNodes.append(OData.CreateNode(XGrids[0], YGrids[i-1]+y, NodeType = 2,
                                                  GridX = 0, GridY = i, GroupId = i))

                if Archetype.CustomSection is not None:
                    if type(Archetype.CustomSection[i-1]) is list:#Archetype.CouplingBeams[i-1] is not None or Archetype.CouplingBeams[i-2] is not None:
                        SplitNodes.append(
                            [OData.CreateNode(XGrids[0] - l_wall / 2. + Archetype.CustomSection[i-1][0].l_w / 2.0,
                                              YGrids[i - 1] + y, NodeType=2, GridX=0, GridY=i, GroupId=i),
                             OData.CreateNode(XGrids[0] + l_wall / 2. - Archetype.CustomSection[i-1][1].l_w / 2.0,
                                              YGrids[i - 1] + y, NodeType=2, GridX=0, GridY=i, GroupId=i),
                            ])
                    else:
                        SplitNodes.append([])
                else:
                    SplitNodes.append([])

        CoreNodes.append(OData.CreateNode(XGrids[0], YGrids[i], GridX=0, GridY=i, GroupId=i))
        PDeltaColumnNode.append(OData.CreateNode(XGrids[1], YGrids[i], GridX=1, GridY=i, GroupId=i))

        if i < len(YGrids) - 1: # if not last story
            if Archetype.CustomSection[i] is not None:
                if type(Archetype.CustomSection[i]) is list:
                    left_be = Archetype.CustomSection[i][0].l_w
                    right_be = Archetype.CustomSection[i][1].l_w
                    SplitNodes.append(
                        [OData.CreateNode(XGrids[0] - l_wall / 2. + left_be / 2., YGrids[i],
                                          GridX=0, GridY=i, GroupId=i),
                         OData.CreateNode(XGrids[0] + l_wall / 2. - right_be / 2., YGrids[i],
                                          GridX=0, GridY=i, GroupId=i),
                         ])
                elif type(Archetype.CustomSection[i-1]) is list:
                    left_be = Archetype.CustomSection[i-1][0].l_w
                    right_be = Archetype.CustomSection[i-1][0].l_w
                    SplitNodes.append(
                        [OData.CreateNode(XGrids[0] - l_wall / 2. + left_be / 2., YGrids[i],
                                          GridX=0, GridY=i, GroupId=i),
                         OData.CreateNode(XGrids[0] + l_wall / 2. - right_be / 2., YGrids[i],
                                          GridX=0, GridY=i, GroupId=i),
                         ])
                else:
                    SplitNodes.append([])
            else:
                SplitNodes.append([])
        elif type(Archetype.CustomSection[i-1]) is list: #if last story
                    left_be = Archetype.CustomSection[i-1][0].l_w
                    right_be = Archetype.CustomSection[i-1][0].l_w
                    SplitNodes.append(
                        [OData.CreateNode(XGrids[0] - l_wall / 2. + left_be / 2., YGrids[i], GridX=0, GridY=i, GroupId=i),
                         OData.CreateNode(XGrids[0] + l_wall / 2. - right_be / 2., YGrids[i], GridX=0, GridY=i, GroupId=i),
                         ])
        else: # No list detected, therefore no split column
            SplitNodes.append([])

    #endregion

    #region ########################## Define Geometric Transformations ##########################
    OData.AddObject(OpenSeesAPI.TCL.CodeTitle('Define Geometric Transformations'))

    PDelta = OData.AddObject(
        OpenSeesAPI.Model.Element.GeomTransf.PDelta(1))

    #endregion

    ##############################################################################
    ### All OpenSEES Objects are adding directly to the Database Beyond This Point
    ##############################################################################

    #region ########################## Define Materials and Sections ##########################
    OData.AddObject(OpenSeesAPI.TCL.CodeTitle('Define Materials and Sections'))

    ElasticRigid = OData.AddObject(
        OpenSeesAPI.Material.UniaxialMaterial.Elastic(
            OData.GetFreeMaterialId(4, 0), 1e16, Notes='This Rigid Material'))

    # Define Core Wall Sections
    CoreWallSections = []

    # Define Core Wall Sections
    for i in range(1, len(YGrids)):
        for j in range(0, NoOfDivisionsPerFloor):
            ind = (i-1) * NoOfDivisionsPerFloor + j + 1

            # Compute Core Wall Fiber Sections
            max_mesh_Size = 6 #Inches
            # NoOfIntPoints = 5

            height = CoreNodes[ind].Y - CoreNodes[ind-1].Y

            cover = 3.0

            bar_size = 10.17

            if Archetype.CustomSection is not None:
                #Check if Split Column Below
                if i >= 2: # Check if second story or more
                    #Check to see if below is a split column and above is not.
                    if type(Archetype.CustomSection[i-1]) is not list and type(Archetype.CustomSection[i-2]) is list:
                        if isinstance(Archetype.CustomSection[i - 2][0], ATCWallObjects.PlanarWallSection): # If planar wall below then combine, else just use section specified
                            # Create a custom sections with the two pier sections joined.
                            Core = ATCWallArchetypeHelpers.ComputeJoinedCustomPlanarWallFiberSection\
                                (OData, Archetype.CustomSection[i - 2][0], Archetype.CustomSection[i - 2][1],
                                 Archetype.CouplingBeamLength[i-2],
                                 cover, height,
                                 NoOfIntPoints,
                                 max_mesh_Size,
                                 Elastic=False,
                                 RegularizeSteel=RegularizeSteel,
                                 RegularizeFracture=RegularizeFracture,
                                 CrushingStrengthRatio=CrushingStrengthRatio,
                                 GfccOGfc=GfccOGfc,
                                 ConcreteMaterialModel=ConcreteMaterialModel,
                                 ConfinementModel=ConfinementModel,
                                 SteelMaterialModel=SteelMaterialModel,
                                 SteelUltimateStrainTension=SteelUltimateStrainTension,
                                 SteelPostYieldStiffness=SteelPostYieldStiffness,
                                 GfcOfpc=GfcOfpc,
                                 UnconfinedBeta=UnconfinedBeta,
                                 Regularized=Regularized,
                                 WallThicknessMultipler = WallThicknessMultipler,
                                 UseForceBased=UseForceBasedElements
                                 )
                            CoreWallSections.append([Core])
                            continue
                        else: # Use section specified instead of joining them (this is the case for I sections)
                            Core = ATCWallArchetypeHelpers.ComputeCustomIWallFiberSection(OData, Archetype.CustomSection[i-1],
                                                                                           cover, height,
                                                                                           NoOfIntPoints, max_mesh_Size,
                                                                                           Elastic=False,
                                                                                           RegularizeSteel=RegularizeSteel,
                                                                                           RegularizeFracture=RegularizeFracture,
                                                                                           CrushingStrengthRatio=CrushingStrengthRatio,
                                                                                           GfccOGfc=GfccOGfc,
                                                                                           ConcreteMaterialModel=ConcreteMaterialModel,
                                                                                           ConfinementModel = ConfinementModel,
                                                                                           SteelMaterialModel=SteelMaterialModel,
                                                                                           SteelUltimateStrainTension=SteelUltimateStrainTension,
                                                                                           SteelPostYieldStiffness=SteelPostYieldStiffness,
                                                                                           GfcOfpc=GfcOfpc,
                                                                                           UnconfinedBeta=UnconfinedBeta,
                                                                                           Regularized=Regularized,
                                                                                           UseForceBased=UseForceBasedElements

                                                                                           )

                            CoreWallSections.append([Core])
                            continue
                #if not split column then continue as usual
                if type(Archetype.CustomSection[i-1]) is not list :
                    if isinstance(Archetype.CustomSection[i-1], ATCWallObjects.PlanarWallSection):  # check if wall is a planar wall section
                        Core = ATCWallArchetypeHelpers.ComputeCustomPlanarWallFiberSection(OData, Archetype.CustomSection[i-1],
                                                                                           cover, height,
                                                                                           NoOfIntPoints, max_mesh_Size,
                                                                                           Elastic=False,
                                                                                           RegularizeSteel=RegularizeSteel,
                                                                                           RegularizeFracture=RegularizeFracture,
                                                                                           CrushingStrengthRatio=CrushingStrengthRatio,
                                                                                           GfccOGfc=GfccOGfc,
                                                                                           ConcreteMaterialModel=ConcreteMaterialModel,
                                                                                           ConfinementModel = ConfinementModel,
                                                                                           SteelMaterialModel=SteelMaterialModel,
                                                                                           SteelUltimateStrainTension=SteelUltimateStrainTension,
                                                                                           SteelPostYieldStiffness=SteelPostYieldStiffness,
                                                                                           GfcOfpc=GfcOfpc,
                                                                                           UseForceBased=UseForceBasedElements,
                                                                                           UnconfinedBeta=UnconfinedBeta, Regularized=Regularized,
                                                                                           WallThicknessMultipler=WallThicknessMultipler,)
                        CoreWallSections.append([Core])

                    elif isinstance(Archetype.CustomSection[i-1], ATCWallObjects.IWallSection): # check if wall is a I wall section
                        Core = ATCWallArchetypeHelpers.ComputeCustomIWallFiberSection(OData, Archetype.CustomSection[i-1],
                                                                                       cover, height,
                                                                                       NoOfIntPoints, max_mesh_Size,
                                                                                       Elastic=False,
                                                                                       RegularizeSteel=RegularizeSteel,
                                                                                       RegularizeFracture=RegularizeFracture,
                                                                                       CrushingStrengthRatio=CrushingStrengthRatio,
                                                                                       GfccOGfc=GfccOGfc,
                                                                                       ConcreteMaterialModel=ConcreteMaterialModel,
                                                                                       ConfinementModel = ConfinementModel,
                                                                                       SteelMaterialModel=SteelMaterialModel,
                                                                                       SteelUltimateStrainTension=SteelUltimateStrainTension,
                                                                                       SteelPostYieldStiffness=SteelPostYieldStiffness,
                                                                                       GfcOfpc=GfcOfpc,
                                                                                           UseForceBased=UseForceBasedElements,
                                                                                      UnconfinedBeta=UnconfinedBeta, Regularized=Regularized,)
                        CoreWallSections.append([Core])
                else: # Split Column Case
                    if isinstance(Archetype.CustomSection[i-1][0], ATCWallObjects.PlanarWallSection): # check if wall is a planar section
                        CoreA = ATCWallArchetypeHelpers.ComputeCustomPlanarWallFiberSection(OData, Archetype.CustomSection[i-1][0],
                                                                                           cover, height,
                                                                                           NoOfIntPoints, max_mesh_Size,
                                                                                           Elastic=False,
                                                                                           RegularizeSteel=RegularizeSteel,
                                                                                           RegularizeFracture=RegularizeFracture,
                                                                                           CrushingStrengthRatio=CrushingStrengthRatio,
                                                                                            GfccOGfc=GfccOGfc,
                                                                                            ConcreteMaterialModel=ConcreteMaterialModel,
                                                                                            ConfinementModel = ConfinementModel,
                                                                                            SteelMaterialModel=SteelMaterialModel,
                                                                                            SteelUltimateStrainTension=SteelUltimateStrainTension,
                                                                                            SteelPostYieldStiffness=SteelPostYieldStiffness,
                                                                                            GfcOfpc=GfcOfpc,
                                                                                           UseForceBased=UseForceBasedElements,
                                                                                            UnconfinedBeta=UnconfinedBeta,Regularized=Regularized,
                                                                                            WallThicknessMultipler=WallThicknessMultipler,)

                        CoreB = ATCWallArchetypeHelpers.ComputeCustomPlanarWallFiberSection(OData, Archetype.CustomSection[i-1][1],
                                                                                           cover, height,
                                                                                           NoOfIntPoints, max_mesh_Size,
                                                                                           Elastic=False,
                                                                                           RegularizeSteel=RegularizeSteel,
                                                                                           RegularizeFracture=RegularizeFracture,
                                                                                           CrushingStrengthRatio=CrushingStrengthRatio,
                                                                                            GfccOGfc=GfccOGfc,
                                                                                            ConcreteMaterialModel=ConcreteMaterialModel,
                                                                                            ConfinementModel = ConfinementModel,
                                                                                            SteelMaterialModel=SteelMaterialModel,
                                                                                            SteelUltimateStrainTension=SteelUltimateStrainTension,
                                                                                            SteelPostYieldStiffness=SteelPostYieldStiffness,
                                                                                            GfcOfpc=GfcOfpc,
                                                                                           UseForceBased=UseForceBasedElements,
                                                                                            UnconfinedBeta=UnconfinedBeta,Regularized=Regularized,
                                                                                            WallThicknessMultipler=WallThicknessMultipler,)

                    elif isinstance(Archetype.CustomSection[i-1][0], ATCWallObjects.TWallSection): # check if wall is a twall section
                        CoreA = ATCWallArchetypeHelpers.ComputeCustomTWallFiberSection(OData,
                                                                                    Archetype.CustomSection[i - 1][0],
                                                                                    cover, height,
                                                                                    NoOfIntPoints,
                                                                                    max_mesh_Size,
                                                                                    Elastic=False,
                                                                                    RegularizeSteel=RegularizeSteel,
                                                                                    RegularizeFracture=RegularizeFracture,
                                                                                    CrushingStrengthRatio=CrushingStrengthRatio,
                                                                                       GfccOGfc=GfccOGfc,
                                                                                       ConcreteMaterialModel=ConcreteMaterialModel,
                                                                                       ConfinementModel = ConfinementModel,
                                                                                       SteelMaterialModel=SteelMaterialModel,
                                                                                       SteelUltimateStrainTension=SteelUltimateStrainTension,
                                                                                       SteelPostYieldStiffness=SteelPostYieldStiffness,
                                                                                       GfcOfpc=GfcOfpc,
                                                                                           UseForceBased=UseForceBasedElements,
                                                                                       UnconfinedBeta=UnconfinedBeta,Regularized=Regularized,
                                                                                       )

                        CoreB = ATCWallArchetypeHelpers.ComputeCustomTWallFiberSection(OData,
                                                                                    Archetype.CustomSection[i - 1][1],
                                                                                    cover, height,
                                                                                    NoOfIntPoints,
                                                                                    max_mesh_Size,
                                                                                    Elastic=False,
                                                                                    RegularizeSteel=RegularizeSteel,
                                                                                    RegularizeFracture=RegularizeFracture,
                                                                                    CrushingStrengthRatio=CrushingStrengthRatio,
                                                                                       GfccOGfc=GfccOGfc,
                                                                                       ConcreteMaterialModel=ConcreteMaterialModel,
                                                                                       ConfinementModel = ConfinementModel,
                                                                                       SteelMaterialModel=SteelMaterialModel,
                                                                                       SteelUltimateStrainTension=SteelUltimateStrainTension,
                                                                                       SteelPostYieldStiffness=SteelPostYieldStiffness,
                                                                                       GfcOfpc=GfcOfpc,
                                                                                           UseForceBased=UseForceBasedElements,
                                                                                       UnconfinedBeta=UnconfinedBeta,Regularized=Regularized,)
                    CoreWallSections.append([CoreA, CoreB])
            else: # use non-custom section type defined in archetype class, this should not be used in the ATC Wall Project
                l_w = Archetype.l_w[i - 1]

                if Archetype.CouplingBeamLength is not None:
                    b_f = (Archetype.b_f[i - 1] - Archetype.CouplingBeamLength[i - 1]) / 2.
                else:
                    b_f = Archetype.b_f[i - 1]

                if Archetype.BoundaryElement is None:
                    Core = ATCWallArchetypeHelpers.ComputeIShapedWallFiberSection(OData, l_w, b_f, Archetype.t[i-1], cover,
                                                      Archetype.fce, Archetype.fye, Archetype.fu,
                                                      bar_size, Archetype.rho[i-1], Archetype.rho_t[i-1],
                                                      height, NoOfIntPoints, max_mesh_Size,
                                                      Elastic=False)
                    CoreWallSections.append([Core])
                else:
                    if Archetype.BoundaryElement[i - 1] == None:
                        bar_size = 4.
                    else:
                        bar_size = 10.17
                    Core = ATCWallArchetypeHelpers.ComputeIShapedWallFiberSection(OData, l_w, b_f, Archetype.t[i - 1],
                                                                                  cover,
                                                                                  Archetype.fce,
                                                                                  Archetype.fye,
                                                                                  Archetype.fu,
                                                                                  bar_size, Archetype.rho[i - 1],
                                                                                  Archetype.rho_t[i - 1],
                                                                                  height, NoOfIntPoints,
                                                                                  max_mesh_Size,
                                                                                  Elastic=False,
                                                                                  boundaryelement=\
                                                                                      Archetype.BoundaryElement[i - 1],
                                                                                  rho_boundaryelement=\
                                                                                    Archetype.RhoBoundaryElement[i-1],
                                                                                  CrushingStrengthRatio=CrushingStrengthRatio,
                                                                                  GfccOGfc=GfccOGfc,
                                                                                  ConcreteMaterialModel=ConcreteMaterialModel,
                                                                                  ConfinementModel = ConfinementModel,
                                                                                  SteelMaterialModel=SteelMaterialModel,
                                                                                  SteelUltimateStrainTension=SteelUltimateStrainTension,
                                                                                  SteelPostYieldStiffness=SteelPostYieldStiffness,
                                                                                  GfcOfpc=GfcOfpc,
                                                                                  UseForceBased=UseForceBasedElements,
                                                                                  UnconfinedBeta=UnconfinedBeta, Regularized=Regularized,
                                                                                  )
                    CoreWallSections.append([Core])

    # Plot Sections
    if PlotSections:
        for i in range(1, len(YGrids)):
            import OSFiberSectionViewer as FiberViewer
            if len(CoreWallSections[i-1]) < 2:
                FiberViewer.ShowFiberSection(CoreWallSections[(i-1) * NoOfDivisionsPerFloor][0][0]._Section._fibers,
                                             'Figures/CoreSection_L%02d-%s.png'%(i,Archetype.Name))
            else:
                FiberViewer.ShowFiberSection(CoreWallSections[(i-1) * NoOfDivisionsPerFloor][0][0]._Section._fibers,
                                             'Figures/CoreSection_L%02d-%s-A.png' % (i, Archetype.Name))
                FiberViewer.ShowFiberSection(CoreWallSections[(i-1) * NoOfDivisionsPerFloor][1][0]._Section._fibers,
                                             'Figures/CoreSection_L%02d-%s-B.png' % (i, Archetype.Name))

    # endregion

    #region ########################## Define Rotational Springs for Plastic Hinge ##########################
    OData.AddObject(OpenSeesAPI.TCL.CodeTitle('Define Rotational Springs for Plastic Hinge'))

    # endregion

    #region ########################## Define Elements ##########################
    OData.AddObject(OpenSeesAPI.TCL.CodeTitle('Define Elements'))

    # Add Core Wall Line Elements
    CoreWallElements = []
    WallShearSprings = []
    for j in range(len(CoreNodes)-1):
        corewallelements = []
        # Check if single or wall with opening
        if len(CoreWallSections[j]) >= 2 : # double wall section
            # Double Wall
            wallsection = CoreWallSections[j][0]
            Sections = ''
            for a in range(NoOfIntPoints):
                Sections += ' %d' % wallsection[a].id

            # Using Force Based Elements
            if UseForceBasedElements:
                FiberElementA = OData.AddElement(
                                OpenSeesAPI.Element.ForceBeamColumnOriginal(OData.GetFreeElementId(5, j+1),
                                                                            SplitNodes[j][0],
                                                                            SplitNodes[j+1][0],
                                                                            NoOfIntPoints,
                                                                            wallsection, #'-sections %s' % Sections,
                                                                            PDelta,
                                                                            Optional='-iter %d %e' % (
                                                                            NoOfIterations, FBE_Tolerance)))
            else:
                bottomNode = SplitNodes[j][0]
                ShearNode = OData.CreateNode(bottomNode.X, bottomNode.Y, NodeType=3)
                FiberElementA = OData.AddElement(
                                OpenSeesAPI.Element.DispBeamColumn(OData.GetFreeElementId(5, j+1),
                                                                    ShearNode,
                                                                    SplitNodes[j+1][0],
                                                                    NoOfIntPoints,
                                                                    wallsection,
                                                                    PDelta, Optional='-integration  Lobatto'
                                                                    ))
                h = SplitNodes[j+1][0].Y - bottomNode.Y
                ShearStiffness = wallsection[0]._MatList[0]._E / h
                ShearSpringMaterial = OData.AddMaterial(OpenSeesAPI.Material.UniaxialMaterial.Elastic(OData.GetFreeMaterialId(2,j),ShearStiffness))
                ShearSpringA = OData.AddElement(OpenSeesAPI.Element.ZeroLength(OData.GetFreeElementId(9, j+1),
                                                                              bottomNode,
                                                                              ShearNode,
                                                                              [ShearSpringMaterial, ElasticRigid, ElasticRigid],
                                                                              [1,2,3]))
            wallsection = CoreWallSections[j][1]
            Sections = ''
            for a in range(NoOfIntPoints):
                Sections += ' %d' % wallsection[a].id

            if UseForceBasedElements:
                FiberElementB = OData.AddElement(
                                OpenSeesAPI.Element.ForceBeamColumnOriginal(OData.GetFreeElementId(5, j+1),
                                                                            SplitNodes[j][1],
                                                                            SplitNodes[j+1][1],
                                                                            NoOfIntPoints,
                                                                            wallsection, #'-sections %s' % Sections,
                                                                            PDelta,
                                                                            Optional='-iter %d %e' % (
                                                                            NoOfIterations, FBE_Tolerance)))
            else:
                # FiberElementB = OData.AddElement(
                #                 OpenSeesAPI.Element.DispBeamColumn(OData.GetFreeElementId(5, j+1),
                #                                                             SplitNodes[j][1],
                #                                                             SplitNodes[j+1][1],
                #                                                             NoOfIntPoints,
                #                                                             wallsection,
                #                                                             PDelta, Optional='-integration  Lobatto'
                #                                                             ))
                bottomNode = SplitNodes[j][1]
                ShearNode = OData.CreateNode(bottomNode.X, bottomNode.Y, NodeType=3)
                FiberElementB = OData.AddElement(
                                OpenSeesAPI.Element.DispBeamColumn(OData.GetFreeElementId(5, j+1),
                                                                    ShearNode,
                                                                    SplitNodes[j+1][1],
                                                                    NoOfIntPoints,
                                                                    wallsection,
                                                                    PDelta, Optional='-integration  Lobatto'
                                                                    ))
                h = SplitNodes[j+1][1].Y - bottomNode.Y
                ShearStiffness = wallsection[0]._MatList[0]._E / h
                ShearSpringMaterial = OData.AddMaterial(OpenSeesAPI.Material.UniaxialMaterial.Elastic(OData.GetFreeMaterialId(2,j),ShearStiffness))
                ShearSpringB = OData.AddElement(OpenSeesAPI.Element.ZeroLength(OData.GetFreeElementId(9, j+1),
                                                                              bottomNode,
                                                                              ShearNode,
                                                                              [ShearSpringMaterial, ElasticRigid, ElasticRigid],
                                                                              [1,2,3]))
                WallShearSprings.append([ShearSpringA, ShearSpringB])
            FiberElement = [FiberElementA, FiberElementB]

        else: #Single wall section
            # Single Core
            wallsection = CoreWallSections[j][0]
            Sections = ''
            for a in range(NoOfIntPoints):
                Sections += ' %d' % wallsection[a].id

            # Using Force Based Elements
            if UseForceBasedElements:
                FiberElement = [OData.AddElement(
                                OpenSeesAPI.Element.ForceBeamColumnOriginal(OData.GetFreeElementId(5, j+1),
                                                                            CoreNodes[j],
                                                                            CoreNodes[j+1],
                                                                            NoOfIntPoints,
                                                                            wallsection, #'-sections %s' % Sections,
                                                                            PDelta,
                                                                            Optional='-iter %d %e' % (
                                                                                NoOfIterations, FBE_Tolerance)))]
            else:
                # FiberElement = [OData.AddElement(
                #                 OpenSeesAPI.Element.DispBeamColumn(OData.GetFreeElementId(5, j+1),
                #                                                             CoreNodes[j],
                #                                                             CoreNodes[j+1],
                #                                                             NoOfIntPoints,
                #                                                             wallsection,
                #                                                             PDelta, Optional='-integration  Lobatto'
                #                                                             ))]
                bottomNode = CoreNodes[j]
                ShearNode = OData.CreateNode(bottomNode.X, bottomNode.Y, NodeType=3)
                FiberElement = [OData.AddElement(
                                OpenSeesAPI.Element.DispBeamColumn(OData.GetFreeElementId(5, j+1),
                                                                    ShearNode,
                                                                    CoreNodes[j + 1],
                                                                    NoOfIntPoints,
                                                                    wallsection,
                                                                    PDelta, Optional='-integration  Lobatto'
                                                                    ))]
                h = CoreNodes[j+1].Y - bottomNode.Y
                ShearStiffness = wallsection[0]._MatList[0]._E / h
                ShearSpringMaterial = OData.AddMaterial(OpenSeesAPI.Material.UniaxialMaterial.Elastic(OData.GetFreeMaterialId(2,j),ShearStiffness))
                ShearSpring = OData.AddElement(OpenSeesAPI.Element.ZeroLength(OData.GetFreeElementId(9, j+1),
                                                                              bottomNode,
                                                                              ShearNode,
                                                                              [ShearSpringMaterial, ElasticRigid, ElasticRigid],
                                                                              [1,2,3]))

                WallShearSprings.append([ShearSpring])
        CoreWallElements.append(FiberElement)

    # Add Coupling Beams and Rigid Beams
    CouplingBeamDivisions = int(np.ceil(DivisionsPerStory / 2.))

    CouplingBeamElements = []
    AllCouplingBeamSections = []
    for i in range(1, len(YGrids)):
        if Archetype.CouplingBeamLength is not None:
            # if last story and wall below is double
            if i == len(YGrids)-1:
                if len(CoreWallSections[(i - 1) * NoOfDivisionsPerFloor]) == 2:
                    # Left Core
                    # Add Rigid Beams
                    NodeI = SplitNodes[i * NoOfDivisionsPerFloor][0]
                    LCNodeJ = OData.CreateNode(NodeI.X + Archetype.CustomSection[i - 1][0].l_w / 2., NodeI.Y, NodeType=2,
                                               GroupId=i)
                    OData.AddElement(OpenSeesAPI.Element.ElasticBeamColumn(OData.GetFreeMaterialId(8, i),
                                                                           NodeI, LCNodeJ, 1e6, 1e6, 1e6,
                                                                           PDelta,
                                                                           _Notes='Left Rigid Beam to Connect Core to Coupling Beam'))

                    # Right Core
                    # Add Rigid Beam
                    NodeI = SplitNodes[i * NoOfDivisionsPerFloor][1]
                    RCNodeJ = OData.CreateNode(NodeI.X - Archetype.CustomSection[i - 1][1].l_w / 2., NodeI.Y, NodeType=2,
                                               GroupId=i)
                    OData.AddElement(OpenSeesAPI.Element.ElasticBeamColumn(OData.GetFreeMaterialId(8, i),
                                                                           RCNodeJ, NodeI, 1e6, 1e6, 1e6,
                                                                           PDelta,
                                                                           _Notes='Right Rigid Beam to Connect Core to Coupling Beam'))

                    # Add Intermediate Nodes for Coupling Beam Meshing
                    CouplingBeamNodes = [LCNodeJ]
                    for k in range(1, CouplingBeamDivisions):
                        dX = (RCNodeJ.X - LCNodeJ.X) / CouplingBeamDivisions
                        CouplingBeamNodes.append(OData.CreateNode(LCNodeJ.X + dX * k, LCNodeJ.Y, NodeType=2,
                                                                  GroupId=i))
                    CouplingBeamNodes.append(RCNodeJ)

                    # Add Fiber Section Beams
                    couplingbeam = Archetype.CouplingBeams[i - 1]
                    Cover = 3.0
                    widthdiagonal = 12
                    heightdiagonal = 8

                    CouplingBeamSections = ATCWallArchetypeHelpers. \
                        DiagonalCouplingBeamSection(OData, couplingbeam.b * WallThicknessMultipler,
                                                    couplingbeam.h,
                                                    Cover,
                                                    couplingbeam.NoOfBarsX,
                                                    couplingbeam.NoOfBarsY,
                                                    widthdiagonal,
                                                    heightdiagonal,
                                                    couplingbeam.BarDia,
                                                    couplingbeam.fce,
                                                    couplingbeam.fye,
                                                    Archetype.fu,
                                                    Archetype.CouplingBeamLength[i - 1],
                                                    tie_spacing=couplingbeam.TieSpacing,
                                                    s_bar_x=(couplingbeam.b - cover) / float(couplingbeam.no_of_ties_x),
                                                    s_bar_y=(couplingbeam.h - cover) / float(couplingbeam.no_of_ties_y),
                                                    tie_size=couplingbeam.TieDia,
                                                    meshsize=3,
                                                    NoOfIntPoints=NoOfIntPoints,
                                                    Elastic=False,
                                                    RegularizeSteel=RegularizeSteel,
                                                    RegularizeFracture=RegularizeFracture,
                                                    CrushingStrengthRatio=CrushingStrengthRatio,
                                                    GfccOGfc=GfccOGfc,
                                                    ConcreteMaterialModel=ConcreteMaterialModel,
                                                    ConfinementModel = ConfinementModel,
                                                    SteelMaterialModel=SteelMaterialModel,
                                                    SteelUltimateStrainTension=SteelUltimateStrainTension,
                                                    SteelPostYieldStiffness=SteelPostYieldStiffness,
                                                    GfcOfpc=GfcOfpc,
                                                    NoOfDivisions=CouplingBeamDivisions,
                                                    UseForceBased=UseForceBasedElements,
                                                    UnconfinedBeta=UnconfinedBeta, Regularized=Regularized,
                                                    )

                    AllCouplingBeamSections.append(CouplingBeamSections)

                    for b in range(CouplingBeamDivisions):
                        Sections = ''
                        for a in range(NoOfIntPoints):
                            Sections += ' %d' % CouplingBeamSections[b][a].id

                        if UseForceBasedElements:
                            FiberElement = OData.AddElement(
                                OpenSeesAPI.Element.ForceBeamColumnOriginal(OData.GetFreeElementId(6, j + 1),
                                                                            CouplingBeamNodes[b],#LCNodeJ,
                                                                            CouplingBeamNodes[b + 1],#RCNodeJ,
                                                                            NoOfIntPoints,
                                                                            CouplingBeamSections[b], #'-sections %s' % Sections,
                                                                            PDelta,
                                                                            Optional='-iter %d %e' % (NoOfIterations, FBE_Tolerance)))
                        else:
                            leftNode = CouplingBeamNodes[b]
                            rightNode = CouplingBeamNodes[b + 1]
                            ShearNode = OData.CreateNode(leftNode.X, leftNode.Y, NodeType=3)
                            FiberElement = OData.AddElement(
                                OpenSeesAPI.Element.DispBeamColumn(OData.GetFreeElementId(5, j + 1),
                                                                   ShearNode,
                                                                   rightNode,
                                                                   NoOfIntPoints,
                                                                   CouplingBeamSections[b],
                                                                   PDelta, Optional='-integration  Lobatto'
                                                                   ))
                            l = rightNode.X - leftNode.X
                            ShearStiffness = CouplingBeamSections[b][0]._MatList[0]._E / l
                            ShearSpringMaterial = OData.AddMaterial(
                                OpenSeesAPI.Material.UniaxialMaterial.Elastic(OData.GetFreeMaterialId(2, j),
                                                                              ShearStiffness))
                            ShearSpring = OData.AddElement(
                                OpenSeesAPI.Element.ZeroLength(OData.GetFreeElementId(9, j + 1),
                                                               leftNode,
                                                               ShearNode,
                                                               [ElasticRigid, ShearSpringMaterial, ElasticRigid],
                                                               [1, 2, 3]))

                        CouplingBeamElements.append(FiberElement)

            # double wall section below and single above or double wall above and single wall below
            elif len(CoreWallSections[(i - 1) * NoOfDivisionsPerFloor]) == 2 \
                    and len(CoreWallSections[i * NoOfDivisionsPerFloor]) == 1 \
                    or len(CoreWallSections[i * NoOfDivisionsPerFloor]) == 2\
                            and len(CoreWallSections[(i-1) * NoOfDivisionsPerFloor]) == 1:
                # Left Core
                # Add Rigid Beams
                NodeI = SplitNodes[i * NoOfDivisionsPerFloor][0]
                LCNodeJ = CoreNodes[i * NoOfDivisionsPerFloor]
                OData.AddElement(OpenSeesAPI.Element.ElasticBeamColumn(OData.GetFreeMaterialId(8, i),
                                                                       NodeI, LCNodeJ, 1e6, 1e6, 1e6,
                                                                       PDelta,
                                                                       _Notes='Left Rigid Beam to Connect Double to Single Core'))

                # Right Core
                # Add Rigid Beam
                NodeI = SplitNodes[i * NoOfDivisionsPerFloor][1]
                RCNodeJ = CoreNodes[i * NoOfDivisionsPerFloor]
                OData.AddElement(OpenSeesAPI.Element.ElasticBeamColumn(OData.GetFreeMaterialId(8, i),
                                                                       RCNodeJ, NodeI, 1e6, 1e6, 1e6,
                                                                       PDelta,
                                                                       _Notes='Right Rigid Beam to Connect Double to Single Core'))

                AllCouplingBeamSections.append([])

            # check if double wall section below and above.
            elif len(CoreWallSections[(i - 1) * NoOfDivisionsPerFloor]) == 2 \
                    and len(CoreWallSections[i * NoOfDivisionsPerFloor]) == 2:
                # Left Core
                # Add Rigid Beams
                NodeI = SplitNodes[i * NoOfDivisionsPerFloor][0]
                LCNodeJ = OData.CreateNode(NodeI.X + Archetype.CustomSection[i-1][0].l_w/2., NodeI.Y, NodeType=2, GroupId=i)
                OData.AddElement(OpenSeesAPI.Element.ElasticBeamColumn(OData.GetFreeMaterialId(8, i),
                                                                       NodeI, LCNodeJ, 1e6, 1e6, 1e6,
                                                                       PDelta,
                                                                       _Notes='Left Rigid Beam to Connect Core to Coupling Beam'))

                # Right Core
                # Add Rigid Beam
                NodeI = SplitNodes[i * NoOfDivisionsPerFloor][1]
                RCNodeJ = OData.CreateNode(NodeI.X - Archetype.CustomSection[i-1][1].l_w/2., NodeI.Y, NodeType=2, GroupId=i)
                OData.AddElement(OpenSeesAPI.Element.ElasticBeamColumn(OData.GetFreeMaterialId(8, i),
                                                                       RCNodeJ, NodeI, 1e6, 1e6, 1e6,
                                                                       PDelta,
                                                                       _Notes='Right Rigid Beam to Connect Core to Coupling Beam'))

                # Add Intermediate Nodes for Coupling Beam Meshing
                CouplingBeamNodes = [LCNodeJ]
                for k in range(1, CouplingBeamDivisions):
                    dX = (RCNodeJ.X - LCNodeJ.X) / CouplingBeamDivisions
                    CouplingBeamNodes.append(OData.CreateNode(LCNodeJ.X + dX * k, LCNodeJ.Y, NodeType=2,
                                                              GroupId=i))
                CouplingBeamNodes.append(RCNodeJ)

                # Add Fiber Section Beams
                couplingbeam = Archetype.CouplingBeams[i - 1]
                Cover = 3.0
                widthdiagonal = 12
                heightdiagonal = 8

                CouplingBeamSections = ATCWallArchetypeHelpers.\
                    DiagonalCouplingBeamSection(OData, couplingbeam.b * WallThicknessMultipler,
                                               couplingbeam.h,
                                               Cover,
                                               couplingbeam.NoOfBarsX,
                                               couplingbeam.NoOfBarsY,
                                               widthdiagonal,
                                               heightdiagonal,
                                               couplingbeam.BarDia,
                                               couplingbeam.fce,
                                               couplingbeam.fye,
                                               Archetype.fu,
                                               Archetype.CouplingBeamLength[i-1],
                                               tie_spacing=couplingbeam.TieSpacing,
                                               s_bar_x=(couplingbeam.b-cover)/float(couplingbeam.no_of_ties_x),
                                               s_bar_y=(couplingbeam.h-cover)/float(couplingbeam.no_of_ties_y),
                                               tie_size=couplingbeam.TieDia,
                                               meshsize=3,
                                               NoOfIntPoints=NoOfIntPoints,
                                               Elastic=False,
                                               RegularizeSteel=RegularizeSteel,
                                               RegularizeFracture=RegularizeFracture,
                                                CrushingStrengthRatio=CrushingStrengthRatio,
                                                GfccOGfc=GfccOGfc,
                                                ConcreteMaterialModel=ConcreteMaterialModel,
                                                ConfinementModel = ConfinementModel,
                                                SteelMaterialModel=SteelMaterialModel,
                                                SteelUltimateStrainTension=SteelUltimateStrainTension,
                                                SteelPostYieldStiffness=SteelPostYieldStiffness,
                                                GfcOfpc=GfcOfpc,
                                                NoOfDivisions=CouplingBeamDivisions,
                                                UseForceBased=UseForceBasedElements,
                                                UnconfinedBeta=UnconfinedBeta, Regularized=Regularized,)

                AllCouplingBeamSections.append(CouplingBeamSections)

                for b in range(CouplingBeamDivisions):
                    Sections = ''
                    for a in range(NoOfIntPoints):
                        Sections += ' %d' % CouplingBeamSections[b][a].id

                    if UseForceBasedElements:
                        FiberElement = OData.AddElement(
                            OpenSeesAPI.Element.ForceBeamColumnOriginal(OData.GetFreeElementId(6, j + 1),
                                                                        CouplingBeamNodes[b],#LCNodeJ,
                                                                        CouplingBeamNodes[b+1],#RCNodeJ,
                                                                        NoOfIntPoints,
                                                                        CouplingBeamSections[b], #'-sections %s' % Sections,
                                                                        PDelta,
                                                                        Optional='-iter %d %e'%(NoOfIterations, FBE_Tolerance)))
                    else:
                        # FiberElement = OData.AddElement(
                        #     OpenSeesAPI.Element.DispBeamColumn(OData.GetFreeElementId(6, j + 1),
                        #                                                 CouplingBeamNodes[b],#LCNodeJ,
                        #                                                 CouplingBeamNodes[b+1],#RCNodeJ,
                        #                                                 NoOfIntPoints,
                        #                                                 CouplingBeamSections[b],
                        #                                                 PDelta, Optional='-integration Lobatto'
                        #                                                 ))
                        leftNode = CouplingBeamNodes[b]
                        rightNode = CouplingBeamNodes[b + 1]
                        ShearNode = OData.CreateNode(leftNode.X, leftNode.Y, NodeType=3)
                        FiberElement = OData.AddElement(
                            OpenSeesAPI.Element.DispBeamColumn(OData.GetFreeElementId(5, j + 1),
                                                               ShearNode,
                                                               rightNode,
                                                               NoOfIntPoints,
                                                               CouplingBeamSections[b],
                                                               PDelta, Optional='-integration  Lobatto'
                                                               ))
                        l = rightNode.X - leftNode.X
                        ShearStiffness = CouplingBeamSections[b][0]._MatList[0]._E / l
                        ShearSpringMaterial = OData.AddMaterial(
                            OpenSeesAPI.Material.UniaxialMaterial.Elastic(OData.GetFreeMaterialId(2, j),
                                                                          ShearStiffness))
                        ShearSpring = OData.AddElement(
                            OpenSeesAPI.Element.ZeroLength(OData.GetFreeElementId(9, j + 1),
                                                           leftNode,
                                                           ShearNode,
                                                           [ElasticRigid, ShearSpringMaterial, ElasticRigid],
                                                           [1, 2, 3]))

                    CouplingBeamElements.append(FiberElement)
            else:
                AllCouplingBeamSections.append([])

    if PlotSections and len(AllCouplingBeamSections) > 0:
        for i in range(1, len(YGrids)):
            if len(AllCouplingBeamSections[i-1]) > 0:
                for k in range(CouplingBeamDivisions):
                    for j in range(NoOfIntPoints):
                        FiberViewer.ShowFiberSection(AllCouplingBeamSections[i-1][k][j]._Section._fibers,
                                                     'Figures/CouplingBeams_L%02d-D%d-I%d_%s.png' % (i, k, j, Archetype.Name))


    # Add Gravity Columns
    if PDeltaColumn:
        Ec = 0.35*57.*(Archetype.fce*1000)**.5
        A = 24*24*14
        I = 24.**4./12
        for i in range(1,len(YGrids)):
            OData.AddElement(OpenSeesAPI.Element.ElasticBeamColumn(OData.GetFreeElementId(9, i),
                                                                   PDeltaColumnNode[i-1], PDeltaColumnNode[i],
                                                                   A*100., Ec, I/1.e5, PDelta,
                                                                   _Notes='PDelta Column at Level %d' % i))

    if hasattr(Archetype, 'BasementProperties'):
        ElasticDia = OData.AddMaterial(OpenSeesAPI.Material.UniaxialMaterial.Elastic(OData.GetFreeMaterialId(9,1), 1.))
        BasementWallStiffness = Archetype.BasementProperties.WallStiffnesses
        BasementFloorStiffness = Archetype.BasementProperties.FloorStiffnesses
        BasementMass = Archetype.BasementProperties.BasementMass
        NoOfBasementLevels = len(BasementMass)

        CoreNodesToConnectTo = []
        for i in range(0,len(YGrids)):
            if len(CoreWallSections[(i - 1) * NoOfDivisionsPerFloor]) == 2: # If coupling beam
                node = SplitNodes[i * NoOfDivisionsPerFloor][1]
            else: # if no coupling beam exists
                node = OData.GetNodesByGrid(0,i)[0]

                CoreNodesToConnectTo.append(node)

        # Create Nodes
        BasementNodes = []
        DiaLength = 120.
        for i in range(0,len(BasementWallStiffness)+1):
            BasementNodes.append(OData.CreateNode(CoreNodesToConnectTo[i].X - DiaLength, CoreNodesToConnectTo[i].Y, NodeType=2))
            if i == 0:
                BasementSupportNode = BasementNodes[0]
                OData.AddConstraint(OpenSeesAPI.Model.Constraint.Fix(BasementNodes[0], [1, 1, 1]))

        # Add Basement Mass
        g = 386.4
        for i in range(1, len(BasementWallStiffness) + 1):
            OData.AddConstraint(OpenSeesAPI.Node.Mass(BasementNodes[i],
                                                      [BasementMass[i - 1] / g, 1.e-6, 1.e-6]))

        # Add Connect Basement Elements
        for i in range(1, len(BasementWallStiffness) + 1):
            OData.AddElement(OpenSeesAPI.Model.Element.Element.ElasticTimoshenkoBeam(
                OData.GetFreeElementId(8, 1), BasementNodes[i - 1], BasementNodes[i],  1e16, 1., 1e16, 1e16, BasementWallStiffness[i-1], PDelta))
            # OData.AddElement(OpenSeesAPI.Model.Element.Element.ElasticBeamColumn(
            #     OData.GetFreeElementId(8, 1), BasementNodes[i], BasementNodes[i - 1], 1e3, 1e3, 1e3, GeoTransfLinear))
            OData.AddElement(OpenSeesAPI.Model.Element.Element.Truss(OData.GetFreeElementId(8, 1), BasementNodes[i], CoreNodesToConnectTo[i], BasementFloorStiffness[i-1] * DiaLength, ElasticDia))
    else:
        NoOfBasementLevels = 0

    # endregion

    #region ########################## Define Restraints/Constraints ##########################
    #Find All Used Nodes and Set them as used
    OData.AssignUsedNodes()

    #Defining Fixity
    SupportZeroLengthElements = []
    SupportNodes = [] # Stores all the nodes with fix supports
    CoreSupportNodes = []
    elementsupportnodes = OData.GetNodesByYCoordinate(0, 1)

    for node in elementsupportnodes:
        if hasattr(node,'Used'):
            if node.Used != True:
                continue
        else:
            continue

        supportnode = OData.CreateNode(node.X,node.Y,NodeType=2,GroupId=0,_Notes='Used to Extract Reactions')
        supportnode.__setattr__('Used',True)
        # if node.X == XGrids[0]:
        SupportZeroLengthElements.append(OData.AddConstraint(
            OpenSeesAPI.Element.ZeroLength(
                OData.GetFreeElementId(9,1),supportnode,node,
                [ElasticRigid, ElasticRigid, ElasticRigid],
                [1,2,3])))
        # else:
        #     SupportZeroLengthElements.append(OData.AddConstraint(
        #         OpenSeesAPI.Element.ZeroLength(
        #             OData.GetFreeElementId(9,1),supportnode,node,
        #             [ElasticPlastic, ElasticRigid, ElasticRigid],
        #             [1,2,3])))

        if node.X == XGrids[-1]:
            OData.AddConstraint(OpenSeesAPI.Model.Constraint.Fix(supportnode, [1, 1, 0]))
        else:
            OData.AddConstraint(OpenSeesAPI.Model.Constraint.Fix(supportnode,[1,1,1]))
            CoreSupportNodes.append(supportnode)
        SupportNodes.append(supportnode)

    # Add Basement Support Node
    if hasattr(Archetype, 'BasementProperties'):
        SupportNodes.append(BasementSupportNode)

    # Create Mass Node
    DiaNode = []
    # DiaNode = PDeltaColumnNode
    for i in range(1, len(YGrids)):
        g = 386.4
        if len(CoreWallSections[(i - 1) * NoOfDivisionsPerFloor]) == 2:  # If coupling beam
            node = SplitNodes[i * NoOfDivisionsPerFloor]
            MassNode = OData.AddConstraint(OpenSeesAPI.Node.Mass(node[0],
                                                             [Archetype.Mass[i-1]/g/2., 1.e-6, 1.e-6]))
            MassNode = OData.AddConstraint(OpenSeesAPI.Node.Mass(node[1],
                                                             [Archetype.Mass[i-1]/g/2., 1.e-6, 1.e-6]))
            DiaNode.append(node[1])
        else:  # if no coupling beam exists
            node = OData.GetNodesByGrid(0, i)[0]
            MassNode = OData.AddConstraint(OpenSeesAPI.Node.Mass(node,
                                                             [Archetype.Mass[i-1]/g, 1.e-6, 1.e-6]))
            DiaNode.append(node)

    # Set Dia Nodes are Used
    list(map(lambda x: x.__setattr__('Used',True),DiaNode))

    #Define Rigid Diaphragm
    if PDeltaColumn:
        RigidMat = OData.AddMaterial(OpenSeesAPI.Material.UniaxialMaterial.Elastic(OData.GetFreeMaterialId(1, 1), 1e8))
        for i in range(1,len(YGrids)):
            if len(CoreWallSections[(i - 1) * NoOfDivisionsPerFloor]) == 2: # If coupling beam
                node = SplitNodes[i * NoOfDivisionsPerFloor][1]
            else: # if no coupling beam exists
                node = OData.GetNodesByGrid(0,i)[0]
            OData.AddElement(OpenSeesAPI.Element.Truss(OData.GetFreeElementId(9,i), node, PDeltaColumnNode[i], 1e8, RigidMat))

    # endregion

    ##############################################################################
    ### Start Writing Elements to the Executible File
    ##############################################################################

    #region ########################## Write Model to TCL File ##########################
    OData.WriteModel()

    # endregion

    #region ########################## Eigenvalue Analysis ##########################
    OData.AddObject(OpenSeesAPI.TCL.CodeTitle('Eigenvalue Analysis'))

    # NoOfModes = min((len(YGrids) - 2 - NoOfBasementLevels),12) # min(len(YGrids)-3,3)#len(YGrids)-1 #<= this was before
    NoOfModes = min((len(YGrids) - 1), 40)
    if not POELFForces:
        if PushOver == True or (T1 == None and T2 == None) or CyclicStatic or ModalDamping:
                OData.AddObject(OpenSeesAPI.Analysis.Eigen(NoOfModes, symmBandLapack=False))
                for mode in range(1, NoOfModes+1):
                    for i in range(1, len(YGrids)):
                        OData.AddObject(OpenSeesAPI.TCL.TCLScript(
                            'set EigenVector%d%d_X [nodeEigenvector %d %d %d]' % (mode, i, DiaNode[i - 1].id, mode, 1)))

                        OData.AddObject(OpenSeesAPI.TCL.TCLScript(
                            'puts \" EigenVector Mode:%d Story:%d $EigenVector%d%d_X\"' % (mode, i, mode, i)))

    #endregion

    #region ########################## Rayleigh Damping ##########################
    OData.AddObject(OpenSeesAPI.TCL.CodeTitle('Rayleigh Damping'))

    # Use TBI Zeta
    if UseTBIZeta:
        H = (YGrids[-1] - YGrids[int(NoOfBasementLevels)]) / 12.
        Zeta = min(0.05,max(0.36/np.sqrt(H), 0.025))

    # Adding Rayleigh Damping to the Mass Matrix Only
    if Zeta != 0.:
        if not(ModalDamping):
            if PushOver == True or (T1 == None and T2 == None):
                # Adding Rayleigh Damping to the Mass Matrix Only
                OData.AddObject(OpenSeesAPI.TCL.TCLScript('set zeta %f'%Zeta))
                OData.AddObject(OpenSeesAPI.TCL.TCLScript('set alpha0 [expr $zeta*$w1*$w2/($w1+$w2)]'))
                OData.AddObject(OpenSeesAPI.TCL.TCLScript('set alpha1 [expr $zeta*2.0/($w1+$w2)]'))
                OData.AddObject(OpenSeesAPI.TCL.TCLScript('rayleigh $alpha0 0 $alpha1 0'))
            else:
                # Adding Rayleigh Damping to the Mass Matrix Only
                OData.AddObject(OpenSeesAPI.TCL.TCLScript('set zeta %f' % Zeta))
                OData.AddObject(OpenSeesAPI.TCL.TCLScript('set w1 [expr 2*3.141592654/%f]'%T1))
                OData.AddObject(OpenSeesAPI.TCL.TCLScript('set w2 [expr 2*3.141592654/%f]'%T2))
                OData.AddObject(OpenSeesAPI.TCL.TCLScript('set alpha0 [expr $zeta*$w1*$w2/($w1+$w2)]'))
                OData.AddObject(OpenSeesAPI.TCL.TCLScript('set alpha1 [expr $zeta*2.0/($w1+$w2)]'))
                OData.AddObject(OpenSeesAPI.TCL.TCLScript('rayleigh $alpha0 0 $alpha1 0'))
        elif not IncludeSupplementalRayleigh:
            OData.AddObject(OpenSeesAPI.TCL.TCLScript('modalDamping %.5f'%Zeta))
        elif IncludeSupplementalRayleigh:
            # Adding Modal Damping Plus Rayleigh Damping
            OData.AddObject(OpenSeesAPI.TCL.TCLScript('set zeta %f' % Zeta))
            ModalDampingCommand = 'modalDamping '
            for i in range(NoOfModes):
                OData.AddObject(OpenSeesAPI.TCL.TCLScript('set zeta_%d [expr $zeta - $zeta*$w%d/($w%d)]' %(i + 1, i + 1, NoOfModes)))
                ModalDampingCommand += '$zeta_%d '%(i+1)
                OData.AddObject(OpenSeesAPI.TCL.TCLScript('puts "Zeta for Mode: %d = $zeta_%d"'%(i+1, i+1)))

            OData.AddObject(OpenSeesAPI.TCL.TCLScript(ModalDampingCommand))
            OData.AddObject(OpenSeesAPI.TCL.TCLScript('set alpha1 [expr 2.*$zeta/($w%d)]' % NoOfModes))

            import itertools
            eletags = "".join(['%s '%x.id for x in list(itertools.chain.from_iterable(CoreWallElements))])
            OData.AddObject(OpenSeesAPI.TCL.TCLScript('region 999999 -ele %s -rayleigh 0 0 $alpha1 0'%eletags))

    # endregion

    #region ########################## Loads ##########################
    OData.AddObject(OpenSeesAPI.TCL.CodeTitle('Loads'))

    AddGravityLoad = True

    # Make Sure Wall Gravity Load is Negative
    Archetype.WallGravityLoad = -1 * np.abs(Archetype.WallGravityLoad)

    # Add Gravity Loads
    # To Core Wall Elements
    Loads = []
    # for i in range(1,len(YGrids)):
    if AddGravityLoad:
        if not(ConstantWallAxialForce):
            for i in range(1, len(YGrids)):
                if len(CoreWallSections[(i - 1) * NoOfDivisionsPerFloor]) == 2:  # If coupling beam
                    node = SplitNodes[i * NoOfDivisionsPerFloor]
                    Loads.append(OpenSeesAPI.Model.Pattern.Load(node[0], [0, Archetype.WallGravityLoad[i - 1] / 2. * WallAxialLoadMultiplier, 0]))
                    Loads.append(OpenSeesAPI.Model.Pattern.Load(node[1], [0, Archetype.WallGravityLoad[i - 1] / 2. * WallAxialLoadMultiplier, 0]))
                else:
                    node = OData.GetNodesByGrid(0, i)[0]
                    Loads.append(OpenSeesAPI.Model.Pattern.Load(node,[0, Archetype.WallGravityLoad[i-1] * WallAxialLoadMultiplier, 0 ]))
        else:
            for i in range(1, len(YGrids)):
                if YGrids[i-1] < heff and YGrids[i] > heff: # Apply Constant Load at Top of Heff
                    if len(CoreWallSections[(i - 1) * NoOfDivisionsPerFloor]) == 2:  # If coupling beam
                        node = SplitNodes[i * NoOfDivisionsPerFloor]
                        Loads.append(OpenSeesAPI.Model.Pattern.Load(node[0], [0, np.sum(Archetype.WallGravityLoad) / 2. * WallAxialLoadMultiplier, 0]))
                        Loads.append(OpenSeesAPI.Model.Pattern.Load(node[1], [0, np.sum(Archetype.WallGravityLoad) / 2. * WallAxialLoadMultiplier, 0]))
                    else:
                        node = OData.GetNodesByGrid(0, i-1)[0]
                        Loads.append(OpenSeesAPI.Model.Pattern.Load(node,[0, np.sum(Archetype.WallGravityLoad) * WallAxialLoadMultiplier, 0 ]))

    # To PDelta Column
    if PDeltaColumn:
        if ApplyPDeltaLoad:
            for i in range(1, len(YGrids)):
                Loads.append(
                    OpenSeesAPI.Model.Pattern.Load(PDeltaColumnNode[i],
                                                   [0, -1*np.abs(Archetype.GravityLoad[i-1]), 0]))

    OData.AddObject(OpenSeesAPI.Model.Pattern.Plain(100,'Linear',Loads))

    # endregion

    #region ########################## Time Series ##########################
    OData.AddObject(OpenSeesAPI.TCL.CodeTitle('Time Series'))

    #Adding Time Series from GMData parameter
    TimeSeries = OpenSeesAPI.Model.TimeSeries.Path(1, Dt, GMData)
    OData.AddObject(TimeSeries)

    # endregion

    #region ########################## Recorders ##########################
    OData.AddObject(OpenSeesAPI.TCL.CodeTitle('Recorder Setup'))

    OutputFolder = ResultDirectory

    if PushOver:
        Optional = ''
    elif TimeHistory:
        Optional = '-dT %f'%Dt
    else:
        Optional = ''

    Displacement_File_Name = '%s-NodeD-%s.dat' % (ModelName, timestamp)
    OData.AddObject(
        OpenSeesAPI.Output.Recorder.Node(OutputFolder + '/' + Displacement_File_Name, [DiaNode[-1]], [1,2,3], 'disp', Optional))

    AllStoriesDisp = '%s-AllNodeD-%s.dat' % (ModelName, timestamp)
    OData.AddObject(
        OpenSeesAPI.Output.Recorder.Node(OutputFolder + '/' + AllStoriesDisp, DiaNode, [1], 'disp', Optional))

    AllStoriesAcceleration = '%s-AllNodeAcceleration-%s.dat' % (ModelName, timestamp)
    OData.AddObject(
        OpenSeesAPI.Output.Recorder.Node(OutputFolder + '/' + AllStoriesAcceleration, DiaNode, [1], 'accel', '-timeSeries %d'%TimeSeries.id))

    CoreDisp = '%s-CoreD-%s.dat' % (ModelName, timestamp)
    OData.AddObject(
        OpenSeesAPI.Output.Recorder.Node(OutputFolder + '/' + CoreDisp, CoreNodes, [1, 2, 3], 'disp', Optional))

    Reaction_File_Name = '%s-NodeBaseShear-%s.dat' % (ModelName, timestamp)
    OData.AddObject(
        OpenSeesAPI.Output.Recorder.Node(OutputFolder + '/' + Reaction_File_Name, SupportNodes, [1], 'reaction', Optional))

    FullReaction_File_Name = '%s-FullNodeReact-%s.dat' % (ModelName, timestamp)
    OData.AddObject(
        OpenSeesAPI.Output.Recorder.Node(OutputFolder + '/' + FullReaction_File_Name, SupportNodes, [1,2,3], 'reaction', Optional))

    CoreSupportNodesReact = '%s-CoreSupportNodesReact-%s.dat' % (ModelName, timestamp)
    OData.AddObject(
        OpenSeesAPI.Output.Recorder.Node(OutputFolder + '/' + CoreSupportNodesReact, CoreSupportNodes, [1,2,3], 'reaction', Optional))

    CoreMoment = '%s-CoreMoment-%s.dat' % (ModelName, timestamp)
    OData.AddObject(
        OpenSeesAPI.Output.Recorder.Node(OutputFolder + '/' + CoreMoment, SupportNodes, [3], 'reaction', Optional))

    AllNodeDispl = '%s-AllNodeDispl-%s.dat' % (ModelName, timestamp)
    AllUsedNodes = [x for x in OData._Nodes if hasattr(x, 'Used')]
    OData.AddObject(
        OpenSeesAPI.Output.Recorder.Node(OutputFolder + '/' + AllNodeDispl, AllUsedNodes, [1, 2, 3], 'disp', Optional))

    StoryDrift = '%s-StoryDrift-%s.dat' % (ModelName, timestamp)
    iNodes = '%s'%(SupportNodes[0].id)
    jNodes = ''
    Count = -1
    for node in DiaNode:
        Count += 1
        if Count < len(DiaNode) - 1:
            iNodes += ' %d' % node.id
            jNodes += ' %d' % node.id
        else:
            jNodes += ' %d' % node.id
    OData.AddObject(OpenSeesAPI.TCL.TCLScript('recorder Drift -file %s -precision 3 -time -iNode %s -jNode %s -dof 1 -perpDirn 2'%(OutputFolder + '/' + StoryDrift, iNodes, jNodes)))

    # Find Out Location of Extreme Fibers
    NoOfSamplePoints = 2
    if Archetype.CustomSection is None:
        SamplePoints = np.linspace(1, -1, NoOfSamplePoints) * Archetype.l_w[
            0] / 2.  # this wont be correct for the top stories
    else:
        if type(Archetype.CustomSection[-1]) is not list:
            SamplePoints = np.linspace(1, -1, NoOfSamplePoints) * \
                           Archetype.CustomSection[-1].l_w / 2.  # this wont be correct for the bottom
        else:
            Length = Archetype.CustomSection[-1][0].l_w # Assumes both Ends have the same section length
            SamplePoints = np.linspace(1, -1, NoOfSamplePoints) * \
                           Length / 2.  # this wont be correct for the bottom

    # Outside Extreme Fiber Strains
    BaseExtremeFiberStrains1 = '%s-BaseExtremeFiberStrains-1-%s.dat' % (ModelName, timestamp)
    BaseExtremeFiberStrains2 = '%s-BaseExtremeFiberStrains-2-%s.dat' % (ModelName, timestamp)
    if len(CoreWallElements[0]) == 1:
        BaseWallElements = CoreWallElements[0]
    else:
        BaseWallElements = CoreWallElements[0]
    OData.AddObject(OpenSeesAPI.TCL.TCLScript(
        'recorder Element -file %s -time %s -ele %s section %d fiber %f %f stressStrain ' % (
            OutputFolder + '/' + BaseExtremeFiberStrains1, Optional,
            ''.join([' %d' % x.id for x in BaseWallElements]),
            1, SamplePoints[0], 0.0)))
    OData.AddObject(OpenSeesAPI.TCL.TCLScript(
        'recorder Element -file %s -time %s -ele %s section %d fiber %f %f stressStrain ' % (
            OutputFolder + '/' + BaseExtremeFiberStrains2, Optional,
            ''.join([' %d' % x.id for x in BaseWallElements]),
            1, SamplePoints[-1], 0.0)))

    def JoinList(x):
        import itertools
        return list(itertools.chain.from_iterable(x))

    def JoinListFilterForPier(x, Single = True):
        new = []
        for i in range(len(x)):
            if Single:
                if type(x[i]) is not list:
                    new.append(x[i])
                elif len(x[i]) == 1:
                    new.append(x[i][0])
            else:
                if type(x[i]) is not list:
                    continue
                elif len(x[i]) > 1:
                    new.extend(x[i])

        return new

    # Extract Bottom Most Core Elements
    CoreWallElementsAtStoryBottom = []
    for i in range(0, len(CoreWallElements), int(DivisionsPerStory)):
        CoreWallElementsAtStoryBottom.append(CoreWallElements[i])

    # Save Shear and Moment
    CoreShearAndMomentAtStoryBottom = '%s-CoreShearMoment-%s.dat' % (ModelName, timestamp)
    if len(CoreWallElementsAtStoryBottom[0]) == 1: # Uncoupled Direction
        OData.AddObject(OpenSeesAPI.TCL.TCLScript(
            'recorder Element -file %s -time %s -ele %s localForce' % (
                OutputFolder + '/' + CoreShearAndMomentAtStoryBottom, Optional,
                ''.join([' %d' % x.id for x in JoinList(CoreWallElementsAtStoryBottom)]))))
    else: # Coupled Direction
        pass

    # Save Strains
    CoreStrainAtStoryBottom1 = '%s-CoreShearAndMomentAtStoryBottom-1-%s.dat' % (ModelName, timestamp)
    CoreStrainAtStoryBottom2 = '%s-CoreShearAndMomentAtStoryBottom-2-%s.dat' % (ModelName, timestamp)
    if len(CoreWallElements[0]) == 1:
        OData.AddObject(OpenSeesAPI.TCL.TCLScript(
            'recorder Element -file %s -time %s -ele %s section %d fiber %f %f stressStrain ' % (
                OutputFolder + '/' + CoreStrainAtStoryBottom1, Optional,
                ''.join([' %d' % x.id for x in JoinList(CoreWallElementsAtStoryBottom)]),
                1, SamplePoints[0], 0.0)))
        OData.AddObject(OpenSeesAPI.TCL.TCLScript(
            'recorder Element -file %s -time %s -ele %s section %d fiber %f %f stressStrain ' % (
                OutputFolder + '/' + CoreStrainAtStoryBottom2, Optional,
                ''.join([' %d' % x.id for x in JoinList(CoreWallElementsAtStoryBottom)]),
                1, SamplePoints[-1], 0.0)))
    else:
        pass

    # Save Coupling Beam Rotation
    if len(CoreWallElementsAtStoryBottom[0]) > 1:  # Coupled Direction
        pass

    OutputFileNames = [Displacement_File_Name, AllStoriesDisp, CoreDisp, Reaction_File_Name,
                       CoreMoment, AllNodeDispl, FullReaction_File_Name, CoreSupportNodesReact,
                       BaseExtremeFiberStrains1, BaseExtremeFiberStrains2]

    # Add Coupling Beam OutputFiles
    if Archetype.CouplingBeams is not None and Archetype.CouplingBeams[0] is not None:
        CouplingExtremeFiberStrainsTopLeft = '%s-CouplingExtremeFiberStrainsTopLeft-%s.dat' % (ModelName, timestamp)
        CouplingExtremeFiberStrainsBottomLeft = '%s-CouplingExtremeFiberStrainsBottomLeft-%s.dat' % (ModelName, timestamp)
        CouplingExtremeFiberStrainsTopRight = '%s-CouplingExtremeFiberStrainsTopRight-%s.dat' % (ModelName, timestamp)
        CouplingExtremeFiberStrainsBottomRight = '%s-CouplingExtremeFiberStrainsBottomRight-%s.dat' % (ModelName, timestamp)

        OutputFileNames.extend([CouplingExtremeFiberStrainsTopLeft, CouplingExtremeFiberStrainsBottomLeft,
                                CouplingExtremeFiberStrainsTopRight, CouplingExtremeFiberStrainsBottomRight])

        LeftBeams = []
        RightBeams = []
        for i in range(0,len(CouplingBeamElements),CouplingBeamDivisions):
            LeftBeams.append(CouplingBeamElements[i])
        for i in range(CouplingBeamDivisions-1,len(CouplingBeamElements),CouplingBeamDivisions):
            RightBeams.append(CouplingBeamElements[i])

        BeamXLoc = Archetype.CouplingBeams[0].b / 2
        BeamYLoc = Archetype.CouplingBeams[0].h
        OData.AddObject(OpenSeesAPI.TCL.TCLScript(
            'recorder Element -file %s -time %s -ele %s section %d fiber %f %f stressStrain ' % (
                OutputFolder + '/' + CouplingExtremeFiberStrainsTopLeft, Optional,
                ''.join([' %d' % x.id for x in LeftBeams]),
                1, BeamYLoc, BeamXLoc )))

        OData.AddObject(OpenSeesAPI.TCL.TCLScript(
            'recorder Element -file %s -time %s -ele %s section %d fiber %f %f stressStrain ' % (
                OutputFolder + '/' + CouplingExtremeFiberStrainsBottomLeft, Optional,
                ''.join([' %d' % x.id for x in LeftBeams]),
                1, 0, BeamXLoc)))

        OData.AddObject(OpenSeesAPI.TCL.TCLScript(
            'recorder Element -file %s -time %s -ele %s section %d fiber %f %f stressStrain ' % (
                OutputFolder + '/' + CouplingExtremeFiberStrainsTopRight, Optional,
                ''.join([' %d' % x.id for x in RightBeams]),
                NoOfIntPoints, BeamYLoc, BeamXLoc )))

        OData.AddObject(OpenSeesAPI.TCL.TCLScript(
            'recorder Element -file %s -time %s -ele %s section %d fiber %f %f stressStrain ' % (
                OutputFolder + '/' + CouplingExtremeFiberStrainsBottomRight, Optional,
                ''.join([' %d' % x.id for x in RightBeams]),
                NoOfIntPoints, 0, BeamXLoc)))

    if not UseForceBasedElements:
        ShearSpringShears = '%s-ShearSpringShears-%s.dat' % (ModelName, timestamp)
        ShearSprings = JoinList(WallShearSprings)
        OData.AddObject(
            OpenSeesAPI.Output.Recorder.Element(OutputFolder + '/' + ShearSpringShears, ShearSprings, 'force', Optional = Optional))

        OutputFileNames.append(ShearSpringShears)

    if EnhancedOutput:
        ### Extract Axial Load From Core Wall Elements
        for j in range(NoOfIntPoints):
            SingleCoreElements = JoinListFilterForPier(CoreWallElements)
            CoreAxialLoad = '%s-CoreAxialLoad-%s-%d.dat'%(ModelName, timestamp, j)
            OData.AddObject(OpenSeesAPI.TCL.TCLScript('recorder Element -file %s -time %s -ele %s section %d force'%(OutputFolder + '/' + CoreAxialLoad, Optional, ''.join([' %d'%x.id for x in SingleCoreElements]), j+1)))

        # Single Pier Elements
        NoOfSamplePoints = 20
        SingleCoreElements = JoinListFilterForPier(CoreWallElements)
        if Archetype.CustomSection is None:
            SamplePoints = np.linspace(1,-1,NoOfSamplePoints)*Archetype.l_w[0]/2. #this wont be correct for the top stories
        else:
            if type(Archetype.CustomSection[-1]) is not list:
                SamplePoints = np.linspace(1, -1, NoOfSamplePoints) * \
                               Archetype.CustomSection[-1].l_w / 2.  # this wont be correct for the bottom
            else:
                Length = Archetype.CustomSection[-1][0].l_w + Archetype.CustomSection[-1][1].l_w  + Archetype.CouplingBeamLength[-1]
                SamplePoints = np.linspace(1, -1, NoOfSamplePoints) * \
                               Length / 2.# this wont be correct for the bottom

        for j in range(NoOfIntPoints):
            for k in range(NoOfSamplePoints):
                CoreAxialStressStrain = '%s-CoreAxialStressStrain-%s-%d-%d.dat'%(ModelName, timestamp, j, k)
                OData.AddObject(OpenSeesAPI.TCL.TCLScript('recorder Element -file %s -time %s -ele %s section %d fiber %f %f stressStrain'%(OutputFolder + '/' + CoreAxialStressStrain, Optional, ''.join([' %d'%x.id for x in SingleCoreElements]), j+1, SamplePoints[k], 0.0)))

        # Double Pier Elements
        NoOfSamplePointsDouble = int(NoOfSamplePoints/2 - 1)
        DoubleCoreElements = JoinListFilterForPier(CoreWallElements, False)
        if len(CoreWallElementsAtStoryBottom[0]) == 1: #Single Wall
            SamplePointsDouble = np.linspace(1, -1, NoOfSamplePointsDouble) * \
                                 JoinListFilterForPier(Archetype.CustomSection)[0].l_w / 2.  # this wont be correct for the bottom
        else:
            SamplePointsDouble = np.linspace(1, -1, NoOfSamplePointsDouble) * \
                                 JoinListFilterForPier(Archetype.CustomSection, False)[0].l_w / 2.  # this wont be correct for the bottom

        for j in range(NoOfIntPoints):
            for k in range(int(NoOfSamplePointsDouble)):
                CoreAxialStressStrain = '%s-CoreAxialStressStrain-Double-%s-%d-%d.dat'%(ModelName, timestamp, j, k)
                OData.AddObject(OpenSeesAPI.TCL.TCLScript('recorder Element -file %s -time %s -ele %s section %d fiber %f %f stressStrain '%(OutputFolder + '/' + CoreAxialStressStrain, Optional, ''.join([' %d'%x.id for x in DoubleCoreElements]), j+1, SamplePointsDouble[k], 0.0)))

        # Extrain Strains and Stresses in Core Wall At Base
        # BaseMaterials = CoreWallSections[0][0][0]._Section._fibers[:4]
        # for j in range(len(BaseMaterials)):
        #     for k in range(NoOfSamplePoints):
        #         BaseStressStrain = '%s-BaseStressStrain-%s-%d-%d.dat' % (ModelName, timestamp, k, j)
        #         OData.AddObject(OpenSeesAPI.TCL.TCLScript(
        #         'recorder Element -file %s -time -ele %s section %d fiber %f %f %d stressStrain' % (
        #         OutputFolder + '/' + BaseStressStrain, CoreWallElements[0].id, 1,
        #         SamplePoints[k], 0.0, BaseMaterials[j]._Mat.id)))

        # Extract Forces From Core Wall Elements
        CoreForces = '%s-CoreWallForces-%s.dat' % (ModelName, timestamp)
        OData.AddObject(
            OpenSeesAPI.Output.Recorder.Element(OutputFolder + '/' + CoreForces, JoinList(CoreWallElements), 'globalForce', Optional=Optional))

        OutputFileNames.append(CoreForces)

        # Extract Node Displacements to all core wall and coupling beam elements
        AllFiberSectionDeformationFiles = []
        for i in range(NoOfIntPoints):
            AllFiberSectionDeformation = '%s-FiberSectionDeformation-Section-%d-%s.dat' % (ModelName, i, timestamp)
            AllFiberSectionDeformationFiles.append(AllFiberSectionDeformation)
            OData.AddObject(OpenSeesAPI.TCL.TCLScript(
                'recorder Element -file %s -time %s -ele %s section %s forceAndDeformation ' % (
                OutputFolder + '/' + AllFiberSectionDeformation, Optional,
                ''.join([' %d' % x.id for x in JoinList(CoreWallElements)]),
                '%d'%(i+1))))

        # Keeping Track of all Output Files
        for j in range(NoOfIntPoints):
            CoreAxialLoad = '%s-CoreAxialLoad-%s-%d.dat' % (ModelName, timestamp, j)
            AllFiberSectionDeformation = '%s-FiberSectionDeformation-Section-%d-%s.dat' % (ModelName, i, timestamp)
            OutputFileNames.append(AllFiberSectionDeformation)
            for k in range(NoOfSamplePoints):
                CoreAxialStressStrain = '%s-CoreAxialStressStrain-%s-%d-%d.dat'%(ModelName, timestamp, j, k)
                OutputFileNames.append(CoreAxialStressStrain)

    # endregion

    #region ########################## Display Results ##########################
    OData.AddObject(OpenSeesAPI.TCL.CodeTitle('Display Results'))

    # endregion

    #region ########################## Gravity Analysis ##########################
    OData.AddObject(OpenSeesAPI.TCL.CodeTitle('Gravity Analysis'))

    if AddGravityLoad:
        NoOfGravitySteps = 100
        OData.AddObject(OpenSeesAPI.Analysis.Constraints.Transformation())
        OData.AddObject(OpenSeesAPI.Analysis.Numberer.RCM())
        OData.AddObject(OpenSeesAPI.Analysis.System.UmfPack())
        # OData.AddObject(OpenSeesAPI.Analysis.System.BandGeneral())
        # OData.AddObject(OpenSeesAPI.Analysis.Test.NormUnbalance(Tolerance, NoOfIterations))
        OData.AddObject(OpenSeesAPI.Analysis.Test.NormDispIncr(Tolerance, NoOfIterations))
        OData.AddObject(OpenSeesAPI.Analysis.Algorithm.KrylovNewton())
        OData.AddObject(OpenSeesAPI.Analysis.Integrator.Static.LoadControl(1.0/NoOfGravitySteps,1,0.2,0.2))
        OData.AddObject(OpenSeesAPI.Analysis.Analysis.Static())
        OData.AddObject(OpenSeesAPI.TCL.TCLScript('set ok [analyze %d]'%NoOfGravitySteps))
        OData.AddObject(OpenSeesAPI.TCL.TCLScript('if {$ok == 0} {puts "Gravity Analysis Success" } else {puts "Gravity Analysis Failed"} '))
        OData.AddObject(OpenSeesAPI.TCL.TCLScript('loadConst -time 0.0'))

    # endregion

    # region ########################## StiffnessCheck ##########################
    OData.AddObject(OpenSeesAPI.TCL.CodeTitle('Stiffness Check'))

    StiffnessCheck = False

    def CheckStiffness(dir):
        OData.AddObject(OpenSeesAPI.TCL.TCLScript('puts "Stiffness Check in dir %d"'%dir))

        # Define Analysis
        OData.AddObject(OpenSeesAPI.Analysis.Constraints.Transformation())
        OData.AddObject(OpenSeesAPI.Analysis.Numberer.RCM())
        # OData.AddObject(OpenSeesAPI.Analysis.System.Mumps(Optional='-ICNTL 50'))
        OData.AddObject(OpenSeesAPI.Analysis.System.UmfPack())
        # OData.AddObject(OpenSeesAPI.Analysis.Test.NormUnbalance(1e-6, 200, 5))
        ControlNode = DiaNode[-1]

        # Load Pattern
        Loads = []
        if dir == 1:
            Nodes = [DiaNode[-1]]
            Load = 100
            for node in Nodes:
                Loads.append(OpenSeesAPI.TCL.TCLScript('load %d %f 0 0 0 0 0' % (node.id, Load)))
            OData.AddObject(OpenSeesAPI.Model.Pattern.Plain(200+dir, 'Linear', Loads))
            StepSize = YGrids[-1] * 0.01

        # Run Analysis
        MaxU = YGrids[-1] * 0.20
        MaxIteration = 1
        OData.AddObject(OpenSeesAPI.TCL.TCLScript('set ok 0;'))
        OData.AddObject(OpenSeesAPI.TCL.TCLScript('set MaxU %f;' % MaxU))
        OData.AddObject(OpenSeesAPI.TCL.TCLScript('set MaxStep %d;' % MaxIteration))
        OData.AddObject(OpenSeesAPI.TCL.TCLScript('set step 0;'))
        OData.AddObject(OpenSeesAPI.TCL.TCLScript('set currentDisp 0;'))

        OData.AddObject(OpenSeesAPI.TCL.TCLScript('set Yielded 0;'))
        OData.AddObject(OpenSeesAPI.TCL.TCLScript('set Stiffness 0.0001;'))
        OData.AddObject(OpenSeesAPI.TCL.TCLScript('set PreviousStiffness 0;'))
        OData.AddObject(OpenSeesAPI.TCL.TCLScript('set PreviousReaction 0;'))
        OData.AddObject(OpenSeesAPI.TCL.TCLScript('set PreviousDisp 0;'))
        OData.AddObject(OpenSeesAPI.TCL.TCLScript('set YieldReaction 0;'))
        OData.AddObject(OpenSeesAPI.TCL.TCLScript('set TotalReaction 0;'))

        OData.AddObject(OpenSeesAPI.TCL.TCLScript('while {$ok == 0 & $step < $MaxStep & $currentDisp < $MaxU} {'))

        OData.AddObject(OpenSeesAPI.Analysis.Test.NormDispIncr(Tolerance, NoOfIterations, 2))
        OData.AddObject(OpenSeesAPI.Analysis.Algorithm.KrylovNewton(MaxDim=10))
        # OData.AddObject(OpenSeesAPI.Analysis.Integrator.Static.DisplacementControl(ControlNode, dir, StepSize))
        OData.AddObject(OpenSeesAPI.Analysis.Integrator.Static.LoadControl(1.0, 1, 0, 0))
        OData.AddObject(OpenSeesAPI.Analysis.Analysis.Static())

        OData.AddObject(OpenSeesAPI.TCL.TCLScript('set PreviousDisp [expr $currentDisp]'))
        OData.AddObject(OpenSeesAPI.TCL.TCLScript('set ok [analyze 1]'))

        import OSAnalysisHelper
        OSAnalysisHelper.PushOverSolutionAlgorithimDispIncr(OData, StepSize / 1.e1, Tolerance, ControlNode)
        OSAnalysisHelper.PushOverSolutionAlgorithimDispIncr(OData, StepSize / 1.e2, Tolerance, ControlNode)
        OSAnalysisHelper.PushOverSolutionAlgorithimDispIncr(OData, StepSize / 1.e3, Tolerance, ControlNode)

        OData.AddObject(OpenSeesAPI.TCL.TCLScript('set currentDisp [nodeDisp %d %d]' % (ControlNode.id, dir)))
        OData.AddObject(OpenSeesAPI.TCL.TCLScript('set RoofDrift [expr $currentDisp/%f]"' % YGrids[-1]))
        OData.AddObject(OpenSeesAPI.TCL.TCLScript('puts "Current Roof Displ: $RoofDrift "'))
        OData.AddObject(OpenSeesAPI.TCL.TCLScript('puts "Running Push Over Step: $step"'))

        OData.AddObject(OpenSeesAPI.TCL.TCLScript('set PreviousStiffness [expr $Stiffness];'))
        OData.AddObject(OpenSeesAPI.TCL.TCLScript('set PreviousReaction $TotalReaction;'))

        # Find Out When The Structure Yields and then Stop analysis at 80% of the Yield Strength
        # GroundFloorColumns=filter(lambda x: SupportNodes.__contains__(x._NodeI),OData._Elements)
        for i in range(len(SupportZeroLengthElements)):
            # OData.AddObject(OpenSeesAPI.TCL.TCLScript('puts [lindex [eleResponse %d forces] 0] ;'%(SupportZeroLengthElements[i].id)))
            OData.AddObject(OpenSeesAPI.TCL.TCLScript('set NodeReaction%d [lindex [eleResponse %d forces] %d];' % (
            i, SupportZeroLengthElements[i].id, dir - 1)))
        OData.AddObject(OpenSeesAPI.TCL.TCLScript('set TotalReaction [expr %s];' % (
        ''.join(['$NodeReaction%d+' % x for x in range(0, len(SupportZeroLengthElements))])[:-1])))

        # NodeReaction is not working
        # for i in range(len(SupportNodes)):
        #     OData.AddObject(OpenSeesAPI.TCL.TCLScript('set NodeReaction%d [nodeReaction %d 1];'%(i,SupportNodes[i].id)))
        # OData.AddObject(OpenSeesAPI.TCL.TCLScript('set TotalReaction [expr %s];'%(''.join(map(lambda x: '$NodeReaction%d+'%x,range(0,len(SupportNodes))))[:-1])))

        # OData.AddObject(OpenSeesAPI.TCL.TCLScript(
        #     'set Stiffness [expr abs(($TotalReaction-$PreviousReaction)/($currentDisp-$PreviousDisp))];'))
        OData.AddObject(OpenSeesAPI.TCL.TCLScript(
                'set Stiffness [expr abs((%f)/($currentDisp-$PreviousDisp))];'%Load))

        OData.AddObject(OpenSeesAPI.TCL.TCLScript('puts "Stiffness (kips per in): $Stiffness";'))
        # Check If 60 Percent Yield Reaction Reached
        OData.AddObject(OpenSeesAPI.TCL.TCLScript(
            'if {[expr 0.6*abs($YieldReaction)] > [expr abs($TotalReaction)] & $Yielded == 1} {break};'))

        # Trigger Yielding
        OData.AddObject(OpenSeesAPI.TCL.TCLScript(
            'if {$step != 0 & [expr abs(($Stiffness-$PreviousStiffness)/$Stiffness)] > 1e-3} {'))
        OData.AddObject(OpenSeesAPI.TCL.TCLScript('set Yielded 1;'))
        OData.AddObject(OpenSeesAPI.TCL.TCLScript('set YieldReaction [expr abs($TotalReaction)]\n};'))

        OData.AddObject(OpenSeesAPI.TCL.TCLScript('set step [expr $step+1]'))

        OData.AddObject(OpenSeesAPI.TCL.TCLScript('}'))
        OData.AddObject(
            OpenSeesAPI.TCL.TCLScript('if {$ok == 0} {puts "Analysis Success"} else { puts "Analysis Failed" }'))

    if StiffnessCheck:
        CheckStiffness(1)
        OData.AddObject(OpenSeesAPI.TCL.TCLScript('reset'))

    # endregion

    #region ########################## Pushover Analysis ##########################
    OData.AddObject(OpenSeesAPI.TCL.CodeTitle('Pushover Analysis'))

    if PushOver and not(StiffnessCheck):
        #Define Analysis
        OData.AddObject(OpenSeesAPI.Analysis.Constraints.Transformation())
        OData.AddObject(OpenSeesAPI.Analysis.Numberer.RCM())
        # OData.AddObject(OpenSeesAPI.Analysis.System.Mumps(Optional='-ICNTL 50'))
        OData.AddObject(OpenSeesAPI.Analysis.System.UmfPack())
        OData.AddObject(OpenSeesAPI.Analysis.Test.EnergyIncr(Tolerance, NoOfIterations, 5))
        ControlNode = DiaNode[-1]
        StepSize = YGrids[-1]*0.0001
        #Load Pattern
        Loads = []

        if POELFForces:
            import ASCEHelper
            Csx = ASCEHelper.ComputeELFForces(Archetype.Mass , YGrids[1:]/12., CuTa)
            for i in range(1, len(YGrids)):
                if len(SplitNodes[i*NoOfDivisionsPerFloor]) == 0:
                    Nodes = [CoreNodes[i*NoOfDivisionsPerFloor]]#OData.GetNodesByYCoordinate(YGrids[i],1)
                else:
                    Nodes = [SplitNodes[i*NoOfDivisionsPerFloor][0]]
                Nodes = list([x for x in Nodes if hasattr(x, 'Used')])  # Filter for used nodes
                for node in Nodes:
                    Loads.append(OpenSeesAPI.TCL.TCLScript('load %d %.2f 0 0' % (node.id, Csx[i-1])))
        elif POModalForces:
            for i in range(1,len(YGrids)):
                if len(SplitNodes[i*NoOfDivisionsPerFloor]) == 0:
                    Nodes = [CoreNodes[i*NoOfDivisionsPerFloor]]#OData.GetNodesByYCoordinate(YGrids[i],1)
                else:
                    Nodes = [SplitNodes[i*NoOfDivisionsPerFloor][0]]
                Nodes = list([x for x in Nodes if hasattr(x,'Used')]) #Filter for used nodes
                for node in Nodes:
                    Loads.append(OpenSeesAPI.TCL.TCLScript('load %d [expr abs($EigenVector%d%d_X)] 0 0'%(node.id, 1,i)))
        else:
            for i in range(1,len(YGrids)-1):
                if YGrids[i-1] < heff and YGrids[i] > heff:
                    if len(SplitNodes[i*NoOfDivisionsPerFloor]) == 0:
                        Nodes = [CoreNodes[i*NoOfDivisionsPerFloor]]  # OData.GetNodesByYCoordinate(YGrids[i],1)
                    else:
                        Nodes = [SplitNodes[i*NoOfDivisionsPerFloor][0], SplitNodes[i*NoOfDivisionsPerFloor][1]]
                    Nodes = list([x for x in Nodes if hasattr(x,'Used')]) #Filter for used nodes
                    for node in Nodes:
                        MomentArm = -1*(heff-YGrids[i-1])
                        Loads.append(OpenSeesAPI.TCL.TCLScript('load %d 1 0 %.2f'%(node.id, MomentArm)))

        OData.AddObject(OpenSeesAPI.Model.Pattern.Plain(200,'Linear', Loads))

        #Run Analysis
        MaxU = YGrids[-1]*MaxPORoofDrift
        MaxIteration = 10000
        OData.AddObject(OpenSeesAPI.TCL.TCLScript('set ok 0;'))
        OData.AddObject(OpenSeesAPI.TCL.TCLScript('set MaxU %f;'%MaxU))
        OData.AddObject(OpenSeesAPI.TCL.TCLScript('set MaxStep %d;'%MaxIteration))
        OData.AddObject(OpenSeesAPI.TCL.TCLScript('set step 0;'))
        OData.AddObject(OpenSeesAPI.TCL.TCLScript('set currentDisp 0;'))

        OData.AddObject(OpenSeesAPI.TCL.TCLScript('set Yielded 0;'))
        OData.AddObject(OpenSeesAPI.TCL.TCLScript('set Stiffness 0.0001;'))
        OData.AddObject(OpenSeesAPI.TCL.TCLScript('set PreviousStiffness 0;'))
        OData.AddObject(OpenSeesAPI.TCL.TCLScript('set PreviousReaction 0;'))
        OData.AddObject(OpenSeesAPI.TCL.TCLScript('set PreviousDisp 0;'))
        OData.AddObject(OpenSeesAPI.TCL.TCLScript('set YieldReaction 0;'))
        OData.AddObject(OpenSeesAPI.TCL.TCLScript('set TotalReaction 0;'))
        OData.AddObject(OpenSeesAPI.TCL.TCLScript('set MaxReaction 0;'))

        OData.AddObject(OpenSeesAPI.TCL.TCLScript('while {$ok == 0 & $step < $MaxStep & $currentDisp < $MaxU} {'))

        OData.AddObject(OpenSeesAPI.Analysis.Test.NormDispIncr(Tolerance, NoOfIterations, 2))
        OData.AddObject(OpenSeesAPI.Analysis.Algorithm.KrylovNewton(MaxDim=10))
        OData.AddObject(OpenSeesAPI.Analysis.Integrator.Static.DisplacementControl(ControlNode, 1, StepSize))
        OData.AddObject(OpenSeesAPI.Analysis.Analysis.Static())

        OData.AddObject(OpenSeesAPI.TCL.TCLScript('set PreviousDisp [expr $currentDisp]'))
        OData.AddObject(OpenSeesAPI.TCL.TCLScript('set ok [analyze 1]'))

        import OSAnalysisHelper
        OSAnalysisHelper.PushOverSolutionAlgorithimConstantAlgorithmDispIncr(OData, StepSize / 1.e1, Tolerance, ControlNode, NoOfIterations)
        OSAnalysisHelper.PushOverSolutionAlgorithimConstantAlgorithmDispIncr(OData, StepSize / 1.e2, Tolerance, ControlNode, NoOfIterations)
        OSAnalysisHelper.PushOverSolutionAlgorithimConstantAlgorithmDispIncr(OData, StepSize / 1.e3, Tolerance, ControlNode, NoOfIterations)
        # OSAnalysisHelper.PushOverSolutionAlgorithimConstantAlgorithm(OData, StepSize / 1.e3, 1.e-6, ControlNode)
        OSAnalysisHelper.PushOverSolutionAlgorithimConstantAlgorithmDispIncr(OData, StepSize / 1.e4, Tolerance*10, ControlNode, NoOfIterations)
        # OSAnalysisHelper.PushOverSolutionAlgorithimConstantAlgorithm(OData, StepSize / 1.e4, 1e-4, ControlNode)

        OData.AddObject(OpenSeesAPI.TCL.TCLScript('set currentDisp [nodeDisp %d 1]'%ControlNode.id))
        OData.AddObject(OpenSeesAPI.TCL.TCLScript('set RoofDrift [expr $currentDisp/%f*100.]"'%YGrids[-1]))
        OData.AddObject(OpenSeesAPI.TCL.TCLScript('puts "Current Roof Drift (per) : $RoofDrift "'))
        OData.AddObject(OpenSeesAPI.TCL.TCLScript('puts "Running Push Over Step: $step"'))

        OData.AddObject(OpenSeesAPI.TCL.TCLScript('set PreviousStiffness [expr $Stiffness];'))
        OData.AddObject(OpenSeesAPI.TCL.TCLScript('set PreviousReaction $TotalReaction;'))

        #Find Out When The Structure Yields and then Stop analysis at 80% of the Yield Strength
        # GroundFloorColumns=filter(lambda x: SupportNodes.__contains__(x._NodeI),OData._Elements)
        for i in range(len(SupportZeroLengthElements)):
            # OData.AddObject(OpenSeesAPI.TCL.TCLScript('puts [lindex [eleResponse %d forces] 0] ;'%(SupportZeroLengthElements[i].id)))
            OData.AddObject(OpenSeesAPI.TCL.TCLScript('set NodeReaction%d [lindex [eleResponse %d forces] 0];'%(i,SupportZeroLengthElements[i].id)))
        OData.AddObject(OpenSeesAPI.TCL.TCLScript('set TotalReaction [expr %s];'%(''.join(['$NodeReaction%d+'%x for x in range(0,len(SupportZeroLengthElements))])[:-1])))

        #NodeReaction is not working
        # for i in range(len(SupportNodes)):
        #     OData.AddObject(OpenSeesAPI.TCL.TCLScript('set NodeReaction%d [nodeReaction %d 1];'%(i,SupportNodes[i].id)))
        # OData.AddObject(OpenSeesAPI.TCL.TCLScript('set TotalReaction [expr %s];'%(''.join(map(lambda x: '$NodeReaction%d+'%x,range(0,len(SupportNodes))))[:-1])))

        OData.AddObject(OpenSeesAPI.TCL.TCLScript('set Stiffness [expr abs(($TotalReaction-$PreviousReaction)/($currentDisp-$PreviousDisp))];'))
        OData.AddObject(OpenSeesAPI.TCL.TCLScript('puts "Stiffness (kips per in): $Stiffness";'))
        #Check If 60 Percent Yield Reaction Reached

        # OData.AddObject(OpenSeesAPI.TCL.TCLScript('if {[expr 0.6*abs($YieldReaction)] > [expr abs($TotalReaction)] & $Yielded == 1} {puts "Reaching 60% of the Strength" ;'))
        # OData.AddObject(OpenSeesAPI.TCL.TCLScript(
        #     'break}'))

        #Trigger Yielding
        OData.AddObject(OpenSeesAPI.TCL.TCLScript('if {$step != 0 & [expr abs(($Stiffness-$PreviousStiffness)/$Stiffness)] > 1e-3} {'))
        OData.AddObject(OpenSeesAPI.TCL.TCLScript('set Yielded 1;'))
        OData.AddObject(OpenSeesAPI.TCL.TCLScript('set YieldReaction [expr abs($TotalReaction)]\n};'))

        OData.AddObject(OpenSeesAPI.TCL.TCLScript('set step [expr $step+1]'))

        # FindMaxReaction
        OData.AddObject(OpenSeesAPI.TCL.TCLScript('if {$step != 0 & $MaxReaction < abs($TotalReaction)} {'))
        OData.AddObject(OpenSeesAPI.TCL.TCLScript('set $MaxReaction [expr abs($TotalReaction)]\n};'))
        OData.AddObject(OpenSeesAPI.TCL.TCLScript('}'))

        # Break if Reaction is lower than 1% of Max.
        OData.AddObject(OpenSeesAPI.TCL.TCLScript('if {$step != 0 & abs($TotalReaction) < [expr 0.01*$MaxReaction]} {'))
        OData.AddObject(OpenSeesAPI.TCL.TCLScript('break;'))
        OData.AddObject(OpenSeesAPI.TCL.TCLScript('}'))

        OData.AddObject(OpenSeesAPI.TCL.TCLScript('}'))

        OData.AddObject(OpenSeesAPI.TCL.TCLScript('if {$ok == 0} {puts "Analysis Success"} else { puts "Analysis Failed" }'))

    # endregion

    #region ########################## Cyclic Static Analysis ##########################
    OData.AddObject(OpenSeesAPI.TCL.CodeTitle('Cyclic Static Analysis'))

    if CyclicStatic:
        #Define Analysis
        OData.AddObject(OpenSeesAPI.Analysis.Constraints.Transformation())
        OData.AddObject(OpenSeesAPI.Analysis.Numberer.RCM())
        OData.AddObject(OpenSeesAPI.Analysis.System.UmfPack())
        # OData.AddObject(OpenSeesAPI.Analysis.System.Mumps(Optional='-ICNTL 50'))
        # OData.AddObject(OpenSeesAPI.Analysis.Test.EnergyIncr(Tolerance, NoOfIterations, 5))
        OData.AddObject(OpenSeesAPI.Analysis.Test.NormDispIncr(Tolerance, NoOfIterations, 0))
        ControlNode = DiaNode[-1]
        #Load Pattern
        Loads = []

        # for i in range(1,len(YGrids)):
        #     Nodes = OData.GetNodesByYCoordinate(YGrids[i],1)
        #     Nodes = list(filter(lambda x: hasattr(x,'Used'),Nodes)) #Filter for used nodes
        #     for node in Nodes:
        #         Loads.append(OpenSeesAPI.TCL.TCLScript('load %d [expr abs($EigenVector%d%d_X)] 0 0'%(node.id, 1,i)))
        # OData.AddObject(OpenSeesAPI.Model.Pattern.Plain(200,'Linear', Loads))

        if POELFForces:
            import ASCEHelper
            Csx = ASCEHelper.ComputeELFForces(Archetype.Mass , np.array(YGrids[1:])/12., CuTa)
            for i in range(1, len(YGrids)):
                if len(SplitNodes[int(i*NoOfDivisionsPerFloor)]) == 0:
                    Nodes = [CoreNodes[int(i*NoOfDivisionsPerFloor)]]#OData.GetNodesByYCoordinate(YGrids[i],1)
                else:
                    Nodes = [SplitNodes[int(i*NoOfDivisionsPerFloor)][0]]
                Nodes = list([x for x in Nodes if hasattr(x, 'Used')])  # Filter for used nodes
                for node in Nodes:
                    Loads.append(OpenSeesAPI.TCL.TCLScript('load %d %.2f 0 0' % (node.id, Csx[i-1])))
        elif POModalForces:
            for i in range(1,len(YGrids)):
                if len(SplitNodes[int(i*NoOfDivisionsPerFloor)]) == 0:
                    Nodes = [CoreNodes[int(i*NoOfDivisionsPerFloor)]]#OData.GetNodesByYCoordinate(YGrids[i],1)
                else:
                    Nodes = [SplitNodes[int(i*NoOfDivisionsPerFloor)][0]]
                Nodes = list([x for x in Nodes if hasattr(x,'Used')]) #Filter for used nodes
                for node in Nodes:
                    Loads.append(OpenSeesAPI.TCL.TCLScript('load %d [expr abs($EigenVector%d%d_X)] 0 0'%(node.id, 1,i)))
        else:
            for i in range(1,len(YGrids)-1):
                if YGrids[i-1] < heff and YGrids[i] > heff:
                    if len(SplitNodes[int(i*NoOfDivisionsPerFloor)]) == 0:
                        Nodes = [CoreNodes[int(i*NoOfDivisionsPerFloor)]]  # OData.GetNodesByYCoordinate(YGrids[i],1)
                    else:
                        Nodes = [SplitNodes[int(i*NoOfDivisionsPerFloor)][0], SplitNodes[int(i*NoOfDivisionsPerFloor)][1]]
                    Nodes = list([x for x in Nodes if hasattr(x,'Used')]) #Filter for used nodes
                    for node in Nodes:
                        MomentArm = -1*(heff-YGrids[i-1])
                        Loads.append(OpenSeesAPI.TCL.TCLScript('load %d 1 0 %.2f'%(node.id, MomentArm)))

        OData.AddObject(OpenSeesAPI.Model.Pattern.Plain(200, 'Linear', Loads))

        # Run Analysis
        OData.AddObject(OpenSeesAPI.TCL.TCLScript('set ok 0;'))
        OData.AddObject(OpenSeesAPI.TCL.TCLScript('set step 0;'))
        OData.AddObject(OpenSeesAPI.TCL.TCLScript('set prevDisp 0'))
        OData.AddObject(
            OpenSeesAPI.TCL.TCLScript('set Drifts [list %s]' % (''.join(['%f \t' % (x) for x in np.array(DriftHistory)*YGrids[-1]]))))

        OData.AddObject(OpenSeesAPI.TCL.TCLScript('foreach targetDisp $Drifts {'))
        OData.AddObject(OpenSeesAPI.Analysis.Algorithm.Newton())

        OData.AddObject(OpenSeesAPI.TCL.TCLScript('set StepSize [expr ($targetDisp-$prevDisp)]'))
        OData.AddObject(
            OpenSeesAPI.TCL.TCLScript('integrator DisplacementControl %d 1 $StepSize' % ControlNode._id))
        OData.AddObject(OpenSeesAPI.TCL.TCLScript('set ok [analyze 1]'))

        # Try Reducing Step Size
        StepReduction = 10.
        OData.AddObject(OpenSeesAPI.TCL.TCLScript('if {$ok != 0} {'))
        OData.AddObject(OpenSeesAPI.TCL.TCLScript('set ok 0'))
        OData.AddObject(
            OpenSeesAPI.TCL.TCLScript('integrator DisplacementControl %d 1 [expr $StepSize/%f]' %(ControlNode._id, StepReduction)))
        OData.AddObject(OpenSeesAPI.TCL.TCLScript('while {$ok == 0 & [expr abs($prevDisp - $targetDisp)] > 0.01 } {'))
        OData.AddObject(OpenSeesAPI.TCL.TCLScript('puts "Running Ministeps: Target: $targetDisp Current: $prevDisp"'))
        OData.AddObject(OpenSeesAPI.Analysis.Algorithm.KrylovNewton(MaxDim=10))
        OData.AddObject(OpenSeesAPI.TCL.TCLScript('set ok [analyze 1 ]'))
        OData.AddObject(OpenSeesAPI.TCL.TCLScript('set prevDisp [nodeDisp %d 1]' % ControlNode._id))
        OData.AddObject(OpenSeesAPI.TCL.TCLScript('}'))
        OData.AddObject(OpenSeesAPI.TCL.TCLScript('}'))

        # Try Reducing Step Size
        StepReduction = 100.
        OData.AddObject(OpenSeesAPI.TCL.TCLScript('if {$ok != 0} {'))
        OData.AddObject(OpenSeesAPI.TCL.TCLScript('set ok 0'))
        OData.AddObject(
            OpenSeesAPI.TCL.TCLScript('integrator DisplacementControl %d 1 [expr $StepSize/%f]' %(ControlNode._id, StepReduction)))
        OData.AddObject(OpenSeesAPI.TCL.TCLScript('while {$ok == 0 & [expr abs($prevDisp - $targetDisp)] > 0.01 } {'))
        OData.AddObject(OpenSeesAPI.TCL.TCLScript('puts "Running MiniMinisteps: Target: $targetDisp Current: $prevDisp"'))
        OData.AddObject(OpenSeesAPI.Analysis.Algorithm.KrylovNewton(MaxDim=10))
        OData.AddObject(OpenSeesAPI.TCL.TCLScript('set ok [analyze 1 ]'))
        OData.AddObject(OpenSeesAPI.TCL.TCLScript('set prevDisp [nodeDisp %d 1]' % ControlNode._id))
        OData.AddObject(OpenSeesAPI.TCL.TCLScript('}'))
        OData.AddObject(OpenSeesAPI.TCL.TCLScript('}'))

        # Try Reducing Step Size
        StepReduction = 1000.
        OData.AddObject(OpenSeesAPI.TCL.TCLScript('if {$ok != 0} {'))
        OData.AddObject(OpenSeesAPI.TCL.TCLScript('set ok 0'))
        OData.AddObject(
            OpenSeesAPI.TCL.TCLScript(
                'integrator DisplacementControl %d 1 [expr $StepSize/%f]' % (ControlNode._id, StepReduction)))
        OData.AddObject(OpenSeesAPI.TCL.TCLScript('while {$ok == 0 & [expr abs($prevDisp - $targetDisp)] > 0.01 } {'))
        OData.AddObject(OpenSeesAPI.TCL.TCLScript('puts "Running MiniMiniMinisteps: Target: $targetDisp Current: $prevDisp"'))
        OData.AddObject(OpenSeesAPI.Analysis.Algorithm.KrylovNewton(MaxDim=10))
        OData.AddObject(OpenSeesAPI.TCL.TCLScript('set ok [analyze 1 ]'))
        OData.AddObject(OpenSeesAPI.TCL.TCLScript('set prevDisp [nodeDisp %d 1]' % ControlNode._id))
        OData.AddObject(OpenSeesAPI.TCL.TCLScript('}'))
        OData.AddObject(OpenSeesAPI.TCL.TCLScript('}'))

        # Try Reducing Step Size
        StepReduction = 10000.
        OData.AddObject(OpenSeesAPI.TCL.TCLScript('if {$ok != 0} {'))
        OData.AddObject(OpenSeesAPI.TCL.TCLScript('set ok 0'))
        OData.AddObject(
            OpenSeesAPI.TCL.TCLScript(
                'integrator DisplacementControl %d 1 [expr $StepSize/%f]' % (ControlNode._id, StepReduction)))
        OData.AddObject(OpenSeesAPI.TCL.TCLScript('while {$ok == 0 & [expr abs($prevDisp - $targetDisp)] > 0.01 } {'))
        OData.AddObject(OpenSeesAPI.TCL.TCLScript('puts "Running MiniMiniMinisteps: Target: $targetDisp Current: $prevDisp"'))
        OData.AddObject(OpenSeesAPI.Analysis.Algorithm.KrylovNewton(MaxDim=10))
        OData.AddObject(OpenSeesAPI.TCL.TCLScript('set ok [analyze 1 ]'))
        OData.AddObject(OpenSeesAPI.TCL.TCLScript('set prevDisp [nodeDisp %d 1]' % ControlNode._id))
        OData.AddObject(OpenSeesAPI.TCL.TCLScript('}'))
        OData.AddObject(OpenSeesAPI.TCL.TCLScript('}'))

        OData.AddObject(OpenSeesAPI.TCL.TCLScript('if {$ok != 0} {break}'))

        OData.AddObject(OpenSeesAPI.TCL.TCLScript('incr step'))
        OData.AddObject(OpenSeesAPI.TCL.TCLScript('set prevDisp [nodeDisp %d 1]'%ControlNode._id))
        OData.AddObject(OpenSeesAPI.TCL.TCLScript('puts "Running Push Over Step: $step, Current: $prevDisp"'))
        # OData.AddObject(OpenSeesAPI.TCL.TCLScript('set prevDisp $targetDisp'))

        OData.AddObject(OpenSeesAPI.TCL.TCLScript('}'))

    # endregion

    #region ########################## Time History Analysis ##########################
    OData.AddObject(OpenSeesAPI.TCL.CodeTitle('Time History Analysis'))

    if TimeHistory:
        # Analysis Options
        OData.AddObject(OpenSeesAPI.Analysis.Constraints.Transformation())
        OData.AddObject(OpenSeesAPI.Analysis.Numberer.RCM())
        OData.AddObject(OpenSeesAPI.Analysis.Test.NormDispIncr(Tolerance, 1000, 0))
        # OData.AddObject(OpenSeesAPI.Analysis.Test.EnergyIncr(Tolerance, 1000))
        if not HHTTransientIntegrator:
            OData.AddObject(OpenSeesAPI.Analysis.Integrator.Transient.Newmark(0.5, 0.25))
        else:
            OData.AddObject(OpenSeesAPI.Analysis.Integrator.Transient.HHT(0.95))

        OData.AddObject(OpenSeesAPI.Analysis.Algorithm.Newton())
        OData.AddObject(OpenSeesAPI.Analysis.Analysis.Transient())

        # Load Pattern
        OData.AddObject(OpenSeesAPI.Model.Pattern.UniformExcitation(400, 1, TimeSeries))

        import OSAnalysisHelper
        AdvancedSolutionAlgorithm = False
        if not(AdvancedSolutionAlgorithm):
            from sys import platform
            if platform == 'darwin' or OpenSeesCommand == 'OpenSees': # Check to see if Mac... OpenSeesSP not available in mac yet
                OData.AddObject(OpenSeesAPI.Analysis.System.UmfPack())
            else:
                # OData.AddObject(OpenSeesAPI.Analysis.System.Mumps('-ICNTL 50'))
                OData.AddObject(OpenSeesAPI.Analysis.System.UmfPack())

            # Run Analysis
            OData.AddObject(OpenSeesAPI.TCL.TCLScript('set ok 0;'))
            OData.AddObject(OpenSeesAPI.TCL.TCLScript('set Nsteps %d;' % len(GMData)))
            OData.AddObject(OpenSeesAPI.TCL.TCLScript('set step 0;'))
            OData.AddObject(OpenSeesAPI.TCL.TCLScript('while {$ok == 0 & $step < [expr $Nsteps +1]} {'))
            OData.AddObject(OpenSeesAPI.TCL.TCLScript('set ok [analyze 1 %f]' % Dt))

            OSAnalysisHelper.SolutionAlgorithimKrylovOnly(OData, Dt / 10., Tolerance, 10, )
            OSAnalysisHelper.SolutionAlgorithimKrylovOnly(OData, Dt / 100., Tolerance, 100, )
            OSAnalysisHelper.SolutionAlgorithimKrylovOnly(OData, Dt / 1000., Tolerance * 10, 1000, )
            OSAnalysisHelper.SolutionAlgorithimKrylovOnly(OData, Dt / 10000., Tolerance * 10, 10000, )
            OSAnalysisHelper.SolutionAlgorithimKrylovOnly(OData, Dt / 1000., Tolerance * 100, 1000, )
            OSAnalysisHelper.SolutionAlgorithimKrylovOnly(OData, Dt / 10000., Tolerance * 100, 10000, )
            OSAnalysisHelper.SolutionAlgorithimKrylovOnly(OData, Dt / 1000., Tolerance * 1000, 1000, MaxDim=20)
            OSAnalysisHelper.SolutionAlgorithimKrylovOnly(OData, Dt / 10000., Tolerance * 1000, 10000, MaxDim=20)

            OData.AddObject(OpenSeesAPI.Analysis.Test.NormDispIncr(Tolerance, 1000, 0))
            OData.AddObject(OpenSeesAPI.Analysis.Algorithm.Newton())

            OData.AddObject(OpenSeesAPI.TCL.TCLScript('puts "Running Time History Step: $step out of %d"' % len(GMData)))
            OData.AddObject(OpenSeesAPI.TCL.TCLScript('set step [expr $step+1]'))
            if TrackPeriod:
                OData.AddObject(OpenSeesAPI.TCL.TCLScript('set lambdaN [eigen %s %d]; \n' % ('', 1)))
                i = 0
                OData.AddObject(OpenSeesAPI.TCL.TCLScript('if {[expr $lambdaN ]> 0} { \n'))
                OData.AddObject(OpenSeesAPI.TCL.TCLScript('set w%d [expr pow($lambdaN,0.5)]; \n' % (i + 1)))
                OData.AddObject(OpenSeesAPI.TCL.TCLScript('set T%d [expr 2.0*$pi/$w%d]; \n'%(i+1,i+1)))
                OData.AddObject(OpenSeesAPI.TCL.TCLScript('puts "InelasticPeriods: [expr $T1]";'))
                OData.AddObject(OpenSeesAPI.TCL.TCLScript('} else {'))
                OData.AddObject(OpenSeesAPI.TCL.TCLScript('puts "InelasticPeriods: 0";'))
                OData.AddObject(OpenSeesAPI.TCL.TCLScript('}'))

            # TrackDrifts = True
            TrackDriftsFileName = os.getcwd().replace('\\','/') + '/TempOutput/ATCWallArchetype/Output_(%s).dat'%OutputTag
            if TrackDrifts:
                NodesToTrack = [SupportNodes[0]] + DiaNode
                OData.AddObject(OpenSeesAPI.TCL.TCLScript('pwd'))
                OData.AddObject(OpenSeesAPI.TCL.TCLScript('set fp [open "%s" a+]'%(TrackDriftsFileName)))
                for i in range(0, len(NodesToTrack)):
                    OData.AddObject(OpenSeesAPI.TCL.TCLScript(
                        'set DiaNodeDisp%d [nodeDisp %s 1]' % (i, NodesToTrack[i].id)))
                for i in range(1, len(NodesToTrack)):
                    dY = NodesToTrack[i].Y - NodesToTrack[i-1].Y
                    OData.AddObject(OpenSeesAPI.TCL.TCLScript(
                        'set DiaDrift%d [expr ($DiaNodeDisp%d - $DiaNodeDisp%d) / %f ]' % (i, i, i-1, dY)))
                    OData.AddObject(OpenSeesAPI.TCL.TCLScript(
                        'puts $fp [format "$step %d Drift'%(i) + ' %.3f"' + ' [expr {$DiaDrift%d}]]'%(i)))

                OData.AddObject(OpenSeesAPI.TCL.TCLScript('close $fp'))

            OData.AddObject(OpenSeesAPI.TCL.TCLScript('}'))
        else:
            OData.AddObject(OpenSeesAPI.Analysis.System.Mumps('-ICNTL 50'))
            OSAnalysisHelper.SenSolutionAlgorithim(OData, Dt, len(GMData), 1.e-6, NoOfIterations=5000)

        OData.AddObject(
            OpenSeesAPI.TCL.TCLScript('if {$ok == 0} {puts "Analysis Success"} else { puts "Analysis Failed" }'))

    # endregion

    #region ########################## Close File ##########################
    OData.AddObject(OpenSeesAPI.TCL.CodeTitle('Close File'))
    OData.AddObject(OpenSeesAPI.TCL.TCLScript('wipe;'))
    OData.AddObject(OpenSeesAPI.TCL.TCLScript('puts "Models Run Complete";'))

    # endregion

    ##############################################################################
    ### Start Running OpenSees File
    ##############################################################################

    #region ########################## Plot Geometry ##########################

    if Viewer:
        import OSViewer
        O = OSViewer.OpenSeesViewerGUI(OData, XGrids, YGrids, TwoDimensional=True)
        O.mainloop()
    #endregion

    ########################## Run OpenSees Script ##########################
    OData.Executable.StartAnalysis(SuppressOutput=SupressOutput)
    TCLFileLines = open(TCLFileDirectory+FileName, 'r').readlines()

    #region ########################## Plot Results ##########################
    OutputFolder = TCLFileDirectory + ResultDirectory

    def ReadFile(File):
        try:
            return np.genfromtxt(File)
        except:
            return np.genfromtxt(File, skip_footer=1)

    Displ = ReadFile(OutputFolder + '/' + Displacement_File_Name)
    Displ = Displ[:,:2]

    AllDispl = ReadFile(OutputFolder + '/' + AllStoriesDisp)

    Reac = ReadFile(OutputFolder + '/' + Reaction_File_Name)

    CoreDisp = ReadFile(OutputFolder + '/' + CoreDisp)

    AllNodeDispl = ReadFile(OutputFolder + '/' + AllNodeDispl)

    AllStoriesAcceleration = ReadFile(OutputFolder + '/' + AllStoriesAcceleration)

    CoreMoment = ReadFile(OutputFolder + '/' + CoreMoment)

    AllReactions = ReadFile(OutputFolder + '/' + FullReaction_File_Name)

    BaseExtremeFiberStrains1 = ReadFile(OutputFolder + '/' + BaseExtremeFiberStrains1)

    BaseExtremeFiberStrains2 = ReadFile(OutputFolder + '/' + BaseExtremeFiberStrains2)

    CoreSupportNodesReact = ReadFile(OutputFolder + '/' + CoreSupportNodesReact)

    StoryDrift = ReadFile(OutputFolder + '/' + StoryDrift)

    if not UseForceBasedElements:
        ShearSpringShears = ReadFile(OutputFolder + '/' + ShearSpringShears)

    # Read Shear and Moment
    try:
        CoreShearAndMomentAtStoryBottom = ReadFile(OutputFolder + '/' + CoreShearAndMomentAtStoryBottom)
    except:
        CoreShearAndMomentAtStoryBottom = None

    # Read Strains
    try:
        CoreStrainAtStoryBottom1 = ReadFile(OutputFolder + '/' + CoreStrainAtStoryBottom1)
        CoreStrainAtStoryBottom2 = ReadFile(OutputFolder + '/' + CoreStrainAtStoryBottom2)
    except:
        CoreStrainAtStoryBottom1 = None
        CoreStrainAtStoryBottom2 = None

    if Archetype.CouplingBeams is not None and Archetype.CouplingBeams[0] is not None:
        CouplingExtremeFiberStrainsTopLeft = ReadFile(OutputFolder + '/' + CouplingExtremeFiberStrainsTopLeft)
        CouplingExtremeFiberStrainsBottomLeft = ReadFile(OutputFolder + '/' + CouplingExtremeFiberStrainsBottomLeft)
        CouplingExtremeFiberStrainsTopRight = ReadFile(OutputFolder + '/' + CouplingExtremeFiberStrainsTopRight)
        CouplingExtremeFiberStrainsBottomRight = ReadFile(OutputFolder + '/' + CouplingExtremeFiberStrainsBottomRight)

    if EnhancedOutput:
        CoreForces = ReadFile(OutputFolder + '/' + CoreForces)

        CoreStress = np.zeros((NoOfSamplePoints,
                               (len(YGrids)-1) * NoOfIntPoints * NoOfDivisionsPerFloor,
                                  len(AllNodeDispl[:])))
        CoreStrain = np.array(CoreStress)
        CoreXLocation = np.zeros((NoOfSamplePoints,
                                  (len(YGrids)-1) * NoOfIntPoints * NoOfDivisionsPerFloor))
        CoreYLocation = np.zeros((NoOfSamplePoints,
                                  (len(YGrids)-1) * NoOfIntPoints * NoOfDivisionsPerFloor))

        CoreCrushingStrain = np.array(CoreXLocation)
        CoreYieldingStrain = np.array(CoreXLocation)
        CoreRuptureStrain = np.array(CoreXLocation)

        AllFiberSectionDeformation = []
        for file in AllFiberSectionDeformationFiles:
            AllFiberSectionDeformation.append(ReadFile(OutputFolder + '/' + file))

        L_ip = ATCWallArchetypeHelpers.GetL_IP(NoOfIntPoints)

        def FindElementIndex(ele, eleList):
            for i in range(len(eleList)):
                if ele.id == eleList[i].id:
                    return i

        for j in range(NoOfIntPoints):
            for k in range(NoOfSamplePoints):
                CoreAxialStressStrain = '%s-CoreAxialStressStrain-%s-%d-%d.dat' % (ModelName, timestamp, j, k)
                temp = ReadFile(OutputFolder + '/' + CoreAxialStressStrain)
                for i in range(1, len(YGrids)):
                    # Single Pier
                    if len(CoreWallElements[i-1]) == 1:
                        for l in range(NoOfDivisionsPerFloor):
                            eleInd = FindElementIndex(CoreWallElements[i-1][0], SingleCoreElements)
                            colNo = (i-1) * NoOfIntPoints * NoOfDivisionsPerFloor + j + NoOfIntPoints * l
                            CoreStress[k][colNo] = temp[:, 1 + ( eleInd * NoOfDivisionsPerFloor + l ) * 2]
                            CoreStrain[k][colNo] = temp[:, 2 + ( eleInd * NoOfDivisionsPerFloor + l ) * 2]

                            if Archetype.CustomSection is None:
                                CoreXLocation[k][colNo] = Archetype.l_w[0]*np.linspace(0,1,NoOfSamplePoints)[k] #this is currently incorrect if l_w changes between stories.
                            else:
                                CoreXLocation[k][colNo] = Archetype.CustomSection[-1].l_w * np.linspace(0, 1, NoOfSamplePoints)[
                                    k]  # this is currently incorrect if l_w changes between stories.
                            ylocI = CoreWallElements[ (i-1) * NoOfDivisionsPerFloor + l ][0]._NodeI.Y
                            ylocJ = CoreWallElements[(i - 1) * NoOfDivisionsPerFloor + l][0]._NodeJ.Y

                            CoreYLocation[k][colNo] = (ylocJ-ylocI) * (np.sum(L_ip[:j]) + L_ip[j]/2.) + ylocI

                            SectionIndex = NoOfDivisionsPerFloor*(i-1) + l
                            CoreCrushingStrain[k][colNo] = ATCWallArchetypeHelpers.GetEpsCU(
                                CoreWallSections[SectionIndex][0][j]._Section._fibers)
                            CoreYieldingStrain[k][colNo] = ATCWallArchetypeHelpers.GetEpsSY(
                                CoreWallSections[SectionIndex][0][j]._Section._fibers)
                            CoreRuptureStrain[k][colNo] = ATCWallArchetypeHelpers.GetEpsSU(
                                CoreWallSections[SectionIndex][0][j]._Section._fibers)

        # Add Double Pier Results
        for j in range(NoOfIntPoints):
            for k in range(int(NoOfSamplePointsDouble)):
                CoreAxialStressStrain = '%s-CoreAxialStressStrain-Double-%s-%d-%d.dat' % (ModelName, timestamp, j, k)
                temp = ReadFile(OutputFolder + '/' + CoreAxialStressStrain)
                for i in range(1, len(YGrids)):
                    if len(CoreWallElements[i - 1]) == 2:
                        for l in range(NoOfDivisionsPerFloor):
                            ElementIndex = (i - 1) * NoOfDivisionsPerFloor  + l
                            eleIndA = FindElementIndex(CoreWallElements[ElementIndex][0], DoubleCoreElements)
                            eleIndB = FindElementIndex(CoreWallElements[ElementIndex][1], DoubleCoreElements)

                            colNo = (i - 1) * NoOfIntPoints * NoOfDivisionsPerFloor + j + NoOfIntPoints * l

                            CoreStress[k][colNo] = temp[:, 1 + (eleIndA  ) * 2]
                            CoreStrain[k][colNo] = temp[:, 2 + (eleIndA  ) * 2]
                            CoreStress[k+NoOfSamplePoints/2 + 1][colNo] = temp[:, 1 + (eleIndB  ) * 2]
                            CoreStrain[k+NoOfSamplePoints/2 + 1][colNo] = temp[:, 2 + (eleIndB  ) * 2]

                            # Find Total Wall Length
                            if type(Archetype.CustomSection[-1]) is list:
                                WallTotalLength = Archetype.CustomSection[-1][0].l_w + Archetype.CustomSection[-1][1].l_w + Archetype.CouplingBeamLength[-1]
                            else:
                                WallTotalLength = Archetype.CustomSection[-1].l_w

                            CoreXLocation[k][colNo] = Archetype.CustomSection[0][0].l_w * \
                                                      np.linspace(0, 1, int(NoOfSamplePointsDouble))[k]  # this is currently incorrect if l_w changes between stories.
                            CoreXLocation[k + NoOfSamplePoints / 2 + 1][colNo] = WallTotalLength - Archetype.CustomSection[0][0].l_w * \
                                                      np.linspace(1, 0, int(NoOfSamplePointsDouble))[k]  # this is currently incorrect if l_w changes between stories.

                            ylocI = CoreWallElements[(i - 1) * NoOfDivisionsPerFloor + l][0]._NodeI.Y
                            ylocJ = CoreWallElements[(i - 1) * NoOfDivisionsPerFloor + l][0]._NodeJ.Y

                            CoreYLocation[k][colNo] = (ylocJ - ylocI) * (np.sum(L_ip[:j]) + L_ip[j] / 2.) + ylocI
                            CoreYLocation[k+NoOfSamplePoints/2  + 1][colNo] = (ylocJ - ylocI) * (np.sum(L_ip[:j]) + L_ip[j] / 2.) + ylocI

                            SectionIndex = NoOfDivisionsPerFloor * (i - 1) + l
                            CoreCrushingStrain[k][colNo] = ATCWallArchetypeHelpers.GetEpsCU(
                                CoreWallSections[SectionIndex][0][j]._Section._fibers)
                            CoreYieldingStrain[k][colNo] = ATCWallArchetypeHelpers.GetEpsSY(
                                CoreWallSections[SectionIndex][0][j]._Section._fibers)
                            CoreRuptureStrain[k][colNo] = ATCWallArchetypeHelpers.GetEpsSU(
                                CoreWallSections[SectionIndex][0][j]._Section._fibers)

                            # Check Index on Section
                            CoreCrushingStrain[k+NoOfSamplePoints/2  + 1][colNo] = ATCWallArchetypeHelpers.GetEpsCU(
                                CoreWallSections[SectionIndex][1][j]._Section._fibers)
                            CoreYieldingStrain[k+NoOfSamplePoints/2  + 1][colNo] = ATCWallArchetypeHelpers.GetEpsSY(
                                CoreWallSections[SectionIndex][1][j]._Section._fibers)
                            CoreRuptureStrain[k+NoOfSamplePoints/2  + 1][colNo] = ATCWallArchetypeHelpers.GetEpsSU(
                                CoreWallSections[SectionIndex][1][j]._Section._fibers)

        for j in range(NoOfIntPoints):
            for k in range(NoOfSamplePointsDouble):
                for i in range(1, len(YGrids)):
                    if len(CoreWallElements[i - 1]) == 2:
                        for l in range(NoOfDivisionsPerFloor):
                            colNo = (i - 1) * NoOfIntPoints * NoOfDivisionsPerFloor + j + NoOfIntPoints * l
                            CoreXLocation[NoOfSamplePoints / 2 - 1,colNo] = CoreXLocation[NoOfSamplePoints / 2 - 2,colNo]
                            CoreXLocation[NoOfSamplePoints / 2 ,colNo] = CoreXLocation[NoOfSamplePoints / 2 + 1,colNo]
                            CoreYLocation[NoOfSamplePoints / 2 - 1,colNo] = CoreYLocation[NoOfSamplePoints / 2 - 2,colNo]
                            CoreYLocation[NoOfSamplePoints / 2 ,colNo] = CoreYLocation[NoOfSamplePoints / 2 + 1,colNo]

        CoreAxialLoads = []
        if len(CoreWallElements[i-1]) == 1:
            for j in range(NoOfIntPoints):
                CoreAxialLoad = '%s-CoreAxialLoad-%s-%d.dat' % (ModelName, timestamp, j)
                temp = ReadFile(OutputFolder + '/' + CoreAxialLoad)
                CoreAxialLoads.append(temp)

            AxialHistory = np.zeros((len(temp[:,0]),NoOfDivisionsPerFloor*NoOfIntPoints*(len(YGrids)-1)))
            MomentHistory = np.zeros((len(temp[:,0]),NoOfDivisionsPerFloor*NoOfIntPoints*(len(YGrids)-1)))
            for j in range(NoOfIntPoints):
                for i in range(NoOfDivisionsPerFloor*(len(YGrids)-1)):
                    AxialHistory[:,i*NoOfIntPoints + j] = CoreAxialLoads[j][:,1 + i * 3]
                    MomentHistory[:, i * NoOfIntPoints + j] = CoreAxialLoads[j][:, 2 + i * 3]
        else:
            AxialHistory = []
            MomentHistory = []

    #### Read Log File

    LogLines = open(TCLFileDirectory + OData.Executable.LogFileName, 'r').readlines()
    # Get Whether Analysis is Successfull
    AnalysisSuccess = LogLines.__contains__('Analysis Success\n')
    GravitySuccess = LogLines.__contains__('Gravity Analysis Success\n')

    TrackedPeriods = []
    T3 = 0.0
    # Get Period
    if T1 is None and T2 is None:
        for line in LogLines:
            if line.startswith('T1 = '):
                T1 = float(line.split()[2])
            if line.startswith('T2 = '):
                T2 = float(line.split()[2])
            if line.startswith('T3 = '):
                T3 = float(line.split()[2])
            if line.startswith('InelasticPeriods:'):
                TrackedPeriods.append(float(line.split()[1]))
    else:
        T3 = None

    ##### Create Output File for Plotter

    GravitySteps = NoOfGravitySteps
    if PushOver:
        t = np.arange(GravitySteps, len(Displ[:, 0]))
    else:
        # This is untested
        DataPoints = len(GMData)
        t = Displ[GravitySteps:DataPoints, 0]

    NodesDispDict = {}
    HingeRotationDict = {}

    for i in range(len(AllUsedNodes)):
        node = AllUsedNodes[i]
        NodesDispDict[node.id] = [AllNodeDispl[:,1+(i)*3],AllNodeDispl[:,2+(i)*3],AllNodeDispl[:,3+(i)*3]]

    class Output:
        def __init__(self):
            self.t = t
            self.NodesDispDict = NodesDispDict
            self.HingeRotationDict = HingeRotationDict

    ##### Compute Interstory Drift
    InterStoryDrifts = []
    InterStoryDriftsWORigid = []
    MaxRoofDrifts = None
    ResidualDrifts = []
    FloorAcceleration = []

    #Compute Max Interstory Drift
    for i in range(1,len(YGrids)):
        if PDeltaColumn:
            node = OData.GetNodesByGrid(1,i,NodeType=1)[0]
            nodebot = OData.GetNodesByGrid(1,i-1,NodeType=1)[0]
            if i >= 2:
                nodebotbot = OData.GetNodesByGrid(1, i - 2, NodeType=1)[0]
        else:
            if len(SplitNodes[i]) > 1:
                node = SplitNodes[i][0]
            else:
                node = CoreNodes[i*NoOfDivisionsPerFloor]
            if len(SplitNodes[i-1]) > 1:
                nodebot = SplitNodes[i-1][0]
            else:
                nodebot = CoreNodes[(i-1)*NoOfDivisionsPerFloor]

            if i >= 2:
                nodebotbot = CoreNodes[(i-2)*NoOfDivisionsPerFloor]

        RelDisp = np.array(NodesDispDict[node.id][0]) - np.array(NodesDispDict[nodebot.id][0])

        InterStoryDrifts.append(max(abs(RelDisp))/(YGrids[i] - YGrids[i-1]))
        ResidualDrifts.append(abs(RelDisp[-1]) / (YGrids[i] - YGrids[i - 1]))
        FloorAcceleration.append(max(abs((AllStoriesAcceleration[:,i])/386.4)))

        if i >= 2:
            RelDispBot = (np.array(NodesDispDict[nodebot.id][0]) - np.array(NodesDispDict[nodebotbot.id][0])) / (
                        YGrids[i- 1] - YGrids[i - 2])
        else:
            RelDispBot = 0

        RelDispWORigidBody = (np.array(NodesDispDict[node.id][0]) - np.array(NodesDispDict[nodebot.id][0])) / (
                    YGrids[i] - YGrids[i - 1]) - RelDispBot

        InterStoryDriftsWORigid.append(max(abs(RelDispWORigidBody)))

    node = OData.GetNodesByGrid(1, len(YGrids)-1, NodeType=1)[0]
    nodebot = OData.GetNodesByGrid(1, NoOfBasementLevels, NodeType=1)[0]
    RelDisp = np.array(NodesDispDict[node.id][0]) - np.array(NodesDispDict[nodebot.id][0])
    MaxRoofDrift = max(abs(RelDisp)) / (YGrids[len(YGrids)-1] - YGrids[NoOfBasementLevels])

    ##### Delete Output files
    if not DebugMode:
        import shutil
        shutil.rmtree(TCLFileDirectory, ignore_errors=True)

    # Compute Moment at Base with Time
    def ComputeMomentAtBase(Moments, AxialLoads, Arm):
        NoOfPiers = len(Moments)
        TotalMoment = np.zeros(len(Moments[0]))
        CouplingMoment = np.zeros(len(Moments[0]))

        if NoOfPiers == 1:
            TotalMoment = Moments[0]
        else:
            CouplingMoment = (-1*AxialLoads[0] * Arm + AxialLoads[1] * Arm)
            TotalMoment = Moments[0] + Moments[1] + CouplingMoment

        return TotalMoment, CouplingMoment

    # Compute Base Moments and Coupling Moments and Degree of Coupling
    if len(SplitNodes[0]) <= 1:
        BasePierMoments = [CoreSupportNodesReact[:, 3 + 3 * 0]]
        BasePierAxialLoads = [CoreSupportNodesReact[:, 2 + 3 * 0]]
        PierMomentArm = 0
    else:
        BasePierMoments = [CoreSupportNodesReact[:, 3 + 3 * 0], CoreSupportNodesReact[:, 3 + 3 * 1]]
        BasePierAxialLoads = [CoreSupportNodesReact[:, 2 + 3 * 0], CoreSupportNodesReact[:, 2 + 3 * 1]]
        PierMomentArm = abs(SplitNodes[0][0].X - SplitNodes[0][1].X)/2.
    TotalMoment, CouplingMoment = ComputeMomentAtBase(BasePierMoments, BasePierAxialLoads, PierMomentArm)

    # Compute Max Strains in the Extreme Fibers for the Walls
    if len(SplitNodes[0]) <= 1:
        MaxExtremeFiberStrain = np.max([BaseExtremeFiberStrains1[:,2 + 2 * 0 ],BaseExtremeFiberStrains2[:,2 + 2 * 0 ]])
        MinExtremeFiberStrain = np.min([BaseExtremeFiberStrains1[:,2 + 2 * 0 ],BaseExtremeFiberStrains2[:,2 + 2 * 0 ]])
    else:
        MaxExtremeFiberStrain = np.max([np.max(BaseExtremeFiberStrains1[:, 2 + 2 * 0 ]),
                                       np.max(BaseExtremeFiberStrains2[:, 2 + 2 * 0 ]),
                                       np.max(BaseExtremeFiberStrains1[:, 2 + 2 * 1 ]),
                                       np.max(BaseExtremeFiberStrains2[:, 2 + 2 * 1 ]),
                                        ])
        MinExtremeFiberStrain = np.min([np.min(BaseExtremeFiberStrains1[:, 2 + 2 * 0 ]),
                                        np.min(BaseExtremeFiberStrains2[:, 2 + 2 * 0 ]),
                                        np.min(BaseExtremeFiberStrains1[:, 2 + 2 * 1 ]),
                                        np.min(BaseExtremeFiberStrains2[:, 2 + 2 * 1 ]),
                                        ])

    # Compute Max Axial Load Per Pier
    if len(SplitNodes[0]) <= 1: # If
        MaxBaseAxialForce = np.max(CoreSupportNodesReact[:, 2 + 3 * 0])
        MinBaseAxialForce = np.min(CoreSupportNodesReact[:, 2 + 3 * 0])
    else:
        MaxBaseAxialForce = np.max([CoreSupportNodesReact[:, 2 + 3 * 0], CoreSupportNodesReact[:, 2 + 3 * 1]])
        MinBaseAxialForce = np.min([CoreSupportNodesReact[:, 2 + 3 * 0], CoreSupportNodesReact[:, 2 + 3 * 1]])

    # Compute Maximum Shear Stress
    MaxShearStress = 0
    if Archetype.CustomSection is not None:
        if len(SplitNodes[0]) <= 1:
            MaxShearStress = np.max(np.abs(CoreSupportNodesReact[:, 1 + 3 * 0])) / Archetype.CustomSection[0].l_w / Archetype.CustomSection[0].t_w
        else:
            MaxShearStress = np.max(np.abs([CoreSupportNodesReact[:, 1 + 3 * 0], CoreSupportNodesReact[:, 1 + 3 * 1]])) / Archetype.CustomSection[0][0].l_w / Archetype.CustomSection[0][0].t_w # Assume Same Size

    # Compute Dissipated Energy
    import scipy.integrate as integrate
    DissipatedEnergy = integrate.trapz(-1*np.sum(Reac[GravitySteps:,1:], axis=1), Displ[GravitySteps:,1])

    # Compute Max Rotation in Coupling Beams
    NoYieldedCouplingBeams = 0
    eps_y = Archetype.fye/29000.
    if Archetype.CouplingBeams is not None and Archetype.CouplingBeams[0] is not None:
        for i in range(0, int(len(CouplingExtremeFiberStrainsTopLeft[0,1:])/2)):
            maxStrain = np.max( [CouplingExtremeFiberStrainsTopLeft[:, 2 + 2*i ],
                                 CouplingExtremeFiberStrainsBottomLeft[:, 2 + 2*i ],
                                 CouplingExtremeFiberStrainsTopRight[:, 2 + 2*i ],
                                 CouplingExtremeFiberStrainsBottomRight[:, 2 + 2*i ]])
            if maxStrain > eps_y:
                NoYieldedCouplingBeams += 1

    ShearAtBot = []
    MomentAtBot = []
    MaxStrainAtBot = []
    MinStrainAtBot = []

    if len(SplitNodes[0]) <= 1: # Uncoupled
        for i in range(len(CoreWallElementsAtStoryBottom)):
            # Compute Max. Shear vs. Story
            ShearAtBot.append(np.max(np.abs(CoreShearAndMomentAtStoryBottom[:, i * 6 + 2])))
            # Compute Max. Moment vs. Story
            MomentAtBot.append(np.max(np.abs(CoreShearAndMomentAtStoryBottom[:, i * 6 + 3])))
            # Compute Min Norm. Strain vs. Story
            MaxStrainAtBot.append(np.max([CoreStrainAtStoryBottom1[:, i * 2 + 2],
                                          CoreStrainAtStoryBottom2[:, i * 2 + 2]]))
            # Compute Max. Norm. Strain vs. Story
            MinStrainAtBot.append(np.min([CoreStrainAtStoryBottom1[:, i * 2 + 2],
                                          CoreStrainAtStoryBottom2[:, i * 2 + 2]]))
    else:
        pass # Coupled Direction

    # Compute Regularized Strains for the Steel and Concrete Materials - this will be used to normalize the strains later
    StoryEpsSu = []
    StoryEpsSy = []
    StoryEpsCu = []
    StoryEpsCuUnconfined = []

    for i in range(len(CoreWallElementsAtStoryBottom)):
    # Compute Crushing and Yield Strains for Each Story
        StoryEpsSu.append(ATCWallArchetypeHelpers.GetEpsSU(
            CoreWallElementsAtStoryBottom[i][0]._Section[0]._Section._fibers))
        StoryEpsSy.append(ATCWallArchetypeHelpers.GetEpsSY(
            CoreWallElementsAtStoryBottom[i][0]._Section[0]._Section._fibers))
        StoryEpsCu.append(ATCWallArchetypeHelpers.GetEpsCU(
            CoreWallElementsAtStoryBottom[i][0]._Section[0]._Section._fibers, True))
        StoryEpsCuUnconfined.append(ATCWallArchetypeHelpers.GetEpsCU(
            CoreWallElementsAtStoryBottom[i][0]._Section[0]._Section._fibers, False))

    class Data:
        def __init__(self):
            self.RoofDisplacements = Displ
            self.Reactions = Reac
            self.BaseShear = np.sum(Reac[:,1:], axis=1)
            self.BaseMoment = TotalMoment
            self.BaseCouplingMoment = CouplingMoment
            self.MaxBaseShear = np.max(np.abs(np.sum(Reac[:,1:], axis=1)))
            self.AllReactions = AllReactions
            self.AllDispl = AllDispl
            self.XGrids = XGrids
            self.YGrids = YGrids

            self.MaxShearStress = MaxShearStress
            self.DissipatedEnergy = DissipatedEnergy
            self.MaxBaseAxialForce = MaxBaseAxialForce
            self.MinBaseAxialForce = MinBaseAxialForce
            self.MaxExtremeFiberStrain = MaxExtremeFiberStrain
            self.MinExtremeFiberStrain = MinExtremeFiberStrain
            self.MaxDegreeOfCoupling = np.max(CouplingMoment/TotalMoment)
            self.NoYieldedCouplingBeams = NoYieldedCouplingBeams
            self.MaxRoofDrift = np.max(np.abs(Displ[:,1]/YGrids[-1]))

            self.CriticalSectionFractureStrain = ATCWallArchetypeHelpers.GetEpsSU(
                CoreWallSections[0][0][0]._Section._fibers)
            self.CriticalSectionYieldingStrain = ATCWallArchetypeHelpers.GetEpsSY(
                CoreWallSections[0][0][0]._Section._fibers)
            self.CriticalSectionCrushingStrain = ATCWallArchetypeHelpers.GetEpsCU(
                CoreWallSections[0][0][0]._Section._fibers, True)

            self.NoOfGravitySteps = NoOfGravitySteps

            self.AllNodeDispl = AllNodeDispl
            self.AllUsedNodes = AllUsedNodes

            self.CoreDisp = CoreDisp
            self.OData = OData

            self.SupportNodes = SupportNodes

            self.CoreWallElements = CoreWallElements

            self.StoryShear = ShearAtBot
            self.StoryMoment = MomentAtBot
            self.StoryMaxStrain = MaxStrainAtBot
            self.StoryMinStrain = MinStrainAtBot

            self.StoryEpsCu = StoryEpsCu
            self.StoryEpsSu = StoryEpsSu
            self.StoryEpsSy = StoryEpsSy
            self.StoryEpsCuUnconfined = StoryEpsCuUnconfined

            self.TCLFileLines = TCLFileLines

            self.CoreShearAndMomentAtStoryBottom = CoreShearAndMomentAtStoryBottom

            self.TrackedPeriods = TrackedPeriods

            if EnhancedOutput:
                self.CoreStress = CoreStress
                self.CoreStrain = CoreStrain
                self.CoreXLocation = CoreXLocation
                self.CoreYLocation = CoreYLocation

                self.CoreYieldingStrain = CoreYieldingStrain
                self.CoreRuptureStrain = CoreRuptureStrain
                self.CoreCrushingStrain = CoreCrushingStrain

                self.CoreMoment = CoreMoment
                self.CoreForces = CoreForces

                self.AllFiberSectionDeformation = AllFiberSectionDeformation

                self.NoOfSamplePoints = NoOfSamplePoints

                self.AxialLoadHistory = AxialHistory
                self.MomentLoadHistory = MomentHistory

                # self.BaseMaterials = BaseMaterials
                # self.BaseStrain = BaseStrain
                # self.BaseStress = BaseStress

            if not UseForceBasedElements:
                self.ShearSpringShears = ShearSpringShears

            self.NoOfIntPoints = NoOfIntPoints
            self.NoOfDivisionsPerFloor = NoOfDivisionsPerFloor

            self.InterStoryDrifts = InterStoryDrifts
            self.MaxInterStoryDrift = np.max(InterStoryDrifts)
            self.MaxRoofDrift = MaxRoofDrift
            self.ResidualDrifts = ResidualDrifts
            self.FloorAcceleration = FloorAcceleration
            self.InterStoryDriftsWORigidBody = InterStoryDriftsWORigid

            self.StoryDriftRaw = StoryDrift

            self.T1 = T1
            self.T2 = T2
            self.T3 = T3
            self.AnalysisSuccess = AnalysisSuccess
            self.GravitySuccess = GravitySuccess

            self.Archetype = Archetype

            self.CoupledDirection = False

            pass

    #region ########################## OpenSees Animation ##########################

    if Animation:
        import OSAnimationHelper
        import matplotlib.pylab as plt

        ani = OSAnimationHelper.StructureTimeHistoryAnimation(Output(), OData, XGrids, YGrids, TwoDimensional=True, PlotText='Concrete Core Wall Archetype', DeflectedShapeScaleFactor=10)

        # ani.save(os.getcwd()+'/Figures/'+'Wall.mp4', fps=100)
        # ani.SaveToMP4andGIF(os.getcwd()+'/Figures/','WALL.mp4')
        ani.setAnimationTime(10)
        ani.SaveToMP4('/Figures/','Wall.mp4')
        # ani.save('Figures/TBI-BRB.gif', writer='imagemagick', fps=4)
        # plt.show()

    #endregion

    return Data()

####################################################################################
#endregion
####################################################################################
