####################################################################################
#region Libraries
####################################################################################

from __future__ import absolute_import
import numpy as np
import os
import OpenSeesAPI
import OSMaterialHelper as OSMat
from six.moves import range

####################################################################################
#endregion
####################################################################################



####################################################################################
#region Defining Helpers
####################################################################################

def ComputeCustomPlanarWallFiberSection(OSDatabase, CustomSection, cover, height, NoOfIntPoints=5, max_mesh_Size=3,
                                        Elastic=True, RegularizeSteel=True, RegularizeFracture=False,
                                        CrushingStrengthRatio=0.2, GfccOGfc = 2.5, ConcreteMaterialModel='Concrete02', ConfinementModel='SaatRazvi',
                                        SteelMaterialModel = 'Steel02', SteelUltimateStrainTension = 0.2,
                                        SteelPostYieldStiffness = 0.00784, GfcOfpc=2., UseForceBased=True,
                                        UnconfinedBeta=0.01, Regularized = True, WallThicknessMultipler = 1.):
    import numpy as np

    bar_size = CustomSection.boundary_bar_size
    fpc = CustomSection.fce
    fy = CustomSection.fye
    fu = CustomSection.fu
    t_wall = CustomSection.t_w * WallThicknessMultipler
    l_wall = CustomSection.l_w

    if bar_size == None:
        A_bar = 0
    else:
        A_bar = np.pi * (bar_size / 8. / 2.) ** 2.

    L_ip = GetL_IP(NoOfIntPoints, ForceBased=UseForceBased)

    ConcConfinedAllIP = []
    SteelConfinedAllIP = []
    ConcUnconfinedAllIP = []
    SteelUnconfinedAllIP = []

    if bar_size is None:
        Abar = 0
        bar_size = 0.01
    else:
        Abar = np.pi * (bar_size / 8. / 2.) ** 2.

    import OSMaterialHelper as OSMat
    for i in range(len(L_ip)):
        Aweb = np.pi * (CustomSection.bar_size_web / 8. / 2.) ** 2.

        if Elastic:
            Ec = 0.5 * 57. * (fpc * 1000) ** .5
            Es = 29000*1e-13
            ElasticConcrete = OSDatabase.AddMaterial(
            OpenSeesAPI.Model.Element.Material.UniaxialMaterial.Elastic(OSDatabase.GetFreeMaterialId(1, 1), Ec,
                                                                        _Notes='Concrete Testing'))
            ElasticSteel = OSDatabase.AddMaterial(
            OpenSeesAPI.Model.Element.Material.UniaxialMaterial.Elastic(OSDatabase.GetFreeMaterialId(1, 1), Es,
                                                                        _Notes='Steel Testing'))
            ConcConfinedAllIP.append(ElasticConcrete)
            SteelConfinedAllIP.append(ElasticSteel)
            ConcUnconfinedAllIP.append(ElasticConcrete)
            SteelUnconfinedAllIP.append(ElasticSteel)

            # Shear
            # Define Elastic Spring for Shear
            Ec = 57. * (fpc * 1000.) ** .5
            Mu = 0.2  # 15
            Gc = Ec / (2. * (1. + Mu))
            Geff = Gc  # 0.1 * Gc  # assume 10% of elastic shear modulus per Oyen (2006)
            ks = 5. / 6.  # reduction factor for shear area per Pugh (2012)

            Ax = (l_wall)*t_wall

            EgX = Geff * ks * Ax

            ShearElasticX = OSDatabase.AddMaterial(
                OpenSeesAPI.Model.Element.Material.UniaxialMaterial.Elastic(OSDatabase.GetFreeMaterialId(1, 1), EgX,
                                                                            _Notes='Shear Spring for Core Wall'))

        else:
            if CustomSection.right_boundary is not None:
                bar_spacing = (CustomSection.right_boundary - 2 * cover) / (len(CustomSection.right_reinf_layout) - 1.)
                tie_spacing = CustomSection.boundary_tie_spacing
            else:
                bar_spacing = 2. * Aweb / CustomSection.web_rho * 100 / t_wall * WallThicknessMultipler # assume two layers
                tie_spacing = CustomSection.boundary_tie_spacing

            tie_size = CustomSection.boundary_tie_bar # pugh assumes #4 only

            Et = 57000. * (fpc * 1000) ** .5
            ftu = 4.0 * (fpc * 1000) ** .5
            Gfrac = 5.07e-4 * 1000.
            Lt = L_ip[i] * height
            DeltaE = 2. * Gfrac / ftu / Lt
            Ets = ftu / DeltaE
            if not RegularizeFracture:
                EtsOEt = 0.05 #Ets/Et
            else:
                EtsOEt = Ets/Et

            # Unonfined Concrete
            ConcUnconfinedAllIP.append(GetPughUnconfinedConcreteMaterial(OSDatabase, fpc,
                                                                         L_ip[i] * height,
                                                                         EtsOEt=EtsOEt,
                                                                         fpuOfpc=UnconfinedBeta,
                                                                         GfcOfpc=GfcOfpc, Regularized = True))

            # ecu_Unconfined = ConcUnconfinedAllIP[-1]._ecu
            ecu_Unconfined = ConcUnconfinedAllIP[-1]._Materials[0]._ecu

            # Assume 4 for Minimum Reinforc.
            ANo4 = np.pi * (CustomSection.bar_size_web / 8. / 2.) ** 2.
            rho_min = 0.0025
            RS_spacing = 2. * ANo4 / t_wall / rho_min * WallThicknessMultipler

            if RegularizeSteel:
                SteelUnconfinedAllIP.append(GetPughSteel02Material(OSDatabase, fy, fu, L_ip[i] * height,
                                                                   ecu_Unconfined, eps_ult_exp=SteelUltimateStrainTension))
            else:
                if SteelMaterialModel == 'Steel02':
                    SteelUnconfinedAllIP.append(
                        GetPughSteel02MaterialUnregularized(OSDatabase, fy, fu, ecu_Unconfined,
                                                            eps_ult_exp=SteelUltimateStrainTension,
                                                            SteelPostYieldStiffness=SteelPostYieldStiffness))
                elif SteelMaterialModel == 'Steel02Fatigue':
                    SteelUnconfinedAllIP.append(
                        GetSteel02FatigueMaterialUnregularized(OSDatabase, fy, fu, ecu_Unconfined,
                                                            eps_ult_exp=SteelUltimateStrainTension,
                                                            SteelPostYieldStiffness=SteelPostYieldStiffness))
                elif SteelMaterialModel == 'ReinforcingSteel':
                        SteelUnconfinedAllIP.append(
                            GetReinforcingSteelMaterialUnregularized(OSDatabase, fy, fu, SteelUltimateStrainTension,
                                                                     RS_spacing, CustomSection.bar_size_web,
                                                                     SteelPostYieldStiffness=SteelPostYieldStiffness))
                elif SteelMaterialModel == 'ReinforcingSteelWithDMBuckling':
                        SteelUnconfinedAllIP.append(
                            GetReinforcingSteelMaterialUnregularized(OSDatabase, fy, fu, SteelUltimateStrainTension,
                                                                     RS_spacing, CustomSection.bar_size_web,
                                                                     SteelPostYieldStiffness=SteelPostYieldStiffness,
                                                                     withBuckling=True))
                elif SteelMaterialModel == 'ReinforcingSteelWithDMBucklingAndFatigue':
                    SteelUnconfinedAllIP.append(
                        GetReinforcingSteelMaterialUnregularized(OSDatabase, fy, fu, SteelUltimateStrainTension,
                                                                 RS_spacing, CustomSection.bar_size_web,
                                                                 withBuckling=True, withFatigue=True,
                                                                 SteelPostYieldStiffness=SteelPostYieldStiffness))
                elif SteelMaterialModel == 'ReinforcingSteelWithFatigue':
                    SteelUnconfinedAllIP.append(
                        GetReinforcingSteelMaterialUnregularized(OSDatabase, fy, fu, SteelUltimateStrainTension,
                                                                 RS_spacing, CustomSection.bar_size_web,
                                                                 withBuckling=False, withFatigue=True,
                                                                 SteelPostYieldStiffness=SteelPostYieldStiffness))
                elif SteelMaterialModel == 'Pinching4':
                    SteelUnconfinedAllIP.append(
                        GetPinching4MaterialUnregularized(OSDatabase, fy, fu, CustomSection.left_boundary,
                                                          SteelUltimateStrainTension, RS_spacing,
                                                          CustomSection.bar_size_web, -1.0,
                                                          SteelPostYieldStiffness=SteelPostYieldStiffness,
                                                          stirrup_dia=CustomSection.bar_size_web,
                                                          no_of_legs_y=1,
                                                          no_of_rows=1,
                                                          cover=cover,
                                                          WebRebar=True,))
                elif SteelMaterialModel == 'Steel02LB':
                    SteelUnconfinedAllIP.append(
                        GetSteel02LBMaterialUnregularized(OSDatabase, fy, fu, CustomSection.left_boundary,
                                                          SteelUltimateStrainTension, RS_spacing,
                                                          CustomSection.bar_size_web, -1.0,
                                                          SteelPostYieldStiffness=SteelPostYieldStiffness,
                                                          stirrup_dia=CustomSection.bar_size_web,
                                                          no_of_legs_y=1,
                                                          no_of_rows=1,
                                                          cover=cover,
                                                          WebRebar=True,))
                elif SteelMaterialModel == 'Steel02WithFatigue':
                    SteelUnconfinedAllIP.append(
                        GetPughSteel02MaterialUnregularizedWithFatigue(OSDatabase, fy, fu, ecu_Unconfined,
                                                                       eps_ult_exp=SteelUltimateStrainTension,
                                                                       SteelPostYieldStiffness=SteelPostYieldStiffness))
                elif SteelMaterialModel == 'Hysteretic':
                    SteelUnconfinedAllIP.append(
                        GetHystereticMaterialUnregularized(OSDatabase, fy, fu, CustomSection.left_boundary,
                                                          SteelUltimateStrainTension, RS_spacing,
                                                          CustomSection.bar_size_web, -1.0,
                                                          SteelPostYieldStiffness=SteelPostYieldStiffness,
                                                          stirrup_dia=CustomSection.bar_size_web,
                                                          no_of_legs_y=1,
                                                          no_of_rows=1,
                                                          cover=cover,
                                                          WebRebar=True,))

            if tie_size == 0 or tie_size is None:
                A_tie = 0.0001
                tie_size = 0.0001
                ConcConfinedAllIP.append(ConcUnconfinedAllIP[-1])
                SteelConfinedAllIP.append(SteelUnconfinedAllIP[-1])

            else:
                A_tie = np.pi * (tie_size / 8. / 2.) ** 2.

                # Confined Concrete
                if ConfinementModel == 'SaatRazvi':
                    ConcConfinedAllIP.append(GetPughConfinedConcreteMaterial(OSDatabase, fpc,
                                                                           CustomSection.left_boundary  - 2 * cover,
                                                                           t_wall-2*cover,
                                                                           bar_spacing, t_wall-2*cover,
                                                                           tie_spacing, tie_spacing,
                                                                           fy, A_tie,
                                                                           CustomSection.boundary_tie_x_no,
                                                                           CustomSection.boundary_tie_y_no,
                                                                           L_ip[i] * height, GfccOGfc=GfccOGfc,
                                                                           EtsOEt=EtsOEt,
                                                                             fpuOfpc=CrushingStrengthRatio,
                                                                             GfcOfpc=GfcOfpc, Regularized = Regularized,
                                                                             ConcreteMaterialModel=ConcreteMaterialModel))
                elif ConfinementModel == 'Mander':
                    if CustomSection.left_boundary != 0:
                       rho_l = np.sum(CustomSection.left_reinf_layout) * A_bar / t_wall / CustomSection.left_boundary
                    else:
                        rho_l = 0.02
                    ConcConfinedAllIP.append(GetManderConfinedConcreteMaterial(OSDatabase, fpc,
                                                                           CustomSection.left_boundary - 2 * cover,
                                                                           t_wall-2*cover, rho_l,
                                                                           bar_spacing, t_wall-2*cover,
                                                                           tie_spacing, tie_spacing,
                                                                           fy, A_tie,
                                                                           CustomSection.boundary_tie_x_no,
                                                                           CustomSection.boundary_tie_y_no,
                                                                           L_ip[i] * height, GfccOGfc=GfccOGfc,
                                                                           EtsOEt=EtsOEt,
                                                                             fpuOfpc=CrushingStrengthRatio,
                                                                               GfcOfpc=GfcOfpc, ConcreteMaterialModel=ConcreteMaterialModel))
                elif ConfinementModel == 'Richart':
                    if CustomSection.left_boundary != 0:
                       rho_l = np.sum(CustomSection.left_reinf_layout) * A_bar / t_wall / CustomSection.left_boundary
                    else:
                        rho_l = 0.02
                    ConcConfinedAllIP.append(GetRichartConfinedConcreteMaterial(OSDatabase, fpc,
                                                                           CustomSection.left_boundary - 2 * cover,
                                                                           t_wall-2*cover, rho_l,
                                                                           bar_spacing, t_wall-2*cover,
                                                                           tie_spacing, tie_spacing,
                                                                           fy, A_tie,
                                                                           CustomSection.boundary_tie_x_no,
                                                                           CustomSection.boundary_tie_y_no,
                                                                           L_ip[i] * height, GfccOGfc=GfccOGfc,
                                                                           EtsOEt=EtsOEt,
                                                                             fpuOfpc=CrushingStrengthRatio,
                                                                               GfcOfpc=GfcOfpc, ConcreteMaterialModel=ConcreteMaterialModel))
                elif ConfinementModel == 'Unconfined':
                    ConcConfinedAllIP.append(ConcUnconfinedAllIP[-1])

                ecu_Confined = ConcConfinedAllIP[-1]._Materials[0]._ecu

                if tie_spacing == None or tie_spacing == 0:
                    tie_spacing = 1e10

                if RegularizeSteel:
                    SteelConfinedAllIP.append(
                        GetPughSteel02Material(OSDatabase, fy, fu, L_ip[i] * height, ecu_Confined,
                                               eps_ult_exp=SteelUltimateStrainTension, SteelPostYieldStiffness=SteelPostYieldStiffness))
                else:
                    if SteelMaterialModel == 'Steel02':
                        SteelConfinedAllIP.append(
                            GetPughSteel02MaterialUnregularized(OSDatabase, fy, fu, ecu_Confined,
                                                                eps_ult_exp=SteelUltimateStrainTension,
                                                                SteelPostYieldStiffness=SteelPostYieldStiffness))
                    elif SteelMaterialModel == 'Steel02Fatigue':
                        SteelConfinedAllIP.append(
                            GetSteel02FatigueMaterialUnregularized(OSDatabase, fy, fu, ecu_Confined,
                                                                   eps_ult_exp=SteelUltimateStrainTension,
                                                                   SteelPostYieldStiffness=SteelPostYieldStiffness))
                    elif SteelMaterialModel == 'ReinforcingSteel':
                        SteelConfinedAllIP.append(
                            GetReinforcingSteelMaterialUnregularized(OSDatabase, fy, fu, SteelUltimateStrainTension,
                                                                     tie_spacing, bar_size,
                                                                     eps_ult_comp=ecu_Confined,
                                                                     SteelPostYieldStiffness=SteelPostYieldStiffness))
                    elif SteelMaterialModel == 'ReinforcingSteelWithDMBuckling':
                        SteelConfinedAllIP.append(
                            GetReinforcingSteelMaterialUnregularized(OSDatabase, fy, fu, SteelUltimateStrainTension,
                                                                     tie_spacing, bar_size,
                                                                     withBuckling=True, SteelPostYieldStiffness=SteelPostYieldStiffness))
                    elif SteelMaterialModel == 'ReinforcingSteelWithDMBucklingAndFatigue':
                        SteelConfinedAllIP.append(
                            GetReinforcingSteelMaterialUnregularized(OSDatabase, fy, fu, SteelUltimateStrainTension, tie_spacing, bar_size,
                                                                     withBuckling=True, withFatigue=True,
                                                                     SteelPostYieldStiffness=SteelPostYieldStiffness))
                    elif SteelMaterialModel == 'ReinforcingSteelWithFatigue':
                        SteelConfinedAllIP.append(
                            GetReinforcingSteelMaterialUnregularized(OSDatabase, fy, fu, SteelUltimateStrainTension, tie_spacing, bar_size,
                                                                     withBuckling=False, withFatigue=True,
                                                                     SteelPostYieldStiffness=SteelPostYieldStiffness))
                    elif SteelMaterialModel == 'Pinching4':
                        SteelConfinedAllIP.append(
                            GetPinching4MaterialUnregularized(OSDatabase, fy, fu, CustomSection.t_w,
                                                              SteelUltimateStrainTension, CustomSection.boundary_tie_spacing,
                                                              CustomSection.boundary_bar_size, -1.0,
                                                              SteelPostYieldStiffness=SteelPostYieldStiffness,
                                                              stirrup_dia=tie_size,
                                                              no_of_legs_y=1,#CustomSection.boundary_tie_y_no,
                                                              no_of_rows=1,#3,
                                                              cover=cover,
                                                              WebRebar=False, ))
                    elif SteelMaterialModel == 'Steel02LB':
                        SteelConfinedAllIP.append(
                            GetSteel02LBMaterialUnregularized(OSDatabase, fy, fu, CustomSection.t_w,
                                                              SteelUltimateStrainTension, CustomSection.boundary_tie_spacing,
                                                              CustomSection.boundary_bar_size, -1.0,
                                                              SteelPostYieldStiffness=SteelPostYieldStiffness,
                                                              stirrup_dia=tie_size,
                                                              no_of_legs_y=1,#CustomSection.boundary_tie_y_no,
                                                              no_of_rows=1,
                                                              cover=cover,
                                                              WebRebar=False, ))
                    elif SteelMaterialModel == 'Steel02WithFatigue':
                        SteelConfinedAllIP.append(
                            GetPughSteel02MaterialUnregularizedWithFatigue(OSDatabase, fy, fu, ecu_Confined,
                                                                           eps_ult_exp=SteelUltimateStrainTension,
                                                                           SteelPostYieldStiffness=SteelPostYieldStiffness))
                    elif SteelMaterialModel == 'Hysteretic':
                        SteelConfinedAllIP.append(
                            GetHystereticMaterialUnregularized(OSDatabase, fy, fu, CustomSection.t_w,
                                                              SteelUltimateStrainTension, CustomSection.boundary_tie_spacing,
                                                              CustomSection.boundary_bar_size, -1.0,
                                                              SteelPostYieldStiffness=SteelPostYieldStiffness,
                                                              stirrup_dia=tie_size,
                                                              no_of_legs_y=1,#CustomSection.boundary_tie_y_no,
                                                              no_of_rows=1,#3,
                                                              cover=cover,
                                                              WebRebar=False, ))

            # Shear
            # Define Elastic Spring for Shear
            Ec = 57. * (fpc * 1000.) ** .5
            Mu = 0.15  # 15
            Gc = Ec / (2. * (1. + Mu))
            Geff = 0.1 * Gc  # assume 10% of elastic shear modulus per Oyen (2006)
            ks = 5. / 6.  # reduction factor for shear area per Pugh (2012)

            Av = l_wall * t_wall

            GcOEc = Gc/Ec

            EgX = 0.4*Ec*Av*ks #Gc * Av

            ShearElasticX = OSDatabase.AddMaterial(
                OpenSeesAPI.Model.Element.Material.UniaxialMaterial.Elastic(OSDatabase.GetFreeMaterialId(1, 1), EgX,
                                                                            _Notes='Shear Spring for Core Wall'))

    Sections = []

    for j in range(len(L_ip)):
        ConcConfined = ConcConfinedAllIP[j]
        SteelConfined = SteelConfinedAllIP[j]
        ConcUnconfined = ConcUnconfinedAllIP[j]
        SteelUnconfined = SteelUnconfinedAllIP[j]

        Fibers = []

        # Confined Parts
        if True:
            left_BE = CustomSection.left_boundary
            right_BE = CustomSection.right_boundary

            l_web = CustomSection.l_w - left_BE - right_BE
            l_wall = CustomSection.l_w

            ## Web
            if l_web != 0.:
                ShiftY = -t_wall / 2.
                ShiftX = -l_wall / 2. + left_BE
                Fibers.append(OpenSeesAPI.Section.FiberSection.Patch.Rect(ConcUnconfined,
                                                                          int(np.ceil((l_wall - left_BE - right_BE ) / max_mesh_Size)),
                                                                          1,#int(np.ceil((t_wall) / max_mesh_Size)), # Making things faster
                                                                          ShiftX,
                                                                          ShiftY,
                                                                          ShiftX + l_wall - left_BE - right_BE,
                                                                          ShiftY + t_wall))

            if CustomSection.left_boundary != 0:
                # Left BE
                ShiftY = -t_wall / 2. + cover
                ShiftX = -l_wall / 2. + cover
                Fibers.append(OpenSeesAPI.Section.FiberSection.Patch.Rect(ConcConfined,
                                                                          int(np.ceil((left_BE - 1 * cover) / max_mesh_Size)),
                                                                          1, #int(np.ceil((t_wall - 2 * cover) / max_mesh_Size)),
                                                                          ShiftX,
                                                                          ShiftY,
                                                                          ShiftX + left_BE - 2 * cover,
                                                                          ShiftY + t_wall - 2 * cover))

            if CustomSection.right_boundary  != 0:
                # Right BE
                ShiftY = -t_wall / 2. + cover
                ShiftX = +l_wall / 2. - right_BE + cover
                Fibers.append(OpenSeesAPI.Section.FiberSection.Patch.Rect(ConcConfined,
                                                                          int(np.ceil((right_BE - cover) / max_mesh_Size)),
                                                                          1, #int(np.ceil((t_wall - 2 * cover) / max_mesh_Size)),
                                                                          ShiftX,
                                                                          ShiftY,
                                                                          ShiftX + right_BE - 2 * cover,
                                                                          ShiftY + t_wall - 2 * cover ))
            # Add Concrete Covers to Boundary Elements
            if CustomSection.left_boundary != 0:
                ShiftY = -t_wall / 2.
                ShiftX = -l_wall / 2.
                Fibers.append(OpenSeesAPI.Section.FiberSection.Patch.Rect(ConcUnconfined,
                                                                          int(np.ceil((left_BE) / max_mesh_Size)),
                                                                          1, #int(np.ceil((cover) / max_mesh_Size)),
                                                                          ShiftX,
                                                                          ShiftY,
                                                                          ShiftX + left_BE,
                                                                          ShiftY + cover))
                ShiftY = -t_wall / 2. + t_wall - cover
                ShiftX = -l_wall / 2.
                Fibers.append(OpenSeesAPI.Section.FiberSection.Patch.Rect(ConcUnconfined,
                                                                          int(np.ceil((left_BE) / max_mesh_Size)),
                                                                          1, #int(np.ceil((cover) / max_mesh_Size)),
                                                                          ShiftX,
                                                                          ShiftY,
                                                                          ShiftX + left_BE,
                                                                          ShiftY + cover))

                ShiftY = -t_wall / 2. + cover
                ShiftX = -l_wall / 2.
                Fibers.append(OpenSeesAPI.Section.FiberSection.Patch.Rect(ConcUnconfined,
                                                                          int(np.ceil((cover) / max_mesh_Size)),
                                                                          1, #int(np.ceil((t_wall - 2 * cover) / max_mesh_Size)),
                                                                          ShiftX,
                                                                          ShiftY,
                                                                          ShiftX + cover,
                                                                          ShiftY + t_wall - 2 * cover))

                ShiftY = -t_wall / 2. + cover
                ShiftX = -l_wall / 2. + left_BE - cover
                Fibers.append(OpenSeesAPI.Section.FiberSection.Patch.Rect(ConcUnconfined,
                                                                          int(np.ceil((cover) / max_mesh_Size)),
                                                                          1, #int(np.ceil((t_wall - 2 * cover) / max_mesh_Size)),
                                                                          ShiftX,
                                                                          ShiftY,
                                                                          ShiftX + cover,
                                                                          ShiftY + t_wall - 2 * cover))
            if CustomSection.right_boundary != 0:
                ShiftY = -t_wall / 2.
                ShiftX = l_wall / 2. - right_BE
                Fibers.append(OpenSeesAPI.Section.FiberSection.Patch.Rect(ConcUnconfined,
                                                                          int(np.ceil((right_BE) / max_mesh_Size)),
                                                                          1, #int(np.ceil((cover) / max_mesh_Size)),
                                                                          ShiftX,
                                                                          ShiftY,
                                                                          ShiftX + right_BE,
                                                                          ShiftY + cover))
                ShiftY = -t_wall / 2. + t_wall - cover
                ShiftX = l_wall / 2. - right_BE
                Fibers.append(OpenSeesAPI.Section.FiberSection.Patch.Rect(ConcUnconfined,
                                                                          int(np.ceil((right_BE) / max_mesh_Size)),
                                                                          1, #int(np.ceil((cover) / max_mesh_Size)),
                                                                          ShiftX,
                                                                          ShiftY,
                                                                          ShiftX + right_BE,
                                                                          ShiftY + cover))

                ShiftY = -t_wall / 2. + cover
                ShiftX = l_wall / 2. - cover
                Fibers.append(OpenSeesAPI.Section.FiberSection.Patch.Rect(ConcUnconfined,
                                                                          int(np.ceil((cover) / max_mesh_Size)),
                                                                          1, #int(np.ceil((t_wall - 2 * cover) / max_mesh_Size)),
                                                                          ShiftX,
                                                                          ShiftY,
                                                                          ShiftX + cover,
                                                                          ShiftY + t_wall - 2 * cover))
                ShiftY = -t_wall / 2. + cover
                ShiftX = l_wall / 2. - right_BE
                Fibers.append(OpenSeesAPI.Section.FiberSection.Patch.Rect(ConcUnconfined,
                                                                          int(np.ceil((cover) / max_mesh_Size)),
                                                                          1, #int(np.ceil((t_wall - 2 * cover) / max_mesh_Size)),
                                                                          ShiftX,
                                                                          ShiftY,
                                                                          ShiftX + cover,
                                                                          ShiftY + t_wall - 2 * cover))
        # Unconfined

        # Define Rebar
        def DefineRebar(l_x, l_y, no_bar_x, no_bar_y, ShiftX, ShiftY, Mat=SteelConfined, A=A_bar):
            for a in range(int(no_bar_y)):
                spacing_y = (l_y - 2.*cover) / float(no_bar_y-1)
                Fibers.append(OpenSeesAPI.Section.FiberSection.Layer.Straight(Mat,
                                                                              no_bar_x,
                                                                              A,
                                                                              ShiftX + cover,
                                                                              ShiftY + cover + spacing_y*(a),
                                                                              ShiftX + l_x - cover,
                                                                              ShiftY + cover + spacing_y*(a)
                                                                              ))

        def DefineRebarUsingLayout(ReinfLayout, l_x, l_y, ShiftX, ShiftY, Mat=SteelConfined, A=A_bar):
            spacing_x = (l_x - 2. * cover) / (len(ReinfLayout) - 1)
            for a in range(len(ReinfLayout)):
                Fibers.append(OpenSeesAPI.Section.FiberSection.Layer.Straight(Mat,
                                                                              ReinfLayout[a],
                                                                              A,
                                                                              ShiftX + cover + spacing_x*(a),
                                                                              ShiftY + cover,
                                                                              ShiftX + cover + spacing_x*(a),
                                                                              ShiftY + l_y - cover
                                                                              ))

        if True: # (Planar Wall with Boundary Element)
            if CustomSection.web_rho != 0.:
                bar_spacing_web = 2. * Aweb / CustomSection.web_rho * 100 / t_wall * WallThicknessMultipler # assume two layers
            else:
                bar_spacing_web = None

            #Web Rebar
            if bar_spacing_web is not None:
                ShiftY = -t_wall / 2.
                ShiftX = -l_wall/2. + left_BE
                DefineRebar(l_wall - right_BE - left_BE, t_wall,
                            np.ceil((l_wall - 2 * cover - right_BE - left_BE) / float(bar_spacing_web)),
                            2., ShiftX, ShiftY, Mat=SteelUnconfined, A=Aweb)

            if CustomSection.left_boundary != 0:
                ### Constant rho in BE and Web, as per Pugh's assumptions
                # Left BE
                ShiftY = -t_wall / 2.
                ShiftX = -l_wall / 2.
                DefineRebarUsingLayout(CustomSection.left_reinf_layout, left_BE, t_wall,
                             ShiftX, ShiftY)

            if CustomSection.right_boundary != 0:
                # Right BE
                ShiftY = -t_wall / 2.
                ShiftX = +l_wall / 2. - right_BE
                DefineRebarUsingLayout(CustomSection.right_reinf_layout, right_BE, t_wall,
                                       ShiftX, ShiftY)

        FiberSection = OpenSeesAPI.Section.FiberSection(OSDatabase.GetFreeMaterialId(1, 1),
                                                                        Fibers, Notes='Core Wall Fiber Section')

        OSDatabase.AddMaterial(FiberSection)

        Section = OpenSeesAPI.Section.Aggregator(OSDatabase.GetFreeMaterialId(1, 1),
                                                                        [ShearElasticX],
                                                                        ['Vy'],
                                                                        FiberSection,
                                                                        _Notes='Core Wall Aggregator For Fiber Section')

        OSDatabase.AddMaterial(Section)

        Sections.append(Section)

    return Sections

def ComputeJoinedCustomPlanarWallFiberSection(OSDatabase, CustomSectionA, CustomSectionB, CouplingBeamLength, cover,
                                              height, NoOfIntPoints=5, max_mesh_Size=3,
                                              Elastic=True, RegularizeSteel=True, RegularizeFracture=False,
                                              CrushingStrengthRatio=0.2, GfccOGfc = 2.5, ConcreteMaterialModel='Concrete02', ConfinementModel='SaatRazvi',
                                              SteelMaterialModel='Steel02', SteelUltimateStrainTension=0.2,
                                              SteelPostYieldStiffness=0.00784, GfcOfpc=2., UseForceBased=True,
                                              UnconfinedBeta=0.01, Regularized = True , WallThicknessMultipler = 1.0):
    import numpy as np

    bar_size = CustomSectionA.boundary_bar_size
    fpc = CustomSectionA.fce
    fy = CustomSectionA.fye
    fu = CustomSectionA.fu
    t_wall = CustomSectionA.t_w * WallThicknessMultipler
    l_wall = CustomSectionA.l_w + CustomSectionB.l_w + CouplingBeamLength

    if bar_size is None:
        A_bar = 0
    else:
        A_bar = np.pi * (bar_size / 8. / 2.) ** 2.

    L_ip = GetL_IP(NoOfIntPoints, ForceBased=UseForceBased)

    ConcConfinedAllIP = []
    SteelConfinedAllIP = []
    ConcUnconfinedAllIP = []
    SteelUnconfinedAllIP = []

    import OSMaterialHelper as OSMat
    for i in range(len(L_ip)):
        if Elastic:
            Ec = 0.5 * 57. * (fpc * 1000) ** .5
            Es = 29000*1e-13
            ElasticConcrete = OSDatabase.AddMaterial(
            OpenSeesAPI.Model.Element.Material.UniaxialMaterial.Elastic(OSDatabase.GetFreeMaterialId(1, 1), Ec,
                                                                        _Notes='Concrete Testing'))
            ElasticSteel = OSDatabase.AddMaterial(
            OpenSeesAPI.Model.Element.Material.UniaxialMaterial.Elastic(OSDatabase.GetFreeMaterialId(1, 1), Es,
                                                                        _Notes='Steel Testing'))
            ConcConfinedAllIP.append(ElasticConcrete)
            SteelConfinedAllIP.append(ElasticSteel)
            ConcUnconfinedAllIP.append(ElasticConcrete)
            SteelUnconfinedAllIP.append(ElasticSteel)

            # Shear
            # Define Elastic Spring for Shear
            Ec = 57. * (fpc * 1000.) ** .5
            Mu = 0.2  # 15
            Gc = Ec / (2. * (1. + Mu))
            Geff = Gc  # 0.1 * Gc  # assume 10% of elastic shear modulus per Oyen (2006)
            ks = 5. / 6.  # reduction factor for shear area per Pugh (2012)

            Ax = (l_wall)*t_wall

            EgX = Geff * ks * Ax

            ShearElasticX = OSDatabase.AddMaterial(
                OpenSeesAPI.Model.Element.Material.UniaxialMaterial.Elastic(OSDatabase.GetFreeMaterialId(1, 1), EgX,
                                                                            _Notes='Shear Spring for Core Wall'))

        else:
            if bar_size is None:
                Abar = 0.
            else:
                Abar = np.pi * (bar_size / 8. / 2.) ** 2.

            Aweb = np.pi * (CustomSectionA.bar_size_web / 8. / 2.) ** 2.

            if CustomSectionA.right_boundary is not None:
                bar_spacing = (CustomSectionA.right_boundary - 2 * cover) / (len(CustomSectionA.right_reinf_layout) - 1.)
                tie_spacing = CustomSectionA.boundary_tie_spacing
            else:
                bar_spacing = 2. * Aweb / CustomSectionA.web_rho * 100 / t_wall * WallThicknessMultipler # assume two layers
                tie_spacing = CustomSectionA.boundary_tie_spacing

            tie_size = CustomSectionA.boundary_tie_bar # pugh assumes #4 only


            Et = 57000. * (fpc * 1000) ** .5
            ftu = 4.0 * (fpc * 1000) ** .5
            Gfrac = 5.07e-4 * 1000.
            Lt = L_ip[i] * height
            DeltaE = 2. * Gfrac / ftu / Lt
            Ets = ftu / DeltaE
            if not RegularizeFracture:
                EtsOEt = 0.05 #Ets/Et
            else:
                EtsOEt = Ets/Et

            # Unonfined Concrete
            ConcUnconfinedAllIP.append(GetPughUnconfinedConcreteMaterial(OSDatabase, fpc,
                                                                         L_ip[i] * height,
                                                                         EtsOEt=EtsOEt, fpuOfpc = UnconfinedBeta,
                                                                         GfcOfpc=GfcOfpc, Regularized = Regularized))

            ecu_Unconfined = ConcUnconfinedAllIP[-1]._Materials[0]._ecu

            # Assume 4 for Minimum Reinforc.
            ANo4 = np.pi * (CustomSectionA.bar_size_web / 8. / 2.) ** 2.
            rho_min = 0.0025
            RS_spacing = 2. * ANo4 / t_wall / rho_min * WallThicknessMultipler

            if RegularizeSteel:
                SteelUnconfinedAllIP.append(GetPughSteel02Material(OSDatabase, fy, fu, L_ip[i] * height, ecu_Unconfined,
                                                                   eps_ult_exp=SteelUltimateStrainTension,
                                                                   SteelPostYieldStiffness=SteelPostYieldStiffness))
            else:
                if SteelMaterialModel == 'Steel02':
                    SteelUnconfinedAllIP.append(
                        GetPughSteel02MaterialUnregularized(OSDatabase, fy, fu, ecu_Unconfined, eps_ult_exp=SteelUltimateStrainTension,
                                                                   SteelPostYieldStiffness=SteelPostYieldStiffness))
                elif SteelMaterialModel == 'Steel02Fatigue':
                    SteelUnconfinedAllIP.append(
                        GetSteel02FatigueMaterialUnregularized(OSDatabase, fy, fu, ecu_Unconfined,
                                                            eps_ult_exp=SteelUltimateStrainTension,
                                                            SteelPostYieldStiffness=SteelPostYieldStiffness))
                elif SteelMaterialModel == 'ReinforcingSteel':
                        SteelUnconfinedAllIP.append(
                            GetReinforcingSteelMaterialUnregularized(OSDatabase, fy, fu, SteelUltimateStrainTension, RS_spacing, CustomSectionA.bar_size_web,
                                                                   SteelPostYieldStiffness=SteelPostYieldStiffness))
                elif SteelMaterialModel == 'ReinforcingSteelWithDMBuckling':
                        SteelUnconfinedAllIP.append(
                            GetReinforcingSteelMaterialUnregularized(OSDatabase, fy, fu, SteelUltimateStrainTension, RS_spacing, CustomSectionA.bar_size_web,
                                                                     withBuckling=True,
                                                                   SteelPostYieldStiffness=SteelPostYieldStiffness))
                elif SteelMaterialModel == 'ReinforcingSteelWithDMBucklingAndFatigue':
                    SteelUnconfinedAllIP.append(
                        GetReinforcingSteelMaterialUnregularized(OSDatabase, fy, fu, SteelUltimateStrainTension, RS_spacing, CustomSectionA.bar_size_web,
                                                                 withBuckling=True, withFatigue=True,
                                                                   SteelPostYieldStiffness=SteelPostYieldStiffness))
                elif SteelMaterialModel == 'ReinforcingSteelWithFatigue':
                    SteelUnconfinedAllIP.append(
                        GetReinforcingSteelMaterialUnregularized(OSDatabase, fy, fu, SteelUltimateStrainTension, RS_spacing, CustomSectionA.bar_size_web,
                                                                 withBuckling=False, withFatigue=True,
                                                                   SteelPostYieldStiffness=SteelPostYieldStiffness))
                elif SteelMaterialModel == 'Pinching4':
                    SteelUnconfinedAllIP.append(
                        GetPinching4MaterialUnregularized(OSDatabase, fy, fu, SteelUltimateStrainTension, RS_spacing, CustomSectionA.bar_size_web, -1.0,
                                                                   SteelPostYieldStiffness=SteelPostYieldStiffness))
                elif SteelMaterialModel == 'Steel02LB':
                    SteelUnconfinedAllIP.append(
                        GetSteel02LBMaterialUnregularized(OSDatabase, fy, fu, SteelUltimateStrainTension, RS_spacing, CustomSectionA.bar_size_web, -1.0,
                                                                   SteelPostYieldStiffness=SteelPostYieldStiffness))
                elif SteelMaterialModel == 'Steel02WithFatigue':
                    SteelUnconfinedAllIP.append(
                        GetPughSteel02MaterialUnregularizedWithFatigue(OSDatabase, fy, fu, ecu_Confined,
                                                                       eps_ult_exp=SteelUltimateStrainTension,
                                                                       SteelPostYieldStiffness=SteelPostYieldStiffness))
                elif SteelMaterialModel == 'Hysteretic':
                    SteelUnconfinedAllIP.append(
                        GetHystereticMaterialUnregularized(OSDatabase, fy, fu, SteelUltimateStrainTension, RS_spacing, CustomSectionA.bar_size_web, -1.0,
                                                                   SteelPostYieldStiffness=SteelPostYieldStiffness))

            if tie_size == 0 or tie_size is None:
                A_tie = 0.0001
                tie_size = 0.0001
                ConcConfinedAllIP.append(ConcUnconfinedAllIP[-1])
                SteelConfinedAllIP.append(SteelUnconfinedAllIP[-1])

            else:
                A_tie = np.pi * (tie_size / 8. / 2.) ** 2.

                # Confined Concrete
                if ConfinementModel == 'SaatRazvi':
                    ConcConfinedAllIP.append(GetPughConfinedConcreteMaterial(OSDatabase, fpc,
                                                                           CustomSectionA.left_boundary  - 2 * cover,
                                                                           t_wall-2*cover,
                                                                           bar_spacing, t_wall-2*cover,
                                                                           tie_spacing, tie_spacing,
                                                                           fy, A_tie,
                                                                           CustomSectionA.boundary_tie_x_no,
                                                                           CustomSectionA.boundary_tie_y_no,
                                                                           L_ip[i] * height, GfccOGfc=GfccOGfc,
                                                                           EtsOEt=EtsOEt,
                                                                             fpuOfpc=CrushingStrengthRatio,
                                                                             GfcOfpc=GfcOfpc, Regularized = Regularized,
                                                                             ConcreteMaterialModel=ConcreteMaterialModel))
                elif ConfinementModel == 'Mander':
                    if CustomSectionA.left_boundary != 0:
                       rho_l = np.sum(CustomSectionA.left_reinf_layout) * A_bar / t_wall / CustomSectionA.left_boundary
                    else:
                        rho_l = 0.02
                    ConcConfinedAllIP.append(GetManderConfinedConcreteMaterial(OSDatabase, fpc,
                                                                           CustomSectionA.left_boundary - 2 * cover,
                                                                           t_wall-2*cover, rho_l,
                                                                           bar_spacing, t_wall-2*cover,
                                                                           tie_spacing, tie_spacing,
                                                                           fy, A_tie,
                                                                           CustomSectionA.boundary_tie_x_no,
                                                                           CustomSectionA.boundary_tie_y_no,
                                                                           L_ip[i] * height, GfccOGfc=GfccOGfc,
                                                                           EtsOEt=EtsOEt,
                                                                             fpuOfpc=CrushingStrengthRatio,
                                                                               GfcOfpc=GfcOfpc, ConcreteMaterialModel=ConcreteMaterialModel))
                elif ConfinementModel == 'Richart':
                    if CustomSectionA.left_boundary != 0:
                       rho_l = np.sum(CustomSectionA.left_reinf_layout) * A_bar / t_wall / CustomSectionA.left_boundary
                    else:
                        rho_l = 0.02
                    ConcConfinedAllIP.append(GetRichartConfinedConcreteMaterial(OSDatabase, fpc,
                                                                           CustomSectionA.left_boundary - 2 * cover,
                                                                           t_wall-2*cover, rho_l,
                                                                           bar_spacing, t_wall-2*cover,
                                                                           tie_spacing, tie_spacing,
                                                                           fy, A_tie,
                                                                           CustomSectionA.boundary_tie_x_no,
                                                                           CustomSectionA.boundary_tie_y_no,
                                                                           L_ip[i] * height, GfccOGfc=GfccOGfc,
                                                                           EtsOEt=EtsOEt,
                                                                             fpuOfpc=CrushingStrengthRatio,
                                                                               GfcOfpc=GfcOfpc, ConcreteMaterialModel=ConcreteMaterialModel))
                elif ConfinementModel == 'Unconfined':
                    ConcConfinedAllIP.append(ConcUnconfinedAllIP[-1])

                ecu_Confined = ConcConfinedAllIP[-1]._Materials[0]._ecu

                if tie_spacing == None or tie_spacing == 0:
                    tie_spacing = 1e10

                if RegularizeSteel:
                    SteelConfinedAllIP.append(
                        GetPughSteel02Material(OSDatabase, fy, fu, L_ip[i] * height, ecu_Confined, eps_ult_exp=SteelUltimateStrainTension,
                                                                   SteelPostYieldStiffness=SteelPostYieldStiffness))
                else:
                    if SteelMaterialModel == 'Steel02':
                        SteelConfinedAllIP.append(
                            GetPughSteel02MaterialUnregularized(OSDatabase, fy, fu, ecu_Confined, eps_ult_exp=SteelUltimateStrainTension,
                                                                   SteelPostYieldStiffness=SteelPostYieldStiffness))
                    elif SteelMaterialModel == 'Steel02Fatigue':
                        SteelConfinedAllIP.append(
                            GetSteel02FatigueMaterialUnregularized(OSDatabase, fy, fu, ecu_Confined,
                                                                   eps_ult_exp=SteelUltimateStrainTension,
                                                                   SteelPostYieldStiffness=SteelPostYieldStiffness))
                    elif SteelMaterialModel == 'ReinforcingSteel':
                        SteelConfinedAllIP.append(
                            GetReinforcingSteelMaterialUnregularized(OSDatabase, fy, fu, SteelUltimateStrainTension, tie_spacing, bar_size,
                                                                     eps_ult_comp=ecu_Confined,
                                                                   SteelPostYieldStiffness=SteelPostYieldStiffness))
                    elif SteelMaterialModel == 'ReinforcingSteelWithDMBuckling':
                        SteelConfinedAllIP.append(
                            GetReinforcingSteelMaterialUnregularized(OSDatabase, fy, fu, SteelUltimateStrainTension, tie_spacing, bar_size,
                                                                     withBuckling=True,
                                                                   SteelPostYieldStiffness=SteelPostYieldStiffness))
                    elif SteelMaterialModel == 'ReinforcingSteelWithDMBucklingAndFatigue':
                        SteelConfinedAllIP.append(
                            GetReinforcingSteelMaterialUnregularized(OSDatabase, fy, fu, SteelUltimateStrainTension, tie_spacing, bar_size,
                                                                     withBuckling=True, withFatigue=True,
                                                                   SteelPostYieldStiffness=SteelPostYieldStiffness))
                    elif SteelMaterialModel == 'ReinforcingSteelWithFatigue':
                        SteelConfinedAllIP.append(
                            GetReinforcingSteelMaterialUnregularized(OSDatabase, fy, fu, SteelUltimateStrainTension, tie_spacing, bar_size,
                                                                     withBuckling=False, withFatigue=True,
                                                                   SteelPostYieldStiffness=SteelPostYieldStiffness))
                    elif SteelMaterialModel == 'Pinching4':
                        SteelConfinedAllIP.append(
                            GetPinching4MaterialUnregularized(OSDatabase, fy, fu, SteelUltimateStrainTension, tie_spacing, bar_size, -1.0,
                                                                   SteelPostYieldStiffness=SteelPostYieldStiffness))
                    elif SteelMaterialModel == 'Steel02LB':
                        SteelConfinedAllIP.append(
                            GetSteel02LBMaterialUnregularized(OSDatabase, fy, fu, SteelUltimateStrainTension, tie_spacing, bar_size, -1.0,
                                                                   SteelPostYieldStiffness=SteelPostYieldStiffness))
                    elif SteelMaterialModel == 'Steel02WithFatigue':
                        SteelConfinedAllIP.append(
                            GetPughSteel02MaterialUnregularizedWithFatigue(OSDatabase, fy, fu, ecu_Confined, eps_ult_exp=SteelUltimateStrainTension,
                                                                   SteelPostYieldStiffness=SteelPostYieldStiffness))
                    elif SteelMaterialModel == 'Hysteretic':
                        SteelConfinedAllIP.append(
                            GetHystereticMaterialUnregularized(OSDatabase, fy, fu, SteelUltimateStrainTension,
                                                               tie_spacing,
                                                               bar_size, -1.0,
                                                               SteelPostYieldStiffness=SteelPostYieldStiffness))

            # Shear
            # Define Elastic Spring for Shear
            Ec = 57. * (fpc * 1000.) ** .5
            Mu = 0.15  # 15
            Gc = Ec / (2. * (1. + Mu))
            Geff = 0.1 * Gc  # assume 10% of elastic shear modulus per Oyen (2006)
            ks = 5. / 6.  # reduction factor for shear area per Pugh (2012)

            Av = l_wall * t_wall

            GcOEc = Gc/Ec

            EgX = 0.4*Ec*Av*ks

            ShearElasticX = OSDatabase.AddMaterial(
                OpenSeesAPI.Model.Element.Material.UniaxialMaterial.Elastic(OSDatabase.GetFreeMaterialId(1, 1), EgX,
                                                                            _Notes='Shear Spring for Core Wall'))

    Sections = []

    for j in range(len(L_ip)):
        ConcConfined = ConcConfinedAllIP[j]
        SteelConfined = SteelConfinedAllIP[j]
        ConcUnconfined = ConcUnconfinedAllIP[j]
        SteelUnconfined = SteelUnconfinedAllIP[j]

        Fibers = []

        # Confined Parts
        if True:
            LW_left_BE = CustomSectionA.left_boundary
            LW_right_BE = CustomSectionA.right_boundary
            LW_web = CustomSectionA.l_w - LW_left_BE - LW_right_BE

            RW_left_BE = CustomSectionB.left_boundary
            RW_right_BE = CustomSectionB.right_boundary
            RW_web = CustomSectionB.l_w - RW_left_BE - RW_right_BE

            MW_web = l_wall - CustomSectionB.l_w - CustomSectionA.l_w

            ### Left Wall
            ## Web
            ShiftY = -t_wall / 2.
            ShiftX = -l_wall / 2. + LW_left_BE
            Fibers.append(OpenSeesAPI.Section.FiberSection.Patch.Rect(ConcUnconfined,
                                                                      int(np.ceil((LW_web) / max_mesh_Size)),
                                                                      1, #int(np.ceil((t_wall) / max_mesh_Size)),
                                                                      ShiftX,
                                                                      ShiftY,
                                                                      ShiftX + LW_web,
                                                                      ShiftY + t_wall ))

            # Left BE
            ShiftY = -t_wall / 2.
            ShiftX = -l_wall / 2.
            if LW_left_BE != 0:
                Fibers.append(OpenSeesAPI.Section.FiberSection.Patch.Rect(ConcConfined,
                                                                          int(np.ceil((LW_left_BE) / max_mesh_Size)),
                                                                          1,  #int(np.ceil((t_wall) / max_mesh_Size)),
                                                                          ShiftX,
                                                                          ShiftY,
                                                                          ShiftX + LW_left_BE,
                                                                          ShiftY + t_wall ))

            # Right BE
            ShiftY = -t_wall / 2.
            ShiftX = -l_wall / 2. + LW_left_BE + LW_web
            if LW_right_BE != 0:
                Fibers.append(OpenSeesAPI.Section.FiberSection.Patch.Rect(ConcConfined,
                                                                          int(np.ceil((LW_right_BE) / max_mesh_Size)),
                                                                          1,  #int(np.ceil((t_wall) / max_mesh_Size)),
                                                                          ShiftX,
                                                                          ShiftY,
                                                                          ShiftX + LW_right_BE,
                                                                          ShiftY + t_wall ))

            ### Middle Wall
            ShiftY = -t_wall / 2.
            ShiftX = -MW_web / 2.
            Fibers.append(OpenSeesAPI.Section.FiberSection.Patch.Rect(ConcUnconfined,
                                                                      int(np.ceil((MW_web) / max_mesh_Size)),
                                                                      1,  #int(np.ceil((t_wall) / max_mesh_Size)),
                                                                      ShiftX,
                                                                      ShiftY,
                                                                      ShiftX + MW_web,
                                                                      ShiftY + t_wall ))

            ### Right Wall
            ## Web
            ShiftY = -t_wall / 2.
            ShiftX = MW_web / 2. + RW_left_BE
            Fibers.append(OpenSeesAPI.Section.FiberSection.Patch.Rect(ConcUnconfined,
                                                                      int(np.ceil((RW_web) / max_mesh_Size)),
                                                                      1,  #int(np.ceil((t_wall) / max_mesh_Size)),
                                                                      ShiftX,
                                                                      ShiftY,
                                                                      ShiftX + RW_web,
                                                                      ShiftY + t_wall))

            # Left BE
            ShiftY = -t_wall / 2.
            ShiftX = MW_web / 2.
            if RW_left_BE != 0:
                Fibers.append(OpenSeesAPI.Section.FiberSection.Patch.Rect(ConcConfined,
                                                                          int(np.ceil((RW_left_BE) / max_mesh_Size)),
                                                                          1,  #int(np.ceil((t_wall) / max_mesh_Size)),
                                                                          ShiftX,
                                                                          ShiftY,
                                                                          ShiftX + RW_left_BE,
                                                                          ShiftY + t_wall))

            # Right BE
            ShiftY = -t_wall / 2.
            ShiftX = l_wall / 2. - RW_right_BE
            if RW_right_BE != 0:
                Fibers.append(OpenSeesAPI.Section.FiberSection.Patch.Rect(ConcConfined,
                                                                          int(np.ceil((RW_right_BE) / max_mesh_Size)),
                                                                          1,  #int(np.ceil((t_wall) / max_mesh_Size)),
                                                                          ShiftX,
                                                                          ShiftY,
                                                                          ShiftX + RW_right_BE,
                                                                          ShiftY + t_wall))

        # Define Rebar
        def DefineRebar(l_x, l_y, no_bar_x, no_bar_y, ShiftX, ShiftY, Mat=SteelConfined, A=A_bar):
            for a in range(int(no_bar_y)):
                spacing_y = (l_y - 2.*cover) / float(no_bar_y-1)
                Fibers.append(OpenSeesAPI.Section.FiberSection.Layer.Straight(Mat,
                                                                              no_bar_x,
                                                                              A,
                                                                              ShiftX + cover,
                                                                              ShiftY + cover + spacing_y*(a),
                                                                              ShiftX + l_x - cover,
                                                                              ShiftY + cover + spacing_y*(a)
                                                                              ))

        def DefineRebarUsingLayout(ReinfLayout, l_x, l_y, ShiftX, ShiftY, Mat=SteelConfined, A=A_bar):
            spacing_x = (l_x - 2. * cover) / (len(ReinfLayout) - 1)
            for a in range(len(ReinfLayout)):
                Fibers.append(OpenSeesAPI.Section.FiberSection.Layer.Straight(Mat,
                                                                              ReinfLayout[a],
                                                                              A,
                                                                              ShiftX + cover + spacing_x*(a),
                                                                              ShiftY + cover,
                                                                              ShiftX + cover + spacing_x*(a),
                                                                              ShiftY + l_y - cover
                                                                              ))

        if True: # (Planar Wall with Boundary Element)
            bar_spacing_web = 2. * Aweb / CustomSectionA.web_rho * 100 / t_wall * WallThicknessMultipler # assume two layers

            ### Left Wall
            #Web Rebar
            ShiftY = -t_wall / 2.
            ShiftX = -l_wall/2. + LW_left_BE
            DefineRebar(LW_web, t_wall,
                        np.ceil((LW_web - 2 * cover) / float(bar_spacing_web)),
                        2., ShiftX, ShiftY, Mat=SteelUnconfined, A=Aweb)

            ### Constant rho in BE and Web, as per Pugh's assumptions
            # Left BE
            ShiftY = -t_wall / 2.
            ShiftX = -l_wall / 2.
            DefineRebarUsingLayout(CustomSectionA.left_reinf_layout, LW_left_BE, t_wall,
                         ShiftX, ShiftY)

            # Right BE
            ShiftY = -t_wall / 2.
            ShiftX = -l_wall / 2. + LW_left_BE + LW_web
            DefineRebarUsingLayout(CustomSectionA.right_reinf_layout, LW_right_BE, t_wall,
                                   ShiftX, ShiftY)

            ### Middle Wall
            ShiftY = -t_wall / 2.
            ShiftX = -MW_web/2.
            DefineRebar(MW_web, t_wall,
                        np.ceil((MW_web - 2 * cover) / float(bar_spacing_web)),
                        2., ShiftX, ShiftY, Mat=SteelUnconfined, A=Aweb)

            bar_spacing_web = 2. * Aweb / CustomSectionA.web_rho * 100 / t_wall  # assume two layers
            ### Right Wall
            #Web Rebar
            ShiftY = -t_wall / 2.
            ShiftX = MW_web / 2. + RW_left_BE
            DefineRebar(RW_web, t_wall,
                        np.ceil((RW_web - 2 * cover) / float(bar_spacing_web)),
                        2., ShiftX, ShiftY, Mat=SteelUnconfined, A=Aweb)

            ### Constant rho in BE and Web, as per Pugh's assumptions
            # Left BE
            ShiftY = -t_wall / 2.
            ShiftX = MW_web / 2.
            DefineRebarUsingLayout(CustomSectionB.left_reinf_layout, RW_left_BE, t_wall,
                         ShiftX, ShiftY)

            # Right BE
            ShiftY = -t_wall / 2.
            ShiftX = l_wall / 2. - RW_right_BE
            DefineRebarUsingLayout(CustomSectionB.right_reinf_layout, RW_right_BE, t_wall,
                                   ShiftX, ShiftY)

        FiberSection = OpenSeesAPI.Section.FiberSection(OSDatabase.GetFreeMaterialId(1, 1),
                                                                        Fibers, Notes='Core Wall Fiber Section')

        OSDatabase.AddMaterial(FiberSection)

        Section = OpenSeesAPI.Section.Aggregator(OSDatabase.GetFreeMaterialId(1, 1),
                                                                        [ShearElasticX],
                                                                        ['Vy'],
                                                                        FiberSection,
                                                                        _Notes='Core Wall Aggregator For Fiber Section')

        OSDatabase.AddMaterial(Section)

        Sections.append(Section)

    return Sections

def ComputeCustomIWallFiberSection(OSDatabase, CustomSection, cover, height, NoOfIntPoints=5, max_mesh_Size=3,
                                        Elastic=True, RegularizeSteel=True, RegularizeFracture=False,
                                   CrushingStrengthRatio=0.2, GfccOGfc = 2.5, ConcreteMaterialModel='Concrete02', ConfinementModel='SaatRazvi',
                                   SteelMaterialModel='Steel02', SteelUltimateStrainTension=0.2,
                                   SteelPostYieldStiffness=0.00784, GfcOfpc=2., UseForceBased=True,
                                   UnconfinedBeta=0.01, Regularized = True):
    import numpy as np

    bar_size = CustomSection.bar_size
    fpc = CustomSection.fce
    fy = CustomSection.fye
    fu = CustomSection.fu

    t_wall = CustomSection.t_w
    l_wall = CustomSection.l_w
    b_wall = CustomSection.b_w

    rho = CustomSection.rho

    A_bar = np.pi * (bar_size / 8. / 2.) ** 2.

    L_ip = GetL_IP(NoOfIntPoints, ForceBased=UseForceBased)

    ConcConfinedAllIP = []
    SteelConfinedAllIP = []
    ConcUnconfinedAllIP = []
    SteelUnconfinedAllIP = []

    for i in range(len(L_ip)):
        Abar = np.pi * (bar_size / 8. / 2.) ** 2.

        if Elastic:
            Ec = 0.5 * 57. * (fpc * 1000) ** .5
            Es = 29000*1e-13
            ElasticConcrete = OSDatabase.AddMaterial(
            OpenSeesAPI.Model.Element.Material.UniaxialMaterial.Elastic(OSDatabase.GetFreeMaterialId(1, 1), Ec,
                                                                        _Notes='Concrete Testing'))
            ElasticSteel = OSDatabase.AddMaterial(
            OpenSeesAPI.Model.Element.Material.UniaxialMaterial.Elastic(OSDatabase.GetFreeMaterialId(1, 1), Es,
                                                                        _Notes='Steel Testing'))
            ConcConfinedAllIP.append(ElasticConcrete)
            SteelConfinedAllIP.append(ElasticSteel)
            ConcUnconfinedAllIP.append(ElasticConcrete)
            SteelUnconfinedAllIP.append(ElasticSteel)

            # Shear
            # Define Elastic Spring for Shear
            Ec = 57. * (fpc * 1000.) ** .5
            Mu = 0.2  # 15
            Gc = Ec / (2. * (1. + Mu))
            Geff = Gc  # 0.1 * Gc  # assume 10% of elastic shear modulus per Oyen (2006)
            ks = 5. / 6.  # reduction factor for shear area per Pugh (2012)

            Ax = (l_wall)*t_wall

            EgX = Geff * ks * Ax

            ShearElasticX = OSDatabase.AddMaterial(
                OpenSeesAPI.Model.Element.Material.UniaxialMaterial.Elastic(OSDatabase.GetFreeMaterialId(1, 1), EgX,
                                                                            _Notes='Shear Spring for Core Wall'))

        else:
            if t_wall - 2 * cover >= 14.:
                rho_mod = (rho * b_wall * t_wall / 100. - A_bar * 2. ) / b_wall / t_wall # Accounts for the fact that 3 bars are at the ends, note this breaks down if the thickness is more than 3*14
                bar_spacing = 2. * A_bar / rho_mod / t_wall
                EdgeBars = 3
            else:
                bar_spacing = 2. * A_bar / rho * 100 / t_wall  # assume two layers
                EdgeBars = 2

            AssumeUnconfined = False
            if CustomSection.tie_vertical_spacing is not None: # If tie spacing not specified
                tie_spacing = CustomSection.tie_vertical_spacing
                tie_size = CustomSection.tie_bar_size # pugh assumes #4 only
            else:
                AssumeUnconfined = True
                tie_spacing = GetTieSpacing(t_wall, bar_size / 8., bar_spacing)
                tie_size = 4 # assume somethi=ng so that it does not crash

            Et = 57000. * (fpc * 1000) ** .5
            ftu = 4.0 * (fpc * 1000) ** .5
            Gfrac = 5.07e-4 * 1000.
            Lt = L_ip[i] * height
            DeltaE = 2. * Gfrac / ftu / Lt
            Ets = ftu / DeltaE
            if not RegularizeFracture:
                EtsOEt = 0.05 #Ets/Et
            else:
                EtsOEt = Ets/Et

            # Unonfined Concrete
            ConcUnconfinedAllIP.append(GetPughUnconfinedConcreteMaterial(OSDatabase, fpc,
                                                                         L_ip[i] * height,
                                                                         EtsOEt=EtsOEt,
                                                                         fpuOfpc = UnconfinedBeta,
                                                                         GfcOfpc=GfcOfpc , Regularized = Regularized))

            # ecu_Unconfined = ConcUnconfinedAllIP[-1]._ecu
            ecu_Unconfined = ConcUnconfinedAllIP[-1]._Materials[0]._ecu

            if RegularizeSteel:
                SteelUnconfinedAllIP.append(GetPughSteel02Material(OSDatabase, fy, fu, L_ip[i] * height, ecu_Unconfined, eps_ult_exp=SteelUltimateStrainTension,
                                                                   SteelPostYieldStiffness=SteelPostYieldStiffness))
            else:
                if SteelMaterialModel == 'Steel02':
                    SteelUnconfinedAllIP.append(
                        GetPughSteel02MaterialUnregularized(OSDatabase, fy, fu, ecu_Unconfined, eps_ult_exp=SteelUltimateStrainTension,
                                                                   SteelPostYieldStiffness=SteelPostYieldStiffness))
                elif SteelMaterialModel == 'Steel02Fatigue':
                    SteelUnconfinedAllIP.append(
                        GetSteel02FatigueMaterialUnregularized(OSDatabase, fy, fu, ecu_Unconfined,
                                                            eps_ult_exp=SteelUltimateStrainTension,
                                                            SteelPostYieldStiffness=SteelPostYieldStiffness))
                elif SteelMaterialModel == 'ReinforcingSteel':
                        SteelUnconfinedAllIP.append(
                            GetReinforcingSteelMaterialUnregularized(OSDatabase, fy, fu, SteelUltimateStrainTension, bar_spacing, bar_size,
                                                                   SteelPostYieldStiffness=SteelPostYieldStiffness))
                elif SteelMaterialModel == 'ReinforcingSteelWithDMBuckling':
                        SteelUnconfinedAllIP.append(
                            GetReinforcingSteelMaterialUnregularized(OSDatabase, fy, fu, SteelUltimateStrainTension, bar_spacing, bar_size,
                                                                     withBuckling=True,
                                                                   SteelPostYieldStiffness=SteelPostYieldStiffness))
                elif SteelMaterialModel == 'ReinforcingSteelWithDMBucklingAndFatigue':
                    SteelUnconfinedAllIP.append(
                        GetReinforcingSteelMaterialUnregularized(OSDatabase, fy, fu, SteelUltimateStrainTension, bar_spacing, bar_size,
                                                                 withBuckling=True, withFatigue=True,
                                                                   SteelPostYieldStiffness=SteelPostYieldStiffness))
                elif SteelMaterialModel == 'ReinforcingSteelWithFatigue':
                    SteelUnconfinedAllIP.append(
                        GetReinforcingSteelMaterialUnregularized(OSDatabase, fy, fu, SteelUltimateStrainTension, bar_spacing, bar_size,
                                                                 withBuckling=False, withFatigue=True,
                                                                   SteelPostYieldStiffness=SteelPostYieldStiffness))
                elif SteelMaterialModel == 'Pinching4':
                    SteelUnconfinedAllIP.append(
                        GetPinching4MaterialUnregularized(OSDatabase, fy, fu, SteelUltimateStrainTension, bar_spacing, bar_size, -1.0,
                                                                   SteelPostYieldStiffness=SteelPostYieldStiffness))
                elif SteelMaterialModel == 'Steel02LB':
                    SteelUnconfinedAllIP.append(
                        GetSteel02LBMaterialUnregularized(OSDatabase, fy, fu, SteelUltimateStrainTension, bar_spacing, bar_size, -1.0,
                                                                   SteelPostYieldStiffness=SteelPostYieldStiffness))
                elif SteelMaterialModel == 'Steel02WithFatigue':
                    SteelUnconfinedAllIP.append(
                        GetPughSteel02MaterialUnregularizedWithFatigue(OSDatabase, fy, fu, ecu_Unconfined,
                                                                       eps_ult_exp=SteelUltimateStrainTension,
                                                                       SteelPostYieldStiffness=SteelPostYieldStiffness))
                elif SteelMaterialModel == 'Hysteretic':
                    SteelUnconfinedAllIP.append(
                        GetHystereticMaterialUnregularized(OSDatabase, fy, fu, SteelUltimateStrainTension, bar_spacing, bar_size, -1.0,
                                                                   SteelPostYieldStiffness=SteelPostYieldStiffness))

            A_tie = np.pi * (tie_size / 8. / 2.) ** 2.

            if AssumeUnconfined:
                ConfinementModel = 'Unconfined'

            # Confined Concrete
            if ConfinementModel == 'SaatRazvi':
                ConcConfinedAllIP.append(GetPughConfinedConcreteMaterial(OSDatabase, fpc,
                                                                         l_wall - 2 * cover,
                                                                         t_wall - 2 * cover,
                                                                         bar_spacing, t_wall - 2 * cover,
                                                                         tie_spacing, tie_spacing,
                                                                         fy, A_tie,
                                                                         np.ceil(l_wall / bar_spacing),
                                                                         EdgeBars,
                                                                         L_ip[i] * height, GfccOGfc=GfccOGfc,
                                                                         EtsOEt=EtsOEt,
                                                                         fpuOfpc=CrushingStrengthRatio,
                                                                         GfcOfpc=GfcOfpc, Regularized = Regularized,
                                                                         ConcreteMaterialModel=ConcreteMaterialModel))
            elif ConfinementModel == 'Mander':
                ConcConfinedAllIP.append(GetManderConfinedConcreteMaterial(OSDatabase, fpc,
                                                                           l_wall - 2 * cover,
                                                                           t_wall - 2 * cover, rho,
                                                                           bar_spacing, t_wall - 2 * cover,
                                                                           tie_spacing, tie_spacing,
                                                                           fy, A_tie,
                                                                           np.ceil(l_wall / bar_spacing),
                                                                           EdgeBars,
                                                                           L_ip[i] * height, GfccOGfc=GfccOGfc,
                                                                           EtsOEt=EtsOEt,
                                                                           fpuOfpc=CrushingStrengthRatio,
                                                                           GfcOfpc=GfcOfpc, ConcreteMaterialModel=ConcreteMaterialModel))
            elif ConfinementModel == 'Richart':
                ConcConfinedAllIP.append(GetRichartConfinedConcreteMaterial(OSDatabase, fpc,
                                                                           l_wall - 2 * cover,
                                                                           t_wall - 2 * cover, rho,
                                                                           bar_spacing, t_wall - 2 * cover,
                                                                           tie_spacing, tie_spacing,
                                                                           fy, A_tie,
                                                                           np.ceil(l_wall / bar_spacing),
                                                                           EdgeBars,
                                                                           L_ip[i] * height, GfccOGfc=GfccOGfc,
                                                                           EtsOEt=EtsOEt,
                                                                           fpuOfpc=CrushingStrengthRatio,
                                                                           GfcOfpc=GfcOfpc, ConcreteMaterialModel=ConcreteMaterialModel))
            elif ConfinementModel == 'Unconfined':
                ConcConfinedAllIP.append(ConcUnconfinedAllIP[-1])

            ecu_Confined = ConcConfinedAllIP[-1]._Materials[0]._ecu

            if tie_spacing == None or tie_spacing == 0:
                tie_spacing = 1e10

            if RegularizeSteel:
                SteelConfinedAllIP.append(
                    GetPughSteel02Material(OSDatabase, fy, fu, L_ip[i] * height, ecu_Confined,
                                           eps_ult_exp=SteelUltimateStrainTension))
            else:
                if SteelMaterialModel == 'Steel02':
                    SteelConfinedAllIP.append(
                        GetPughSteel02MaterialUnregularized(OSDatabase, fy, fu, ecu_Confined,
                                                            eps_ult_exp=SteelUltimateStrainTension))
                elif SteelMaterialModel == 'Steel02Fatigue':
                    SteelConfinedAllIP.append(
                        GetSteel02FatigueMaterialUnregularized(OSDatabase, fy, fu, ecu_Confined,
                                                            eps_ult_exp=SteelUltimateStrainTension,
                                                            SteelPostYieldStiffness=SteelPostYieldStiffness))
                elif SteelMaterialModel == 'ReinforcingSteel':
                    SteelConfinedAllIP.append(
                        GetReinforcingSteelMaterialUnregularized(OSDatabase, fy, fu, SteelUltimateStrainTension,
                                                                 tie_spacing, bar_size,
                                                                 eps_ult_comp=ecu_Confined))
                elif SteelMaterialModel == 'ReinforcingSteelWithDMBuckling':
                    SteelConfinedAllIP.append(
                        GetReinforcingSteelMaterialUnregularized(OSDatabase, fy, fu, SteelUltimateStrainTension,
                                                                 tie_spacing, bar_size,
                                                                 withBuckling=True))
                elif SteelMaterialModel == 'ReinforcingSteelWithDMBucklingAndFatigue':
                    SteelConfinedAllIP.append(
                        GetReinforcingSteelMaterialUnregularized(OSDatabase, fy, fu, SteelUltimateStrainTension,
                                                                 tie_spacing, bar_size,
                                                                 withBuckling=True, withFatigue=True))
                elif SteelMaterialModel == 'ReinforcingSteelWithFatigue':
                    SteelConfinedAllIP.append(
                        GetReinforcingSteelMaterialUnregularized(OSDatabase, fy, fu, SteelUltimateStrainTension,
                                                                 tie_spacing, bar_size,
                                                                 withBuckling=False, withFatigue=True))
                elif SteelMaterialModel == 'Pinching4':
                    SteelConfinedAllIP.append(
                        GetPinching4MaterialUnregularized(OSDatabase, fy, fu, SteelUltimateStrainTension, tie_spacing,
                                                          bar_size, -1.0))
                elif SteelMaterialModel == 'Steel02LB':
                    SteelConfinedAllIP.append(
                        GetSteel02LBMaterialUnregularized(OSDatabase, fy, fu, SteelUltimateStrainTension, tie_spacing,
                                                          bar_size, -1.0))
                elif SteelMaterialModel == 'Steel02WithFatigue':
                    SteelConfinedAllIP.append(
                        GetPughSteel02MaterialUnregularizedWithFatigue(OSDatabase, fy, fu, ecu_Confined,
                                                                       eps_ult_exp=SteelUltimateStrainTension,
                                                                       SteelPostYieldStiffness=SteelPostYieldStiffness))
                elif SteelMaterialModel == 'Hysteretic':
                    SteelConfinedAllIP.append(
                        GetHystereticMaterialUnregularized(OSDatabase, fy, fu, SteelUltimateStrainTension, tie_spacing,
                                                           bar_size, -1.0,
                                                           SteelPostYieldStiffness=SteelPostYieldStiffness))

            # Shear
            # Define Elastic Spring for Shear
            Ec = 57. * (fpc * 1000.) ** .5
            Mu = 0.15  # 15
            Gc = Ec / (2. * (1. + Mu))
            Geff = 0.1 * Gc  # assume 10% of elastic shear modulus per Oyen (2006)
            ks = 5. / 6.  # reduction factor for shear area per Pugh (2012)

            Av = l_wall * t_wall

            GcOEc = Gc/Ec

            EgX = 0.4*Ec*Av*ks

            ShearElasticX = OSDatabase.AddMaterial(
                OpenSeesAPI.Model.Element.Material.UniaxialMaterial.Elastic(OSDatabase.GetFreeMaterialId(1, 1), EgX,
                                                                            _Notes='Shear Spring for Core Wall'))

    Sections = []

    for j in range(len(L_ip)):
        ConcConfined = ConcConfinedAllIP[j]
        SteelConfined = SteelConfinedAllIP[j]
        ConcUnconfined = ConcUnconfinedAllIP[j]
        SteelUnconfined = SteelUnconfinedAllIP[j]

        Fibers = []

        # Confined Parts
        if True:
            l_web = CustomSection.l_w - 2. * CustomSection.t_w
            l_wall = CustomSection.l_w
            t_wall = CustomSection.t_w
            b_wall = CustomSection.b_w

            ## Web
            ShiftY = -t_wall / 2.
            ShiftX = -l_web / 2.
            Fibers.append(OpenSeesAPI.Section.FiberSection.Patch.Rect(ConcConfined,
                                                                      int(np.ceil((l_web) / max_mesh_Size)),
                                                                      int(np.ceil((t_wall) / max_mesh_Size)),
                                                                      ShiftX,
                                                                      ShiftY,
                                                                      ShiftX + l_web,
                                                                      ShiftY + t_wall))

            # Left Flange
            ShiftY = -b_wall / 2.
            ShiftX = -l_wall / 2.
            Fibers.append(OpenSeesAPI.Section.FiberSection.Patch.Rect(ConcConfined,
                                                                      int(np.ceil((t_wall) / max_mesh_Size)),
                                                                      int(np.ceil((b_wall) / max_mesh_Size)),
                                                                      ShiftX,
                                                                      ShiftY,
                                                                      ShiftX + t_wall,
                                                                      ShiftY + b_wall))

            # Right Flange
            ShiftY = -b_wall / 2.
            ShiftX = l_wall / 2. - t_wall
            Fibers.append(OpenSeesAPI.Section.FiberSection.Patch.Rect(ConcConfined,
                                                                      int(np.ceil((t_wall) / max_mesh_Size)),
                                                                      int(np.ceil((b_wall) / max_mesh_Size)),
                                                                      ShiftX,
                                                                      ShiftY,
                                                                      ShiftX + t_wall,
                                                                      ShiftY + b_wall))

        # Unconfined
        # Save this for cover...


        # Define Rebar
        def DefineRebar(l_x, l_y, no_bar_x, no_bar_y, ShiftX, ShiftY, Mat=SteelConfined, A=A_bar):
            for a in range(int(no_bar_y)):
                spacing_y = (l_y - 2.*cover) / float(no_bar_y-1)
                Fibers.append(OpenSeesAPI.Section.FiberSection.Layer.Straight(Mat,
                                                                              no_bar_x,
                                                                              A,
                                                                              ShiftX + cover,
                                                                              ShiftY + cover + spacing_y*(a),
                                                                              ShiftX + l_x - cover,
                                                                              ShiftY + cover + spacing_y*(a)
                                                                              ))

        def DefineRebarUsingLayout(ReinfLayout, l_x, l_y, ShiftX, ShiftY, Mat=SteelConfined, A=A_bar):
            spacing_x = (l_x - 2. * cover) / (len(ReinfLayout) - 1)
            for a in range(len(ReinfLayout)):
                Fibers.append(OpenSeesAPI.Section.FiberSection.Layer.Straight(Mat,
                                                                              ReinfLayout[a],
                                                                              A,
                                                                              ShiftX + cover + spacing_x*(a),
                                                                              ShiftY + cover,
                                                                              ShiftX + cover + spacing_x*(a),
                                                                              ShiftY + l_y - cover
                                                                              ))

        if True: # (Planar Wall with Boundary Element)
            # Web Rebar
            ShiftX = -l_web / 2.
            ShiftY = -t_wall / 2.
            DefineRebar(l_web , t_wall,
                        np.ceil((l_wall - 2 * cover) / float(bar_spacing)),
                        2., ShiftX, ShiftY, Mat=SteelConfined, A=Abar)

            # Left Flange Rebar
            ShiftX = -l_wall / 2.
            ShiftY = -b_wall / 2.
            DefineRebar(t_wall , b_wall ,
                        2., np.ceil((b_wall - 2 * cover) / float(bar_spacing)),
                        ShiftX, ShiftY, Mat=SteelConfined, A=Abar)

            # Right Flange Rebar
            ShiftX = +l_web / 2.
            ShiftY = -b_wall / 2.
            DefineRebar(t_wall , b_wall ,
                        2., np.ceil((b_wall - 2 * cover) / float(bar_spacing)),
                        ShiftX, ShiftY, Mat=SteelConfined, A=Abar)

            if EdgeBars == 3:
                ShiftX = -l_wall / 2.
                ShiftY = -b_wall / 2.
                DefineRebar(t_wall, b_wall,
                            1., 2,
                            ShiftX + t_wall / 2. - cover, ShiftY, Mat=SteelConfined, A=Abar)

                ShiftX = +l_web / 2.
                ShiftY = -b_wall / 2.
                DefineRebar(t_wall, b_wall,
                            1., 2,
                            ShiftX + t_wall / 2. - cover, ShiftY, Mat=SteelConfined, A=Abar)

        FiberSection = OpenSeesAPI.Section.FiberSection(OSDatabase.GetFreeMaterialId(1, 1),
                                                                        Fibers, Notes='Core Wall Fiber Section')

        OSDatabase.AddMaterial(FiberSection)

        Section = OpenSeesAPI.Section.Aggregator(OSDatabase.GetFreeMaterialId(1, 1),
                                                                        [ShearElasticX],
                                                                        ['Vy'],
                                                                        FiberSection,
                                                                        _Notes='Core Wall Aggregator For Fiber Section')

        OSDatabase.AddMaterial(Section)
        Sections.append(Section)

    return Sections

def ComputeCustomTWallFiberSection(OSDatabase, CustomSection, cover, height, NoOfIntPoints=5, max_mesh_Size=3,
                                   Elastic=True, RegularizeSteel=True, RegularizeFracture=False,
                                   CrushingStrengthRatio = 0.2, GfccOGfc = 2.5, ConcreteMaterialModel='Concrete02', ConfinementModel='SaatRazvi',
                                   SteelMaterialModel='Steel02', SteelUltimateStrainTension=0.2,
                                   SteelPostYieldStiffness=0.00784, GfcOfpc=2., UseForceBased=True,
                                   UnconfinedBeta=0.01, Regularized = True):
    import numpy as np

    bar_size = CustomSection.bar_size
    fpc = CustomSection.fce
    fy = CustomSection.fye
    fu = CustomSection.fu

    t_wall = CustomSection.t_w
    l_wall = CustomSection.l_w
    b_wall = CustomSection.b_w

    rho = CustomSection.rho

    A_bar = np.pi * (bar_size / 8. / 2.) ** 2.

    L_ip = GetL_IP(NoOfIntPoints, ForceBased=UseForceBased)

    ConcConfinedAllIP = []
    SteelConfinedAllIP = []
    ConcUnconfinedAllIP = []
    SteelUnconfinedAllIP = []

    for i in range(len(L_ip)):
        Abar = np.pi * (bar_size / 8. / 2.) ** 2.

        if Elastic:
            Ec = 0.5 * 57. * (fpc * 1000) ** .5
            Es = 29000 * 1e-13
            ElasticConcrete = OSDatabase.AddMaterial(
                OpenSeesAPI.Model.Element.Material.UniaxialMaterial.Elastic(OSDatabase.GetFreeMaterialId(1, 1), Ec,
                                                                            _Notes='Concrete Testing'))
            ElasticSteel = OSDatabase.AddMaterial(
                OpenSeesAPI.Model.Element.Material.UniaxialMaterial.Elastic(OSDatabase.GetFreeMaterialId(1, 1), Es,
                                                                            _Notes='Steel Testing'))
            ConcConfinedAllIP.append(ElasticConcrete)
            SteelConfinedAllIP.append(ElasticSteel)
            ConcUnconfinedAllIP.append(ElasticConcrete)
            SteelUnconfinedAllIP.append(ElasticSteel)

            # Shear
            # Define Elastic Spring for Shear
            Ec = 57. * (fpc * 1000.) ** .5
            Mu = 0.2  # 15
            Gc = Ec / (2. * (1. + Mu))
            Geff = Gc  # 0.1 * Gc  # assume 10% of elastic shear modulus per Oyen (2006)
            ks = 5. / 6.  # reduction factor for shear area per Pugh (2012)

            Ax = (l_wall) * t_wall

            EgX = Geff * ks * Ax

            ShearElasticX = OSDatabase.AddMaterial(
                OpenSeesAPI.Model.Element.Material.UniaxialMaterial.Elastic(OSDatabase.GetFreeMaterialId(1, 1), EgX,
                                                                            _Notes='Shear Spring for Core Wall'))

        else:
            bar_spacing = 2. * A_bar / rho * 100 / t_wall  # assume two layers

            if CustomSection.tie_vertical_spacing is not None: # If tie spacing not specified
                tie_spacing = CustomSection.tie_vertical_spacing
                tie_size = CustomSection.tie_bar_size # pugh assumes #4 only
            else:
                tie_spacing = GetTieSpacing(t_wall, bar_size / 8., bar_spacing)
                tie_size = 4 # pugh assumes #4 only

            Et = 57000. * (fpc * 1000) ** .5
            ftu = 4.0 * (fpc * 1000) ** .5
            Gfrac = 5.07e-4 * 1000.
            Lt = L_ip[i] * height
            DeltaE = 2. * Gfrac / ftu / Lt
            Ets = ftu / DeltaE
            if not RegularizeFracture:
                EtsOEt = 0.05  # Ets/Et
            else:
                EtsOEt = Ets / Et

            # Unonfined Concrete
            ConcUnconfinedAllIP.append(GetPughUnconfinedConcreteMaterial(OSDatabase, fpc,
                                                                         L_ip[i] * height,
                                                                         EtsOEt=EtsOEt,
                                                                         fpuOfpc = UnconfinedBeta,
                                                                         GfcOfpc=GfcOfpc , Regularized = Regularized))

            # ecu_Unconfined = ConcUnconfinedAllIP[-1]._ecu
            ecu_Unconfined = ConcUnconfinedAllIP[-1]._Materials[0]._ecu

            if RegularizeSteel:
                SteelUnconfinedAllIP.append(GetPughSteel02Material(OSDatabase, fy, fu, L_ip[i] * height, ecu_Unconfined, eps_ult_exp=SteelUltimateStrainTension,
                                                                   SteelPostYieldStiffness=SteelPostYieldStiffness))
            else:
                if SteelMaterialModel == 'Steel02':
                    SteelUnconfinedAllIP.append(
                        GetPughSteel02MaterialUnregularized(OSDatabase, fy, fu, ecu_Unconfined, eps_ult_exp=SteelUltimateStrainTension,
                                                                   SteelPostYieldStiffness=SteelPostYieldStiffness))
                elif SteelMaterialModel == 'Steel02Fatigue':
                    SteelUnconfinedAllIP.append(
                        GetSteel02FatigueMaterialUnregularized(OSDatabase, fy, fu, ecu_Unconfined,
                                                            eps_ult_exp=SteelUltimateStrainTension,
                                                            SteelPostYieldStiffness=SteelPostYieldStiffness))
                elif SteelMaterialModel == 'ReinforcingSteel':
                        SteelUnconfinedAllIP.append(
                            GetReinforcingSteelMaterialUnregularized(OSDatabase, fy, fu, SteelUltimateStrainTension, bar_spacing, bar_size,
                                                                   SteelPostYieldStiffness=SteelPostYieldStiffness))
                elif SteelMaterialModel == 'ReinforcingSteelWithDMBuckling':
                        SteelUnconfinedAllIP.append(
                            GetReinforcingSteelMaterialUnregularized(OSDatabase, fy, fu, SteelUltimateStrainTension, bar_spacing, bar_size,
                                                                     withBuckling=True,
                                                                   SteelPostYieldStiffness=SteelPostYieldStiffness))
                elif SteelMaterialModel == 'ReinforcingSteelWithDMBucklingAndFatigue':
                    SteelUnconfinedAllIP.append(
                        GetReinforcingSteelMaterialUnregularized(OSDatabase, fy, fu, SteelUltimateStrainTension, bar_spacing, bar_size,
                                                                 withBuckling=True, withFatigue=True,
                                                                   SteelPostYieldStiffness=SteelPostYieldStiffness))
                elif SteelMaterialModel == 'ReinforcingSteelWithFatigue':
                    SteelUnconfinedAllIP.append(
                        GetReinforcingSteelMaterialUnregularized(OSDatabase, fy, fu, SteelUltimateStrainTension, bar_spacing, bar_size,
                                                                 withBuckling=False, withFatigue=True,
                                                                   SteelPostYieldStiffness=SteelPostYieldStiffness))
                elif SteelMaterialModel == 'Pinching4':
                    SteelUnconfinedAllIP.append(
                        GetPinching4MaterialUnregularized(OSDatabase, fy, fu, SteelUltimateStrainTension, bar_spacing, bar_size, -1.0,
                                                                   SteelPostYieldStiffness=SteelPostYieldStiffness))
                elif SteelMaterialModel == 'Steel02LB':
                    SteelUnconfinedAllIP.append(
                        GetSteel02LBMaterialUnregularized(OSDatabase, fy, fu, SteelUltimateStrainTension, bar_spacing, bar_size, -1.0,
                                                                   SteelPostYieldStiffness=SteelPostYieldStiffness))
                elif SteelMaterialModel == 'Steel02WithFatigue':
                    SteelUnconfinedAllIP.append(
                        GetPughSteel02MaterialUnregularizedWithFatigue(OSDatabase, fy, fu, ecu_Unconfined,
                                                                       eps_ult_exp=SteelUltimateStrainTension,
                                                                       SteelPostYieldStiffness=SteelPostYieldStiffness))
                elif SteelMaterialModel == 'Hysteretic':
                    SteelConfinedAllIP.append(
                        GetHystereticMaterialUnregularized(OSDatabase, fy, fu, SteelUltimateStrainTension, bar_spacing, bar_size, -1.0,
                                                                   SteelPostYieldStiffness=SteelPostYieldStiffness))

            A_tie = np.pi * (tie_size / 8. / 2.) ** 2.

            # Confined Concrete
            if ConfinementModel == 'SaatRazvi':
                ConcConfinedAllIP.append(GetPughConfinedConcreteMaterial(OSDatabase, fpc,
                                                                         l_wall - 2 * cover,
                                                                         t_wall - 2 * cover,
                                                                         bar_spacing, t_wall - 2 * cover,
                                                                         tie_spacing, tie_spacing,
                                                                         fy, A_tie,
                                                                         np.ceil(l_wall / bar_spacing),
                                                                         np.ceil(t_wall / bar_spacing),
                                                                         L_ip[i] * height, GfccOGfc=GfccOGfc,
                                                                         EtsOEt=EtsOEt,
                                                                         fpuOfpc=CrushingStrengthRatio,
                                                                         GfcOfpc=GfcOfpc, Regularized = Regularized,
                                                                         ConcreteMaterialModel=ConcreteMaterialModel))
            elif ConfinementModel == 'Mander':
                if CustomSection.left_boundary != 0:
                    rho_l = np.sum(CustomSection.left_reinf_layout) * A_bar / t_wall / CustomSection.left_boundary
                else:
                    rho_l = 0.02
                ConcConfinedAllIP.append(GetManderConfinedConcreteMaterial(OSDatabase, fpc,
                                                                           l_wall - 2 * cover,
                                                                           t_wall - 2 * cover, rho_l,
                                                                           bar_spacing, t_wall - 2 * cover,
                                                                           tie_spacing, tie_spacing,
                                                                           fy, A_tie,
                                                                           np.ceil(l_wall / bar_spacing),
                                                                           np.ceil(t_wall / bar_spacing),
                                                                           L_ip[i] * height, GfccOGfc=GfccOGfc,
                                                                           EtsOEt=EtsOEt,
                                                                           fpuOfpc=CrushingStrengthRatio,
                                                                           GfcOfpc=GfcOfpc,
                                                                           ConcreteMaterialModel=ConcreteMaterialModel))
            elif ConfinementModel == 'Richart':
                if CustomSection.left_boundary != 0:
                    rho_l = np.sum(CustomSection.left_reinf_layout) * A_bar / t_wall / CustomSection.left_boundary
                else:
                    rho_l = 0.02
                ConcConfinedAllIP.append(GetRichartConfinedConcreteMaterial(OSDatabase, fpc,
                                                                           l_wall - 2 * cover,
                                                                           t_wall - 2 * cover, rho_l,
                                                                           bar_spacing, t_wall - 2 * cover,
                                                                           tie_spacing, tie_spacing,
                                                                           fy, A_tie,
                                                                           np.ceil(l_wall / bar_spacing),
                                                                           np.ceil(t_wall / bar_spacing),
                                                                           L_ip[i] * height, GfccOGfc=GfccOGfc,
                                                                           EtsOEt=EtsOEt,
                                                                           fpuOfpc=CrushingStrengthRatio,
                                                                           GfcOfpc=GfcOfpc, ConcreteMaterialModel=ConcreteMaterialModel))
            elif ConfinementModel == 'Unconfined':
                ConcConfinedAllIP.append(ConcUnconfinedAllIP[-1])

            ecu_Confined = ConcConfinedAllIP[-1]._Materials[0]._ecu

            if tie_spacing == None or tie_spacing == 0:
                tie_spacing = 1e10

            if RegularizeSteel:
                SteelConfinedAllIP.append(
                    GetPughSteel02Material(OSDatabase, fy, fu, L_ip[i] * height, ecu_Confined,
                                           eps_ult_exp=SteelUltimateStrainTension,
                                                                   SteelPostYieldStiffness=SteelPostYieldStiffness))
            else:
                if SteelMaterialModel == 'Steel02':
                    SteelConfinedAllIP.append(
                        GetPughSteel02MaterialUnregularized(OSDatabase, fy, fu, ecu_Confined,
                                                            eps_ult_exp=SteelUltimateStrainTension,
                                                                   SteelPostYieldStiffness=SteelPostYieldStiffness))
                elif SteelMaterialModel == 'Steel02Fatigue':
                    SteelConfinedAllIP.append(
                        GetSteel02FatigueMaterialUnregularized(OSDatabase, fy, fu, ecu_Confined,
                                                            eps_ult_exp=SteelUltimateStrainTension,
                                                            SteelPostYieldStiffness=SteelPostYieldStiffness))
                elif SteelMaterialModel == 'ReinforcingSteel':
                    SteelConfinedAllIP.append(
                        GetReinforcingSteelMaterialUnregularized(OSDatabase, fy, fu, SteelUltimateStrainTension,
                                                                 tie_spacing, bar_size,
                                                                 eps_ult_comp=ecu_Confined,
                                                                   SteelPostYieldStiffness=SteelPostYieldStiffness))
                elif SteelMaterialModel == 'ReinforcingSteelWithDMBuckling':
                    SteelConfinedAllIP.append(
                        GetReinforcingSteelMaterialUnregularized(OSDatabase, fy, fu, SteelUltimateStrainTension,
                                                                 tie_spacing, bar_size,
                                                                 withBuckling=True,
                                                                   SteelPostYieldStiffness=SteelPostYieldStiffness))
                elif SteelMaterialModel == 'ReinforcingSteelWithDMBucklingAndFatigue':
                    SteelConfinedAllIP.append(
                        GetReinforcingSteelMaterialUnregularized(OSDatabase, fy, fu, SteelUltimateStrainTension,
                                                                 tie_spacing, bar_size,
                                                                 withBuckling=True, withFatigue=True,
                                                                   SteelPostYieldStiffness=SteelPostYieldStiffness))
                elif SteelMaterialModel == 'ReinforcingSteelWithFatigue':
                    SteelConfinedAllIP.append(
                        GetReinforcingSteelMaterialUnregularized(OSDatabase, fy, fu, SteelUltimateStrainTension,
                                                                 tie_spacing, bar_size,
                                                                 withBuckling=False, withFatigue=True,
                                                                   SteelPostYieldStiffness=SteelPostYieldStiffness))
                elif SteelMaterialModel == 'Pinching4':
                    SteelConfinedAllIP.append(
                        GetPinching4MaterialUnregularized(OSDatabase, fy, fu, SteelUltimateStrainTension, tie_spacing,
                                                          bar_size, -1.0,
                                                                   SteelPostYieldStiffness=SteelPostYieldStiffness))
                elif SteelMaterialModel == 'Steel02LB':
                    SteelConfinedAllIP.append(
                        GetSteel02LBMaterialUnregularized(OSDatabase, fy, fu, SteelUltimateStrainTension, tie_spacing,
                                                          bar_size, -1.0,
                                                                   SteelPostYieldStiffness=SteelPostYieldStiffness))
                elif SteelMaterialModel == 'Steel02WithFatigue':
                    SteelConfinedAllIP.append(
                        GetPughSteel02MaterialUnregularizedWithFatigue(OSDatabase, fy, fu, ecu_Confined,
                                                                       eps_ult_exp=SteelUltimateStrainTension,
                                                                       SteelPostYieldStiffness=SteelPostYieldStiffness))
                elif SteelMaterialModel == 'Hysteretic':
                    SteelConfinedAllIP.append(
                        GetHystereticMaterialUnregularized(OSDatabase, fy, fu, SteelUltimateStrainTension, tie_spacing,
                                                           bar_size, -1.0,
                                                           SteelPostYieldStiffness=SteelPostYieldStiffness))

            # Shear
            # Define Elastic Spring for Shear
            Ec = 57. * (fpc * 1000.) ** .5
            Mu = 0.15  # 15
            Gc = Ec / (2. * (1. + Mu))
            Geff = 0.1 * Gc  # assume 10% of elastic shear modulus per Oyen (2006)
            ks = 5. / 6.  # reduction factor for shear area per Pugh (2012)

            Av = l_wall * t_wall

            GcOEc = Gc / Ec

            EgX = 0.4*Ec*Av*ks

            ShearElasticX = OSDatabase.AddMaterial(
                OpenSeesAPI.Model.Element.Material.UniaxialMaterial.Elastic(OSDatabase.GetFreeMaterialId(1, 1), EgX,
                                                                            _Notes='Shear Spring for Core Wall'))

    Sections = []

    for j in range(len(L_ip)):
        ConcConfined = ConcConfinedAllIP[j]
        SteelConfined = SteelConfinedAllIP[j]
        ConcUnconfined = ConcUnconfinedAllIP[j]
        SteelUnconfined = SteelUnconfinedAllIP[j]

        Fibers = []

        # Confined Parts
        if True:
            l_web = CustomSection.l_w - CustomSection.t_w
            l_wall = CustomSection.l_w
            t_wall = CustomSection.t_w
            b_wall = CustomSection.b_w

            ## Web
            ShiftY = -t_wall / 2.
            ShiftX = -l_web / 2.
            Fibers.append(OpenSeesAPI.Section.FiberSection.Patch.Rect(ConcConfined,
                                                                      int(np.ceil((l_web) / max_mesh_Size)),
                                                                      1, #int(np.ceil((t_wall) / max_mesh_Size)),
                                                                      ShiftX,
                                                                      ShiftY,
                                                                      ShiftX + l_web,
                                                                      ShiftY + t_wall))

            if CustomSection.flange_is_left:
                # Left Flange
                ShiftY = -b_wall / 2.
                ShiftX = -l_web / 2. - t_wall
                Fibers.append(OpenSeesAPI.Section.FiberSection.Patch.Rect(ConcConfined,
                                                                          int(np.ceil((t_wall) / max_mesh_Size)),
                                                                          1,  #int(np.ceil((b_wall) / max_mesh_Size)),
                                                                          ShiftX,
                                                                          ShiftY,
                                                                          ShiftX + t_wall,
                                                                          ShiftY + b_wall))
            else:
                # Right Flange
                ShiftY = -b_wall / 2.
                ShiftX = l_web / 2.
                Fibers.append(OpenSeesAPI.Section.FiberSection.Patch.Rect(ConcConfined,
                                                                          int(np.ceil((t_wall) / max_mesh_Size)),
                                                                          1,  #int(np.ceil((b_wall) / max_mesh_Size)),
                                                                          ShiftX,
                                                                          ShiftY,
                                                                          ShiftX + t_wall,
                                                                          ShiftY + b_wall))

        # Unconfined
        # Save this for cover...


        # Define Rebar
        def DefineRebar(l_x, l_y, no_bar_x, no_bar_y, ShiftX, ShiftY, Mat=SteelConfined, A=A_bar):
            for a in range(int(no_bar_y)):
                spacing_y = (l_y - 2. * cover) / float(no_bar_y - 1)
                Fibers.append(OpenSeesAPI.Section.FiberSection.Layer.Straight(Mat,
                                                                              no_bar_x,
                                                                              A,
                                                                              ShiftX + cover,
                                                                              ShiftY + cover + spacing_y * (a),
                                                                              ShiftX + l_x - cover,
                                                                              ShiftY + cover + spacing_y * (a)
                                                                              ))

        def DefineRebarUsingLayout(ReinfLayout, l_x, l_y, ShiftX, ShiftY, Mat=SteelConfined, A=A_bar):
            spacing_x = (l_x - 2. * cover) / (len(ReinfLayout) - 1)
            for a in range(len(ReinfLayout)):
                Fibers.append(OpenSeesAPI.Section.FiberSection.Layer.Straight(Mat,
                                                                              ReinfLayout[a],
                                                                              A,
                                                                              ShiftX + cover + spacing_x * (a),
                                                                              ShiftY + cover,
                                                                              ShiftX + cover + spacing_x * (a),
                                                                              ShiftY + l_y - cover
                                                                              ))

        if True:  # (Planar Wall with Boundary Element)
            # Web Rebar
            ShiftX = -l_web / 2.
            ShiftY = -t_wall / 2.
            DefineRebar(l_web , t_wall ,
                        np.ceil((l_web - 2 * cover) / float(bar_spacing)),
                        2., ShiftX, ShiftY, Mat=SteelConfined, A=Abar)

            if CustomSection.flange_is_left:
                # Flange Rebar
                ShiftX = -l_web / 2. - t_wall
                ShiftY = -b_wall / 2.
                DefineRebar(t_wall, b_wall,
                            2., np.ceil((b_wall - 2 * cover) / float(bar_spacing)),
                            ShiftX, ShiftY, Mat=SteelConfined, A=Abar)
            else:
                # Flange Rebar
                ShiftX = +l_web / 2.
                ShiftY = -b_wall / 2.
                DefineRebar(t_wall, b_wall,
                            2., np.ceil((b_wall - 2 * cover) / float(bar_spacing)),
                            ShiftX, ShiftY, Mat=SteelConfined, A=Abar)

        FiberSection = OpenSeesAPI.Section.FiberSection(OSDatabase.GetFreeMaterialId(1, 1),
                                                        Fibers, Notes='Core Wall Fiber Section')

        OSDatabase.AddMaterial(FiberSection)

        Section = OpenSeesAPI.Section.Aggregator(OSDatabase.GetFreeMaterialId(1, 1),
                                                 [ShearElasticX],
                                                 ['Vy'],
                                                 FiberSection,
                                                 _Notes='Core Wall Aggregator For Fiber Section')

        OSDatabase.AddMaterial(Section)
        Sections.append(Section)

    return Sections

def ComputeIShapedWallFiberSection(OSDatabase, l_w, b_f, t_wall, cover, fpc, fy, fu, bar_size, rho, rho_t,
                                   height, NoOfIntPoints=5, max_mesh_Size=3, Elastic=True,
                                   boundaryelement = None, rho_boundaryelement = None,
                                   CrushingStrengthRatio=0.2, GfccOGfc=2.5, ConcreteMaterialModel='Concrete02', ConfinementModel='SaatRazvi',
                                   SteelMaterialModel='Steel02', SteelUltimateStrainTension=0.2,
                                   SteelPostYieldStiffness=0.00784, GfcOfpc=2., UseForceBased=True,
                                   UnconfinedBeta=0.01, Regularized = True):
    import numpy as np

    A_bar = np.pi * (bar_size / 8. / 2.) ** 2.

    L_ip = GetL_IP(NoOfIntPoints, ForceBased=UseForceBased)

    ConcConfinedAllIP = []
    SteelConfinedAllIP = []
    ConcUnconfinedAllIP = []
    SteelUnconfinedAllIP = []

    import OSMaterialHelper as OSMat
    for i in range(len(L_ip)):
        if Elastic:
            Ec = 57. * (fpc * 1000) ** .5
            Es = 29000*1e-13
            ElasticConcrete = OSDatabase.AddMaterial(
            OpenSeesAPI.Model.Element.Material.UniaxialMaterial.Elastic(OSDatabase.GetFreeMaterialId(1, 1), Ec,
                                                                        _Notes='Concrete Testing'))
            ElasticSteel = OSDatabase.AddMaterial(
            OpenSeesAPI.Model.Element.Material.UniaxialMaterial.Elastic(OSDatabase.GetFreeMaterialId(1, 1), Es,
                                                                        _Notes='Steel Testing'))
            ConcConfinedAllIP.append(ElasticConcrete)
            SteelConfinedAllIP.append(ElasticSteel)
            ConcUnconfinedAllIP.append(ElasticConcrete)
            SteelUnconfinedAllIP.append(ElasticSteel)

            # Shear
            # Define Elastic Spring for Shear
            Ec = 57. * (fpc * 1000.) ** .5
            Mu = 0.2  # 15
            Gc = Ec / (2. * (1. + Mu))
            Geff = Gc  # 0.1 * Gc  # assume 10% of elastic shear modulus per Oyen (2006)
            ks = 5. / 6.  # reduction factor for shear area per Pugh (2012)

            Ax = b_f*4.*t_wall + (l_w - 2.*t_wall)*t_wall

            EgX = Geff * ks * Ax

            ShearElasticX = OSDatabase.AddMaterial(
                OpenSeesAPI.Model.Element.Material.UniaxialMaterial.Elastic(OSDatabase.GetFreeMaterialId(1, 1), EgX,
                                                                            _Notes='Shear Spring for Core Wall'))

        else:
            Abar = np.pi * (bar_size / 8. / 2.) ** 2.
            Aweb = np.pi * (4. / 8. / 2.) ** 2.

            if rho_boundaryelement is None:
                bar_spacing = 2.*Abar/rho*100/t_wall #assume two layers
            else:
                bar_spacing = 2. * Abar / rho_boundaryelement * 100 / t_wall  # assume two layers

            tie_size = 4. # pugh assumes #4 only

            A_tie = np.pi * (tie_size / 8. / 2.) ** 2.

            # This is how pugh computes transverse shear spacing, may not match ACI 318
            # s1 = 1./3.*t_wall
            # s2 = 6.*(bar_size / 8.)
            # s3 = max(4,min(4+(14-bar_spacing),6))
            # tie_spacing = float(min(s1,s2,s3))

            tie_spacing = GetTieSpacing(t_wall, bar_size/8., bar_spacing)

            # Confined Concrete
            if ConfinementModel == 'SaatRazvi':
                ConcConfinedAllIP.append(GetPughConfinedConcreteMaterial(OSDatabase, fpc,
                                                                         CustomSection.left_boundary - 2 * cover,
                                                                         t_wall - 2 * cover,
                                                                         bar_spacing, t_wall - 2 * cover,
                                                                         tie_spacing, tie_spacing,
                                                                         fy, A_tie,
                                                                         CustomSection.boundary_tie_x_no,
                                                                         CustomSection.boundary_tie_y_no,
                                                                         L_ip[i] * height, GfccOGfc=GfccOGfc,
                                                                         EtsOEt=EtsOEt,
                                                                         fpuOfpc=CrushingStrengthRatio,
                                                                         GfcOfpc=GfcOfpc, Regularized = Regularized,
                                                                         ConcreteMaterialModel=ConcreteMaterialModel))
            elif ConfinementModel == 'Mander':
                if CustomSection.left_boundary != 0:
                    rho_l = np.sum(CustomSection.left_reinf_layout) * A_bar / t_wall / CustomSection.left_boundary
                else:
                    rho_l = 0.02
                ConcConfinedAllIP.append(GetManderConfinedConcreteMaterial(OSDatabase, fpc,
                                                                           CustomSection.left_boundary - 2 * cover,
                                                                           t_wall - 2 * cover, rho_l,
                                                                           bar_spacing, t_wall - 2 * cover,
                                                                           tie_spacing, tie_spacing,
                                                                           fy, A_tie,
                                                                           CustomSection.boundary_tie_x_no,
                                                                           CustomSection.boundary_tie_y_no,
                                                                           L_ip[i] * height, GfccOGfc=GfccOGfc,
                                                                           EtsOEt=EtsOEt,
                                                                           fpuOfpc=CrushingStrengthRatio,
                                                                           GfcOfpc=GfcOfpc, ConcreteMaterialModel=ConcreteMaterialModel))
            elif ConfinementModel == 'Richart':
                if CustomSection.left_boundary != 0:
                    rho_l = np.sum(CustomSection.left_reinf_layout) * A_bar / t_wall / CustomSection.left_boundary
                else:
                    rho_l = 0.02
                ConcConfinedAllIP.append(GetRichartConfinedConcreteMaterial(OSDatabase, fpc,
                                                                           CustomSection.left_boundary - 2 * cover,
                                                                           t_wall - 2 * cover, rho_l,
                                                                           bar_spacing, t_wall - 2 * cover,
                                                                           tie_spacing, tie_spacing,
                                                                           fy, A_tie,
                                                                           CustomSection.boundary_tie_x_no,
                                                                           CustomSection.boundary_tie_y_no,
                                                                           L_ip[i] * height, GfccOGfc=GfccOGfc,
                                                                           EtsOEt=EtsOEt,
                                                                           fpuOfpc=CrushingStrengthRatio,
                                                                           GfcOfpc=GfcOfpc, ConcreteMaterialModel=ConcreteMaterialModel))
            elif ConfinementModel == 'Unconfined':
                ConcConfinedAllIP.append(ConcUnconfinedAllIP[-1])

            ecu_Confined = ConcConfinedAllIP[-1]._ecu

            SteelConfinedAllIP.append(OSMat.GetPughSteel02Material(OSDatabase, fy, fu, L_ip[i] * height, ecu_Confined, eps_ult_exp=0.2))

            # Unonfined Concrete
            ConcUnconfinedAllIP.append(OSMat.GetPughUnconfinedConcreteMaterial(OSDatabase, fpc, L_ip[i] * height, Regularized = Regularized))#, fpuOfpc = 0.01))

            ecu_Unconfined = ConcUnconfinedAllIP[-1]._ecu

            SteelUnconfinedAllIP.append(OSMat.GetPughSteel02Material(OSDatabase, fy, fu, L_ip[i] * height, ecu_Unconfined, eps_ult_exp=0.2))

            # Shear
            # Define Elastic Spring for Shear
            Ec = 57. * (fpc * 1000.) ** .5
            Mu = 0.15  # 15
            Gc = Ec / (2. * (1. + Mu))
            Geff = 0.1 * Gc  # assume 10% of elastic shear modulus per Oyen (2006)
            ks = 5. / 6.  # reduction factor for shear area per Pugh (2012)

            # Ax = b_f * 4. * t_wall + (l_w - 2. * t_wall) * t_wall * 2.0

            Av = l_w * t_wall * 2.0

            GcOEc = Gc/Ec

            EgX = 0.4*Ec*Av*ks

            ShearElasticX = OSDatabase.AddMaterial(
                OpenSeesAPI.Model.Element.Material.UniaxialMaterial.Elastic(OSDatabase.GetFreeMaterialId(1, 1), EgX,
                                                                            _Notes='Shear Spring for Core Wall'))

    Sections = []

    for j in range(len(L_ip)):
        ConcConfined = ConcConfinedAllIP[j]
        SteelConfined = SteelConfinedAllIP[j]
        ConcUnconfined = ConcUnconfinedAllIP[j]
        SteelUnconfined = SteelUnconfinedAllIP[j]

        Fibers = []

        # Confined Parts

        l_flange = 2.*b_f
        if t_wall == b_f:
            l_web = l_w
        else:
            l_web = l_w - 2.*t_wall

        if not(t_wall == b_f):
            ## Bottom Flange
            ShiftY = -l_flange/2.
            ShiftX = -l_web/2. - t_wall
            Fibers.append(OpenSeesAPI.Section.FiberSection.Patch.Rect(ConcConfined,
                                                                      int(np.ceil((t_wall) / max_mesh_Size)),
                                                                      int(np.ceil((l_flange) / max_mesh_Size)),
                                                                      ShiftX,
                                                                      ShiftY,
                                                                      ShiftX + t_wall,
                                                                      ShiftY + l_flange))

            ## Top Flange
            ShiftY = -l_flange/2.
            ShiftX = +l_web/2.
            Fibers.append(OpenSeesAPI.Section.FiberSection.Patch.Rect(ConcConfined,
                                                                      int(np.ceil((t_wall) / max_mesh_Size)),
                                                                      int(np.ceil((l_flange) / max_mesh_Size)),
                                                                      ShiftX,
                                                                      ShiftY,
                                                                      ShiftX + t_wall,
                                                                      ShiftY + l_flange))

        if t_wall == b_f and boundaryelement != None:
            BE = boundaryelement
            ## Web
            ShiftY = -t_wall / 2.
            ShiftX = -l_web / 2. + BE
            Fibers.append(OpenSeesAPI.Section.FiberSection.Patch.Rect(ConcUnconfined,
                                                                      int(np.ceil((l_web - BE*2.) / max_mesh_Size)),
                                                                      int(np.ceil((t_wall) / max_mesh_Size)),
                                                                      ShiftX,
                                                                      ShiftY,
                                                                      ShiftX + l_web - 2.*BE,
                                                                      ShiftY + t_wall ))
            # Left BE
            ShiftY = -t_wall / 2.
            ShiftX = -l_web / 2.
            Fibers.append(OpenSeesAPI.Section.FiberSection.Patch.Rect(ConcConfined,
                                                                      int(np.ceil((BE) / max_mesh_Size)),
                                                                      int(np.ceil((t_wall) / max_mesh_Size)),
                                                                      ShiftX,
                                                                      ShiftY,
                                                                      ShiftX + BE,
                                                                      ShiftY + t_wall ))
            # Right BE
            ShiftY = -t_wall / 2.
            ShiftX = +l_web / 2. - BE
            Fibers.append(OpenSeesAPI.Section.FiberSection.Patch.Rect(ConcConfined,
                                                                      int(np.ceil((BE) / max_mesh_Size)),
                                                                      int(np.ceil((t_wall) / max_mesh_Size)),
                                                                      ShiftX,
                                                                      ShiftY,
                                                                      ShiftX + BE,
                                                                      ShiftY + t_wall ))

        else:
            ## Web
            ShiftY = -t_wall / 2.
            ShiftX = -l_web/2.
            Fibers.append(OpenSeesAPI.Section.FiberSection.Patch.Rect(ConcConfined,
                                                                      int(np.ceil((l_web) / max_mesh_Size)),
                                                                      int(np.ceil((t_wall) / max_mesh_Size)),
                                                                      ShiftX,
                                                                      ShiftY,
                                                                      ShiftX + l_web,
                                                                      ShiftY + t_wall))

        # Unconfined

        # Define Rebar
        def DefineRebar(l_x, l_y, no_bar_x, no_bar_y, ShiftX, ShiftY, Mat=SteelConfined, A=A_bar):
            for a in range(int(no_bar_y)):
                spacing_y = (l_y - 2.*cover) / float(no_bar_y-1)
                Fibers.append(OpenSeesAPI.Section.FiberSection.Layer.Straight(Mat,
                                                                              no_bar_x,
                                                                              A,
                                                                              ShiftX + cover,
                                                                              ShiftY + cover + spacing_y*(a),
                                                                              ShiftX + l_x - cover,
                                                                              ShiftY + cover + spacing_y*(a)
                                                                              ))

        if not (t_wall == b_f): # Flange Rebar (I Section)
            #Bottom Flange Rebar
            ShiftY = -l_flange/2.
            ShiftX = -l_web/2. - t_wall
            DefineRebar(t_wall, l_flange,
                        2, np.ceil((l_flange-2*cover)/float(bar_spacing)),
                        ShiftX, ShiftY)

            #Top Flange Rebar
            ShiftY = -l_flange/2.
            ShiftX = +l_web/2.
            DefineRebar(t_wall, l_flange,
                        2, np.ceil((l_flange-2*cover)/float(bar_spacing)),
                        ShiftX, ShiftY)

        if t_wall == b_f and boundaryelement != None: # (Planar Wall with Boundary Element)
            # #Web Rebar
            # ShiftY = -t_wall
            # ShiftX = -l_web/2.
            # DefineRebar(l_web, t_wall*2.,
            #             np.ceil((l_web-2*cover) / float(bar_spacing)), 4., ShiftX, ShiftY)

            bar_spacing_web = 2. * Aweb / rho * 100 / t_wall  # assume two layers

            #Web Rebar
            ShiftY = -t_wall / 2.
            ShiftX = -l_web/2. + BE
            DefineRebar(l_web - 2*BE, t_wall,
                        np.ceil((l_web - 2 * cover - 2.*BE) / float(bar_spacing_web)),
                        2., ShiftX, ShiftY, Mat=SteelUnconfined, A=Aweb)

            ### Constant rho in BE and Web, as per Pugh's assumptions
            # Left BE
            ShiftY = -t_wall / 2.
            ShiftX = -l_web / 2.
            DefineRebar(BE, t_wall,
                        np.ceil((BE) / float(bar_spacing)), 2., ShiftX, ShiftY)

            # Right BE
            ShiftY = -t_wall / 2.
            ShiftX = +l_web / 2. - BE
            DefineRebar(BE, t_wall,
                        np.ceil((BE) / float(bar_spacing)), 2., ShiftX, ShiftY)
        else:
            #Web Rebar
            ShiftY = -t_wall / 2.
            ShiftX = -l_web/2.
            DefineRebar(l_web, t_wall,
                        np.ceil((l_web-2*cover) / float(bar_spacing))*0.25/rho, 2., ShiftX, ShiftY)

        FiberSection = OpenSeesAPI.Section.FiberSection(OSDatabase.GetFreeMaterialId(1, 1),
                                                                        Fibers, Notes='Core Wall Fiber Section')

        OSDatabase.AddMaterial(FiberSection)

        Section = OpenSeesAPI.Section.Aggregator(OSDatabase.GetFreeMaterialId(1, 1),
                                                                        [ShearElasticX],
                                                                        ['Vy'],
                                                                        FiberSection,
                                                                        _Notes='Core Wall Aggregator For Fiber Section')

        OSDatabase.AddMaterial(Section)

        Sections.append(Section)

    return Sections

def ComputeTShapedWallFiberSection(OSDatabase, l_w, b_f, t_wall, cover, fpc, fy, fu, bar_size,
                                   rho, rho_t, height, NoOfIntPoints=5, max_mesh_Size=3,
                                   Elastic=True, LeftFlange=True, CrushingStrengthRatio=0.2, GfccOGfc=2.5,
                                   ConcreteMaterialModel='SaatRazvi',
                                   SteelMaterialModel='Steel02', SteelUltimateStrainTension=0.2,
                                   SteelPostYieldStiffness=0.00784, GfcOfpc=2., UseForceBased=True,
                                   UnconfinedBeta=0.01, Regularized = True,):
    import numpy as np

    A_bar = np.pi * (bar_size / 8. / 2.) ** 2.

    L_ip = GetL_IP(NoOfIntPoints, ForceBased=UseForceBased)

    ConcConfinedAllIP = []
    SteelConfinedAllIP = []

    import OSMaterialHelper as OSMat
    for i in range(len(L_ip)):
        if Elastic:
            Ec = 57. * (fpc * 1000) ** .5
            Es = 29000*1e-13
            ElasticConcrete = OSDatabase.AddMaterial(
            OpenSeesAPI.Model.Element.Material.UniaxialMaterial.Elastic(OSDatabase.GetFreeMaterialId(1, 1), Ec,
                                                                        _Notes='Concrete Testing'))
            ElasticSteel = OSDatabase.AddMaterial(
            OpenSeesAPI.Model.Element.Material.UniaxialMaterial.Elastic(OSDatabase.GetFreeMaterialId(1, 1), Es,
                                                                        _Notes='Steel Testing'))
            ConcConfinedAllIP.append(ElasticConcrete)
            SteelConfinedAllIP.append(ElasticSteel)

            # Shear
            # Define Elastic Spring for Shear
            Ec = 57. * (fpc * 1000.) ** .5
            Mu = 0.2  # 15
            Gc = Ec / (2. * (1. + Mu))
            Geff = Gc  # 0.1 * Gc  # assume 10% of elastic shear modulus per Oyen (2006)
            ks = 5. / 6.  # reduction factor for shear area per Pugh (2012)

            Ax = b_f*4.*t_wall + (l_w - 2.*t_wall)*t_wall

            EgX = Geff * ks * Ax

            ShearElasticX = OSDatabase.AddMaterial(
                OpenSeesAPI.Model.Element.Material.UniaxialMaterial.Elastic(OSDatabase.GetFreeMaterialId(1, 1), EgX,
                                                                            _Notes='Shear Spring for Core Wall'))

        else:
            Abar = np.pi * (bar_size / 8. / 2.) ** 2.
            bar_spacing = 2.*Abar/rho*100/t_wall #assume two layers

            tie_size = 4. # pugh assumes #4 only

            A_tie = np.pi * (tie_size / 8. / 2.) ** 2.

            # This is how pugh computes transverse shear spacing, may not match ACI 318
            # s1 = 1./3.*t_wall
            # s2 = 6.*(bar_size / 8.)
            # s3 = max(4,min(4+(14-bar_spacing),6))
            # tie_spacing = float(min(s1,s2,s3))

            tie_spacing = GetTieSpacing(t_wall, bar_size/8., bar_spacing)

            # Confined Concrete
            ConcConfinedAllIP.append(OSMat.GetPughConfinedConcreteMaterial(OSDatabase, fpc, bar_spacing, t_wall-2*cover,
                                                                           bar_spacing, t_wall-2*cover,
                                                                           tie_spacing, tie_spacing,
                                                                           fy, A_tie, 1., 2.,
                                                                           L_ip[i] * height, GfccOGfc=GfccOGfc,
                                                                           GfcOfpc=GfcOfpc, Regularized = Regularized))

            ecu_Confined = ConcConfinedAllIP[-1]._ecu

            SteelConfinedAllIP.append(OSMat.GetPughSteel02Material(OSDatabase, fy, fu, L_ip[i] * height, ecu_Confined, eps_ult_exp=0.2))

            # Shear
            # Define Elastic Spring for Shear
            Ec = 57. * (fpc * 1000.) ** .5
            Mu = 0.15  # 15
            Gc = Ec / (2. * (1. + Mu))
            Geff = 0.1 * Gc  # assume 10% of elastic shear modulus per Oyen (2006)
            ks = 5. / 6.  # reduction factor for shear area per Pugh (2012)

            # Ax = b_f * 4. * t_wall + (l_w - 2. * t_wall) * t_wall * 2.0

            Av = l_w * t_wall * 2.0

            GcOEc = Gc/Ec

            EgX = 0.4*Ec*Av*ks

            ShearElasticX = OSDatabase.AddMaterial(
                OpenSeesAPI.Model.Element.Material.UniaxialMaterial.Elastic(OSDatabase.GetFreeMaterialId(1, 1), EgX,
                                                                            _Notes='Shear Spring for Core Wall'))

    Sections = []

    for j in range(len(L_ip)):
        ConcConfined = ConcConfinedAllIP[j]
        SteelConfined = SteelConfinedAllIP[j]

        Fibers = []

        # Confined Parts

        l_flange = b_f
        l_web = l_w - t_wall

        if LeftFlange:
            ## Bottom Flange
            ShiftY = -l_flange/2.
            ShiftX = -l_web/2. - t_wall
            Fibers.append(OpenSeesAPI.Section.FiberSection.Patch.Rect(ConcConfined,
                                                                      int(np.ceil((t_wall) / max_mesh_Size)),
                                                                      int(np.ceil((l_flange) / max_mesh_Size)),
                                                                      ShiftX,
                                                                      ShiftY,
                                                                      ShiftX + t_wall,
                                                                      ShiftY + l_flange))

        if not(LeftFlange):
            ## Top Flange
            ShiftY = -l_flange/2.
            ShiftX = +l_web/2.
            Fibers.append(OpenSeesAPI.Section.FiberSection.Patch.Rect(ConcConfined,
                                                                      int(np.ceil((t_wall) / max_mesh_Size)),
                                                                      int(np.ceil((l_flange) / max_mesh_Size)),
                                                                      ShiftX,
                                                                      ShiftY,
                                                                      ShiftX + t_wall,
                                                                      ShiftY + l_flange))

        ## Web
        ShiftY = -t_wall/2.
        ShiftX = -l_web/2.
        Fibers.append(OpenSeesAPI.Section.FiberSection.Patch.Rect(ConcConfined,
                                                                  int(np.ceil((l_web) / max_mesh_Size)),
                                                                  int(np.ceil((t_wall) / max_mesh_Size)),
                                                                  ShiftX,
                                                                  ShiftY,
                                                                  ShiftX + l_web,
                                                                  ShiftY + t_wall))

        # Unconfined

        # Define Rebar
        def DefineRebar(l_x, l_y, no_bar_x, no_bar_y, ShiftX, ShiftY):
            # SteelConfined
            # Reinf. In the X
            # Fibers.append(OpenSeesAPI.Section.FiberSection.Layer.Straight(SteelConfined,
            #                                                               no_bar_x,
            #                                                               A_bar,
            #                                                               ShiftX + cover,
            #                                                               ShiftY + cover,
            #                                                               ShiftX + l_x - cover,
            #                                                               ShiftY + cover
            #                                                               ))
            #
            # Fibers.append(OpenSeesAPI.Section.FiberSection.Layer.Straight(SteelConfined,
            #                                                               no_bar_x,
            #                                                               A_bar,
            #                                                               ShiftX + cover,
            #                                                               ShiftY + l_y - cover,
            #                                                               ShiftX + l_x - cover,
            #                                                               ShiftY + l_y - cover
            #                                                               ))
            #
            # # Reinf. In the Y
            # if no_bar_y > 2:
            #     for a in range(int(no_bar_y-1)):
            #         spacing_y = (l_y - 2.*cover) / float(no_bar_y-1)
            #         Fibers.append(OpenSeesAPI.Section.FiberSection.Layer.Straight(SteelConfined,
            #                                                                       no_bar_x,
            #                                                                       A_bar,
            #                                                                       ShiftX + cover,
            #                                                                       ShiftY + cover + spacing_y*(a+1),
            #                                                                       ShiftX + l_x - cover,
            #                                                                       ShiftY + cover + spacing_y*(a+1)
            #                                                                       ))
            #
            # # Reinf. In the Y
            for a in range(int(no_bar_y)):
                spacing_y = (l_y - 2.*cover) / float(no_bar_y-1)
                Fibers.append(OpenSeesAPI.Section.FiberSection.Layer.Straight(SteelConfined,
                                                                              no_bar_x,
                                                                              A_bar,
                                                                              ShiftX + cover,
                                                                              ShiftY + cover + spacing_y*(a),
                                                                              ShiftX + l_x - cover,
                                                                              ShiftY + cover + spacing_y*(a)
                                                                              ))

        if LeftFlange:
            #Bottom Flange Rebar
            ShiftY = -l_flange/2.
            ShiftX = -l_web/2. - t_wall
            DefineRebar(t_wall, l_flange,
                        2, np.ceil((l_flange-2*cover)/float(bar_spacing)),
                        ShiftX, ShiftY)

        if not(LeftFlange):
            #Top Flange Rebar
            ShiftY = -l_flange/2.
            ShiftX = +l_web/2.
            DefineRebar(t_wall, l_flange,
                        2, np.ceil((l_flange-2*cover)/float(bar_spacing)),
                        ShiftX, ShiftY)

        #Web Rebar
        ShiftY = -t_wall/2.
        ShiftX = -l_web/2.
        DefineRebar(l_web, t_wall,
                    np.ceil((l_web-2*cover) / float(bar_spacing)), 2., ShiftX, ShiftY)

        FiberSection = OpenSeesAPI.Section.FiberSection(OSDatabase.GetFreeMaterialId(1, 1),
                                                                        Fibers, Notes='Core Wall Fiber Section')

        OSDatabase.AddMaterial(FiberSection)

        Section = OpenSeesAPI.Section.Aggregator(OSDatabase.GetFreeMaterialId(1, 1),
                                                                        [ShearElasticX],
                                                                        ['Vy'],
                                                                        FiberSection,
                                                                        _Notes='Core Wall Aggregator For Fiber Section')

        OSDatabase.AddMaterial(Section)

        Sections.append(Section)

    return Sections

def DiagonalCouplingBeamSection(OSDatabase, l_x, l_y, cover, NoDiagonalBarsX, NoDiagonalBarsY, WidthDiagonal,
                                DepthDiagonal, BarSize, fpc, fy, fu, BeamLength, tie_spacing = 4,
                                s_bar_x = 6, s_bar_y = 6, tie_size=5, meshsize = 6,
                                NoOfIntPoints=5, Elastic = True, RegularizeSteel = True, RegularizeFracture=False,
                                CrushingStrengthRatio=0.2, GfccOGfc=2.5, ConcreteMaterialModel='Concrete02', ConfinementModel='SaatRazvi',
                                SteelMaterialModel='Steel02', SteelUltimateStrainTension=0.2,
                                SteelPostYieldStiffness=0.00784, GfcOfpc=2., UseForceBased=True,
                                NoOfDivisions = 1,
                                UnconfinedBeta=0.01, Regularized=True,
                                ):
    import OSMaterialHelper as OSMat

    # cover = 0. #0.75 + float(tie_size)/8./2. # To Centriod of Stirrup

    A_tie = np.pi * (tie_size / 8. / 2.) ** 2.

    L_ip = GetL_IP(NoOfIntPoints, ForceBased=UseForceBased)

    AllSections = []
    # Preliminary Concrete Models
    for k in range(NoOfDivisions):
        BeamSegmentLength = BeamLength/NoOfDivisions
        ConcConfinedAllIP = []
        SteelConfinedAllIP = []
        for i in range(len(L_ip)):

            no_ties_x = np.ceil(l_x / tie_spacing)
            no_ties_y = np.ceil(l_y / tie_spacing)

            if Elastic:
                Ec = 0.5 * 57. * (fpc * 1000) ** .5
                Es = 29000 * 1e-10
                ElasticConcrete = OSDatabase.AddMaterial(
                    OpenSeesAPI.Model.Element.Material.UniaxialMaterial.Elastic(OSDatabase.GetFreeMaterialId(1, 1), Ec,
                                                                                _Notes='Concrete Testing'))
                ElasticSteel = OSDatabase.AddMaterial(
                    OpenSeesAPI.Model.Element.Material.UniaxialMaterial.Elastic(OSDatabase.GetFreeMaterialId(1, 1), Es,
                                                                                _Notes='Steel Testing'))
                ConcConfinedAllIP.append(ElasticConcrete)
                SteelConfinedAllIP.append(ElasticSteel)

                # Shear
                # Define Elastic Spring for Shear
                Ec = 57. * (fpc * 1000) ** .5
                Mu = 0.2  # 15
                Gc = Ec / (2. * (1. + Mu))
                Geff = Gc  # 0.1 * Gc  # assume 10% of elastic shear modulus per Oyen (2006)
                ks = 5. / 6.  # reduction factor for shear area per Pugh (2012)
                Eg = Geff * ks * l_x * l_y

                ShearElastic = OSDatabase.AddMaterial(
                    OpenSeesAPI.Model.Element.Material.UniaxialMaterial.Elastic(OSDatabase.GetFreeMaterialId(1, 1), Eg,
                                                                                _Notes='Shear Spring of Coupling Beam '))

            else:
                Et = 57000. * (fpc * 1000) ** .5
                ftu = 4.0 * (fpc * 1000)**.5
                Gfrac = 5.07e-4 * 1000 #kip/in
                Lt = L_ip[i] * BeamLength
                DeltaE = 2. * Gfrac / ftu / Lt
                Ets = ftu / DeltaE
                if not RegularizeFracture:
                    EtsOEt = 0.05 #Ets/Et
                else:
                    EtsOEt = Ets/Et

                # Confined Concrete
                if ConfinementModel == 'SaatRazvi':
                    ConcConfinedAllIP.append(GetPughConfinedConcreteMaterial(OSDatabase, fpc, l_x, l_y, s_bar_x, s_bar_y,
                                                                           tie_spacing, tie_spacing, fy, A_tie,
                                                                           no_ties_x,
                                                                           no_ties_y, L_ip[i] * BeamLength,
                                                                             GfccOGfc=GfccOGfc, EtsOEt=EtsOEt,
                                                                             fpuOfpc=CrushingStrengthRatio,
                                                                             GfcOfpc=GfcOfpc, Regularized = Regularized,
                                                                             ConcreteMaterialModel=ConcreteMaterialModel))
                elif ConfinementModel == 'Mander':
                    rho_l = np.sum([int(l_x / s_bar_x), int(l_y / s_bar_y)]) * 2. * A_tie / l_x / l_y
                    ConcConfinedAllIP.append(GetManderConfinedConcreteMaterial(OSDatabase, fpc,
                                                                               l_x, l_y, rho_l,
                                                                               s_bar_x, s_bar_y,
                                                                               tie_spacing, tie_spacing,
                                                                               fy, A_tie,
                                                                               no_ties_x, no_ties_y,
                                                                               L_ip[i] * BeamLength, GfccOGfc=GfccOGfc,
                                                                               EtsOEt=EtsOEt,
                                                                               fpuOfpc=CrushingStrengthRatio,
                                                                               GfcOfpc=GfcOfpc,
                                                                               ConcreteMaterialModel=ConcreteMaterialModel))
                elif ConfinementModel == 'Richart':
                    rho_l = np.sum([int(l_x / s_bar_x), int(l_y / s_bar_y)]) * 2. * A_tie / l_x / l_y
                    ConcConfinedAllIP.append(GetRichartConfinedConcreteMaterial(OSDatabase, fpc,
                                                                               l_x, l_y, rho_l,
                                                                               s_bar_x, s_bar_y,
                                                                               tie_spacing, tie_spacing,
                                                                               fy, A_tie,
                                                                               no_ties_x, no_ties_y,
                                                                               L_ip[i] * BeamLength, GfccOGfc=GfccOGfc,
                                                                               EtsOEt=EtsOEt,
                                                                               fpuOfpc=CrushingStrengthRatio,
                                                                               GfcOfpc=GfcOfpc, ConcreteMaterialModel=ConcreteMaterialModel))
                elif ConfinementModel == 'Unconfined':
                    ConcConfinedAllIP.append(GetPughConfinedConcreteMaterial(OSDatabase, fpc, l_x, l_y, s_bar_x, s_bar_y,
                                                                           tie_spacing, tie_spacing, fy, A_tie,
                                                                           no_ties_x,
                                                                           no_ties_y, L_ip[i] * BeamLength,
                                                                             GfccOGfc=GfccOGfc, EtsOEt=EtsOEt,
                                                                             fpuOfpc=CrushingStrengthRatio,
                                                                             GfcOfpc=GfcOfpc, Regularized = Regularized,
                                                                             ConcreteMaterialModel= ConcreteMaterialModel))

                ecu_Confined = ConcConfinedAllIP[-1]._Materials[0]._ecu

                # SteelConfinedAllIP.append(GetPughSteel02Material(OSDatabase, fy, fu, L_ip[i] * BeamLength, ecu_Confined))

                if tie_spacing == None or tie_spacing == 0:
                    tie_spacing = 1e10

                if RegularizeSteel:
                    SteelConfinedAllIP.append(
                        GetPughSteel02Material(OSDatabase, fy, fu, L_ip[i] * BeamLength, ecu_Confined,
                                               eps_ult_exp=SteelUltimateStrainTension,
                                                                       SteelPostYieldStiffness=SteelPostYieldStiffness))
                else:
                    if SteelMaterialModel == 'Steel02':
                        SteelConfinedAllIP.append(
                            GetPughSteel02MaterialUnregularized(OSDatabase, fy, fu, ecu_Confined,
                                                                eps_ult_exp=SteelUltimateStrainTension,
                                                                       SteelPostYieldStiffness=SteelPostYieldStiffness))
                    elif SteelMaterialModel == 'Steel02Fatigue':
                        SteelConfinedAllIP.append(
                            GetSteel02FatigueMaterialUnregularized(OSDatabase, fy, fu, ecu_Confined,
                                                                   eps_ult_exp=SteelUltimateStrainTension,
                                                                   SteelPostYieldStiffness=SteelPostYieldStiffness))
                    elif SteelMaterialModel == 'ReinforcingSteel':
                        SteelConfinedAllIP.append(
                            GetReinforcingSteelMaterialUnregularized(OSDatabase, fy, fu, SteelUltimateStrainTension,
                                                                     tie_spacing, BarSize,
                                                                     eps_ult_comp=ecu_Confined,
                                                                       SteelPostYieldStiffness=SteelPostYieldStiffness))
                    elif SteelMaterialModel == 'ReinforcingSteelWithDMBuckling':
                        SteelConfinedAllIP.append(
                            GetReinforcingSteelMaterialUnregularized(OSDatabase, fy, fu, SteelUltimateStrainTension,
                                                                     tie_spacing, BarSize,
                                                                     withBuckling=True,
                                                                       SteelPostYieldStiffness=SteelPostYieldStiffness))
                    elif SteelMaterialModel == 'ReinforcingSteelWithDMBucklingAndFatigue':
                        SteelConfinedAllIP.append(
                            GetReinforcingSteelMaterialUnregularized(OSDatabase, fy, fu, SteelUltimateStrainTension,
                                                                     tie_spacing, BarSize,
                                                                     withBuckling=True, withFatigue=True,
                                                                       SteelPostYieldStiffness=SteelPostYieldStiffness))
                    elif SteelMaterialModel == 'ReinforcingSteelWithFatigue':
                        SteelConfinedAllIP.append(
                            GetReinforcingSteelMaterialUnregularized(OSDatabase, fy, fu, SteelUltimateStrainTension,
                                                                     tie_spacing, BarSize,
                                                                     withBuckling=False, withFatigue=True,
                                                                       SteelPostYieldStiffness=SteelPostYieldStiffness))
                    elif SteelMaterialModel == 'Pinching4':
                        SteelConfinedAllIP.append(
                            GetPinching4MaterialUnregularized(OSDatabase, fy, fu, SteelUltimateStrainTension, tie_spacing,
                                                              BarSize, -1.0,
                                                                       SteelPostYieldStiffness=SteelPostYieldStiffness))
                    elif SteelMaterialModel == 'Steel02LB':
                        SteelConfinedAllIP.append(
                            GetSteel02LBMaterialUnregularized(OSDatabase, fy, fu, SteelUltimateStrainTension, tie_spacing,
                                                              BarSize, -1.0,
                                                                       SteelPostYieldStiffness=SteelPostYieldStiffness))
                    elif SteelMaterialModel == 'Steel02WithFatigue':
                        SteelConfinedAllIP.append(
                            GetPughSteel02MaterialUnregularizedWithFatigue(OSDatabase, fy, fu, ecu_Confined,
                                                                           eps_ult_exp=SteelUltimateStrainTension,
                                                                           SteelPostYieldStiffness=SteelPostYieldStiffness))
                    elif SteelMaterialModel == 'Hysteretic':
                        SteelConfinedAllIP.append(
                            GetHystereticMaterialUnregularized(OSDatabase, fy, fu, SteelUltimateStrainTension, tie_spacing,
                                                              BarSize, -1.0,
                                                                       SteelPostYieldStiffness=SteelPostYieldStiffness))
                # Shear
                # Define Elastic Spring for Shear
                Ec = 57. * (fpc * 1000) ** .5
                Mu = 0.15
                Gc = Ec / (2. * (1. + Mu))
                Geff = 0.1 * Gc  # assume 10% of elastic shear modulus per Oyen (2006)
                ks = 5. / 6.  # reduction factor for shear area per Pugh (2012)
                Eg = Geff * ks * l_x * l_y

                ShearElastic = OSDatabase.AddMaterial(
                    OpenSeesAPI.Model.Element.Material.UniaxialMaterial.Elastic(OSDatabase.GetFreeMaterialId(1, 1), Eg,
                                                                                _Notes='Shear Spring for Core Wall'))

        Alpha = np.arctan((l_y-2*cover-DepthDiagonal)/BeamLength) #Bar Angle

        Abar = np.pi*(BarSize/8.)**2./4.*np.cos(Alpha) #Steel Area Rebar in Horizontal Direction

        Sections = []

        for i in range(NoOfIntPoints):

            InclinationShift = (BeamSegmentLength * k + (np.sum(L_ip[:i])/np.sum(L_ip[:])+L_ip[i]/np.sum(L_ip[:])/2.)*BeamSegmentLength)*np.tan(Alpha)

            Fibers = []
            #Define Concrete
            Fibers.append(OpenSeesAPI.Material.Section.FiberSection.Patch.Rect(ConcConfinedAllIP[i],
                                                                               np.ceil(l_x / meshsize),
                                                                               np.ceil(l_y / meshsize),
                                                                               0., 0.,
                                                                               l_x - 0., l_y - 0.,
                                                                               Notes='Coupling Beam'))

            #Define Steel
            SpacingY = DepthDiagonal/(NoDiagonalBarsY-1.)
            for j in range(NoDiagonalBarsY):
                Fibers.append(
                    OpenSeesAPI.Material.Section.FiberSection.Layer.Straight(SteelConfinedAllIP[i], NoDiagonalBarsX,
                                                                             Abar, cover,
                                                                             l_y - cover - SpacingY*j - InclinationShift,
                                                                             cover + WidthDiagonal,
                                                                             l_y - cover - SpacingY*j - InclinationShift))

                Fibers.append(
                    OpenSeesAPI.Material.Section.FiberSection.Layer.Straight(SteelConfinedAllIP[i], NoDiagonalBarsX,
                                                                             Abar, l_x - cover,
                                                                             cover + SpacingY*j + InclinationShift,
                                                                             l_x - cover - WidthDiagonal,
                                                                             cover + SpacingY*j + InclinationShift))

            FiberSection = OSDatabase.AddMaterial(
                OpenSeesAPI.Material.Section.FiberSection(OSDatabase.GetFreeMaterialId(1,1),
                                                          Fibers))

            Sections.append(OSDatabase.AddMaterial(OpenSeesAPI.Section.Aggregator(OSDatabase.GetFreeMaterialId(1, 1),
                                                     [ShearElastic],
                                                     ['Vy'],
                                                     FiberSection)))

        AllSections.append(Sections)
    return AllSections

def GetTieSpacing(thickness, l_bar_dia, l_bar_spacing):
    s1 = 1./4.*thickness
    s2 = 6.*l_bar_dia
    s3 = max(4,min(4+(14.-l_bar_spacing)/3.,6))

    return min(s1, s2, s3)

def GetEpsCU(fibers, CheckIfConfined = False):
    for fiber in fibers:
        if isinstance(fiber._Mat._Materials[0], OpenSeesAPI.Material.UniaxialMaterial.Concrete02):
            if not(CheckIfConfined):
                return fiber._Mat._Materials[0]._ecu
            else:
                if fiber._Mat._Materials[0]._Notes == 'Confined Conc.':
                    return fiber._Mat._Materials[0]._ecu

def GetEpsSU(fibers):
    for fiber in fibers:
        if isinstance(fiber._Mat._Materials[0], OpenSeesAPI.Material.UniaxialMaterial.MinMax):
            return fiber._Mat._Materials[0]._maxStrain

def GetEpsSY(fibers):
    for fiber in fibers:
        if isinstance(fiber._Mat._Materials[0],OpenSeesAPI.Material.UniaxialMaterial.MinMax):
            try:
                fy = fiber._Mat._Materials[0]._OtherMaterial._Fy
                E0 = fiber._Mat._Materials[0]._OtherMaterial._E0
            except:
                try:
                    fy = fiber._Mat._Materials[0]._OtherMaterial._ePf1
                    E0 = fiber._Mat._Materials[0]._OtherMaterial._ePf1/fiber._Mat._Materials[0]._OtherMaterial._ePd1
                except:
                    try:
                        fy = fiber._Mat._Materials[0]._OtherMaterial._s1p
                        eps_y = fiber._Mat._Materials[0]._OtherMaterial._e1p
                        E0 = fy/eps_y
                    except:
                        try:
                            fy = fiber._Mat._Materials[0]._OtherMaterial._OtherMaterial._Fy
                            E0 = fiber._Mat._Materials[0]._OtherMaterial._OtherMaterial._E0
                        except:
                            fy = fiber._Mat._Materials[0]._OtherMaterial._fy
                            E0 = fiber._Mat._Materials[0]._OtherMaterial._Es
            return fy/float(E0)

def GetPughConfinedConcreteMaterial(OSDatabase, fpc, b_x, b_y, s_bars_x, s_bars_y, s_ties_x, s_ties_y, f_y,
                                    A_tie, no_tie_x, no_tie_y, L_ip, GfccOGfc=1.7, fpuOfpc=0.2, EtsOEt = 0.05,
                                    GfcOfpc=2., Regularized = True, eps_u_steel = 0.2, ConcreteMaterialModel='Concrete02'):
    import OpenSeesAPI
    import math

    # Saatcioglu and Razvi 1992
    ksi_to_mpa = 6.8947
    f_l_x = A_tie * f_y * no_tie_x / s_ties_x / b_x * ksi_to_mpa  # Uniform Confining Pressure
    k2_x = min(0.26 * (b_x / s_ties_x * b_x / s_bars_x * 1. / f_l_x) ** .5, 1.)
    f_1e_x = k2_x * f_l_x

    f_l_y = A_tie * f_y * no_tie_y / s_ties_y / b_y * ksi_to_mpa  # Uniform Confining Pressure
    k2_y = min(0.26 * (b_y / s_ties_y * b_y / s_bars_y * 1. / f_l_y) ** .5, 1.)
    f_1e_y = k2_y * f_l_y

    f_1e = (f_1e_x * b_x + f_1e_y * b_y) / (b_x + b_y)
    k_1 = 6.7 * f_1e ** (-0.17)

    fpc_confined = fpc + (k_1 * f_1e) / ksi_to_mpa

    Kc = fpc_confined / fpc

    f20 = fpuOfpc * Kc * fpc

    rho_x = A_tie * math.ceil(b_x / s_bars_x + 1.) / s_ties_x / b_x
    rho_y = A_tie * math.ceil(b_y / s_bars_y + 1.) / s_ties_y / b_y
    rho = rho_x + rho_y

    e085 = 0.0038

    e01 = 2. * fpc / (57. * (fpc * 1000.) ** .5)  # 0.002

    ec = e01 * (1. + 5. * k_1 * f_1e / fpc / ksi_to_mpa) # This is the original SAAT and Razvi Equation
    # ec = 2. * Kc * fpc / 57. / (Kc * fpc * 1000) ** .5  # Pugh Dissertation Equation 2.4

    # e1 = e01 * (1. + 5. * Kc)

    # e20 = 260. * rho * e1 + e085  # Equation 14/15 in Saat and Razvi

    f_t = 4.0 * (fpc * 1000) ** .5 / 1000

    Ec = 57. * (fpc * 1000) ** .5
    Et = Ec
    Ets = EtsOEt * Et

    lam = 0.1  # Not Sure where this is mentioned in Pugh Thesis

    Gfc = GfcOfpc * fpc * 0.0393701
    Gfcc = GfccOGfc * Gfc

    fpeak = Kc * fpc

    # e20u_reg = Gfcc / 0.6 / fpeak / L_ip - 0.8 * fpeak / Ec + ec # Pugh's Equations
    beta = fpuOfpc
    if Regularized:
        e20u_reg = 2. * Gfcc / L_ip / fpeak / (1. + beta) + ec * (beta + 1.0) / 2.  # Revised
    else:
        e20u_reg = 0.004 + 1.4 * rho * f_y * 0.10 / fpeak # Revised

    e20u_reg = max(ec + 0.0001, e20u_reg)

    if ConcreteMaterialModel == 'Concrete02':
        Conc = OpenSeesAPI.Model.Element.Material.UniaxialMaterial.Concrete02(OSDatabase.GetFreeMaterialId(1, 1),
                                                                             -1 * fpc_confined, -1 * ec, -1 * f20,
                                                                             -1 * e20u_reg, lam, f_t, Ets,
                                                                             _Notes='Confined Conc.')
    elif ConcreteMaterialModel == 'Concrete02IS':
        Conc = OpenSeesAPI.Model.Element.Material.UniaxialMaterial.Concrete02IS(OSDatabase.GetFreeMaterialId(1, 1), Ec,
                                                                             -1 * fpc_confined, -1 * ec, -1 * f20,
                                                                             -1 * e20u_reg, lam, f_t, Ets,
                                                                             _Notes='Confined Conc.')
    elif ConcreteMaterialModel == 'ConcreteCM':
        Conc = OpenSeesAPI.Model.Element.Material.UniaxialMaterial.ConcreteCM(OSDatabase.GetFreeMaterialId(1, 1),
                                                                             -1 * fpc_confined, -1 * ec, Ec, -1 * e20u_reg,
                                                                              Ec, f_t, f_t / Ets,
                                                                             _Notes='Confined Conc.')
    elif ConcreteMaterialModel == 'Concrete04':
        Conc = OpenSeesAPI.Model.Element.Material.UniaxialMaterial.Concrete04(OSDatabase.GetFreeMaterialId(1, 1),
                                                                             -1 * fpc_confined, -1 * ec, -1 * e20u_reg, Ec,
                                                                              f_t, f_t / Ets, 0.1,
                                                                             _Notes='Confined Conc.')
    elif ConcreteMaterialModel == 'Concrete02Pugh':
        ec = 2. * fpc / (57. * (fpc * 1000.) ** .5)
        e20u_reg = 2. * Gfcc / L_ip / fpeak / (1. + beta) + ec * (beta + 1.0) / 2.
        Conc = OpenSeesAPI.Model.Element.Material.UniaxialMaterial.Concrete02(OSDatabase.GetFreeMaterialId(1, 1),
                                                                             -1 * fpc_confined, -1 * ec, -1 * f20,
                                                                             -1 * e20u_reg, lam, f_t, Ets,
                                                                             _Notes='Confined Conc.')

    OSDatabase.AddMaterial(Conc)
    Elastic = OpenSeesAPI.Model.Element.Material.UniaxialMaterial.Elastic(OSDatabase.GetFreeMaterialId(1, 1), 1.e-2)
    OSDatabase.AddMaterial(Elastic)
    Mat = OpenSeesAPI.Model.Element.Material.UniaxialMaterial.Parallel(OSDatabase.GetFreeMaterialId(1, 1),
                                                                       [Conc, Elastic])

    # Mat = OpenSeesAPI.Model.Element.Material.UniaxialMaterial.ConcreteCM(OSDatabase.GetFreeMaterialId(1, 1),
    #                                                                      -1 * fpc_confined, -1 * ec, Ec,
    #                                                                      7, e20u_reg, f_t, Ets, 1.2, 1000,
    #                                                                      _Notes='Confined Conc.', _ecu = e20u_reg)

    OSDatabase.AddMaterial(Mat)

    return Mat

def GetPughUnconfinedConcreteMaterial(OSDatabase, fpc, L_ip, fpuOfpc=0.2, EtsOEt = 0.05,
                                      GfcOfpc=2., Regularized = True):

    import OpenSeesAPI

    e0 = 2.*fpc/57./(fpc*1000)**.5              # Pugh Dissertation Equation 2.3

    fpuOfpc = fpuOfpc #Try this for now.

    f20 = fpuOfpc*fpc

    e20 = 0.008

    f_t = 4.0 *(fpc*1000)**.5 / 1000.

    Ec = 57.*(fpc*1000)**.5
    Et = Ec
    Ets = EtsOEt * Et

    lam = 0.1                       # Not Sure where this is mentioned in Pugh Thesis

    Gfc = GfcOfpc*fpc*0.0393701

    fpeak = fpc

    # e20u_reg = Gfc / 0.6 / fpeak / L_ip - 0.8 * fpeak / Ec + e0

    beta = fpuOfpc
    if Regularized:
        e20u_reg = 2. * Gfc / L_ip / fpeak / (1. + beta) + e0 * (beta + 1.0) / 2.  # Revised
    else:
        e20u_reg = 0.004  # Revised

    Conc = OpenSeesAPI.Model.Element.Material.UniaxialMaterial.Concrete02(OSDatabase.GetFreeMaterialId(1,1), -1*fpc,  -1*e0, -1*f20,  -1*e20u_reg, lam, f_t, Ets, _Notes='Unconfined Conc.')
    OSDatabase.AddMaterial(Conc)
    Elastic = OpenSeesAPI.Model.Element.Material.UniaxialMaterial.Elastic(OSDatabase.GetFreeMaterialId(1, 1), 1.e-2)
    OSDatabase.AddMaterial(Elastic)
    Mat = OpenSeesAPI.Model.Element.Material.UniaxialMaterial.Parallel(OSDatabase.GetFreeMaterialId(1, 1),
                                                                       [Conc, Elastic])

    # Mat = OpenSeesAPI.Model.Element.Material.UniaxialMaterial.ConcreteCM(OSDatabase.GetFreeMaterialId(1, 1),
    #                                                                      -1 * fpc, -1 * eu, Ec,
    #                                                                      7, e20u_reg, f_t, Ets, 1.2, 1000,
    #                                                                      _Notes='Confined Conc.', _ecu = e20u_reg)

    OSDatabase.AddMaterial(Mat)
    return Mat

def GetManderConfinedConcreteMaterial(OSDatabase, fpc, b_x, b_y, rho_l, s_bars_x, s_bars_y, s_ties_x, s_ties_y, f_y,
                                      A_tie, no_tie_x, no_tie_y, L_ip, GfccOGfc=1.7, fpuOfpc=0.2, EtsOEt = 0.05,
                                      GfcOfpc=2., ConcreteMaterialModel='Concrete02'):

    import OpenSeesAPI
    import math

    ksi_to_mpa = 6.8947

    # Mander et al 1988
    rho_cc = rho_l # ratio of area of longitudinal reinforcement to area of core of section
    b_c = (b_x)
    d_c = (b_y)
    A_c = b_c * d_c  #area of core of section enclosed by the center lines of the perimeter spiral or hoop

    # Compute the number of ineffectual parabolas
    A_i = (no_tie_x - 1.) * (s_bars_x ** 2.0) / 6.0 + (no_tie_y - 1.) * (s_bars_y ** 2.0) / 6.0
    d_tie = 2. * ( A_tie / np.pi ) ** .5
    sp = max(s_ties_x, s_ties_y) - d_tie
    A_e = (A_c - A_i) * ( 1 - sp / 2. / b_c) * ( 1 - sp / 2 / d_c)

    k_e = (1. - A_i / b_c / d_c) * ( 1 - sp / 2. / b_c ) * ( 1 - sp / 2. / d_c ) / ( 1. - rho_cc )

    A_sx = A_tie * no_tie_x
    A_sy = A_tie * no_tie_y

    rho_x = A_sx / s_ties_x / b_c
    rho_y = A_sy / s_ties_y / d_c

    f_yh = f_y
    f_lx = rho_x * f_yh
    f_ly = rho_y * f_yh

    fp_lx = k_e * rho_x * f_yh
    fp_ly = k_e * rho_y * f_yh

    # fp_cc = fpc * ( -1.254 + 2.254 * np.sqrt(1. + 7.94 * fp_l / fpc) - 2. * fp_l / fpc ) # For Circular Columns.. Incorrect

    def FindManderFccOFcRectCol(LargeConfiningStressRatio, SmallConfiningStressRatio):
        data = np.loadtxt('References/Mander/Mander1988Fig4DataPoints.dat')
        smallratiokeys = data[:, 0]

        if LargeConfiningStressRatio == SmallConfiningStressRatio:
            return np.interp(LargeConfiningStressRatio, data[-2, :], data[-1, :])

        def RemoveZeros(points):
            for i in range(len(points)):
                if points[i] == 0:
                    return points[:i]
            return points

        for i in range(0, len(data) - 4, 2):
            if SmallConfiningStressRatio >= smallratiokeys[i] and SmallConfiningStressRatio <= smallratiokeys[i + 2]:
                x1 = np.interp(LargeConfiningStressRatio, RemoveZeros(data[i + 0, 1:]), RemoveZeros(data[i + 1, 1:]))
                x2 = np.interp(LargeConfiningStressRatio, RemoveZeros(data[i + 2, 1:]), RemoveZeros(data[i + 3, 1:]))
                return np.interp(SmallConfiningStressRatio, [smallratiokeys[i], smallratiokeys[i + 2]], [x1, x2])

        return None

    def FindManderFccOFcRectColUsingEquations(LargeConfiningStressRatio, SmallConfiningStressRatio):
        if LargeConfiningStressRatio != SmallConfiningStressRatio:
            xbar = (LargeConfiningStressRatio + SmallConfiningStressRatio) * 0.5
            r = SmallConfiningStressRatio / LargeConfiningStressRatio
            A = 6.8886 - ( 0.6069 + 17.275 * r) * np.exp( -4.989 * r )
            B = 4.5 / (5. / A * (0.9849 - 0.6306 * np.exp(-3.8939 * r) ) - 0.1) - 5.

            K = 1. + A * xbar * (0.1 + 0.9 / (1. +  B * xbar))
        else:
            xbar = (LargeConfiningStressRatio + SmallConfiningStressRatio) * 0.5
            K = -1.254 + 2.254 * np.sqrt( 1. + 7.94 * xbar) - 2. * xbar

        return K

    # fp_cc = fpc * FindManderFccOFcRectCol(max(fp_lx, fp_ly)/fpc, min(fp_lx, fp_ly)/fpc)
    fp_cc = fpc * FindManderFccOFcRectColUsingEquations(max(fp_lx, fp_ly)/fpc, min(fp_lx, fp_ly)/fpc)

    e_co = 2. * fpc / (57. * (fpc * 1000.) ** .5) # Use this instead of 0.002 by Mander 1988
    xbar = (fp_lx + fp_ly) / 2.0 / fpc
    if 1. + 7.94 * xbar < 0:
        pass

    k_1 = ( 2.254 * ( np.sqrt( 1. + 7.94 * xbar) - 1.) - 2. * xbar ) / xbar
    k_2 = 5.0 * k_1 #5. * k_1
    e_cc = e_co * (1. + k_2 * xbar)

    fpc_confined = fp_cc

    Kc = fpc_confined / fpc

    ec = e_cc

    f20 = fpuOfpc * Kc * fpc

    # Compute Tension
    f_t = 4.0 * (fpc * 1000) ** .5 / 1000

    Ec = 57. * (fpc * 1000) ** .5
    Et = Ec
    Ets = EtsOEt * Et

    lam = 0.1  # Not Sure where this is mentioned in Pugh Thesis

    Gfc = GfcOfpc * fpc * 0.0393701
    Gfcc = GfccOGfc * Gfc

    fpeak = Kc * fpc

    # e20u_reg = Gfcc / 0.6 / fpeak / L_ip - 0.8 * fpeak / Ec + ec # old equation
    beta = fpuOfpc
    e20u_reg = 2. * Gfcc / L_ip / fpeak / (1. + beta) + ec * (beta + 1.0) /  2.   # Revised

    e20u_reg = max(ec + 0.0001, e20u_reg)

    if ConcreteMaterialModel == 'Concrete02':
        Conc = OpenSeesAPI.Model.Element.Material.UniaxialMaterial.Concrete02(OSDatabase.GetFreeMaterialId(1, 1),
                                                                             -1 * fpc_confined, -1 * ec, -1 * f20,
                                                                             -1 * e20u_reg, lam, f_t, Ets,
                                                                             _Notes='Confined Conc.')
    elif ConcreteMaterialModel == 'Concrete02IS':
        Conc = OpenSeesAPI.Model.Element.Material.UniaxialMaterial.Concrete02IS(OSDatabase.GetFreeMaterialId(1, 1), Ec,
                                                                             -1 * fpc_confined, -1 * ec, -1 * f20,
                                                                             -1 * e20u_reg, lam, f_t, Ets,
                                                                             _Notes='Confined Conc.')
    elif ConcreteMaterialModel == 'ConcreteCM':
        Conc = OpenSeesAPI.Model.Element.Material.UniaxialMaterial.ConcreteCM(OSDatabase.GetFreeMaterialId(1, 1),
                                                                             -1 * fpc_confined, -1 * ec, Ec, -1 * e20u_reg,
                                                                              Ec, f_t, f_t / Ets,
                                                                             _Notes='Confined Conc.')
    elif ConcreteMaterialModel == 'Concrete04':
        Conc = OpenSeesAPI.Model.Element.Material.UniaxialMaterial.Concrete04(OSDatabase.GetFreeMaterialId(1, 1),
                                                                             -1 * fpc_confined, -1 * ec, -1 * e20u_reg, Ec,
                                                                              f_t, f_t / Ets, 0.1,
                                                                             _Notes='Confined Conc.')

    OSDatabase.AddMaterial(Conc)
    Elastic = OpenSeesAPI.Model.Element.Material.UniaxialMaterial.Elastic(OSDatabase.GetFreeMaterialId(1, 1), 1.e-2)
    OSDatabase.AddMaterial(Elastic)
    Mat = OpenSeesAPI.Model.Element.Material.UniaxialMaterial.Parallel(OSDatabase.GetFreeMaterialId(1, 1),
                                                                       [Conc, Elastic])

    # Mat = OpenSeesAPI.Model.Element.Material.UniaxialMaterial.ConcreteCM(OSDatabase.GetFreeMaterialId(1, 1),
    #                                                                      -1 * fpc_confined, -1 * ec, Ec,
    #                                                                      7, e20u_reg, f_t, Ets, 1.2, 1000,
    #                                                                      _Notes='Confined Conc.', _ecu = e20u_reg)

    OSDatabase.AddMaterial(Mat)

    return Mat

def GetRichartConfinedConcreteMaterial(OSDatabase, fpc, b_x, b_y, rho_l, s_bars_x, s_bars_y, s_ties_x, s_ties_y, f_y,
                                      A_tie, no_tie_x, no_tie_y, L_ip, GfccOGfc=1.7, fpuOfpc=0.2, EtsOEt = 0.05,
                                      GfcOfpc=2., ConcreteMaterialModel='Concrete02'):

    import OpenSeesAPI
    import math

    ksi_to_mpa = 6.8947

    # Mander et al 1988
    rho_cc = rho_l # ratio of area of longitudinal reinforcement to area of core of section
    b_c = (b_x)
    d_c = (b_y)
    A_c = b_c * d_c  #area of core of section enclosed by the center lines of the perimeter spiral or hoop

    # Compute the number of ineffectual parabolas
    A_i = (no_tie_x - 1.) * (s_bars_x ** 2.0) / 6.0 + (no_tie_y - 1.) * (s_bars_y ** 2.0) / 6.0
    d_tie = 2. * ( A_tie / np.pi ) ** .5
    sp = max(s_ties_x, s_ties_y) - d_tie
    A_e = (A_c - A_i) * ( 1 - sp / 2. / b_c) * ( 1 - sp / 2 / d_c)

    k_e = (1. - A_i / b_c / d_c) * ( 1 - sp / 2. / b_c ) * ( 1 - sp / 2. / d_c ) / ( 1. - rho_cc )

    A_sx = A_tie * no_tie_x
    A_sy = A_tie * no_tie_y

    rho_x = A_sx / s_ties_x / b_c
    rho_y = A_sy / s_ties_y / d_c

    f_yh = f_y
    f_lx = rho_x * f_yh
    f_ly = rho_y * f_yh

    fp_lx = k_e * rho_x * f_yh
    fp_ly = k_e * rho_y * f_yh

    fp_cc = fpc + (fp_lx + fp_ly) / 2. * 4.1

    e_co = 2. * fpc / (57. * (fpc * 1000.) ** .5) # Use this instead of 0.002 by Mander 1988
    xbar = (fp_lx + fp_ly) / 2.0 / fpc

    k_1 = 4.1
    k_2 = 5.0 * k_1
    e_cc = e_co * (1. + k_2 * xbar)

    fpc_confined = fp_cc

    Kc = fpc_confined / fpc

    ec = e_cc

    f20 = fpuOfpc * Kc * fpc

    f_t = 4.0 * (fpc * 1000) ** .5 / 1000

    Ec = 57. * (fpc * 1000) ** .5
    Et = Ec
    Ets = EtsOEt * Et

    lam = 0.1  # Not Sure where this is mentioned in Pugh Thesis

    Gfc = GfcOfpc * fpc * 0.0393701
    Gfcc = GfccOGfc * Gfc

    fpeak = Kc * fpc

    # e20u_reg = Gfcc / 0.6 / fpeak / L_ip - 0.8 * fpeak / Ec + ec
    beta = fpuOfpc
    e20u_reg = 2. * Gfcc / L_ip / fpeak / (1. + beta) + ec * (beta + 1.0) /  2.   # Revised

    e20u_reg = max(ec + 0.0001, e20u_reg)


    if ConcreteMaterialModel == 'Concrete02':
        Conc = OpenSeesAPI.Model.Element.Material.UniaxialMaterial.Concrete02(OSDatabase.GetFreeMaterialId(1, 1),
                                                                             -1 * fpc_confined, -1 * ec, -1 * f20,
                                                                             -1 * e20u_reg, lam, f_t, Ets,
                                                                             _Notes='Confined Conc.')
    elif ConcreteMaterialModel == 'Concrete02IS':
        Conc = OpenSeesAPI.Model.Element.Material.UniaxialMaterial.Concrete02IS(OSDatabase.GetFreeMaterialId(1, 1), Ec,
                                                                             -1 * fpc_confined, -1 * ec, -1 * f20,
                                                                             -1 * e20u_reg, lam, f_t, Ets,
                                                                             _Notes='Confined Conc.')
    elif ConcreteMaterialModel == 'ConcreteCM':
        Conc = OpenSeesAPI.Model.Element.Material.UniaxialMaterial.ConcreteCM(OSDatabase.GetFreeMaterialId(1, 1),
                                                                             -1 * fpc_confined, -1 * ec, Ec, -1 * e20u_reg,
                                                                              Ec, f_t, f_t / Ets,
                                                                             _Notes='Confined Conc.')
    elif ConcreteMaterialModel == 'Concrete04':
        Conc = OpenSeesAPI.Model.Element.Material.UniaxialMaterial.Concrete04(OSDatabase.GetFreeMaterialId(1, 1),
                                                                             -1 * fpc_confined, -1 * ec, -1 * e20u_reg, Ec,
                                                                              f_t, f_t / Ets, 0.1,
                                                                             _Notes='Confined Conc.')

    OSDatabase.AddMaterial(Conc)
    Elastic = OpenSeesAPI.Model.Element.Material.UniaxialMaterial.Elastic(OSDatabase.GetFreeMaterialId(1, 1), 1.e-2)
    OSDatabase.AddMaterial(Elastic)
    Mat = OpenSeesAPI.Model.Element.Material.UniaxialMaterial.Parallel(OSDatabase.GetFreeMaterialId(1, 1),
                                                                       [Conc, Elastic])

    # Mat = OpenSeesAPI.Model.Element.Material.UniaxialMaterial.ConcreteCM(OSDatabase.GetFreeMaterialId(1, 1),
    #                                                                      -1 * fpc_confined, -1 * ec, Ec,
    #                                                                      7, e20u_reg, f_t, Ets, 1.2, 1000,
    #                                                                      _Notes='Confined Conc.', _ecu = e20u_reg)

    OSDatabase.AddMaterial(Mat)

    return Mat

def GetPughSteel02Material(OSDatabase, fy, fu, L_ip, eps_ult_comp, eps_ult_exp = .12, L_gage = 8,
                           SteelPostYieldStiffness=0.00784,):
    import OpenSeesAPI
    #L_gage = 8 # As per Pugh et al 2015 Paper, Use this when the gage length is not available
    RegularizeHardening = False
    if RegularizeHardening:
        Es = 29000.
        eps_y = fy/Es
        eps_ult = eps_y + L_gage/L_ip * (eps_ult_exp - eps_y)

        b_reg = (fu-fy)/(Es*(eps_ult - eps_y))

        Mat1 = OpenSeesAPI.Model.Element.Material.UniaxialMaterial.Steel02(OSDatabase.GetFreeMaterialId(1,1), fy, Es, b_reg, 20.0, 0.925, 0.15)
        OSDatabase.AddMaterial(Mat1)
        Mat2 = OpenSeesAPI.Model.Element.Material.UniaxialMaterial.MinMax(OSDatabase.GetFreeMaterialId(1,1), Mat1, eps_ult_comp, eps_ult)
        OSDatabase.AddMaterial(Mat2)
        Elastic = OpenSeesAPI.Model.Element.Material.UniaxialMaterial.Elastic(OSDatabase.GetFreeMaterialId(1, 1), 1.e-2)
        OSDatabase.AddMaterial(Elastic)
        Mat3 = OpenSeesAPI.Model.Element.Material.UniaxialMaterial.Parallel(OSDatabase.GetFreeMaterialId(1, 1),
                                                                           [Mat2, Elastic])
        OSDatabase.AddMaterial(Mat3)
    else:
        Es = 29000.
        eps_y = fy / Es

        # b = (fu - fy) / (eps_ult_exp - eps_y) / Es
        # def minimizer(x):
        #      return (fu + fy) / 2. * (eps_ult_exp - eps_y) * L_gage / L_ip - fy * (x[0] - eps_y) - b * Es / 2. * (x[0] - eps_y) ** 2.
        #
        # from scipy.optimize import minimize
        # res = minimize(minimizer, [0.2], method='nelder-mead', options = {'xtol': 1e-8, 'disp': True})

        eps_ult = eps_y + L_gage/L_ip * (eps_ult_exp - eps_y)#res.x[0]
        b_reg = SteelPostYieldStiffness

        Mat1 = OpenSeesAPI.Model.Element.Material.UniaxialMaterial.Steel02(OSDatabase.GetFreeMaterialId(1, 1), fy, Es,
                                                                           b_reg, 20.0, 0.925, 0.15)
        OSDatabase.AddMaterial(Mat1)
        Mat2 = OpenSeesAPI.Model.Element.Material.UniaxialMaterial.MinMax(OSDatabase.GetFreeMaterialId(1, 1), Mat1,
                                                                          eps_ult_comp, eps_ult)
        OSDatabase.AddMaterial(Mat2)
        Elastic = OpenSeesAPI.Model.Element.Material.UniaxialMaterial.Elastic(OSDatabase.GetFreeMaterialId(1, 1), 1.e-2)
        OSDatabase.AddMaterial(Elastic)
        Mat3 = OpenSeesAPI.Model.Element.Material.UniaxialMaterial.Parallel(OSDatabase.GetFreeMaterialId(1, 1),
                                                                            [Mat2, Elastic])
        OSDatabase.AddMaterial(Mat3)

    return Mat3

def GetPughSteel02MaterialUnregularized(OSDatabase, fy, fu, eps_ult_comp, eps_ult_exp = .12,
                                        SteelPostYieldStiffness=0.00784,):
    import OpenSeesAPI
    #L_gage = 8 # As per Pugh et al 2015 Paper, Use this when the gage length is not available

    Es = 29000.
    eps_y = fy/Es
    eps_ult = eps_ult_exp

    b_reg = SteelPostYieldStiffness#(fu-fy)/(Es*(eps_ult - eps_y))

    Mat1 = OpenSeesAPI.Model.Element.Material.UniaxialMaterial.Steel02(OSDatabase.GetFreeMaterialId(1,1), fy, Es, b_reg, 20.0, 0.925, 0.15)
    OSDatabase.AddMaterial(Mat1)
    Mat2 = OpenSeesAPI.Model.Element.Material.UniaxialMaterial.MinMax(OSDatabase.GetFreeMaterialId(1,1), Mat1, eps_ult_comp, eps_ult)
    OSDatabase.AddMaterial(Mat2)
    Elastic = OpenSeesAPI.Model.Element.Material.UniaxialMaterial.Elastic(OSDatabase.GetFreeMaterialId(1, 1), 1.e-2)
    OSDatabase.AddMaterial(Elastic)
    Mat3 = OpenSeesAPI.Model.Element.Material.UniaxialMaterial.Parallel(OSDatabase.GetFreeMaterialId(1, 1),
                                                                       [Mat2, Elastic])
    OSDatabase.AddMaterial(Mat3)

    return Mat3

def GetSteel02FatigueMaterialUnregularized(OSDatabase, fy, fu, eps_ult_comp, eps_ult_exp = .12,
                                        SteelPostYieldStiffness=0.00784,):
    if 'Steel02FatigueCd' in OSDatabase._OtherOptions:
        Cd = OSDatabase._OtherOptions['Steel02FatigueCd']
    else:
        Cd = 0.2
    if 'Steel02FatigueCf' in OSDatabase._OtherOptions:
        Cf = OSDatabase._OtherOptions['Steel02FatigueCf']
    else:
        Cf = 0.12
    if 'Steel02FatigueAlpha' in OSDatabase._OtherOptions:
        Alpha = OSDatabase._OtherOptions['Steel02FatigueAlpha']
    else:
        Alpha = 0.44
    if 'Steel02FatigueBeta' in OSDatabase._OtherOptions:
        Beta = OSDatabase._OtherOptions['Steel02FatigueBeta']
    else:
        Beta = 0.45

    import OpenSeesAPI
    #L_gage = 8 # As per Pugh et al 2015 Paper, Use this when the gage length is not available

    Es = 29000.
    eps_y = fy/Es
    eps_ult = eps_ult_exp

    b_reg = SteelPostYieldStiffness #(fu-fy)/(Es*(eps_ult - eps_y))

    # Cd = 0.20
    # Cf = 0.12
    # Alpha = 0.44
    # Beta = 0.45

    Mat1 = OpenSeesAPI.Model.Element.Material.UniaxialMaterial.Steel02Fatigue(OSDatabase.GetFreeMaterialId(1,1), fy, Es, b_reg, Cd, Cf, Alpha, Beta, -1e6, 1e6)
    OSDatabase.AddMaterial(Mat1)
    Mat2 = OpenSeesAPI.Model.Element.Material.UniaxialMaterial.MinMax(OSDatabase.GetFreeMaterialId(1,1), Mat1, eps_ult_comp, eps_ult)
    OSDatabase.AddMaterial(Mat2)
    Elastic = OpenSeesAPI.Model.Element.Material.UniaxialMaterial.Elastic(OSDatabase.GetFreeMaterialId(1, 1), 1.e-2)
    OSDatabase.AddMaterial(Elastic)
    Mat3 = OpenSeesAPI.Model.Element.Material.UniaxialMaterial.Parallel(OSDatabase.GetFreeMaterialId(1, 1),
                                                                       [Mat2, Elastic])
    OSDatabase.AddMaterial(Mat3)

    return Mat3

def GetPughSteel02MaterialUnregularizedWithFatigue(OSDatabase, fy, fu, eps_ult_comp, eps_ult_exp = .12,
                                        SteelPostYieldStiffness=0.00784,):
    import OpenSeesAPI
    #L_gage = 8 # As per Pugh et al 2015 Paper, Use this when the gage length is not available

    Es = 29000.
    eps_y = fy/Es
    eps_ult = eps_ult_exp

    b_reg = SteelPostYieldStiffness#(fu-fy)/(Es*(eps_ult - eps_y))

    Mat1 = OpenSeesAPI.Model.Element.Material.UniaxialMaterial.Steel02(OSDatabase.GetFreeMaterialId(1,1), fy, Es, b_reg, 20.0, 0.925, 0.15)
    OSDatabase.AddMaterial(Mat1)
    Mat1b = OpenSeesAPI.Model.Element.Material.UniaxialMaterial.Fatigue(OSDatabase.GetFreeMaterialId(1, 1), Mat1, 0.12, -0.44)
    OSDatabase.AddMaterial(Mat1b)
    Mat2 = OpenSeesAPI.Model.Element.Material.UniaxialMaterial.MinMax(OSDatabase.GetFreeMaterialId(1,1), Mat1b, eps_ult_comp, eps_ult)
    OSDatabase.AddMaterial(Mat2)
    Elastic = OpenSeesAPI.Model.Element.Material.UniaxialMaterial.Elastic(OSDatabase.GetFreeMaterialId(1, 1), 1.e-2)
    OSDatabase.AddMaterial(Elastic)
    Mat3 = OpenSeesAPI.Model.Element.Material.UniaxialMaterial.Parallel(OSDatabase.GetFreeMaterialId(1, 1),
                                                                       [Mat2, Elastic])
    OSDatabase.AddMaterial(Mat3)

    return Mat3

def GetReinforcingSteelMaterialUnregularized(OSDatabase, fy, fu, eps_ult_exp = .12,
                                             stirrup_spacing = 6, long_bar_size = 4,
                                             Cf = 0.26, a = .506, Cd = 0.389, # Cf = 0.12, a = .44, Cd = 0.35,
                                             withBuckling = False, withFatigue = False,
                                             eps_ult_comp = -0.12,
                                             SteelPostYieldStiffness=0.00784,):
    import OpenSeesAPI
    #L_gage = 8 # As per Pugh et al 2015 Paper, Use this when the gage length is not available

    Es = 29000.
    Esh = 0.2 * Es
    esh = 0.02
    eult = 0.12
    Isr = float(stirrup_spacing) / ( long_bar_size / 8. )

    alpha = 1.0

    if withBuckling == False and withFatigue == False:
        Mat1 = OpenSeesAPI.Model.Element.Material.UniaxialMaterial.ReinforcingSteel.ReinforcingSteel(OSDatabase.GetFreeMaterialId(1, 1), fy,
                                                                                            fu, Es, Esh, esh, eult
                                                                                            )
        Mat2 = OpenSeesAPI.Model.Element.Material.UniaxialMaterial.MinMax(OSDatabase.GetFreeMaterialId(1, 1), Mat1,
                                                                          eps_ult_comp, eps_ult_exp)
    elif withFatigue == False and withBuckling == True:
        Mat1 = OpenSeesAPI.Model.Element.Material.UniaxialMaterial.ReinforcingSteel.DMBuck(OSDatabase.GetFreeMaterialId(1, 1), fy,
                                                                                            fu, Es, Esh, esh, eult, Isr, alpha,
                                                                                            Optional='')
        Mat2 = OpenSeesAPI.Model.Element.Material.UniaxialMaterial.MinMax(OSDatabase.GetFreeMaterialId(1, 1), Mat1,
                                                                          -1 * eps_ult_exp, eps_ult_exp)
    elif withFatigue == True and withBuckling == True:
        Mat1 = OpenSeesAPI.Model.Element.Material.UniaxialMaterial.ReinforcingSteel.DMBuck(OSDatabase.GetFreeMaterialId(1,1),
                                                                                    fy, fu, Es, Esh, esh, eult, Isr, alpha,
                                                                                    Cf = Cf, a = a, Cd = Cd,
                                                                                    Optional= '')
        Mat2 = OpenSeesAPI.Model.Element.Material.UniaxialMaterial.MinMax(OSDatabase.GetFreeMaterialId(1, 1), Mat1,
                                                                          -1 * eps_ult_exp, eps_ult_exp)
    elif withFatigue == True and withBuckling == False:
        Mat1 = OpenSeesAPI.Model.Element.Material.UniaxialMaterial.ReinforcingSteel.ReinforcingSteel(OSDatabase.GetFreeMaterialId(1,1),
                                                                                    fy, fu, Es, Esh, esh, eult,
                                                                                    Optional= ' -CMFatigue %f %f %f '%(Cf, a, Cd))
        Mat2 = OpenSeesAPI.Model.Element.Material.UniaxialMaterial.MinMax(OSDatabase.GetFreeMaterialId(1, 1), Mat1,
                                                                          -1 * eps_ult_exp, eps_ult_exp)

    OSDatabase.AddMaterial(Mat1)
    OSDatabase.AddMaterial(Mat2)

    Elastic = OpenSeesAPI.Model.Element.Material.UniaxialMaterial.Elastic(OSDatabase.GetFreeMaterialId(1, 1), 1.e-2)
    OSDatabase.AddMaterial(Elastic)
    Mat3 = OpenSeesAPI.Model.Element.Material.UniaxialMaterial.Parallel(OSDatabase.GetFreeMaterialId(1, 1),
                                                                       [Mat2, Elastic])
    OSDatabase.AddMaterial(Mat3)

    return Mat3

def GetPinching4MaterialUnregularized(OSDatabase, fy, fu, l_be, eps_ult_exp = .12,
                                             stirrup_spacing = 6, long_bar_size = 8.,
                                             eps_ult_comp = -0.12,
                                      SteelPostYieldStiffness=0.00784, stirrup_dia = 4.,
                                      no_of_legs_y = 3, no_of_rows = 3,
                                      cover = 3, WebRebar = False):
    import OpenSeesAPI
    #Define Rotational Spring
    Es = 29000.
    eps_y = fy / Es
    eps_ult = eps_ult_exp
    Esh = (fu - fy) * ( eps_ult - eps_y)
    b = SteelPostYieldStiffness#(fu - fy) / (eps_ult - eps_y) / Es

    L = stirrup_spacing
    D = long_bar_size/8.

    if not(WebRebar):
        stress, strain = GetMaekawaInelasticBarBuckling(fy, fu, eps_ult_exp, stirrup_spacing, long_bar_size,
                                                            stirrup_dia, l_be,
                                       no_of_legs_y, no_of_rows, cover = cover)
    else:
        stress, strain = GetMaekawaInelasticBarBuckling(fy, fu, eps_ult_exp, stirrup_spacing, long_bar_size,
                                                            stirrup_dia, l_be,
                                                            no_of_legs_y, no_of_rows, cover=cover, UseMax=True)
    # stress = [0.0, fy, sigma_star, sigma_res]
    # strain = [0.0, eps_y, eps_star, eps_res]

    Ks = 0.02 * Es

    ePf1 = fy
    ePd1 = eps_y
    ePf2 = b * (eps_ult - eps_y) * Es + fy
    ePd2 = eps_ult
    ePf3 = 0.01 * fy
    ePd3 = eps_ult + 0.001
    ePf4 = 0.01 * fy
    ePd4 = eps_ult + 0.1

    eNf1 = -stress[1]
    eNd1 = -strain[1]
    eNf2 = -stress[2]
    eNd2 = -strain[2]
    eNf3 = -stress[3]
    eNd3 = -strain[3]
    eNf4 = -stress[3]
    eNd4 = -strain[3] - 1.0

    rDispP = 1.0
    rForceP = 1.0
    uForceP = 0.0
    gK1 = 1.
    gK2 = 1.
    gK3 = 0.
    gK4 = 0.
    gKLim = 0.001
    gD1 = 1.
    gD2 = 1.
    gD3 = 0.
    gD4 = 0.
    gDLim = 0.001
    gF1 = 1.
    gF2 = 1.
    gF3 = 0.
    gF4 = 0.
    gFLim = 0.001
    gE = 1.0
    dmgType = 'cycle'

    rDispN = 1.0
    rForceN = 1.0
    uForceN = 0.0

    Mat1 = OpenSeesAPI.Model.Element.Material.UniaxialMaterial.Pinching4(OSDatabase.GetFreeMaterialId(1, 1), ePf1, ePd1, ePf2, ePd2, ePf3, ePd3, ePf4, ePd4, rDispP, rForceP, uForceP, gK1, gK2, gK3, gK4, gKLim, gD1, gD2, gD3, gD4, gDLim, gF1, gF2, gF3, gF4, gFLim, gE, dmgType, eNf1, eNd1, eNf2, eNd2, eNf3, eNd3, eNf4, eNd4, rDispN , rForceN , uForceN )
    Mat2 = OpenSeesAPI.Model.Element.Material.UniaxialMaterial.MinMax(OSDatabase.GetFreeMaterialId(1, 1), Mat1,
                                                                      eps_ult_comp, eps_ult_exp)

    OSDatabase.AddMaterial(Mat1)
    OSDatabase.AddMaterial(Mat2)

    Elastic = OpenSeesAPI.Model.Element.Material.UniaxialMaterial.Elastic(OSDatabase.GetFreeMaterialId(1, 1), 1.e-2)
    OSDatabase.AddMaterial(Elastic)
    Mat3 = OpenSeesAPI.Model.Element.Material.UniaxialMaterial.Parallel(OSDatabase.GetFreeMaterialId(1, 1),
                                                                       [Mat2, Elastic])
    OSDatabase.AddMaterial(Mat3)

    return Mat3

def GetHystereticMaterialUnregularized(OSDatabase, fy, fu, l_be, eps_ult_exp = .12,
                                             stirrup_spacing = 6, long_bar_size = 8.,
                                             eps_ult_comp = -0.12,
                                      SteelPostYieldStiffness=0.00784, stirrup_dia = 4.,
                                      no_of_legs_y = 3, no_of_rows = 3,
                                      cover = 3, WebRebar = False):
    import OpenSeesAPI
    #Define Rotational Spring
    Es = 29000.
    eps_y = fy / Es
    eps_ult = eps_ult_exp
    Esh = (fu - fy) * ( eps_ult - eps_y)
    b = SteelPostYieldStiffness#(fu - fy) / (eps_ult - eps_y) / Es

    L = stirrup_spacing
    D = long_bar_size/8.

    if not(WebRebar):
        stress, strain = GetMaekawaInelasticBarBuckling(fy, fu, eps_ult_exp, stirrup_spacing, long_bar_size,
                                                            stirrup_dia, l_be,
                                       no_of_legs_y, no_of_rows, cover = cover)
    else:
        stress, strain = GetMaekawaInelasticBarBuckling(fy, fu, eps_ult_exp, stirrup_spacing, long_bar_size,
                                                            stirrup_dia, l_be,
                                                            no_of_legs_y, no_of_rows, cover=cover, UseMax=True)
    # stress = [0.0, fy, sigma_star, sigma_res]
    # strain = [0.0, eps_y, eps_star, eps_res]

    Ks = 0.02 * Es

    ePf1 = fy
    ePd1 = eps_y
    ePf2 = b * (eps_ult - eps_y) * Es + fy
    ePd2 = eps_ult
    ePf3 = 0.01 * fy
    ePd3 = eps_ult + 0.001
    ePf4 = 0.01 * fy
    ePd4 = eps_ult + 0.1

    eNf1 = -stress[1]
    eNd1 = -strain[1]
    eNf2 = -stress[2]
    eNd2 = -strain[2]
    eNf3 = -stress[3]
    eNd3 = -strain[3]
    eNf4 = -stress[3]
    eNd4 = -strain[3] - 1.0

    Mat1 = OpenSeesAPI.Model.Element.Material.UniaxialMaterial.Hysteretic(OSDatabase.GetFreeMaterialId(1, 1), ePf1, ePd1, ePf2, ePd2, eNf1, eNd1, eNf3, eNd3,
                                                                          1., 1., 0.0, 0.0, ePf3, ePd3, eNf4, eNd4
                                                                          )
    Mat2 = OpenSeesAPI.Model.Element.Material.UniaxialMaterial.MinMax(OSDatabase.GetFreeMaterialId(1, 1), Mat1,
                                                                      eps_ult_comp, eps_ult_exp)

    OSDatabase.AddMaterial(Mat1)
    OSDatabase.AddMaterial(Mat2)

    Elastic = OpenSeesAPI.Model.Element.Material.UniaxialMaterial.Elastic(OSDatabase.GetFreeMaterialId(1, 1), 1.e-2)
    OSDatabase.AddMaterial(Elastic)
    Mat3 = OpenSeesAPI.Model.Element.Material.UniaxialMaterial.Parallel(OSDatabase.GetFreeMaterialId(1, 1),
                                                                       [Mat2, Elastic])
    OSDatabase.AddMaterial(Mat3)

    return Mat3

def GetSteel02LBMaterialUnregularized(OSDatabase, fy, fu, l_be, eps_ult_exp = .12,
                                     stirrup_spacing = 6, long_bar_size = 8.,
                                     eps_ult_comp = -0.12,
                                      SteelPostYieldStiffness=0.00784, stirrup_dia = 4.,
                                      no_of_legs_y = 3, no_of_rows = 3,
                                      cover = 3, WebRebar = False):
    import OpenSeesAPI
    #Define Rotational Spring
    Es = 29000.
    eps_y = fy / Es
    eps_ult = eps_ult_exp
    Esh = (fu - fy) * ( eps_ult - eps_y)
    b = SteelPostYieldStiffness#(fu - fy) / (eps_ult - eps_y) / Es

    L = stirrup_spacing
    D = long_bar_size/8.

    if not(WebRebar):
        stress, strain = GetMaekawaInelasticBarBuckling(fy, fu, eps_ult_exp, stirrup_spacing, long_bar_size,
                                                            stirrup_dia, l_be,
                                       no_of_legs_y, no_of_rows, cover = cover)
    else:
        stress, strain = GetMaekawaInelasticBarBuckling(fy, fu, eps_ult_exp, stirrup_spacing, long_bar_size,
                                                            stirrup_dia, l_be,
                                                            no_of_legs_y, no_of_rows, cover=cover, UseMax=True)
    # stress = [0.0, fy, sigma_star, sigma_res]
    # strain = [0.0, eps_y, eps_star, eps_res]

    eps_star = strain[2]
    sigma_star = stress[2]

    eps_res = strain[3]
    sigma_res = stress[3]

    Mat1 = OpenSeesAPI.Model.Element.Material.UniaxialMaterial.Steel02LB(OSDatabase.GetFreeMaterialId(1, 1), fy, Es, b,
                                                                         0.02 * Es,  strain[1], sigma_res)
    Mat2 = OpenSeesAPI.Model.Element.Material.UniaxialMaterial.MinMax(OSDatabase.GetFreeMaterialId(1, 1), Mat1,
                                                                      eps_ult_comp, eps_ult_exp)

    OSDatabase.AddMaterial(Mat1)
    OSDatabase.AddMaterial(Mat2)

    Elastic = OpenSeesAPI.Model.Element.Material.UniaxialMaterial.Elastic(OSDatabase.GetFreeMaterialId(1, 1), 1.e-2)
    OSDatabase.AddMaterial(Elastic)
    Mat3 = OpenSeesAPI.Model.Element.Material.UniaxialMaterial.Parallel(OSDatabase.GetFreeMaterialId(1, 1),
                                                                       [Mat2, Elastic])
    OSDatabase.AddMaterial(Mat3)

    return Mat3

def GetL_IP(NoOfIntPoints, IntScheme = 'Lobatto', ForceBased=True):
    if ForceBased:
        if IntScheme == 'Lobatto':
            if NoOfIntPoints == 3:
                L_ip = [0.1666667, 0.6666667, 0.1666667]
            elif NoOfIntPoints == 5:
                L_ip = [0.05000000, 0.27222222, 0.35555556, 0.27222222, 0.05000000]
            else:
                L_ip = [0.023810, 0.138413, 0.215873, 0.243810, 0.215873, 0.138413, 0.023810]
        else:
            if NoOfIntPoints == 3:
                L_ip = [0.27778, 0.4444445, 0.27778]
            elif NoOfIntPoints == 5:
                L_ip = [0.1184635, 0.2393145, 0.2844445, 0.2393145, 0.1184635]
            else:
                L_ip = np.array([0.1294849662, 0.2797053915, 0.3818300505, 0.4179591837, 0.3818300505, 0.2797053915, 0.1294849662,])/2. # Find Points
    else:
        if NoOfIntPoints == 3:
            L_ip = np.ones(NoOfIntPoints)
        elif NoOfIntPoints == 5:
            L_ip = np.ones(NoOfIntPoints)
        else:
            L_ip = np.ones(NoOfIntPoints)

    return L_ip

def GetMaekawaInelasticBarBuckling(fy=60., fu=105., eps_ult_exp=0.12,
                                   stirrup_spacing=3., bar_dia=8.,
                                   stirrup_dia=4., l_be=30, no_of_legs=2,
                                   no_of_be_rows=2, cover=3,
                                   UseMax=False
                                   ):
    # Based on the Kashani et al. 2016 Paper, Implemented by Kamal Ahmed,

    import numpy as np
    ksi_to_mpa = 4.44 / 25.4 ** 2. * 1000
    in_to_mm = 25.4

    fy = fy * ksi_to_mpa  # 543.0 MPa
    fu = fu * ksi_to_mpa  # 677.0  # MPa
    # eps_ult_exp = .124
    stirrup_spacing = stirrup_spacing * in_to_mm  # mm

    db = bar_dia / 8. * in_to_mm  # 16.0   # dia of long bars [mm]
    dt = stirrup_dia / 8. * in_to_mm  # 6.0    # dia of confining bars [mm]

    Ab = np.pi * db ** 2 / 4.0
    At = np.pi * dt ** 2 / 4.0

    Es = 29000. * ksi_to_mpa  # 181000.    # MPa
    Et = 29000. * ksi_to_mpa  # 231176.    # MPa

    eps_y = fy/Es#0.003
    eps_ult = eps_ult_exp
    Esh = (fu - fy) * (eps_ult - eps_y)
    # b = SteelPostYieldStiffness#(fu - fy) / (eps_ult - eps_y) / Es

    le = (l_be - 2. * cover) * in_to_mm  # mm
    nl = no_of_legs  # no. of legs
    nb = no_of_be_rows  # no. of bars

    I = np.pi * db ** 4 / 64.0
    EI = 0.5 * Es * I * np.sqrt(fy / 400.0)  # flexural rigidity

    Kt = (Et * At / le) * (nl / nb)
    ci = 1.0
    n = 1
    while n < 50:
        ter2, ter2b = 0.0, 0.0
        ter1 = 2 * np.pi ** 4 * EI / (n ** 3 * stirrup_spacing ** 3)
        ter1b = 2 * np.pi ** 4 * EI / ((n + 1) ** 3 * stirrup_spacing ** 3)

        for j in range(1, n + 1):
            ter2 += ci / 4 * (1 - np.cos(2 * j * np.pi / n))
            ter2b += ci / 4 * (1 - np.cos(2 * j * np.pi / (n + 1)))

        ter3 = -np.pi ** 2 / (2 * n * stirrup_spacing)
        ter3b = -np.pi ** 2 / (2 * (n + 1) * stirrup_spacing)

        sol = np.linalg.solve(np.array([[ter2, ter3], [ter2b, ter3b]]), np.array([-ter1, -ter1b]))
        Kn = sol[0]
        # print('n = %d: Kt = %0.1f, Kn = %0.1f' % (n, Kt, Kn))

        if Kn <= Kt:
            break
        else:
            n += 1
            continue

    L = stirrup_spacing * n  # effective length
    D = db
    LOD = L / float(D)

    if UseMax == True:
        LOD = 15.

    eps_star = eps_y * (55. - 2.3 * np.sqrt(fy / 100.) * LOD)
    eps_star = max(eps_star, 7. * eps_y)

    sigma_star_l = fy + Esh * (eps_star - eps_y)
    alpha = 0.75
    sigma_star = sigma_star_l * alpha * (1.1 - 0.016 * np.sqrt(fy / 100) * LOD)
    sigma_res = 0.2 * fy

    eps_res = (sigma_star - sigma_res) / 0.02 / Es + eps_star

    stress = [0.0, fy, sigma_star, sigma_res]
    stress_ksi = np.array(stress) / ksi_to_mpa
    strain = [0.0, eps_y, eps_star, eps_res]

    return stress_ksi, strain

def GetAxialLoadRatio(Archetype):
    import ATCWallArchetypeObjects as ATCWallObjects
    Section = Archetype.CustomSection
    WallAxialLoadRatio = []
    if isinstance(Section[0], ATCWallObjects.PlanarWallSection):
        for i in range(len(Section)):
            fpcAg = Archetype.fpc * Section[i].l_w * Section[i].t_w
            WallAxialLoadRatio.append(np.sum(Archetype.WallGravityLoad[i:]) /  fpcAg)

    elif isinstance(Section[0], ATCWallObjects.IWallSection):
        for i in range(len(Section)):
            fpcAg = Archetype.fpc * ( Section[i].l_w * Section[i].t_w +
                                      Section[i].b_w * Section[i].t_w * 2 -
                                      Section[i].t_w * Section[i].t_w * 2. )
            WallAxialLoadRatio.append(np.sum(Archetype.WallGravityLoad[i:]) / fpcAg)

    return WallAxialLoadRatio

####################################################################################
#endregion
####################################################################################

