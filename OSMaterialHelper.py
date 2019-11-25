from __future__ import absolute_import
from __future__ import print_function
from six.moves import range
__author__ = 'marafi'

def GetIbarraPeakOrientedModel(id,b,h,dp,fy,fpc,rho,rhop,rho_v,N,s,LsOH,L,Ig,a_sl=1,Imperial=True,t_slab=None, rho_slab=None, eff_wdith=None,_Notes=None):
    import OpenSeesAPI

    StiffnessFactor1 = 11
    StiffnessFactor2 = 1.1

    #a_sl=0 if bond slip of long. bars is prevented

    #LsOH #Column Height to Depth
    fpc_exp = fpc*1.25

    d = h - dp

    v = float(N)/(fpc*b*h) #Checked with Haselton Spreadsheet

    if Imperial:
        c_units = 6.9 # For Imperial
    else:
        c_units = 1 # For Metric

    dp_col = 1.25
    s_n = (float(s)/dp_col)*(fy*c_units/100.)**0.5

    rho_sh = rho_v
    theta_p = 0.12*(1.+0.55*a_sl)*(0.16)**v*(0.02+40*rho_sh)**0.43*0.54**(0.01*c_units*fpc)*0.66**(0.1*s_n)*2.27**(10.0*(rho+rhop)) #Checked with Haselton Spreadsheet
    # theta_p = 0.13*(1.+0.55*a_sl)*(0.13)**v*(0.02+40*rho_sh)**0.65*0.57**(0.01*c_units*fpc)

    if rhop != rho:
        theta_p_pos = (max(0.01,rho*fy/fpc)/max(0.01,rhop*fy/fpc))**0.225*theta_p #Checked with Haselton Spreadsheet
        theta_p_neg = (max(0.01,rhop*fy/fpc)/max(0.01,rho*fy/fpc))**0.225*theta_p #Checked with Haselton Spreadsheet
    else:
        theta_p_pos = theta_p
        theta_p_neg = theta_p

    theta_pc = min(0.76*0.031**v*(0.02+40.*rho_sh)*1.02,0.10) #Checked with Haselton Spreadsheet

    # Lambda = 30*0.03**v
    Lambda = 170.7*0.27**v*(0.1)**(s/d) #Checked with Haselton Spreadsheet
    Lambda = Lambda*StiffnessFactor1*0.35/0.2 #To Adjust for Stiffness (Replicating Haselton Spreadsheet)

    My = ComputeMyPanagiotakosFardis2001(b,d,dp,fy,fpc,rho,rhop,0.0,N, Imperial)
    MyNeg = ComputeMyPanagiotakosFardis2001(b,d,dp,fy,fpc,rhop,rho,0.0,N, Imperial)*-1

    Hardening = 1.25*0.89**v*0.91**(0.01*c_units*fpc)-1.0 #Checked with Haselton Spreadsheet
    # Hardening = 0.13

    if Imperial:
        Es = 29000.
        Ec = 57000.*(fpc*1000.)**.5/1000
    else:
        Es = 207.e3
        Ec = 4700.*(fpc)**.5
    n = Es / Ec

    if Imperial:
        EIg = 57000.*(fpc*1000.)**.5/1000.*Ig
    else:
        EIg = 4750.*(fpc)**.5*Ig

    # alpha_eff = max(min(0.6, (-0.07+0.59*v+0.07*LsOH)), 0.2) #Use when element will see yielding
    alpha_eff = max(0.35, min((-0.02+0.98*v+0.09*LsOH), 0.8)) #Use when element will not likely see yielding
    # If Beam Set at the Minimum (This is what Haselton Does in his Spreadsheet)
    if rhop != rho:
        alpha_eff = 0.35
    Ky = alpha_eff*EIg

    # theta_y_pos = My/(6.0*(Ky*StiffnessFactor1)/L)
    # theta_y_neg = -1.0*MyNeg/(6.0*(Ky*StiffnessFactor1)/L)

    theta_p_pos = theta_p_pos#-theta_y_pos
    theta_p_neg = theta_p_neg#-theta_y_neg

    c=1

	# set a_mem [expr ($n+1.0)*($Mycol_12*($McMy-1.0)) / ($Ks_col_1*$th_pP)];	# strain hardening ratio of spring
	# set b [expr ($a_mem)/(1.0+$n*(1.0-$a_mem))];

    elstk = StiffnessFactor1*6.*Ky/L
    AlphaY_pos = ((1+Hardening)*My-My)/theta_p_pos/elstk
    AlphaY_pos = AlphaY_pos*(-1*StiffnessFactor2*AlphaY_pos)/(AlphaY_pos*(AlphaY_pos-StiffnessFactor2))

    AlphaY_neg = -1*((1+Hardening)*MyNeg-MyNeg)/theta_p_neg/elstk
    AlphaY_neg = AlphaY_neg*(-1*StiffnessFactor2*AlphaY_neg)/(AlphaY_neg*(AlphaY_neg-StiffnessFactor2))

    Mat = OpenSeesAPI.Model.Element.Material.UniaxialMaterial.ModIMKPeakOriented(id, elstk, AlphaY_pos, AlphaY_neg, My, MyNeg, Lambda, Lambda, Lambda, Lambda, c,c,c,c,theta_p_pos,theta_p_neg,theta_pc,theta_pc,0.01,0.01,(theta_p+theta_pc)*2,(theta_p+theta_pc)*2,1,1,_Notes=_Notes)

    # AlphaCap = -1*(Hardening+1.0)*My/theta_pc/elstk
    # AlphaCap = AlphaCap*-1*StiffnessFactor2*AlphaCap/(AlphaCap*(AlphaCap-StiffnessFactor2))
    #
    # Mat = OpenSeesAPI.Model.Element.Material.UniaxialMaterial.Clough(id, elstk, My, MyNeg, AlphaY_pos, 0.01, AlphaCap, theta_p_pos, -1*theta_p_neg, Lambda, 0,0,Lambda,c,c,c,c,_Notes =_Notes)

    return Mat

def ComputeMyPanagiotakosFardis2001(b,d,dp,fy,fpc,rho,rhop,rho_v,N, Imperial=True):
    if Imperial:
        Es=29000.
        Ec=57000.*(fpc*1000.)**.5/1000
    else:
        Es=199.5*1e3
        Ec=4750.*(fpc)**.5
    n=Es/Ec

    #dp #the distance of the center of the compression reinforcement from the extreme compression fibers

    deltap = float(dp)/d

    ## Determine Whether Tension or Compression Controlled

    # Tension Steel Controlled
    A = rho + rhop + rho_v+N/b/d/fy
    B = rho + rhop*deltap + 0.5*rho_v*(1.0+deltap)+N/b/d/fy

    ky_t = (n**2.*A**2.+2*n*B)**.5-n*A
    Phi_y_t = fy/(Es*(1.0-ky_t)*d)

    # Compression Zone Controlled
    A = rho + rhop + rho_v - N/1.8/n/b/d/fpc
    B = rho + rhop*deltap + 0.5*rho_v*(1.0+deltap)

    ky_c = (n**2*A**2+2*n*B)**.5-n*A
    Phi_y_c = 1.8*fpc/Ec/ky_c/d

    if Phi_y_t<Phi_y_c:
        Phi_y = Phi_y_t
        ky = ky_t
    else:
        Phi_y = Phi_y_c
        ky = ky_c

    My = b*d**3.*Phi_y*(Ec*ky**2./2.*(0.5*(1.+deltap)-ky/3.)+Es/2.*((1-ky)*rho+(ky-deltap)*rhop+rho_v/6.*(1.-deltap))*(1.-deltap))

    return My

def ComputeIg(h,b,d,dp,rho,rho_p,n):
    TopSteelArea = h*b*rho_p*(n-1.)
    BotSteelArea = h*b*rho*(n-1.)
    ConcArea = float(h)*b
    SumAyTop = ConcArea*h/2.+TopSteelArea*dp+BotSteelArea*d
    y_top = SumAyTop/ConcArea
    I = float(b)*h**3./12.+TopSteelArea*(dp-y_top)**2.+BotSteelArea*(d-y_top)**2.+ConcArea*(h/2.-y_top)**2
    return I

def ComputeIgWithSlab(h,b,d,dp,rho,rho_p,n,t_slab,rho_slab,effective_width):
    TopSteelArea = h*b*rho_p*(n-1.)
    BotSteelArea = h*b*rho*(n-1.)
    ConcArea = float(h)*b
    SlabArea = effective_width*t_slab-b*t_slab
    SlabReinf = SlabArea*rho_slab*(n-1.)
    SumAyTop = ConcArea*h/2.+TopSteelArea*dp+BotSteelArea*d+SlabArea*t_slab/2.+SlabReinf*t_slab/2.
    y_top = SumAyTop/(ConcArea+SlabArea)
    I = float(b)*h**3./12.+(effective_width-b)*float(t_slab)**3/12.+TopSteelArea*(dp-y_top)**2.+BotSteelArea*(d-y_top)**2.+ConcArea*(h/2.-y_top)**2.+SlabArea*(t_slab/2.-y_top)**2.+SlabReinf*(t_slab/2.-y_top)**2
    return I

def ComputeJointStiffness(fpc, ColumnWidth, ColumnDepth, BeamHeight):
    # Stiffness = 22.0*(fpc*1000)**.5*ColumnWidth*ColumnDepth*BeamHeight/0.004/1000
    Stiffness = 0.11*57000.0*(fpc*1000)**.5*ColumnWidth*ColumnDepth*BeamHeight/1000 # This has been checked with the Haselton Spreadsheet
    return Stiffness

def GetCloughModel(id,b,h,dp,fy,fpc,rho,rhop,rho_v,N,s,LsOH,L,Ig,a_sl=1,Imperial=True,t_slab=None, rho_slab=None, eff_wdith=None,_Notes=None):
    import OpenSeesAPI

    StiffnessFactor1 = 11
    StiffnessFactor2 = 1.1

    #a_sl=0 if bond slip of long. bars is prevented

    #LsOH #Column Height to Depth
    d = h - dp

    v = float(N)/(fpc*b*h)

    s_n = (float(s)/d)*(fy/100.)**0.5

    if Imperial:
        c_units = 6.9 # For Imperial
    else:
        c_units = 1 # For Metric

    rho_sh = rho_v
    # theta_p = 0.12*(1.+0.55*a_sl)*(0.16)**v*(0.02+40*rho_sh)**0.43*0.54**(0.01*c_units*fpc)*0.66**(0.1*s_n)*2.27**(10.0*(rho+rhop))
    theta_p = 0.13*(1.+0.55*a_sl)*(0.13)**v*(0.02+40*rho_sh)**0.65*0.57**(0.01*c_units*fpc)

    if rhop != rho:
        theta_p_pos = (max(0.01,rho*fy/fpc)/max(0.1,rhop*fy/fpc))**0.175*theta_p
        theta_p_neg = (max(0.01,rhop*fy/fpc)/max(0.1,rho*fy/fpc))**0.175*theta_p
    else:
        theta_p_pos = theta_p
        theta_p_neg = theta_p

    # From Spreadsheet
    theta_pc = min(0.76*0.031**v*(0.02+40.*rho_sh)*1.02,0.10)

    # From Spreadsheet
    Lambda = 170.7*0.27**v*(0.1)**(s/d)*StiffnessFactor1

    My = ComputeMyPanagiotakosFardis2001(b,d,dp,fy,fpc,rho,rhop,0.0,N, Imperial)
    MyNeg = ComputeMyPanagiotakosFardis2001(b,d,dp,fy,fpc,rhop,rho,0.0,N, Imperial)*-1

    # From Spreadsheet
    Hardening = 1.25*0.89**v*0.91**(0.01*c_units*fpc)-1.0
    # Hardening = 0.13

    if Imperial:
        Es=30000.
        Ec=57000.*(fpc*1000.)**.5/1000
    else:
        Es=207.e3
        Ec=4700.*(fpc)**.5
    n=Es/Ec

    # Ig = ComputeIg(h,b,d,dp,rho,rhop,n)
    if Imperial:
        EIg = 57000.*(fpc*1000.)**.5/1000.*Ig
    else:
        EIg = 4700.*(fpc)**.5*Ig

    alpha_eff = GetEIeffOEIg(v, LsOH)
    Ky = alpha_eff*EIg*StiffnessFactor1

    theta_y_pos = My/(6.0*(Ky)/L)
    theta_y_neg = -1.0*MyNeg/(6.0*(Ky)/L)

    theta_p_pos = theta_p_pos#-theta_y_pos
    theta_p_neg = theta_p_neg#-theta_y_neg

    c=1

    AlphaY_pos = StiffnessFactor1*(Hardening*My)/theta_p_pos/(6.*(Ky/StiffnessFactor2)/L)
    AlphaY_pos = AlphaY_pos/(1.0+StiffnessFactor1/StiffnessFactor2*(1.0-AlphaY_pos))

    AlphaY_neg = StiffnessFactor1*(Hardening*-1.*MyNeg)/theta_p_neg/(6.*(Ky/StiffnessFactor2)/L)
    AlphaY_neg = AlphaY_neg/(1.0+StiffnessFactor1/StiffnessFactor2*(1.0-AlphaY_neg))

    AlphaCap = -1*(Hardening+1.0)*My/theta_pc/(6.0*Ky/L)
    AlphaCap = AlphaCap*-1*StiffnessFactor2*AlphaCap/(AlphaCap*(AlphaCap-StiffnessFactor2))

    Mat = OpenSeesAPI.Model.Element.Material.UniaxialMaterial.Clough(id, 6.*Ky/L, My, MyNeg, AlphaY_pos, 0.01, AlphaCap, theta_p_pos, -1*theta_p_neg, Lambda, 0,0,Lambda,c,c,c,c,_Notes =_Notes)

    return Mat

def GetEIeffOEIg(v,LsOH):
    # alpha_eff = max(min(0.6,(-0.07+0.59*v+0.07*LsOH)),0.2) #Use when element will see yielding
    alpha_eff = max(0.35,min((-0.02+0.98*v+0.09*LsOH),0.8)) #Use when element will not likely see yielding
    return alpha_eff

def GetTrussMaterialModel(id, As, E, Fy, Ry, w, Beta):
    Ko = E

    s1p = Ry*Fy
    e1p = s1p/Ko
    s2p = w * Ry * Fy
    e2p = 10*e1p
    s3p = s2p + (e1p*10)*0.0125*Ko
    e3p = e2p * 2

    s1n = -1*Ry*Fy
    e1n = -1*e1p
    s2n = -1*Beta*w*Ry*Fy
    e2n = -1*e2p
    s3n = -1*(s2p + 0.0125*Ko*(e1p*10))
    e3n = -1*e3p

    pinchx = 1.0
    pinchy = 1.0
    damage1 = 0.001
    damage2 = 0.001
    degradation = 0.001

    import OpenSeesAPI
    Material = OpenSeesAPI.Model.Element.Material.UniaxialMaterial.Hysteretic(id, s1p, e1p,s2p, e2p, s3p, e3p, s1n,e1n, s2n, e2n, s3n, e3n, pinchx, pinchy, damage1, damage2, degradation)

    return Material

def GetConfinedConcreteMaterial(id, fpc):
    #fpc is in ksi
    ec = 0.004
    fpcu = fpc - 0.15*fpc
    ecu = 0.014

    import OpenSeesAPI
    return OpenSeesAPI.Model.Element.Material.UniaxialMaterial.Concrete01(id, -1*fpc,  -1*ec, -1*fpcu,  -1*ecu)

def GetUnconfinedConcreteMaterial(id, fpc):
    #fpc is in ksi
    ec = 0.002
    fpcu = 0
    ecu = 0.006

    import OpenSeesAPI
    return OpenSeesAPI.Model.Element.Material.UniaxialMaterial.Concrete01(id, -1*fpc,  -1*ec, -1*fpcu,  -1*ecu)

def GetSteelMaterial(id, fy):
    #fpc is in ksi
    Ec = 29000
    b = 0.01

    import OpenSeesAPI
    return OpenSeesAPI.Model.Element.Material.UniaxialMaterial.Steel02(id, fy, Ec, b)

def GetSteelMaterialWithMinMax(ODatabase, fy, fu, epsilon_comp, b=0.01):
    #fpc is in ksi
    Es = 29000

    import OpenSeesAPI
    Steel = ODatabase.AddMaterial(OpenSeesAPI.Material.UniaxialMaterial.Steel02(ODatabase.GetFreeMaterialId(1,1), fy, Es, b))
    minStrain = epsilon_comp
    maxStrain = fy/float(Es)+(fu-fy)/float(Es)/float(b)
    MinMax = ODatabase.AddMaterial(OpenSeesAPI.Material.UniaxialMaterial.MinMax(ODatabase.GetFreeMaterialId(1,1), Steel, minStrain, maxStrain))

    return MinMax

def GetPughConfinedConcreteMaterial(OSDatabase, fpc, b_x, b_y, s_bars_x, s_bars_y, s_ties_x, s_ties_y, f_y, A_tie, no_tie_x, no_tie_y, L_ip,
                                    GfccOGfc = 1.7, GfcOfpc = 2.0, fy_hoop=70.2):

    import OpenSeesAPI
    import math
    # Saatcioglu and Razvi 1992
    ksi_to_mpa = 6.8947
    f_l_x = A_tie*fy_hoop*no_tie_x/s_ties_x/b_x*ksi_to_mpa    #Uniform Confining Pressure
    k2_x = min(0.26*(b_x/s_ties_x*b_x/s_bars_x*1./f_l_x)**.5,1.)
    f_1e_x = k2_x*f_l_x

    f_l_y = A_tie*fy_hoop*no_tie_y/s_ties_y/b_y*ksi_to_mpa    #Uniform Confining Pressure
    k2_y = min(0.26*(b_y/s_ties_y*b_y/s_bars_y*1./f_l_y)**.5,1.)
    f_1e_y = k2_y*f_l_y

    f_1e = (f_1e_x*b_x + f_1e_y*b_y)/(b_x+b_y)
    k_1 = 6.7*f_1e**(-0.17)

    fpc_confined = fpc + (k_1*f_1e)/ksi_to_mpa

    Kc = fpc_confined/fpc

    e01 = 2. * fpc / (57. * (fpc * 1000.) ** .5)  # 0.002

    ec = e01 * (1. + 5. * k_1 * f_1e / fpc / ksi_to_mpa) # This is the original SAAT and Razvi Equation
    # ec = 2. * Kc * fpc / 57. / (Kc * fpc * 1000) ** .5  # Pugh Dissertation Equation 2.4

    f20 = 0.2*Kc*fpc

    A_s = A_tie*math.ceil(b_x/s_ties_x+1.) + A_tie*math.ceil(b_y/s_ties_y+1.)

    s_bar_avg = (s_bars_x+s_bars_y)/2.

    rho = A_s/s_bar_avg / (b_x + b_y)

    e085 = 0.0038

    e01 = 2.*fpc/(57.*(fpc*1000.)**.5)#0.002

    e1 = e01*(1.+5.*Kc)

    e20 = 260.*rho*e1 + e085        # Equation 14/15 in Saat and Razvi

    f_t = 4.0*(fpc*1000)**.5/1000

    Ec = 57.*(fpc*1000)**.5
    Et = Ec
    Ets = 0.05 * Et

    lam = 0.1                       # Not Sure where this is mentioned in Pugh Thesis

    Gfc = GfcOfpc*fpc*0.0393701
    Gfcc = GfccOGfc*Gfc

    fpeak = Kc*fpc

    e20u_reg = Gfcc/0.6/fpeak/L_ip - 0.8 * fpeak / Ec + ec

    Mat =  OpenSeesAPI.Model.Element.Material.UniaxialMaterial.Concrete02(OSDatabase.GetFreeMaterialId(1,1), -1*fpc_confined,  -1*ec, -1*f20,  -1*e20u_reg, lam, f_t, Ets, _Notes='Confined Conc.')

    OSDatabase.AddMaterial(Mat)

    return Mat

def GetPughConfinedConcreteCMMaterial(OSDatabase, fpc, b_x, b_y, s_bars_x, s_bars_y, s_ties_x, s_ties_y, f_y, A_tie, no_tie_x, no_tie_y, L_ip):

    import OpenSeesAPI
    import math
    # Saatcioglu and Razvi 1992
    ksi_to_mpa = 6.8947
    f_l_x = A_tie*f_y*no_tie_x/s_ties_x/b_x*ksi_to_mpa    #Uniform Confining Pressure
    k2_x = min(0.26*(b_x/s_ties_x*b_x/s_bars_x*1./f_l_x)**.5,1.)
    f_1e_x = k2_x*f_l_x

    f_l_y = A_tie*f_y*no_tie_y/s_ties_y/b_y*ksi_to_mpa    #Uniform Confining Pressure
    k2_y = min(0.26*(b_y/s_ties_y*b_y/s_bars_y*1./f_l_y)**.5,1.)
    f_1e_y = k2_y*f_l_y

    f_1e = (f_1e_x*b_x + f_1e_y*b_y)/(b_x+b_y)
    k_1 = 6.7*f_1e**(-0.17)

    fpc_confined = fpc + (k_1*f_1e)/ksi_to_mpa

    Kc = fpc_confined/fpc

    ec = 2.*Kc*fpc/57./(Kc*fpc*1000)**.5        # Pugh Dissertation Equation 2.4

    f20 = 0.2*Kc*fpc

    A_s = A_tie*math.ceil(b_x/s_ties_x+1.) + A_tie*math.ceil(b_y/s_ties_y+1.)
    rho = A_s/max(s_bars_x,s_bars_y) / (b_x + b_y)

    e085 = 0.0038

    e01 = 0.002

    e1 = e01*(1.+5.*Kc)

    e20 = 260.*rho*e1 + e085        # Equation 14/15 in Saat and Razvi

    f_t = 4.0*(fpc*1000)**.5/1000

    Ec = 57.*(fpc*1000)**.5
    Et = Ec
    Ets = 0.05 * Et

    lam = 0.1                       # Not Sure where this is mentioned in Pugh Thesis

    Gfc = 2.*fpc*0.0393701
    Gfcc = 1.7*Gfc

    fpeak = Kc*fpc

    e20u_reg = Gfcc/0.6/fpeak/L_ip - 0.8 * fpeak / Ec + ec

    rc = 7 # Shape Factor, See Tsai's equation
    rt = 1.2 # Shape Factor, See Tsai's equation
    et = f_t/Et
    xcrn = 1.035 # Non0dimensional critical strain
    xcrp = 10000

    Mat =  OpenSeesAPI.Material.UniaxialMaterial.ConcreteCM(OSDatabase.GetFreeMaterialId(1,1),
                                                            -1*fpc_confined,  -1*ec, Ec, rc, xcrn, f_t, et, rt, xcrp)
    OSDatabase.AddMaterial(Mat)
    return Mat

def GetPughConfinedConcreteMaterialUnregularized(OSDatabase, fpc, b_x, b_y, s_bars_x, s_bars_y, s_ties_x, s_ties_y, f_y, A_tie, no_tie_x, no_tie_y, L_ip):

    import OpenSeesAPI
    import math
    # Saatcioglu and Razvi 1992
    ksi_to_mpa = 6.8947
    f_l_x = A_tie*f_y*no_tie_x/s_ties_x/b_x*ksi_to_mpa    #Uniform Confining Pressure
    k2_x = min(0.26*(b_x/s_ties_x*b_x/s_bars_x*1./f_l_x)**.5,1.)
    f_1e_x = k2_x*f_l_x

    f_l_y = A_tie*f_y*no_tie_y/s_ties_y/b_y*ksi_to_mpa    #Uniform Confining Pressure
    k2_y = min(0.26*(b_y/s_ties_y*b_y/s_bars_y*1./f_l_y)**.5,1.)
    f_1e_y = k2_y*f_l_y

    f_1e = (f_1e_x*b_x + f_1e_y*b_y)/(b_x+b_y)
    k_1 = 6.7*f_1e**(-0.17)

    fpc_confined = fpc + (k_1*f_1e)/ksi_to_mpa

    Kc = 2.3#fpc_confined/fpc

    ec = 2.*Kc*fpc/57./(Kc*fpc*1000)**.5        # Pugh Dissertation Equation 2.4

    f20 = 0.2*Kc*fpc

    A_s = A_tie*math.ceil(b_x/s_ties_x+1.) + A_tie*math.ceil(b_y/s_ties_y+1.)
    rho = A_s/max(s_bars_x,s_bars_y) / (b_x + b_y)

    e085 = 0.0038

    e01 = 0.002

    e1 = e01*(1.+5.*Kc)

    e20 = 260.*rho*e1 + e085        # Equation 14/15 in Saat and Razvi

    f_t = 4.0*(fpc*1000)**.5/1000

    Ec = 57.*(fpc*1000)**.5
    Et = Ec
    Ets = 0.05 * Et

    lam = 0.1                       # Not Sure where this is mentioned in Pugh Thesis

    Gfc = 2.*fpc*0.0393701

    Kc = 2.3 # Pugh Recommends 1.7, Laura says to use something higher, 2.3 maybe

    Gfcc = Kc*Gfc

    fpeak = Kc*fpc

    #Using Saat and Razvi prediction
    e20u_reg = e20#Gfcc/0.6/fpeak/L_ip + 0.8 * fpeak / Ec + ec

    Mat =  OpenSeesAPI.Model.Element.Material.UniaxialMaterial.Concrete02(OSDatabase.GetFreeMaterialId(1,1), -1*fpc_confined,  -1*ec, -1*f20,  -1*e20u_reg, lam, f_t, Ets)

    OSDatabase.AddMaterial(Mat)
    return Mat

def GetPughUnconfinedConcreteMaterial(OSDatabase, fpc, L_ip, fpuOfpc=0.2, GfcOfpc = 2.0):

    import OpenSeesAPI

    eu = 2.*fpc/57./(fpc*1000)**.5              # Pugh Dissertation Equation 2.3

    f20 = fpuOfpc*fpc

    e20 = 0.008

    f_t = 4.0*(fpc*1000)**.5/1000.

    Ec = 57.*(fpc*1000)**.5
    Et = Ec
    Ets = 0.05 * Et

    lam = 0.1                       # Not Sure where this is mentioned in Pugh Thesis

    Gfc = GfcOfpc*fpc*0.0393701

    fpeak = fpc

    e20u_reg = Gfc / 0.6 / fpeak / L_ip - 0.8 * fpeak / Ec + eu

    Mat = OpenSeesAPI.Model.Element.Material.UniaxialMaterial.Concrete02(OSDatabase.GetFreeMaterialId(1,1), -1*fpc,  -1*eu, -1*f20,  -1*e20u_reg, lam, f_t, Ets, _Notes='Unconfined Conc.')

    # Mat = OpenSeesAPI.Model.Element.Material.UniaxialMaterial.ConcreteCM(OSDatabase.GetFreeMaterialId(1, 1),
    #                                                                      -1 * fpc, -1 * eu, Ec,
    #                                                                      7, e20u_reg, f_t, Ets, 1.2, 1000,
    #                                                                      _Notes='Confined Conc.', _ecu = e20u_reg)

    OSDatabase.AddMaterial(Mat)
    return Mat

def GetPughSteel02Material(OSDatabase, fy, fu, L_ip, eps_ult_comp, eps_ult_exp = .12, L_gage = 8, Regularize=True):
    import OpenSeesAPI
    #L_gage = 8 # As per Pugh et al 2015 Paper, Use this when the gage length is not available

    Es = 29000
    eps_y = fy/Es
    if Regularize:
        eps_ult = eps_y + L_gage/L_ip * (eps_ult_exp - eps_y)
        b_reg = (fu-fy)/(Es*(eps_ult_exp - eps_y))
    else:
        eps_ult = eps_ult_exp
        b_reg = (fu-fy)/(Es*(eps_ult_exp - eps_y))

    Mat1 = OpenSeesAPI.Model.Element.Material.UniaxialMaterial.Steel02(OSDatabase.GetFreeMaterialId(1,1), fy, Es, b_reg, 20.0, 0.925, 0.15)
    Mat2 = OpenSeesAPI.Model.Element.Material.UniaxialMaterial.MinMax(OSDatabase.GetFreeMaterialId(1,1), Mat1, eps_ult_comp, eps_ult)

    OSDatabase.AddMaterial(Mat1)
    OSDatabase.AddMaterial(Mat2)

    return Mat2

def GetPughSteel02MaterialUnregularized(OSDatabase, fy, fu, L_ip, eps_ult_comp, eps_ult_exp = .12, L_gage = 8, ):
    import OpenSeesAPI
    #L_gage = 8 # As per Pugh et al 2015 Paper, Use this when the gage length is not available

    Es = 29000.
    eps_y = fy/Es
    eps_ult = eps_ult_exp

    b_reg = (fu-fy)/(Es*(eps_ult - eps_y))

    Mat1 = OpenSeesAPI.Model.Element.Material.UniaxialMaterial.Steel02(OSDatabase.GetFreeMaterialId(1,1), fy, Es, b_reg, 20.0, 0.925, 0.15)
    Mat2 = OpenSeesAPI.Model.Element.Material.UniaxialMaterial.MinMax(OSDatabase.GetFreeMaterialId(1,1), Mat1, eps_ult_comp, eps_ult)

    OSDatabase.AddMaterial(Mat1)
    OSDatabase.AddMaterial(Mat2)

    return Mat2

def ComputePlanarShearWallFiberSection(OSDatabase, l_w, t_w, l_boundary, cover, bound_reinf_layout = [3,3,3,3,3,3,3], bound_reinf_size = 4 , web_reinf_spacing = 6, web_reinf_size = 2, fpc = 6, fy = 60, fu = 90, tie_spacing = 2, tie_size = 2., no_ties_x = 7, no_ties_y = 2, max_mesh_Size=3, NoOfIntPoints = 5, h = 48):

    import OpenSeesAPI
    import numpy as np

    if NoOfIntPoints == 3:
        # L_ip = [0.6666667, 0.1666667, 0.6666667]
        L_ip = [0.1666667, 0.6666667, 0.6666667]
    elif NoOfIntPoints == 5:
        # L_ip = [0.35555556, 0.27222222, 0.05000000, 0.27222222, 0.35555556]
        L_ip = [0.05000000, 0.27222222, 0.35555556, 0.27222222, 0.05000000]
    else:
        # L_ip = [0.243810, 0.215873, 0.138413, 0.023810, 0.138413, 0.215873, 0.243810]
        L_ip = [0.023810, 0.138413, 0.215873, 0.243810, 0.215873,  0.138413, 0.023810]

    # Axial
    # ConcConfined = GetConfinedConcreteMaterial(OSDatabase.GetFreeMaterialId(1,1) , fpc)
    # ConcUncofined = GetUnconfinedConcreteMaterial(OSDatabase.GetFreeMaterialId(1,1) , fpc)
    # Steel = GetSteelMaterial(OSDatabase.GetFreeMaterialId(1,1) , fy)
    # SteelMinMax = OpenSeesAPI.Model.Element.Material.UniaxialMaterial.MinMax(OSDatabase.GetFreeMaterialId(1,1) , Steel, -0.0038, 0.16)

    # Axial according to Pugh
    s_bar_x = (l_boundary-2.*cover)/len(bound_reinf_layout)
    s_bar_y = (t_w-2.*cover)/np.mean(bound_reinf_layout)

    A_tie = np.pi*(tie_size/8./2.)**2.

    ConcConfinedAllIP = []
    UnconfinedConcAllIP = []
    SteelConfinedAllIP = []
    SteelUnconfinedAllIP = []
    for i in range(len(L_ip)):
        ConcConfinedAllIP.append(GetPughConfinedConcreteMaterial(OSDatabase, fpc, l_boundary, t_w, s_bar_x, s_bar_y, tie_spacing, tie_spacing, fy, A_tie, no_ties_x, no_ties_y, L_ip[i]*h))
        UnconfinedConcAllIP.append(GetPughUnconfinedConcreteMaterial(OSDatabase, fpc, L_ip[i]*h))
        ecu_Confined = ConcConfinedAllIP[-1]._ecu
        ecu_Unconfined = UnconfinedConcAllIP[-1]._ecu
        SteelConfinedAllIP.append(GetPughSteel02Material(OSDatabase, fy, fu, L_ip[i]*h, ecu_Confined))
        SteelUnconfinedAllIP.append(GetPughSteel02Material(OSDatabase, fy, fu, L_ip[i]*h, ecu_Unconfined))

    # Shear
    # Define Elastic Spring for Shear
    Ec = 57. * (fpc * 1000)**.5
    Mu = 0.15
    Gc = Ec / (2. * (1. + Mu))
    Geff = 0.1 * Gc  # assume 10% of elastic shear modulus per Oyen (2006)
    ks = 5. / 6.  # reduction factor for shear area per Pugh (2012)
    Eg = Geff * ks * l_w * t_w

    ShearElastic = OpenSeesAPI.Model.Element.Material.UniaxialMaterial.Elastic(OSDatabase.GetFreeMaterialId(1, 1), Eg)

    OSDatabase.AddMaterial(ShearElastic)

    Sections = []

    for j in range(len(L_ip)):
        ConcConfined = ConcConfinedAllIP[j]
        ConcUncofined = UnconfinedConcAllIP[j]
        SteelConfined = SteelConfinedAllIP[j]
        SteelUncofined = SteelUnconfinedAllIP[j]

        Fibers = []

        #Confined
        Fibers.append(OpenSeesAPI.Model.Element.Material.Section.FiberSection.Patch.Rect(ConcConfined,
                               int((l_boundary - 2 * cover)/max_mesh_Size + 1),
                               int((t_w-2*cover)/max_mesh_Size + 1),
                               cover, cover, l_boundary-cover, t_w-cover
                               ))

        XDistTo2ndBoundary = l_w-l_boundary
        Fibers.append(OpenSeesAPI.Model.Element.Material.Section.FiberSection.Patch.Rect(ConcConfined,
                               int((l_boundary - 2 * cover) / max_mesh_Size + 1),
                               int((t_w - 2 * cover) / max_mesh_Size + 1),
                               cover + XDistTo2ndBoundary, cover, XDistTo2ndBoundary + l_boundary - cover, t_w - cover
                               ))

        #Unconfined

        #Left and Right Covers
        Fibers.append(OpenSeesAPI.Model.Element.Material.Section.FiberSection.Patch.Rect(ConcUncofined,
                               int(cover / max_mesh_Size + 1),
                               int(t_w / max_mesh_Size + 1),
                               0, 0, cover, t_w
                               ))

        Fibers.append(OpenSeesAPI.Model.Element.Material.Section.FiberSection.Patch.Rect(ConcUncofined,
                               int(cover / max_mesh_Size + 1),
                               int(t_w  / max_mesh_Size + 1),
                               l_boundary - cover, 0, l_boundary, t_w
                               ))

        Fibers.append(OpenSeesAPI.Model.Element.Material.Section.FiberSection.Patch.Rect(ConcUncofined,
                                 int(cover / max_mesh_Size + 1),
                                 int(t_w / max_mesh_Size + 1),
                                 XDistTo2ndBoundary,
                                 0,
                                 XDistTo2ndBoundary+cover,
                                 t_w
                                 ))

        Fibers.append(OpenSeesAPI.Model.Element.Material.Section.FiberSection.Patch.Rect(ConcUncofined,
                                 int(cover / max_mesh_Size + 1),
                                 int(t_w / max_mesh_Size + 1),
                                 XDistTo2ndBoundary + l_boundary - cover,
                                 0,
                                 XDistTo2ndBoundary + l_boundary,
                                 t_w
                                 ))

        # Top and Bottom
        Fibers.append(OpenSeesAPI.Model.Element.Material.Section.FiberSection.Patch.Rect(ConcUncofined,
                                 int((l_boundary - 2 * cover) / max_mesh_Size + 1),
                                 int(cover / max_mesh_Size + 1),
                                 cover, 0, l_boundary - cover, cover
                                 ))

        Fibers.append(OpenSeesAPI.Model.Element.Material.Section.FiberSection.Patch.Rect(ConcUncofined,
                                 int((l_boundary - 2 * cover) / max_mesh_Size + 1),
                                 int(cover / max_mesh_Size + 1),
                                 cover, t_w - cover, l_boundary - cover, t_w
                                 ))

        Fibers.append(OpenSeesAPI.Model.Element.Material.Section.FiberSection.Patch.Rect(ConcUncofined,
                                 int((l_boundary - 2 * cover) / max_mesh_Size + 1),
                                 int(cover / max_mesh_Size + 1),
                                 XDistTo2ndBoundary + cover, 0, XDistTo2ndBoundary + l_boundary - cover, cover
                                 ))

        Fibers.append(OpenSeesAPI.Model.Element.Material.Section.FiberSection.Patch.Rect(ConcUncofined,
                                 int((l_boundary - 2 * cover) / max_mesh_Size + 1),
                                 int(cover / max_mesh_Size + 1),
                                 XDistTo2ndBoundary + cover, t_w - cover,
                                 XDistTo2ndBoundary + l_boundary - cover,
                                 t_w
                                 ))

        #Web
        Fibers.append(OpenSeesAPI.Model.Element.Material.Section.FiberSection.Patch.Rect(ConcUncofined,
                                 int((l_w - 2. * l_boundary) / max_mesh_Size + 1),
                                 int(t_w / max_mesh_Size + 1),
                                 l_boundary,
                                 0,
                                 XDistTo2ndBoundary,
                                 t_w
                                 ))

        # Boundary Reinf.
        import numpy as np
        spacing = np.linspace(0,l_boundary-2.*cover, len(bound_reinf_layout))

        for i in range(len(bound_reinf_layout)):
            Fibers.append(OpenSeesAPI.Model.Element.Material.Section.FiberSection.Layer.Straight(SteelConfined,
                                         bound_reinf_layout[i],
                                         np.pi*(bound_reinf_size/8./2.)**2.,
                                         cover + spacing[i],
                                         cover,
                                         cover + spacing[i],
                                         t_w - cover
                                         ))

            Fibers.append(OpenSeesAPI.Model.Element.Material.Section.FiberSection.Layer.Straight(SteelConfined,
                                         bound_reinf_layout[i],
                                         np.pi * (
                                         bound_reinf_size / 8. / 2.) ** 2.,
                                         XDistTo2ndBoundary + cover  + spacing[i],
                                         cover,
                                         XDistTo2ndBoundary + cover  + spacing[i],
                                         t_w - cover
                                         ))

        no_web_bars = int((l_w - 2.*l_boundary)/web_reinf_spacing)+1-2
        Fibers.append(OpenSeesAPI.Model.Element.Material.Section.FiberSection.Layer.Straight(SteelUncofined,
                                     no_web_bars,
                                     np.pi * (
                                         web_reinf_size / 8. / 2.) ** 2.,
                                     l_boundary + web_reinf_spacing,
                                     cover,
                                     XDistTo2ndBoundary - web_reinf_spacing,
                                     cover
                                     ))

        Fibers.append(OpenSeesAPI.Model.Element.Material.Section.FiberSection.Layer.Straight(SteelUncofined,
                                     no_web_bars,
                                     np.pi * (
                                         web_reinf_size / 8. / 2.) ** 2.,
                                     l_boundary + web_reinf_spacing,
                                     t_w - cover,
                                     XDistTo2ndBoundary - web_reinf_spacing,
                                     t_w - cover
                                     ))

        FiberSection = OpenSeesAPI.Model.Element.Material.Section.FiberSection(OSDatabase.GetFreeMaterialId(1,1), Fibers)
        OSDatabase.AddMaterial(FiberSection)

        #Add Shear
        Section = OpenSeesAPI.Model.Element.Material.Section.Aggregator(OSDatabase.GetFreeMaterialId(1,1), [ShearElastic], ['Vy'], FiberSection)
        OSDatabase.AddMaterial(Section)

        Sections.append(Section)

    return Sections

def ComputePlanarShearWallFiberSectionDispBased(OSDatabase, l_w, t_w, l_boundary, cover,
                                                bound_reinf_layout = [3,3,3,3,3,3,3],
                                                bound_reinf_size = 4 , web_reinf_spacing = 6, web_reinf_size = 2,
                                                fpc = 6, fy = 60, fu = 90, tie_spacing = 2, tie_size = 2.,
                                                no_ties_x = 7, no_ties_y = 2,
                                                max_mesh_Size=3, NoOfIntPoints = 5, h = 48, GfcOfpc=0.56,
                                                GfccOGfc=2.0, fy_hoop=70.2, eps_ult_be=0.12,
                                                web_bars_no = None, web_tie_no = None,
                                                web_tie_spacing = None,
                                                web_bar_near_boundary=None):

    import OpenSeesAPI
    import numpy as np

    L_ip = np.ones(NoOfIntPoints)

    # Axial according to Pugh
    s_bar_x = (l_boundary-2.*cover-tie_spacing/8.)/(len(bound_reinf_layout)-1)
    s_bar_y = (t_w-2.*cover-tie_spacing/8.)/np.mean(bound_reinf_layout)

    A_tie = np.pi*(tie_size/8./2.)**2.

    ConcConfinedAllIP = []
    UnconfinedConcAllIP = []
    WebConfinedConcAllIP = []
    SteelConfinedAllIP = []
    SteelUnconfinedAllIP = []
    SteelWebConfinedAllIP = []
    for i in range(len(L_ip)):
        UnconfinedConcAllIP.append(GetPughUnconfinedConcreteMaterial(OSDatabase, fpc, L_ip[i]*h, GfcOfpc=GfcOfpc))
        if web_tie_no != None:
            b_x = l_w - l_boundary*2.
            b_y = t_w - 2. * cover - tie_size / 8.
            WebConfinedConcAllIP.append(
                GetPughConfinedConcreteMaterial(OSDatabase, fpc, b_x, b_y,
                                                b_x / (web_tie_no + 1.), b_y,
                                                web_tie_spacing, web_tie_spacing, fy, A_tie, web_tie_no,
                                                2., L_ip[i] * h, GfcOfpc=GfcOfpc, GfccOGfc=GfccOGfc,
                                                fy_hoop=fy_hoop))
        else:
            WebConfinedConcAllIP.append(UnconfinedConcAllIP[-1])


        if no_ties_x == 0:
            ConcConfinedAllIP.append(UnconfinedConcAllIP[-1])
        else:
            ConcConfinedAllIP.append(GetPughConfinedConcreteMaterial(OSDatabase, fpc, l_boundary - 2.*cover - tie_size/8., t_w - 2.*cover - tie_size/8., s_bar_x, s_bar_y,
                                                                 tie_spacing, tie_spacing, fy, A_tie, no_ties_x,
                                                                 no_ties_y, L_ip[i]*h, GfcOfpc=GfcOfpc, GfccOGfc=GfccOGfc, fy_hoop=fy_hoop))
        ecu_Confined = ConcConfinedAllIP[-1]._ecu
        ecu_Unconfined = UnconfinedConcAllIP[-1]._ecu
        SteelConfinedAllIP.append(GetPughSteel02Material(OSDatabase, fy, fu,
                                                         L_ip[i]*h,
                                                         ecu_Confined,
                                                         eps_ult_exp=eps_ult_be,
                                                         Regularize=False))
        SteelUnconfinedAllIP.append(GetPughSteel02Material(OSDatabase, fy, fu,
                                                           L_ip[i]*h,
                                                           ecu_Unconfined,
                                                           eps_ult_exp=eps_ult_be,
                                                           Regularize=False))
        if web_tie_no != None:
            SteelWebConfinedAllIP.append(GetPughSteel02Material(OSDatabase, fy, fu,
                                                               L_ip[i] * h,
                                                               WebConfinedConcAllIP[-1]._ecu,
                                                               eps_ult_exp=eps_ult_be,
                                                               Regularize=False))
        else:
            SteelWebConfinedAllIP.append(SteelUnconfinedAllIP[-1])

    # Shear
    # Define Elastic Spring for Shear
    Ec = 57. * (fpc * 1000)**.5
    Mu = 0.15
    Gc = Ec / (2. * (1. + Mu))
    Geff = 0.1 * Gc  # assume 10% of elastic shear modulus per Oyen (2006)
    ks = 5. / 6.  # reduction factor for shear area per Pugh (2012)
    Eg = Geff * ks * l_w * t_w

    ShearElastic = OpenSeesAPI.Model.Element.Material.UniaxialMaterial.Elastic(OSDatabase.GetFreeMaterialId(1, 1), Eg)

    OSDatabase.AddMaterial(ShearElastic)

    Sections = []

    for j in range(len(L_ip)):
        if web_tie_no != None:
            ConcConfined = ConcConfinedAllIP[j]
            ConcUncofined = WebConfinedConcAllIP[j]
            SteelConfined = SteelConfinedAllIP[j]
            SteelUncofined = SteelWebConfinedAllIP[j]
        else:
            ConcConfined = ConcConfinedAllIP[j]
            ConcUncofined = UnconfinedConcAllIP[j]
            SteelConfined = SteelConfinedAllIP[j]
            SteelUncofined = SteelUnconfinedAllIP[j]

        Fibers = []

        #Confined
        Fibers.append(OpenSeesAPI.Model.Element.Material.Section.FiberSection.Patch.Rect(ConcConfined,
                               int((l_boundary - 2 * cover)/max_mesh_Size + 1),
                               int((t_w-2*cover)/max_mesh_Size + 1),
                               cover, cover, l_boundary-cover, t_w-cover
                               ))

        XDistTo2ndBoundary = l_w-l_boundary
        Fibers.append(OpenSeesAPI.Model.Element.Material.Section.FiberSection.Patch.Rect(ConcConfined,
                               int((l_boundary - 2 * cover) / max_mesh_Size + 1),
                               int((t_w - 2 * cover) / max_mesh_Size + 1),
                               cover + XDistTo2ndBoundary, cover, XDistTo2ndBoundary + l_boundary - cover, t_w - cover
                               ))

        #Unconfined

        #Left and Right Covers
        Fibers.append(OpenSeesAPI.Model.Element.Material.Section.FiberSection.Patch.Rect(ConcUncofined,
                               int(cover / max_mesh_Size + 1),
                               int(t_w / max_mesh_Size + 1),
                               0, 0, cover, t_w
                               ))

        Fibers.append(OpenSeesAPI.Model.Element.Material.Section.FiberSection.Patch.Rect(ConcUncofined,
                               int(cover / max_mesh_Size + 1),
                               int(t_w  / max_mesh_Size + 1),
                               l_boundary - cover, 0, l_boundary, t_w
                               ))

        Fibers.append(OpenSeesAPI.Model.Element.Material.Section.FiberSection.Patch.Rect(ConcUncofined,
                                 int(cover / max_mesh_Size + 1),
                                 int(t_w / max_mesh_Size + 1),
                                 XDistTo2ndBoundary,
                                 0,
                                 XDistTo2ndBoundary+cover,
                                 t_w
                                 ))

        Fibers.append(OpenSeesAPI.Model.Element.Material.Section.FiberSection.Patch.Rect(ConcUncofined,
                                 int(cover / max_mesh_Size + 1),
                                 int(t_w / max_mesh_Size + 1),
                                 XDistTo2ndBoundary + l_boundary - cover,
                                 0,
                                 XDistTo2ndBoundary + l_boundary,
                                 t_w
                                 ))

        # Top and Bottom
        Fibers.append(OpenSeesAPI.Model.Element.Material.Section.FiberSection.Patch.Rect(ConcUncofined,
                                 int((l_boundary - 2 * cover) / max_mesh_Size + 1),
                                 int(cover / max_mesh_Size + 1),
                                 cover, 0, l_boundary - cover, cover
                                 ))

        Fibers.append(OpenSeesAPI.Model.Element.Material.Section.FiberSection.Patch.Rect(ConcUncofined,
                                 int((l_boundary - 2 * cover) / max_mesh_Size + 1),
                                 int(cover / max_mesh_Size + 1),
                                 cover, t_w - cover, l_boundary - cover, t_w
                                 ))

        Fibers.append(OpenSeesAPI.Model.Element.Material.Section.FiberSection.Patch.Rect(ConcUncofined,
                                 int((l_boundary - 2 * cover) / max_mesh_Size + 1),
                                 int(cover / max_mesh_Size + 1),
                                 XDistTo2ndBoundary + cover, 0, XDistTo2ndBoundary + l_boundary - cover, cover
                                 ))

        Fibers.append(OpenSeesAPI.Model.Element.Material.Section.FiberSection.Patch.Rect(ConcUncofined,
                                 int((l_boundary - 2 * cover) / max_mesh_Size + 1),
                                 int(cover / max_mesh_Size + 1),
                                 XDistTo2ndBoundary + cover, t_w - cover,
                                 XDistTo2ndBoundary + l_boundary - cover,
                                 t_w
                                 ))

        #Web
        Fibers.append(OpenSeesAPI.Model.Element.Material.Section.FiberSection.Patch.Rect(ConcUncofined,
                                 int((l_w - 2. * l_boundary) / max_mesh_Size + 1),
                                 int(t_w / max_mesh_Size + 1),
                                 l_boundary,
                                 0,
                                 XDistTo2ndBoundary,
                                 t_w
                                 ))

        # Boundary Reinf.
        import numpy as np
        spacing = np.linspace(0,l_boundary-2.*cover, len(bound_reinf_layout))

        for i in range(len(bound_reinf_layout)):
            Fibers.append(OpenSeesAPI.Model.Element.Material.Section.FiberSection.Layer.Straight(SteelConfined,
                                         bound_reinf_layout[i],
                                         np.pi*(bound_reinf_size/8./2.)**2.,
                                         cover + spacing[i],
                                         cover,
                                         cover + spacing[i],
                                         t_w - cover
                                         ))

            Fibers.append(OpenSeesAPI.Model.Element.Material.Section.FiberSection.Layer.Straight(SteelConfined,
                                         bound_reinf_layout[i],
                                         np.pi * (
                                         bound_reinf_size / 8. / 2.) ** 2.,
                                         XDistTo2ndBoundary + cover  + spacing[i],
                                         cover,
                                         XDistTo2ndBoundary + cover  + spacing[i],
                                         t_w - cover
                                         ))
        if web_bar_near_boundary == None:
            web_bar_near_boundary = False

        if web_bars_no == None:
            no_web_bars = np.round((l_w - 2.*l_boundary)/web_reinf_spacing)-1
        else:
            no_web_bars = web_bars_no
            if web_bar_near_boundary:
                web_reinf_spacing = (l_w - 2.*l_boundary + 2* cover) / (no_web_bars - 1.)
            else:
                web_reinf_spacing = (l_w - 2. * l_boundary + 2 * cover) / (no_web_bars + 1.)
        if web_bar_near_boundary:
            Fibers.append(OpenSeesAPI.Model.Element.Material.Section.FiberSection.Layer.Straight(SteelUncofined,
                                         no_web_bars,
                                         np.pi * (
                                             web_reinf_size / 8. / 2.) ** 2.,
                                         l_boundary,
                                         cover,
                                         XDistTo2ndBoundary,
                                         cover
                                         ))

            Fibers.append(OpenSeesAPI.Model.Element.Material.Section.FiberSection.Layer.Straight(SteelUncofined,
                                         no_web_bars,
                                         np.pi * (
                                             web_reinf_size / 8. / 2.) ** 2.,
                                         l_boundary,
                                         t_w - cover,
                                         XDistTo2ndBoundary,
                                         t_w - cover
                                         ))
        else:
            Fibers.append(OpenSeesAPI.Model.Element.Material.Section.FiberSection.Layer.Straight(SteelUncofined,
                                                                                                 no_web_bars,
                                                                                                 np.pi * (
                                                                                                     web_reinf_size / 8. / 2.) ** 2.,
                                                                                                 l_boundary + web_reinf_spacing,
                                                                                                 cover,
                                                                                                 XDistTo2ndBoundary - web_reinf_spacing,
                                                                                                 cover
                                                                                                 ))

            Fibers.append(OpenSeesAPI.Model.Element.Material.Section.FiberSection.Layer.Straight(SteelUncofined,
                                                                                                 no_web_bars,
                                                                                                 np.pi * (
                                                                                                     web_reinf_size / 8. / 2.) ** 2.,
                                                                                                 l_boundary + web_reinf_spacing,
                                                                                                 t_w - cover,
                                                                                                 XDistTo2ndBoundary - web_reinf_spacing,
                                                                                                 t_w - cover
                                                                                                 ))

        FiberSection = OpenSeesAPI.Model.Element.Material.Section.FiberSection(OSDatabase.GetFreeMaterialId(1,1), Fibers)
        OSDatabase.AddMaterial(FiberSection)

        #Add Shear
        Section = OpenSeesAPI.Model.Element.Material.Section.Aggregator(OSDatabase.GetFreeMaterialId(1,1), [ShearElastic], ['Vy'], FiberSection)
        OSDatabase.AddMaterial(Section)

        Sections.append(Section)

    return Sections


def ComputePlanarShearWallFiberSection3D(OSDatabase, l_w, t_w, l_boundary, cover, bound_reinf_layout = [3,3,3,3,3,3,3], bound_reinf_size = 4 , web_reinf_spacing = 6, web_reinf_size = 2, fpc = 6, fy = 60, fu = 90, tie_spacing = 2, tie_size = 2., no_ties_x = 7, no_ties_y = 2, max_mesh_Size=3, NoOfIntPoints = 5, h = 48):

    import OpenSeesAPI
    import numpy as np

    if NoOfIntPoints == 3:
        # L_ip = [0.6666667, 0.1666667, 0.6666667]
        L_ip = [0.1666667, 0.6666667, 0.6666667]
    elif NoOfIntPoints == 5:
        # L_ip = [0.35555556, 0.27222222, 0.05000000, 0.27222222, 0.35555556]
        L_ip = [0.05000000, 0.27222222, 0.35555556, 0.27222222, 0.05000000]
    else:
        # L_ip = [0.243810, 0.215873, 0.138413, 0.023810, 0.138413, 0.215873, 0.243810]
        L_ip = [0.023810, 0.138413, 0.215873, 0.243810, 0.215873,  0.138413, 0.023810]

    # Axial
    # ConcConfined = GetConfinedConcreteMaterial(OSDatabase.GetFreeMaterialId(1,1) , fpc)
    # ConcUncofined = GetUnconfinedConcreteMaterial(OSDatabase.GetFreeMaterialId(1,1) , fpc)
    # Steel = GetSteelMaterial(OSDatabase.GetFreeMaterialId(1,1) , fy)
    # SteelMinMax = OpenSeesAPI.Model.Element.Material.UniaxialMaterial.MinMax(OSDatabase.GetFreeMaterialId(1,1) , Steel, -0.0038, 0.16)

    # Axial according to Pugh
    s_bar_x = (l_boundary-2.*cover)/len(bound_reinf_layout)
    s_bar_y = (t_w-2.*cover)/np.mean(bound_reinf_layout)

    A_tie = np.pi*(tie_size/8./2.)**2.

    ConcConfinedAllIP = []
    UnconfinedConcAllIP = []
    SteelConfinedAllIP = []
    SteelUnconfinedAllIP = []
    for i in range(len(L_ip)):
        ConcConfinedAllIP.append(GetPughConfinedConcreteMaterial(OSDatabase, fpc, l_boundary, t_w, s_bar_x, s_bar_y, tie_spacing, tie_spacing, fy, A_tie, no_ties_x, no_ties_y, L_ip[i]*h))
        UnconfinedConcAllIP.append(GetPughUnconfinedConcreteMaterial(OSDatabase, fpc, L_ip[i]*h))
        ecu_Confined = ConcConfinedAllIP[-1]._ecu
        ecu_Unconfined = UnconfinedConcAllIP[-1]._ecu
        SteelConfinedAllIP.append(GetPughSteel02Material(OSDatabase, fy, fu, L_ip[i]*h, ecu_Confined))
        SteelUnconfinedAllIP.append(GetPughSteel02Material(OSDatabase, fy, fu, L_ip[i]*h, ecu_Unconfined))

    # Shear
    # Define Elastic Spring for Shear
    Ec = 57. * (fpc * 1000)**.5
    Mu = 0.15
    Gc = Ec / (2. * (1. + Mu))
    Geff = Gc#0.1 * Gc  # assume 10% of elastic shear modulus per Oyen (2006)
    ks = 5. / 6.  # reduction factor for shear area per Pugh (2012)
    Eg = Geff * ks * l_w * t_w

    a = max(l_w, t_w)
    b = min(l_w, t_w)
    J = a * b ** 3 * (1. / 3. - 0.21 * b / float(a) * (1. - b ** 4 / 12. / (float(a) ** 4)))

    ShearElastic = OSDatabase.AddMaterial(OpenSeesAPI.Model.Element.Material.UniaxialMaterial.Elastic(OSDatabase.GetFreeMaterialId(1, 1), Eg))

    TorsionalElastic = OSDatabase.AddMaterial(OpenSeesAPI.Model.Element.Material.UniaxialMaterial.Elastic(OSDatabase.GetFreeMaterialId(1, 1), Geff*J))


    Sections = []

    for j in range(len(L_ip)):
        ConcConfined = ConcConfinedAllIP[j]
        ConcUncofined = UnconfinedConcAllIP[j]
        SteelConfined = SteelConfinedAllIP[j]
        SteelUncofined = SteelUnconfinedAllIP[j]

        Fibers = []

        #Confined
        Fibers.append(OpenSeesAPI.Model.Element.Material.Section.FiberSection.Patch.Rect(ConcConfined,
                               int((l_boundary - 2 * cover)/max_mesh_Size + 1),
                               int((t_w-2*cover)/max_mesh_Size + 1),
                               cover, cover, l_boundary-cover, t_w-cover
                               ))

        XDistTo2ndBoundary = l_w-l_boundary
        Fibers.append(OpenSeesAPI.Model.Element.Material.Section.FiberSection.Patch.Rect(ConcConfined,
                               int((l_boundary - 2 * cover) / max_mesh_Size + 1),
                               int((t_w - 2 * cover) / max_mesh_Size + 1),
                               cover + XDistTo2ndBoundary, cover, XDistTo2ndBoundary + l_boundary - cover, t_w - cover
                               ))

        #Unconfined

        #Left and Right Covers
        Fibers.append(OpenSeesAPI.Model.Element.Material.Section.FiberSection.Patch.Rect(ConcUncofined,
                               int(cover / max_mesh_Size + 1),
                               int(t_w / max_mesh_Size + 1),
                               0, 0, cover, t_w
                               ))

        Fibers.append(OpenSeesAPI.Model.Element.Material.Section.FiberSection.Patch.Rect(ConcUncofined,
                               int(cover / max_mesh_Size + 1),
                               int(t_w  / max_mesh_Size + 1),
                               l_boundary - cover, 0, l_boundary, t_w
                               ))

        Fibers.append(OpenSeesAPI.Model.Element.Material.Section.FiberSection.Patch.Rect(ConcUncofined,
                                 int(cover / max_mesh_Size + 1),
                                 int(t_w / max_mesh_Size + 1),
                                 XDistTo2ndBoundary,
                                 0,
                                 XDistTo2ndBoundary+cover,
                                 t_w
                                 ))

        Fibers.append(OpenSeesAPI.Model.Element.Material.Section.FiberSection.Patch.Rect(ConcUncofined,
                                 int(cover / max_mesh_Size + 1),
                                 int(t_w / max_mesh_Size + 1),
                                 XDistTo2ndBoundary + l_boundary - cover,
                                 0,
                                 XDistTo2ndBoundary + l_boundary,
                                 t_w
                                 ))

        # Top and Bottom
        Fibers.append(OpenSeesAPI.Model.Element.Material.Section.FiberSection.Patch.Rect(ConcUncofined,
                                 int((l_boundary - 2 * cover) / max_mesh_Size + 1),
                                 int(cover / max_mesh_Size + 1),
                                 cover, 0, l_boundary - cover, cover
                                 ))

        Fibers.append(OpenSeesAPI.Model.Element.Material.Section.FiberSection.Patch.Rect(ConcUncofined,
                                 int((l_boundary - 2 * cover) / max_mesh_Size + 1),
                                 int(cover / max_mesh_Size + 1),
                                 cover, t_w - cover, l_boundary - cover, t_w
                                 ))

        Fibers.append(OpenSeesAPI.Model.Element.Material.Section.FiberSection.Patch.Rect(ConcUncofined,
                                 int((l_boundary - 2 * cover) / max_mesh_Size + 1),
                                 int(cover / max_mesh_Size + 1),
                                 XDistTo2ndBoundary + cover, 0, XDistTo2ndBoundary + l_boundary - cover, cover
                                 ))

        Fibers.append(OpenSeesAPI.Model.Element.Material.Section.FiberSection.Patch.Rect(ConcUncofined,
                                 int((l_boundary - 2 * cover) / max_mesh_Size + 1),
                                 int(cover / max_mesh_Size + 1),
                                 XDistTo2ndBoundary + cover, t_w - cover,
                                 XDistTo2ndBoundary + l_boundary - cover,
                                 t_w
                                 ))

        #Web
        Fibers.append(OpenSeesAPI.Model.Element.Material.Section.FiberSection.Patch.Rect(ConcUncofined,
                                 int((l_w - 2. * l_boundary) / max_mesh_Size + 1),
                                 int(t_w / max_mesh_Size + 1),
                                 l_boundary,
                                 0,
                                 XDistTo2ndBoundary,
                                 t_w
                                 ))

        # Boundary Reinf.
        import numpy as np
        spacing = np.linspace(0,l_boundary-2.*cover, len(bound_reinf_layout))

        for i in range(len(bound_reinf_layout)):
            Fibers.append(OpenSeesAPI.Model.Element.Material.Section.FiberSection.Layer.Straight(SteelConfined,
                                         bound_reinf_layout[i],
                                         np.pi*(bound_reinf_size/8./2.)**2.,
                                         cover + spacing[i],
                                         cover,
                                         cover + spacing[i],
                                         t_w - cover
                                         ))

            Fibers.append(OpenSeesAPI.Model.Element.Material.Section.FiberSection.Layer.Straight(SteelConfined,
                                         bound_reinf_layout[i],
                                         np.pi * (
                                         bound_reinf_size / 8. / 2.) ** 2.,
                                         XDistTo2ndBoundary + cover  + spacing[i],
                                         cover,
                                         XDistTo2ndBoundary + cover  + spacing[i],
                                         t_w - cover
                                         ))

        no_web_bars = int((l_w - 2.*l_boundary)/web_reinf_spacing)+1-2
        Fibers.append(OpenSeesAPI.Model.Element.Material.Section.FiberSection.Layer.Straight(SteelUncofined,
                                     no_web_bars,
                                     np.pi * (
                                         web_reinf_size / 8. / 2.) ** 2.,
                                     l_boundary + web_reinf_spacing,
                                     cover,
                                     XDistTo2ndBoundary - web_reinf_spacing,
                                     cover
                                     ))

        Fibers.append(OpenSeesAPI.Model.Element.Material.Section.FiberSection.Layer.Straight(SteelUncofined,
                                     no_web_bars,
                                     np.pi * (
                                         web_reinf_size / 8. / 2.) ** 2.,
                                     l_boundary + web_reinf_spacing,
                                     t_w - cover,
                                     XDistTo2ndBoundary - web_reinf_spacing,
                                     t_w - cover
                                     ))

        FiberSection = OpenSeesAPI.Model.Element.Material.Section.FiberSection(OSDatabase.GetFreeMaterialId(1,1), Fibers, GJ=Geff*J)
        OSDatabase.AddMaterial(FiberSection)

        #Add Shear
        Section = OpenSeesAPI.Model.Element.Material.Section.Aggregator(OSDatabase.GetFreeMaterialId(1,1), [ShearElastic, ShearElastic, TorsionalElastic], ['Vy', 'Vz', 'T'], FiberSection)
        OSDatabase.AddMaterial(Section)

        Sections.append(Section)

    return Sections

def ComputeOnlyConcWall3D(OSDatabase, l_w, t_w, l_boundary, cover, bound_reinf_layout = [3,3,3,3,3,3,3], bound_reinf_size = 4 , web_reinf_spacing = 6, web_reinf_size = 2, fpc = 6, fy = 60, fu = 90, tie_spacing = 2, tie_size = 2., no_ties_x = 7, no_ties_y = 2, max_mesh_Size=3, NoOfIntPoints = 5, h = 48):

    import OpenSeesAPI
    import numpy as np

    if NoOfIntPoints == 3:
        # L_ip = [0.6666667, 0.1666667, 0.6666667]
        L_ip = [0.1666667, 0.6666667, 0.6666667]
    elif NoOfIntPoints == 5:
        # L_ip = [0.35555556, 0.27222222, 0.05000000, 0.27222222, 0.35555556]
        L_ip = [0.05000000, 0.27222222, 0.35555556, 0.27222222, 0.05000000]
    else:
        # L_ip = [0.243810, 0.215873, 0.138413, 0.023810, 0.138413, 0.215873, 0.243810]
        L_ip = [0.023810, 0.138413, 0.215873, 0.243810, 0.215873,  0.138413, 0.023810]

    # Axial
    # ConcConfined = GetConfinedConcreteMaterial(OSDatabase.GetFreeMaterialId(1,1) , fpc)
    # ConcUncofined = GetUnconfinedConcreteMaterial(OSDatabase.GetFreeMaterialId(1,1) , fpc)
    # Steel = GetSteelMaterial(OSDatabase.GetFreeMaterialId(1,1) , fy)
    # SteelMinMax = OpenSeesAPI.Model.Element.Material.UniaxialMaterial.MinMax(OSDatabase.GetFreeMaterialId(1,1) , Steel, -0.0038, 0.16)

    # Axial according to Pugh
    s_bar_x = (l_boundary-2.*cover)/len(bound_reinf_layout)
    s_bar_y = (t_w-2.*cover)/np.mean(bound_reinf_layout)

    A_tie = np.pi*(tie_size/8./2.)**2.

    ConcConfinedAllIP = []
    UnconfinedConcAllIP = []
    SteelConfinedAllIP = []
    SteelUnconfinedAllIP = []
    for i in range(len(L_ip)):
        ConcConfinedAllIP.append(GetPughConfinedConcreteMaterial(OSDatabase, fpc, l_boundary, t_w, s_bar_x, s_bar_y, tie_spacing, tie_spacing, fy, A_tie, no_ties_x, no_ties_y, L_ip[i]*h))
        UnconfinedConcAllIP.append(GetPughUnconfinedConcreteMaterial(OSDatabase, fpc, L_ip[i]*h))
        ecu_Confined = ConcConfinedAllIP[-1]._ecu
        ecu_Unconfined = UnconfinedConcAllIP[-1]._ecu
        SteelConfinedAllIP.append(GetPughSteel02Material(OSDatabase, fy, fu, L_ip[i]*h, ecu_Confined))
        SteelUnconfinedAllIP.append(GetPughSteel02Material(OSDatabase, fy, fu, L_ip[i]*h, ecu_Unconfined))

    # Shear
    # Define Elastic Spring for Shear
    Ec = 57. * (fpc * 1000)**.5
    Mu = 0.15
    Gc = Ec / (2. * (1. + Mu))
    Geff = Gc#0.1 * Gc  # assume 10% of elastic shear modulus per Oyen (2006)
    ks = 5. / 6.  # reduction factor for shear area per Pugh (2012)
    Eg = Geff * ks * l_w * t_w

    a = max(l_w, t_w)
    b = min(l_w, t_w)
    J = a * b ** 3 * (1. / 3. - 0.21 * b / float(a) * (1. - b ** 4 / 12. / (float(a) ** 4)))

    ShearElastic = OSDatabase.AddMaterial(OpenSeesAPI.Model.Element.Material.UniaxialMaterial.Elastic(OSDatabase.GetFreeMaterialId(1, 1), Eg))

    TorsionalElastic = OSDatabase.AddMaterial(OpenSeesAPI.Model.Element.Material.UniaxialMaterial.Elastic(OSDatabase.GetFreeMaterialId(1, 1), Geff*J))

    Sections = []

    for j in range(len(L_ip)):
        ConcConfined = ConcConfinedAllIP[j]
        ConcUncofined = UnconfinedConcAllIP[j]
        SteelConfined = SteelConfinedAllIP[j]
        SteelUncofined = SteelUnconfinedAllIP[j]

        Fibers = []

        #Confined
        Fibers.append(OpenSeesAPI.Model.Element.Material.Section.FiberSection.Patch.Rect(ConcConfined,
                               int((l_boundary - 2 * cover)/max_mesh_Size + 1),
                               int((t_w-2*cover)/max_mesh_Size + 1),
                               cover, cover, l_boundary-cover, t_w-cover
                               ))

        XDistTo2ndBoundary = l_w-l_boundary
        Fibers.append(OpenSeesAPI.Model.Element.Material.Section.FiberSection.Patch.Rect(ConcConfined,
                               int((l_boundary - 2 * cover) / max_mesh_Size + 1),
                               int((t_w - 2 * cover) / max_mesh_Size + 1),
                               cover + XDistTo2ndBoundary, cover, XDistTo2ndBoundary + l_boundary - cover, t_w - cover
                               ))

        #Unconfined

        #Left and Right Covers
        Fibers.append(OpenSeesAPI.Model.Element.Material.Section.FiberSection.Patch.Rect(ConcUncofined,
                               int(cover / max_mesh_Size + 1),
                               int(t_w / max_mesh_Size + 1),
                               0, 0, cover, t_w
                               ))

        Fibers.append(OpenSeesAPI.Model.Element.Material.Section.FiberSection.Patch.Rect(ConcUncofined,
                               int(cover / max_mesh_Size + 1),
                               int(t_w  / max_mesh_Size + 1),
                               l_boundary - cover, 0, l_boundary, t_w
                               ))

        Fibers.append(OpenSeesAPI.Model.Element.Material.Section.FiberSection.Patch.Rect(ConcUncofined,
                                 int(cover / max_mesh_Size + 1),
                                 int(t_w / max_mesh_Size + 1),
                                 XDistTo2ndBoundary,
                                 0,
                                 XDistTo2ndBoundary+cover,
                                 t_w
                                 ))

        Fibers.append(OpenSeesAPI.Model.Element.Material.Section.FiberSection.Patch.Rect(ConcUncofined,
                                 int(cover / max_mesh_Size + 1),
                                 int(t_w / max_mesh_Size + 1),
                                 XDistTo2ndBoundary + l_boundary - cover,
                                 0,
                                 XDistTo2ndBoundary + l_boundary,
                                 t_w
                                 ))

        # Top and Bottom
        Fibers.append(OpenSeesAPI.Model.Element.Material.Section.FiberSection.Patch.Rect(ConcUncofined,
                                 int((l_boundary - 2 * cover) / max_mesh_Size + 1),
                                 int(cover / max_mesh_Size + 1),
                                 cover, 0, l_boundary - cover, cover
                                 ))

        Fibers.append(OpenSeesAPI.Model.Element.Material.Section.FiberSection.Patch.Rect(ConcUncofined,
                                 int((l_boundary - 2 * cover) / max_mesh_Size + 1),
                                 int(cover / max_mesh_Size + 1),
                                 cover, t_w - cover, l_boundary - cover, t_w
                                 ))

        Fibers.append(OpenSeesAPI.Model.Element.Material.Section.FiberSection.Patch.Rect(ConcUncofined,
                                 int((l_boundary - 2 * cover) / max_mesh_Size + 1),
                                 int(cover / max_mesh_Size + 1),
                                 XDistTo2ndBoundary + cover, 0, XDistTo2ndBoundary + l_boundary - cover, cover
                                 ))

        Fibers.append(OpenSeesAPI.Model.Element.Material.Section.FiberSection.Patch.Rect(ConcUncofined,
                                 int((l_boundary - 2 * cover) / max_mesh_Size + 1),
                                 int(cover / max_mesh_Size + 1),
                                 XDistTo2ndBoundary + cover, t_w - cover,
                                 XDistTo2ndBoundary + l_boundary - cover,
                                 t_w
                                 ))

        #Web
        Fibers.append(OpenSeesAPI.Model.Element.Material.Section.FiberSection.Patch.Rect(ConcUncofined,
                                 int((l_w - 2. * l_boundary) / max_mesh_Size + 1),
                                 int(t_w / max_mesh_Size + 1),
                                 l_boundary,
                                 0,
                                 XDistTo2ndBoundary,
                                 t_w
                                 ))

        FiberSection = OpenSeesAPI.Model.Element.Material.Section.FiberSection(OSDatabase.GetFreeMaterialId(1,1), Fibers, GJ=Geff*J)
        OSDatabase.AddMaterial(FiberSection)

        #Add Shear
        Section = OpenSeesAPI.Model.Element.Material.Section.Aggregator(OSDatabase.GetFreeMaterialId(1,1), [ShearElastic, ShearElastic, TorsionalElastic], ['Vy', 'Vz', 'T'], FiberSection)
        OSDatabase.AddMaterial(Section)

        Sections.append(Section)

    return Sections

def GetManderConfinedConcreteMaterial01(id, fpc, ecu, eccu, fyh, cover, LBE, tBE, BENumberVerticalBarLayers, AbBE, BENumberVerticalBars, AbhBE, BETransverseBarSpacing):
    import math

    s = BETransverseBarSpacing # spacing of spiral
    dc = LBE-cover
    bc = tBE-2.0*cover
    wi = (LBE-cover)/BENumberVerticalBarLayers
    sp = s-float('10')
    rhocc = (AbBE*BENumberVerticalBars)/(bc*dc)
    Asx = 2*AbhBE
    Asy = BENumberVerticalBarLayers*AbhBE
    rhox = Asx/(s*dc)
    rhoy = Asy/(s*bc)
    ke = ((1-((BENumberVerticalBarLayers-1)*(wi**2/(6*bc*dc))))*(1-sp/(2*bc))*(1-sp/(2*dc)))/(1-rhocc)
    flpx = ke*rhox*fyh
    flpy = ke*rhoy*fyh
    flp = (flpx+flpy)/2
    fpcc = fpc*(-1.254+2.254*math.sqrt(1+(7.94*flp)/fpc)-(2*flp)/fpc)
    ecc = ecu*(1+5*(fpcc/fpc-1))
    fpccu  = 0.5*fpcc
    import OpenSeesAPI
    return OpenSeesAPI.Model.Element.Material.UniaxialMaterial.Concrete01(id, -1*fpcc,  -1*ecc, -1*fpccu,  -1*eccu)

def GetManderConfinedConcreteMaterial04(id, fpc, ecu, Ec, fct, et, beta, fyh, cover, LBE, tBE, BENumberVerticalBarLayers, AbBE, BENumberVerticalBars, AbhBE, BETransverseBarSpacing):
    import math

    """
    fpc = concrete strength
    ecu = compression strain at max strength
    Ecu = modulus of elasticty of concrete
    fct = tensile strength
    et = tensile strain capacity
    beta = constant to define logarithmic decay in tension
    cover = concrete cover to BE reinf.
    LBE = length of BE
    tBE = thickness of BE
    BENumberVerticalBayLayers = Number of layers of vertical reinf. in boundary element
    AbBE = area of single vertical bar in BE
    BENumberVerticalBars = Number of vertical bars in the BE
    AbhBE = Area of transverse reinforcing
    BETransverseBarSpacing = Spacing of transverse reinforcing
    """

    s = BETransverseBarSpacing # spacing of spiral
    dc = LBE-cover
    bc = tBE-2.0*cover
    wi = (LBE-cover)/BENumberVerticalBarLayers
    sp = s-float('10')
    rhocc = (AbBE*BENumberVerticalBars)/(bc*dc)
    Asx = 2*AbhBE
    Asy = BENumberVerticalBarLayers*AbhBE
    rhox = Asx/(s*dc)
    rhoy = Asy/(s*bc)
    ke = ((1-((BENumberVerticalBarLayers-1)*(wi**2/(6*bc*dc))))*(1-sp/(2*bc))*(1-sp/(2*dc)))/(1-rhocc)
    flpx = ke*rhox*fyh
    flpy = ke*rhoy*fyh
    flp = (flpx+flpy)/2
    fpcc = fpc*(-1.254+2.254*math.sqrt(1+(7.94*flp)/fpc)-(2*flp)/fpc)
    ecc = ecu*(1+5*(fpcc/fpc-1))
    eccu = 2.0*ecc
    fpccu  = 0.5*fpcc
    import OpenSeesAPI
    return OpenSeesAPI.Model.Element.Material.UniaxialMaterial.Concrete04(id, -1*fpcc,  -1*ecc, -1*eccu, Ec, fct, et, beta)

def GetReinforcingSteelMaterial(id, fy, fu, Es, barSize, BETransverseBarSpacing):
    import math

    # Material Curve Stuff
    Esh = 0.05 * Es
    esh = 0.005
    eult = 0.50

    # Buckling Stuff
    Lsr = BETransverseBarSpacing/float(barSize)
    beta = 1.0
    r = 0.5
    gamma = 0.9

    Cf = 0.5  # ductility constant used to admust the number of cycles until failure. A higher value will result in lower damage for each cycle
    alpha = 0.45  # used to calibrate damage from one strain range to an equaveltn damage at aother strain range
    Cd = 0.27  # strength reduction constant. A lager value will result in a lower strength reduction for each cycle

    # alpha = 0.506
    # Cf = 0.26
    # Cd = 0.389



    import OpenSeesAPI
    return OpenSeesAPI.Model.Element.Material.UniaxialMaterial.ReinforcingSteel.GABuck(id,fy,fu,Es,Esh,esh,eult,Lsr,beta,r,gamma,Cf,alpha,Cd)

def MultiLayeredShellElement(OData, fc, ft, fcu, eps0, epsu, esptu, stc, ShearModulus, fy, fu, strainhardening):

    import OpenSeesAPI
    OpenSeesAPI.Database.Collector.AddMaterial()
    ConcreteMaterial = OData.AddMaterial(OpenSeesAPI.Material.NDMaterial.PlaneStressUserMaterial(OData.GetFreeMaterialId(1,1), fc, ft, fcu, eps0, epsu, esptu, stc))

    PlateShearStress = OData.AddMaterial(OpenSeesAPI.Material.NDMaterial.PlateFromPlaneStress(OData.GetFreeMaterialId(1,1), ConcreteMaterial, ShearModulus))

    E = 29000
    ConfinedSteel = OData.AddMaterial(OpenSeesAPI.Material.UniaxialMaterial.Steel02(OData.GetFreeMaterialId(1,1), fy, E, strainhardening))

    UnconfinedSteel = OData.AddMaterial(OpenSeesAPI.Material.UniaxialMaterial.Steel02(OData.GetFreeMaterialId(1,1), fy, E, strainhardening))

    LongitudinalSteelConfined = OData.AddMaterial(OpenSeesAPI.Material.NDMaterial.PlateRebar(OData.GetFreeMaterialId(1,1), ConfinedSteel, 90))
    LongitudinalSteelUnconfined = OData.AddMaterial(OpenSeesAPI.Material.NDMaterial.PlateRebar(OData.GetFreeMaterialId(1,1), UnconfinedSteel, 90))

    TransverseShear = OData.AddMaterial(OpenSeesAPI.Material.NDMaterial.PlateRebar(OData.GetFreeMaterialId(1,1), UnconfinedSteel, 0))

def SaatRazviConfinement(fc, fy, Abe, sBe, t, sHoop, rhoBars):
    import numpy as np

    Abe = Abe * 25.4 ^ 2
    sBe = sBe * 25.4
    t = t * 25.4
    cover = 2.5 * 25.4
    sHoop = sHoop * 25.4

    eps0 = 2 * fc / (57000 * np.sqrt(fc))
    eps085 = np.interp([0.2 * fc, fc], [0.008, eps0], 0.85 * fc)

    class data():
        pass

    xDir = data()
    yDir = data()

    xDir.barsA = np.array([Abe, Abe]) * 0.5
    xDir.barsFy = np.array([fy, fy]) / 1000 * 6.895
    xDir.b = sBe
    xDir.s1 = sBe
    xDir.sHoop = sHoop
    yDir.barsA = Abe * np.ones(rhoBars)
    yDir.barsFy = fy / 1000 * 6.895 * np.ones(rhoBars)
    yDir.b = t - 2 * cover
    yDir.s1 = (t - 2 * cover) / (rhoBars - 1)
    yDir.sHoop = sHoop

    # strength
    flx = sum(xDir.barsA*xDir.barsFy) / (xDir.b * xDir.sHoop)
    k2x = min(0.26 * np.sqrt(xDir.b ** 2 / xDir.s1 / xDir.sHoop / flx), 1)
    flex = k2x * flx
    fly = sum(yDir.barsA*yDir.barsFy) / (yDir.b * yDir.sHoop)
    k2y = min(0.26 * np.sqrt(yDir.b ** 2 / yDir.s1 / yDir.sHoop / fly), 1)
    fley = k2y * fly
    fle = (flex * xDir.b + fley * yDir.b) / (xDir.b + yDir.b)
    k1 = 6.7 * fle ** -0.17
    fcc = fc / 1000 * 6.895 + k1 * fle
    kc = (fcc / 6.895 * 1000) / fc

    # ductility
    K = k1 * fle / (fc / 1000 * 6.895)
    eps1 = eps0 * (1 + 5 * K)
    sAve = (xDir.sHoop + yDir.sHoop) / 2
    rho = (sum(xDir.barsA) + sum(yDir.barsA)) / sAve / (xDir.b + yDir.b)

    eps85 = eps085 + 260 * rho * eps1
    eps20 = np.interp([kc * fc, 0.85 * kc * fc], [eps1, eps85], 0.2 * kc * fc)
    if eps20 <= 0:
        print('strain determination went south')

    eps1 = 2 * kc * fc / (57000 * np.sqrt(fc))

    return eps1, kc, eps20
