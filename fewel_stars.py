import numpy as np
import math
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import matplotlib as mpl

#import cPickle as pickle

def read(tstep):
#for tstep in range(0,150):
    cell_no = 100
    cell_size = 25
    max_t = 10000 
    max_iso = 280
    
    H_m_abund    = 7.057e-1
    He_m_abund   = 2.750e-1
    C_m_abund    = 2.873e-3
    C12_m_abund  = 2.87e-3
    C13_m_abund  = 3.22e-5
    N_m_abund    = 1.105e-3 # 8.300e-4
    N14_m_abund  = 1.105e-3
    N15_m_abund  = 4.363e-6
    O_m_abund    = 6.063e-3 # 5.161e-3 Asplund05 ; 6.063e-3 <- standard (assumes A=8.73 e.g. GS01) ; 8.421e-3 ; 9.618e-3 ; 9.39e-3 gratton/yong/ngc6752paper
    O16_m_abund  = 6.063e-3
    O17_m_abund  = 2.456e-6
    O18_m_abund  = 1.340e-5
    Mg_m_abund   = 6.6e-4 # Chiappini use ~ 5.13e-4
    Si_m_abund   = 7.109e-4
    Si28_m_abund = 6.53e-4
    Si29_m_abund = 3.426e-5
    Si30_m_abund = 2.352e-5
    Fe_m_abund   = 1.273e-3
    Fe56_m_abund = 1.169e-3
    Fe57_m_abund = 2.855e-5
    Fe58_m_abund = 3.687e-6
    Na_m_abund   = 3.441e-5
    Al_m_abund   = 5.8e-5
    Al26_m_abund = 0.0
    Ne_m_abund   = 1.749E-03
    Ne20_m_abund = 1.619E-03
    Ne21_m_abund = 4.127e-6
    Ne22_m_abund = 1.302e-4
    Mg25_m_abund = 6.766e-5 # 7.76837E-05
    Mg26_m_abund = 7.76e-5  # 8.91307E-05
    Mg24_m_abund = 5.148e-4 # 5.91799E-04
    Ga_m_abund   = 6.592e-8 # was 3.69e-8 (ref pg 10 Rauscher et al ApJ 576)
    Ge_m_abund   = 2.187e-7 #; was 1.31e-7 (ref pg 10 Rauscher et al ApJ 576)
    S_m_abund    = 3.83e-4 # (ref G&S01)
    Zn_m_abund   = 2.087e-6 # (ref pg 10 Rauscher et al ApJ 576)
    P_m_abund    = 6.309e-6 # G&S
    K_m_abund    = 3.71e-6
    Ca_m_abund   = 6.32e-5
    Ti_m_abund   = 3.09e-6
    Mn_m_abund   = 1.315e-5
    Co_m_abund   = 3.38e-6
    Ni_m_abund   = 7.31e-5
    Cu_m_abund   = 8.73e-7
    V_m_abund    = 3.77e-7
    Cr_m_abund   = 1.8e-5
    As_m_abund   = 1.24e-8
    Se_m_abund   = 1.43e-7
    Br_m_abund   = 2.41e-8
    Kr_m_abund   = 1.21e-7
    B_m_abund   = 4.79e-9
    F_m_abund   = 4.049e-7
    Cl_m_abund   = 4.7e-6 # 2.6e-6 ;  4.7e-6
    Ar_m_abund   = 7.09e-5
    Sc_m_abund   = 3.998e-8
    Pb_m_abund   = 1.64e-8
    Sr_m_abund   = 5.143e-8
    Y_m_abund   = 1.090e-8
    Zr_m_abund   = 2.563e-8
    Sn_m_abund   = 1.156e-8
    Ba_m_abund   = 1.557e-8
    La_m_abund   = 1.322e-9
    Eu_m_abund   = 3.471e-10
    Dy_m_abund   = 1.583e-9
    Sr84_m_abund	=	2.800e-10
    Sr86_m_abund	=	5.050e-9
    Sr87_m_abund	=	3.320e-9
    Sr88_m_abund	=	4.320e-8
    Zr90_m_abund	=	1.340e-8
    Zr91_m_abund	=	2.950e-9
    Zr92_m_abund	=	4.560e-9
    Zr94_m_abund	=	4.710e-9
    Zr96_m_abund	=	7.740e-10
    Ba130_m_abund	=	1.540e-11
    Ba132_m_abund	=	1.510e-11
    Ba134_m_abund	=	3.690e-10
    Ba135_m_abund	=	1.010e-9
    Ba136_m_abund	=	1.210e-9
    Ba137_m_abund	=	1.750e-9
    Ba138_m_abund	=	1.120e-8
    fname = 'stars_t' + str(tstep) + '.dat'
    data = None
    file = open(fname, 'r')
    data  = np.loadtxt(fname, dtype=float)

    #;;;;;;;; reading in Isotopic values
    #radius = np.zeros(max_r+1,dtype=float)
    #time = np.zeros(max_t+1,dtype=float)
    #isotope = np.zeros((max_iso+1,max_r+1,max_t+1),dtype=float)
    #SFR = np.zeros((max_r+1,max_t+1),dtype=float)
    #
    #
    #R_SNII = np.zeros((max_r+1,max_t+1),dtype=float)
    #R_SNIa = np.zeros((max_r+1,max_t+1),dtype=float)
    #R_IMS = np.zeros((max_r+1,max_t+1),dtype=float)
    #gas_dens = np.zeros((max_r+1,max_t+1),dtype=float)
    #rem_dens = np.zeros((max_r+1,max_t+1),dtype=float)
    #tot_dens = np.zeros((max_r+1,max_t+1),dtype=float)
    #star_dens = np.zeros((max_r+1,max_t+1),dtype=float)
    #gas_frac = np.zeros((max_r+1,max_t+1),dtype=float)
    H = np.zeros(data.shape[0],dtype=float)
    #H2 = np.zeros((max_r+1,max_t+1),dtype=float)
    He = np.zeros(data.shape[0],dtype=float)
    #Pb = np.zeros((max_r+1,max_t+1),dtype=float)
    #Sr = np.zeros((max_r+1,max_t+1),dtype=float)
    #Y  = np.zeros((max_r+1,max_t+1),dtype=float)
    #Zr = np.zeros((max_r+1,max_t+1),dtype=float)
    #Sn = np.zeros((max_r+1,max_t+1),dtype=float)
    #Ba = np.zeros((max_r+1,max_t+1),dtype=float)
    #La = np.zeros((max_r+1,max_t+1),dtype=float)
    #Eu = np.zeros((max_r+1,max_t+1),dtype=float)
    #Dy = np.zeros((max_r+1,max_t+1),dtype=float)
    Carbon = np.zeros(data.shape[0],dtype=float)
    #C12 = np.zeros((max_r+1,max_t+1),dtype=float)
    #C13 = np.zeros((max_r+1,max_t+1),dtype=float)
    #Al26 = np.zeros((max_r+1,max_t+1),dtype=float)
    N = np.zeros(data.shape[0],dtype=float)
    #N14 = np.zeros((max_r+1,max_t+1),dtype=float)
    #N15 = np.zeros((max_r+1,max_t+1),dtype=float)
    O = np.zeros(data.shape[0],dtype=float)
    #O16 = np.zeros((max_r+1,max_t+1),dtype=float)
    #O17 = np.zeros((max_r+1,max_t+1),dtype=float)
    #O18 = np.zeros((max_r+1,max_t+1),dtype=float)
    #Sulf = np.zeros((max_r+1,max_t+1),dtype=float)
    #Zn = np.zeros((max_r+1,max_t+1),dtype=float)
    Mg = np.zeros(data.shape[0],dtype=float)
    #Mg26 = np.zeros((max_r+1,max_t+1),dtype=float)
    #Mg25 = np.zeros((max_r+1,max_t+1),dtype=float)
    #Mg24 = np.zeros((max_r+1,max_t+1),dtype=float)
    #Gal = np.zeros((max_r+1,max_t+1),dtype=float)
    #Ger = np.zeros((max_r+1,max_t+1),dtype=float)
    #Neon = np.zeros((max_r+1,max_t+1),dtype=float)
    #Ne20 = np.zeros((max_r+1,max_t+1),dtype=float)
    #Ne21 = np.zeros((max_r+1,max_t+1),dtype=float)
    #Ne22 = np.zeros((max_r+1,max_t+1),dtype=float)
    Si = np.zeros(data.shape[0],dtype=float)
    #Si28 = np.zeros((max_r+1,max_t+1),dtype=float)
    #Si29 = np.zeros((max_r+1,max_t+1),dtype=float)
    #Si30 = np.zeros((max_r+1,max_t+1),dtype=float)
    Fe = np.zeros(data.shape[0],dtype=float)
    #Fe56 = np.zeros((max_r+1,max_t+1),dtype=float)
    #Fe57 = np.zeros((max_r+1,max_t+1),dtype=float)
    #Fe58 = np.zeros((max_r+1,max_t+1),dtype=float)
    #Na = np.zeros((max_r+1,max_t+1),dtype=float)
    #Al = np.zeros((max_r+1,max_t+1),dtype=float)
    #Phos= np.zeros((max_r+1,max_t+1),dtype=float)
    #Kpot= np.zeros((max_r+1,max_t+1),dtype=float)
    #Cal= np.zeros((max_r+1,max_t+1),dtype=float)
    #Ti= np.zeros((max_r+1,max_t+1),dtype=float)
    #Mn= np.zeros((max_r+1,max_t+1),dtype=float)
    #Co= np.zeros((max_r+1,max_t+1),dtype=float)
    #Ni= np.zeros((max_r+1,max_t+1),dtype=float)
    #Cu= np.zeros((max_r+1,max_t+1),dtype=float)
    #Fluo= np.zeros((max_r+1,max_t+1),dtype=float)
    #Cl= np.zeros((max_r+1,max_t+1),dtype=float)
    #Ar= np.zeros((max_r+1,max_t+1),dtype=float)
    #Sc= np.zeros((max_r+1,max_t+1),dtype=float)
    #V= np.zeros((max_r+1,max_t+1),dtype=float)
    #Cr= np.zeros((max_r+1,max_t+1),dtype=float)
    #As= np.zeros((max_r+1,max_t+1),dtype=float)
    #Se= np.zeros((max_r+1,max_t+1),dtype=float)
    #Br= np.zeros((max_r+1,max_t+1),dtype=float)
    #Kr= np.zeros((max_r+1,max_t+1),dtype=float)
    #Boron= np.zeros((max_r+1,max_t+1),dtype=float)
    #Z = np.zeros((max_r+1,max_t+1),dtype=float)
    metallicity = np.zeros(data.shape[0],dtype=float)
    #ConO = np.zeros((max_r+1,max_t+1),dtype=float)
    #NonO = np.zeros((max_r+1,max_t+1),dtype=float)
    #MgonFe = np.zeros((max_r+1,max_t+1),dtype=float)
    #Mg25onFe = np.zeros((max_r+1,max_t+1),dtype=float)
    #Mg26onFe = np.zeros((max_r+1,max_t+1),dtype=float)
    #Mg24onFe = np.zeros((max_r+1,max_t+1),dtype=float)
    #GaonFe = np.zeros((max_r+1,max_t+1),dtype=float)
    #GeonFe = np.zeros((max_r+1,max_t+1),dtype=float)
    #NeonFe = np.zeros((max_r+1,max_t+1),dtype=float)
    #MnonFe = np.zeros((max_r+1,max_t+1),dtype=float)
    #PbonFe = np.zeros((max_r+1,max_t+1),dtype=float)
    #SronFe = np.zeros((max_r+1,max_t+1),dtype=float)
    #YonFe  = np.zeros((max_r+1,max_t+1),dtype=float)
    #ZronFe = np.zeros((max_r+1,max_t+1),dtype=float)
    #SnonFe = np.zeros((max_r+1,max_t+1),dtype=float)
    #BaonFe = np.zeros((max_r+1,max_t+1),dtype=float)
    #BaonEu = np.zeros((max_r+1,max_t+1),dtype=float)
    #SionEu = np.zeros((max_r+1,max_t+1),dtype=float)
    #BaonY = np.zeros((max_r+1,max_t+1),dtype=float)
    #LaonFe = np.zeros((max_r+1,max_t+1),dtype=float)
    #EuonFe = np.zeros((max_r+1,max_t+1),dtype=float)
    #DyonFe = np.zeros((max_r+1,max_t+1),dtype=float)
    #KonFe = np.zeros((max_r+1,max_t+1),dtype=float)
    #PonFe = np.zeros((max_r+1,max_t+1),dtype=float)
    SionO = np.zeros(data.shape[0],dtype=float)
    MgonFe = np.zeros(data.shape[0],dtype=float)

    O_t = np.zeros((cell_no,cell_no,max_t),dtype=float)

    #Si28onFe = np.zeros((max_r+1,max_t+1),dtype=float)
    #Si29onFe = np.zeros((max_r+1,max_t+1),dtype=float)
    #Si30onFe = np.zeros((max_r+1,max_t+1),dtype=float)
    #SonFe = np.zeros((max_r+1,max_t+1),dtype=float)
    #SonH = np.zeros((max_r+1,max_t+1),dtype=float)
    FeonH = np.zeros(data.shape[0],dtype=float)
    #MgonH = np.zeros((max_r+1,max_t+1),dtype=float)
    #FeonH2 = np.zeros((max_r+1,max_t+1),dtype=float)
    #OonH = np.zeros((max_r+1,max_t+1),dtype=float)
    #OonH_sq = np.zeros((max_r+1,max_t+1),dtype=float)
    #OonFe = np.zeros((max_r+1,max_t+1),dtype=float)
    #HeonFe = np.zeros((max_r+1,max_t+1),dtype=float)
    #ZnonFe = np.zeros((max_r+1,max_t+1),dtype=float)
    #ConFe = np.zeros((max_r+1,max_t+1),dtype=float)
    #CoonFe = np.zeros((max_r+1,max_t+1),dtype=float)
    #TionFe = np.zeros((max_r+1,max_t+1),dtype=float)
    #NionFe = np.zeros((max_r+1,max_t+1),dtype=float)
    #ClonFe = np.zeros((max_r+1,max_t+1),dtype=float)
    #CaonFe = np.zeros((max_r+1,max_t+1),dtype=float)
    #CuonFe = np.zeros((max_r+1,max_t+1),dtype=float)
    #VonFe = np.zeros((max_r+1,max_t+1),dtype=float)
    #CronFe = np.zeros((max_r+1,max_t+1),dtype=float)
    #Ca = np.zeros((max_r+1,max_t+1),dtype=float)
    #Ga = np.zeros((max_r+1,max_t+1),dtype=float)
    #Ir = np.zeros((max_r+1,max_t+1),dtype=float)
    #Os = np.zeros((max_r+1,max_t+1),dtype=float)
    #Re = np.zeros((max_r+1,max_t+1),dtype=float)
    #W = np.zeros((max_r+1,max_t+1),dtype=float)
    #Yb = np.zeros((max_r+1,max_t+1),dtype=float)
    #Tm = np.zeros((max_r+1,max_t+1),dtype=float)
    #Er = np.zeros((max_r+1,max_t+1),dtype=float)
    #Ho = np.zeros((max_r+1,max_t+1),dtype=float)
    #Dy = np.zeros((max_r+1,max_t+1),dtype=float)
    #Tb = np.zeros((max_r+1,max_t+1),dtype=float)
    #Gd = np.zeros((max_r+1,max_t+1),dtype=float)
    #Eu = np.zeros((max_r+1,max_t+1),dtype=float)
    #Sm = np.zeros((max_r+1,max_t+1),dtype=float)
    #Nd = np.zeros((max_r+1,max_t+1),dtype=float)
    #Pr = np.zeros((max_r+1,max_t+1),dtype=float)
    #Ce = np.zeros((max_r+1,max_t+1),dtype=float)
    #La = np.zeros((max_r+1,max_t+1),dtype=float)
    #Ba = np.zeros((max_r+1,max_t+1),dtype=float)
    #Cs = np.zeros((max_r+1,max_t+1),dtype=float)
    #Germanium = np.zeros((max_r+1,max_t+1),dtype=float)
    #Iodine = np.zeros((max_r+1,max_t+1),dtype=float)
    #Te = np.zeros((max_r+1,max_t+1),dtype=float)
    #I = np.zeros((max_r+1,max_t+1),dtype=float)
    #Sn = np.zeros((max_r+1,max_t+1),dtype=float)
    #Sb = np.zeros((max_r+1,max_t+1),dtype=float)
    #Pd = np.zeros((max_r+1,max_t+1),dtype=float)
    #Ru = np.zeros((max_r+1,max_t+1),dtype=float)
    #Rh = np.zeros((max_r+1,max_t+1),dtype=float)
    #Mo = np.zeros((max_r+1,max_t+1),dtype=float)
    #In = np.zeros((max_r+1,max_t+1),dtype=float)
    #Pt = np.zeros((max_r+1,max_t+1),dtype=float)
    #Au = np.zeros((max_r+1,max_t+1),dtype=float)
    #Ta = np.zeros((max_r+1,max_t+1),dtype=float)
    #Hf = np.zeros((max_r+1,max_t+1),dtype=float)
    #Lu = np.zeros((max_r+1,max_t+1),dtype=float)
    #Cadmium = np.zeros((max_r+1,max_t+1),dtype=float)
    #Ag = np.zeros((max_r+1,max_t+1),dtype=float)
    #Nb = np.zeros((max_r+1,max_t+1),dtype=float)
    #Rb = np.zeros((max_r+1,max_t+1),dtype=float)
    #Y = np.zeros((max_r+1,max_t+1),dtype=float)
    #Xe = np.zeros((max_r+1,max_t+1),dtype=float)
    #AsonFe = np.zeros((max_r+1,max_t+1),dtype=float)
    #SeonFe = np.zeros((max_r+1,max_t+1),dtype=float)
    #BronFe = np.zeros((max_r+1,max_t+1),dtype=float)
    #KronFe = np.zeros((max_r+1,max_t+1),dtype=float)
    #BonFe = np.zeros((max_r+1,max_t+1),dtype=float)
    #FonFe = np.zeros((max_r+1,max_t+1),dtype=float)
    #AronFe = np.zeros((max_r+1,max_t+1),dtype=float)
    #SconFe = np.zeros((max_r+1,max_t+1),dtype=float)
    #NonFe = np.zeros((max_r+1,max_t+1),dtype=float)
    #FeonO = np.zeros((max_r+1,max_t+1),dtype=float)
    #NaonO = np.zeros((max_r+1,max_t+1),dtype=float)
    #FonO = np.zeros((max_r+1,max_t+1),dtype=float)
    #NaonFe = np.zeros((max_r+1,max_t+1),dtype=float)
    #NionFe = np.zeros((max_r+1,max_t+1),dtype=float)
    #AlonFe = np.zeros((max_r+1,max_t+1),dtype=float)
    #AlonMg = np.zeros((max_r+1,max_t+1),dtype=float)
    #OonMg = np.zeros((max_r+1,max_t+1),dtype=float)
    #SionMg = np.zeros((max_r+1,max_t+1),dtype=float)
    #OonHgrad = np.zeros((max_t+1),dtype=float)
    #infallrate = np.zeros((max_r+1,max_t+1),dtype=float)
    #infallrate_pure = np.zeros((max_r+1,max_t+1),dtype=float)
    #Sr84	= np.zeros((max_r+1,max_t+1),dtype=float)
    #Sr86	= np.zeros((max_r+1,max_t+1),dtype=float)
    #Sr87	= np.zeros((max_r+1,max_t+1),dtype=float)
    #Sr88	= np.zeros((max_r+1,max_t+1),dtype=float)
    #Zr90	= np.zeros((max_r+1,max_t+1),dtype=float)
    #Zr91	= np.zeros((max_r+1,max_t+1),dtype=float)
    #Zr92	= np.zeros((max_r+1,max_t+1),dtype=float)
    #Zr94	= np.zeros((max_r+1,max_t+1),dtype=float)
    #Zr96	= np.zeros((max_r+1,max_t+1),dtype=float)
    #Ba130	= np.zeros((max_r+1,max_t+1),dtype=float)
    #Ba132	= np.zeros((max_r+1,max_t+1),dtype=float)
    #Ba134	= np.zeros((max_r+1,max_t+1),dtype=float)
    #Ba135	= np.zeros((max_r+1,max_t+1),dtype=float)
    #Ba136	= np.zeros((max_r+1,max_t+1),dtype=float)
    #Ba137	= np.zeros((max_r+1,max_t+1),dtype=float)
    #Ba138	= np.zeros((max_r+1,max_t+1),dtype=float)
    #GaonMg = np.zeros((max_r+1,max_t+1),dtype=float)
    #GeonMg = np.zeros((max_r+1,max_t+1),dtype=float)
    #AsonMg = np.zeros((max_r+1,max_t+1),dtype=float)
    #SeonMg = np.zeros((max_r+1,max_t+1),dtype=float)
    #KronMg = np.zeros((max_r+1,max_t+1),dtype=float)
    #SronMg = np.zeros((max_r+1,max_t+1),dtype=float)
    #YonMg = np.zeros((max_r+1,max_t+1),dtype=float)
    #ZronMg = np.zeros((max_r+1,max_t+1),dtype=float)
    #SnonMg = np.zeros((max_r+1,max_t+1),dtype=float)
    #BaonMg = np.zeros((max_r+1,max_t+1),dtype=float)
    #LaonMg = np.zeros((max_r+1,max_t+1),dtype=float)
    #EuonMg = np.zeros((max_r+1,max_t+1),dtype=float)
    #DyonMg = np.zeros((max_r+1,max_t+1),dtype=float)
    #PbonMg = np.zeros((max_r+1,max_t+1),dtype=float)
    #MgonSi = np.zeros((max_r+1,max_t+1),dtype=float)
    #Li = np.zeros((max_r+1,max_t+1),dtype=float)
    #Hg = np.zeros((max_r+1,max_t+1),dtype=float)
    #Tl = np.zeros((max_r+1,max_t+1),dtype=float)
    #Bi = np.zeros((max_r+1,max_t+1),dtype=float)


    avgmass = np.zeros(data.shape[0],dtype=float)
    age = np.zeros(data.shape[0],dtype=float)

    logmass = np.zeros(data.shape[0],dtype=float)
    avglogmass = np.zeros(data.shape[0],dtype=float)
    count = 0
#    #nugrid iso version
#        age[i] = data[count][0]
#        H[i] = data[count][1] + data[count][2]
#        #                He[i][j][k] = data[count][1]
#        #                Carbon[i][j][k] = data[count][2]
#        #                N[i][j][k] = data[count][3]
#        O[i] = data[count][10] + data[count][11] + data[count][12]
#        Mg[i] = data[count][19] + data[count][20] + data[count][21]
#        Si[i] = data[count][22] + data[count][23] + data[count][24]
#        Fe[i] = data[count][60] + data[count][61] + data[count][62] + data[count][63]
#
#
    for i in range(0,data.shape[0]):
        age[i] = data[count][0]
        H[i] = data[count][1]
        #                He[i][j][k] = data[count][1]
        #                Carbon[i][j][k] = data[count][2]
        #                N[i][j][k] = data[count][3]
        O[i] = data[count][4]
        Mg[i] = data[count][5]
        Si[i] = data[count][6]
        Fe[i] = data[count][7]
        metallicity[i] = np.sum(data[count][2:])
        count+=1


    lmetallicity = np.zeros(metallicity.shape[0],dtype=float)
    FeonH_avg = 0.0
    metal_avg = 0.0
    for i in range(0,data.shape[0]):
        FeonH[i] =   math.log10(((Fe[i]+1e-20)/(H[i]+1e-20))) - math.log10( Fe_m_abund/H_m_abund )
        MgonFe[i] =   math.log10(((Mg[i]+1e-20)/(Fe[i]+1e-20))) - math.log10( Mg_m_abund/Fe_m_abund )
        SionO[i] =   math.log10(((Si[i]+1e-20)/(O[i]+1e-20))) - math.log10( Si_m_abund/O_m_abund )

        FeonH_avg += FeonH[i]
        metal_avg += metallicity[i]
        lmetallicity[i] = math.log10(metallicity[i]+1e-20)
    print("Average [Fe/H] of ", FeonH_avg/(data.shape[0]))
    print("Number of stars alive ", data.shape[0])

    weights = np.ones_like(FeonH)/float(len(FeonH))

    plt.hist(FeonH, weights=weights,bins=20)
    plt.ylabel("Fraction of Stars")
    plt.xlabel("Fe/H]")
#plt.yscale('log')
#    plt.xlim([-5.5,1])
    plt.show()
    plt.clf()
    plt.hist2d(age/1000,FeonH, bins=(data.shape[0]/500), norm=LogNorm())
    plt.xlabel("Age (Gyr)")
    plt.ylabel("[Fe/H]")
#plt.ylim([-5.5,0])
    plt.savefig('Age_metal.eps')
    plt.show()
    plt.clf()
    plt.hist2d(FeonH, MgonFe, bins=(data.shape[0]/200), norm=LogNorm())
#    plt.xlim([-5.5,-0.5])
#    plt.ylim([-1.5,1.2])
    plt.ylabel("[Mg/Fe]")
    plt.xlabel("[Fe/H]")
    plt.show()
    #
    #for i in range(0,342999):
    #    x[i] = data[4+i*4]
    #for j in range(0,342999):
    #    y[j] = data[5+j*4]
    #for k in range(0,342999):
    #    z[k] = data[6+k*4]
    #for l in range(0,342999):
    #    mass[l] = data[7+l*4]
    #    logmass[l] = math.log10(mass[l])
    #
    ### Creating a 3D plot...
    #zmass = np.zeros(342999,dtype=float)
    #zx = np.zeros(342999,dtype=float)
    #zy = np.zeros(342999,dtype=float)
    #nzmass = np.zeros(342999,dtype=float)
    #nzx = np.zeros(342999,dtype=float)
    #nzy = np.zeros(342999,dtype=float)
    #
    #
    #for i in range(0,324999):
    #    if (z[i] == 35*28.5):
    #        zmass[i] = mass[i]
    #        zx[i] = x[i]
    #        zy[i] = y[i]
    #    else:
    #        i+=1
    #for i in range(0,324999):
    #    if (zmass[i] == 0):
    #        nzmass = np.delete(zmass,i)
    #        nzx = np.delete(zx,i)
    #        nzy = np.delete(zy,i)



    #                z[k] = cell_size*k
    #        count += 1
    #                if (Fe[i][j][k] <=0 ):
    #                    print("Fe", Fe[i][j][k])
    #                    print("mass", mass[i][j][k])

    #    for i in range(0,cell_no):
    #        for j in range(0,cell_no):
    #            for k in range(0,cell_no):
    #                avgmass[i][j] += mass[i][j][k]

    #            SionO_t[i][j][tstep] = Si[i][j][30]/O[i][j][30]
    #            O_t[i][j][tstep] = O[i][j][30]
    #            plt.plot(O_t[i][j][tstep],SionO_t[i][j][tstep], marker='.')

    #            avgmass[i][j] /= cell_no
    #            avglogmass[i][j] = math.log10(avgmass[i][j])

    #    print np.max(mass)
    #    print np.min(mass)
    # plt.plot(O_t[i][j][tstep],SionO_t[i][j][tstep], marker='.')



