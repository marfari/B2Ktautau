import ROOT 
import sys
import numpy as np
from array import array
import sympy as sp
from pygsl  import multiroots, errno
import pygsl._numobj as numx
import matplotlib.pyplot as plt
from numba import njit
import time

np_load_old = np.load
np.load = lambda *a,**k: np_load_old(*a, allow_pickle=True, **k)

np.set_printoptions(suppress=True)

# VARIABLES AND CONSTANTS
dimM = 22
dimX = 23
dimC = 24

BV_offline = ROOT.Math.XYZPoint(0,0,0)

X = []
X_ERR = []
F = []
D = []
for i in range(dimM+dimX):
    X.append(array('d', [0]))
    X_ERR.append(array('d', [0]))
for i in range(dimM+dimX+dimC):
    F.append(array('d', [0]))
for i in range(21):
    D.append(array('d', [0]))

STATUS = array('i', [0])
NITER = array('i', [0])
tolerance = array('d', [0])
chi2 = array('d', [0])
MB = array('d', [0])
dMB = array('d', [0])
init = array('i', [0])

MB0 = array('d', [0])
MB1 = array('d', [0])
MB2 = array('d', [0])
MB3 = array('d', [0])
MB4 = array('d', [0])

mtau = 1776.86
mkaon = 493.677
mnu = 0

x = []
f = []
status = -1
nIter = 0
chi2_value = 0
tol_value = 100
x_true = []

lambdify_f = 0
lambdify_df = 0
lambdify_dfm = 0
lambdify_chi2 = 0
lambdify_bH = 0

computeDerivatives = False
change_L = 0

N_local_minima = array('i', [0])
M_local_minima = []
Chi2_local_minima = []
kMaxN = 20
M_local = array('d', [0]*kMaxN)
Chi2_local = array('d', [0]*kMaxN)
theTime = array('d', [0])

m_symbols = sp.IndexedBase('m_symbols')
W_symbols = sp.IndexedBase('W_symbols')
RPz_symbol = sp.Symbol('RPz_symbol')
xm_symbols = sp.symbols('PVx PVy PVz DV1x DV1y DV1z p3pi1x p3pi1y p3pi1z E3pi1 DV2x DV2y DV2z p3pi2x p3pi2y p3pi2z E3pi2 RPx RPy pKx pKy pKz')
xu_symbols = sp.symbols('BVx BVy BVz pBx pBy pBz mB_squared ptau1x ptau1y ptau1z Etau1 pnu1x pnu1y pnu1z Enu1 ptau2x ptau2y ptau2z Etau2 pnu2x pnu2y pnu2z Enu2')
lb_symbols = sp.symbols('lb0 lb1 lb2 lb3 lb4 lb5 lb6 lb7 lb8 lb9 lb10 lb11 lb12 lb13 lb14 lb15 lb16 lb17 lb18 lb19 lb20 lb21 lb22 lb23')
x_symbols = sp.symbols('PVx PVy PVz DV1x DV1y DV1z p3pi1x p3pi1y p3pi1z E3pi1 DV2x DV2y DV2z p3pi2x p3pi2y p3pi2z E3pi2 RPx RPy pKx pKy pKz BVx BVy BVz pBx pBy pBz mB_squared ptau1x ptau1y ptau1z Etau1 pnu1x pnu1y pnu1z Enu1 ptau2x ptau2y ptau2z Etau2 pnu2x pnu2y pnu2z Enu2 lb0 lb1 lb2 lb3 lb4 lb5 lb6 lb7 lb8 lb9 lb10 lb11 lb12 lb13 lb14 lb15 lb16 lb17 lb18 lb19 lb20 lb21 lb22 lb23')
x_symbols_2nd_test = sp.symbols('lb0 lb1 lb2 lb3 lb4 lb5 lb6 lb7 lb8 lb9 lb10 lb11 lb12 lb13 lb14 lb15 lb16 lb17 lb18 lb19 lb20 lb21 lb22 lb23 PVx PVy PVz DV1x DV1y DV1z p3pi1x p3pi1y p3pi1z E3pi1 DV2x DV2y DV2z p3pi2x p3pi2y p3pi2z E3pi2 RPx RPy pKx pKy pKz BVx BVy BVz pBx pBy pBz mB_squared ptau1x ptau1y ptau1z Etau1 pnu1x pnu1y pnu1z Enu1 ptau2x ptau2y ptau2z Etau2 pnu2x pnu2y pnu2z Enu2')

# FUNCTIONS:
def chisquare(params):
    m = params[0]
    W = params[1]

    chi2 = 0
    for i in range(dimM):
        for j in range(dimM):
            chi2 += ( m[i] - xm_symbols[i] )*W[i,j]*( m[j] - xm_symbols[j] )
    return chi2

def chisquare_value(x, params):
    m = params[0]
    W = params[1]
    RPz = params[2]

    return lambdify_chi2(x,m,W,RPz)

def lagrangian(params):
    RPz = params[2]

    def EB(a, b, c, d):
        return (a + b**2 + c**2 + d**2)**(0.5)

    def EK(e, f, g):
        return (mkaon**2 + e**2 + f**2 + g**2)**(0.5)

    # xm
    PVx = xm_symbols[0]
    PVy = xm_symbols[1]
    PVz = xm_symbols[2]
    DV1x = xm_symbols[3]
    DV1y = xm_symbols[4]
    DV1z = xm_symbols[5]
    p3pi1x = xm_symbols[6]
    p3pi1y = xm_symbols[7]
    p3pi1z = xm_symbols[8]
    E3pi1 = xm_symbols[9]
    DV2x = xm_symbols[10]
    DV2y = xm_symbols[11]
    DV2z = xm_symbols[12]
    p3pi2x = xm_symbols[13]
    p3pi2y = xm_symbols[14]
    p3pi2z = xm_symbols[15]
    E3pi2 = xm_symbols[16]
    RPx = xm_symbols[17]
    RPy = xm_symbols[18]
    pKx = xm_symbols[19]
    pKy = xm_symbols[20]
    pKz = xm_symbols[21]

    # xu
    BVx = xu_symbols[0]
    BVy = xu_symbols[1]
    BVz = xu_symbols[2]
    pBx = xu_symbols[3]
    pBy = xu_symbols[4]
    pBz = xu_symbols[5]
    mB_squared = xu_symbols[6]
    ptau1x = xu_symbols[7]
    ptau1y = xu_symbols[8]
    ptau1z = xu_symbols[9]
    Etau1 = xu_symbols[10]
    pnu1x = xu_symbols[11]
    pnu1y = xu_symbols[12]
    pnu1z = xu_symbols[13]
    Enu1 = xu_symbols[14]
    ptau2x = xu_symbols[15]
    ptau2y = xu_symbols[16]
    ptau2z = xu_symbols[17]
    Etau2 = xu_symbols[18]
    pnu2x = xu_symbols[19]
    pnu2y = xu_symbols[20]
    pnu2z = xu_symbols[21]
    Enu2 = xu_symbols[22]

    g = []
    # BV must be in K+ trajectory
    g.append( pKz*( BVx - RPx ) - pKx*( BVz - RPz ) ) 
    g.append( pKz*( BVy - RPy ) - pKy*( BVz - RPz ) ) 
    # ptau1 must point back to BV
    g.append( ptau1x*( DV1z - BVz ) - ptau1z*( DV1x - BVx ) )
    g.append( ptau1y*( DV1z - BVz ) - ptau1z*( DV1y - BVy ) )
    # ptau2 must point back to BV
    g.append( ptau2x*( DV2z - BVz ) - ptau2z*( DV2x - BVx ) )
    g.append( ptau2y*( DV2z - BVz ) - ptau2z*( DV2y - BVy ) ) 
    # pB must point back to the PV
    g.append( pBx*( BVz - PVz ) - pBz*( BVx - PVx ) ) 
    g.append( pBy*( BVz - PVz ) - pBz*( BVy - PVy ) ) 
    # 4-momentum conservation in DV1
    g.append( ptau1x - p3pi1x - pnu1x ) 
    g.append( ptau1y - p3pi1y - pnu1y ) 
    g.append( ptau1z - p3pi1z - pnu1z ) 
    g.append( Etau1 - E3pi1 - Enu1 ) 
    # tau+ and anti-nu mass constraints
    g.append( Etau1 - (mtau**2 + ptau1x**2 + ptau1y**2 + ptau1z**2)**(0.5) ) 
    g.append( Enu1 - (mnu**2 + pnu1x**2 + pnu1y**2 + pnu1z**2)**(0.5) ) 
    # 4-momentum conservation in DV2
    g.append( ptau2x - p3pi2x - pnu2x ) 
    g.append( ptau2y - p3pi2y - pnu2y ) 
    g.append( ptau2z - p3pi2z - pnu2z ) 
    g.append( Etau2 - E3pi2 - Enu2 ) 
    # tau- and nu mass constraints
    g.append( Etau2 - (mtau**2 + ptau2x**2 + ptau2y**2 + ptau2z**2)**(0.5) ) 
    g.append( Enu2 - (mnu**2 + pnu2x**2 + pnu2y**2 + pnu2z**2)**(0.5) ) 
    # 4-momentum conservation in BV
    g.append( pBx - ptau1x - ptau2x - pKx ) 
    g.append( pBy - ptau1y - ptau2y - pKy ) 
    g.append( pBz - ptau1z - ptau2z - pKz ) 
    g.append(EB(mB_squared,pBx,pBy,pBz) - Etau1 - Etau2 - EK(pKx,pKy,pKz))

    L = chisquare(params)
    for i in range(dimC):
        L += lb_symbols[i]*g[i]
    return L

##njit
def equations_f_exp(params):
    function = lagrangian(params)
    f = np.array( [ function.diff(x_symbols[i]) for i in range(dimM+dimX+dimC) ] )
    return f

##njit
def equations_df_exp(params):
    # This is the matrix A
    f = equations_f_exp(params)
    df = np.array( [ [f[i].diff(x_symbols[j]) for j in range(dimM+dimX+dimC)] for i in range(dimM+dimX+dimC)] )
    return df

##njit
def equations_dfm_exp(params):
    # This is the matrix B
    f = equations_f_exp(params)
    dfm = np.array( [ [f[i].diff(m_symbols[j]) for j in range(dimM)] for i in range(dimM+dimX+dimC)] )
    return dfm

def bordered_Hessian(params):
    function = lagrangian(params)
    f = np.array( [ function.diff(x_symbols_2nd_test[i]) for i in range(dimM+dimX+dimC) ] )
    df = np.array( [ [f[i].diff(x_symbols_2nd_test[j]) for j in range(dimM+dimX+dimC)] for i in range(dimM+dimX+dimC)] )
    return df

def get_principal_minor(H, i):
    A = np.delete(H,np.s_[i:dimM+dimX+dimC],axis=1)
    B = np.delete(A,np.s_[i:dimM+dimX+dimC],axis=0)

    return B

def second_derivative_test(x,params):
    m = params[0]
    W = params[1]
    RPz = params[2]

    # Bordered Hessian
    xprime = np.zeros(dimM+dimX+dimC)
    for i in range(dimC):
        xprime[i] = x[dimM+dimX+i]
    for i in range(dimM+dimX):
        xprime[dimC+i] = x[i]
    
    H = lambdify_bH(xprime,m,W,RPz)

    passes_test = True

    global D

    for i in range(49,69+1):
        s = np.sign( np.linalg.det(get_principal_minor(H,i)) ) 
        D[i-49][0] = s

        if( s < 0 ):
            # print(np.linalg.det(get_principal_minor(H,i)))
            passes_test = False

    return passes_test

##njit
def equations_f(x,params):
    m = params[0]
    W = params[1]
    RPz = params[2]

    return lambdify_f(x,m,W,RPz)

##njit
def equations_df(x,params): 
    m = params[0]
    W = params[1]
    RPz = params[2]

    df = lambdify_df(x,m,W,RPz)

    return df

##njit
def equations_fdf(x, params):
    f = equations_f(x, params)
    df =  equations_df(x, params)
    return f, df

def U_cov(x,params,V):
    m = params[0]
    W = params[1]
    RPz = params[2]

    A = lambdify_df(x,m,W,RPz)
    B = lambdify_dfm(x,m,W,RPz)

    A_inverse = np.linalg.inv(A)
    C = -np.dot(A_inverse,B)

    C_transpose = np.transpose(C)

    U = np.dot(V,C_transpose)
    U = np.dot(C,U)

    return U

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< Initialisations <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#njit
def makeTransformation_point( pK, RP, point, flag ):
    deltaX = RP.x()
    deltaY = RP.y()
    deltaZ = RP.z()

    myShift = ROOT.Math.Translation3D( -deltaX, -deltaY, -deltaZ )

    if(pK.y() != 0):
        phi = -1*np.arctan(pK.x()/pK.y())
    else:
        phi = 0
    if((pK.z() != 0) and (pK.y() != 0)):
        theta = -1*np.arctan(pK.Rho() * (pK.y()/np.abs(pK.y()))/pK.z())
    else:
        theta = 0

    myRotation = ROOT.Math.EulerAngles( phi, theta, 0 )

    myTransformation = ROOT.Math.Transform3D( myRotation )*ROOT.Math.Transform3D( myShift )
    myTransformation_inv = myTransformation.Inverse()

    if(flag):
        point_t = myTransformation_inv(point)
    else:
        point_t = myTransformation(point)
    return point_t

#njit
def makeTransformation_vec( pK, RP, vector, flag ):
    deltaX = RP.x()
    deltaY = RP.y()
    deltaZ = RP.z()

    myShift = ROOT.Math.Translation3D( -deltaX, -deltaY, -deltaZ )

    if(pK.y() != 0):
        phi = -1*np.arctan(pK.x()/pK.y())
    else:
        phi = 0
    if( (pK.z() != 0) and (pK.y() != 0) ):
        theta = -1*np.arctan(pK.Rho() * (pK.y()/np.abs(pK.y()))/pK.z())
    else:
        theta = 0

    myRotation = ROOT.Math.EulerAngles( phi, theta, 0 )

    myTransformation = ROOT.Math.Transform3D( myRotation )*ROOT.Math.Transform3D( myShift )
    myTransformation_inv = myTransformation.Inverse()

    if(flag):
        vector_t = myTransformation_inv(vector)
    else:
        vector_t = myTransformation(vector)
    return vector_t

def x_initial_estimate(init, BV, species, params):
    m = params[0]
    W = params[1]
    RPz = params[2]
    eventNumber = params[3]

    x0 = np.zeros(dimM+dimX+dimC)

    # Initial values for xm (=m)
    for i in range(dimM):
        x0[i] = m[i]

    PV = ROOT.Math.XYZPoint( m[0], m[1], m[2] )
    DV1 = ROOT.Math.XYZPoint( m[3], m[4], m[5] )
    p3pi1 = ROOT.Math.XYZVector( m[6], m[7], m[8] )
    E3pi1 = m[9]
    DV2 = ROOT.Math.XYZPoint( m[10], m[11], m[12] )
    p3pi2 = ROOT.Math.XYZVector( m[13], m[14], m[15] )
    E3pi2 = m[16]
    RP = ROOT.Math.XYZPoint( m[17], m[18], RPz )
    pK = ROOT.Math.XYZVector( m[19], m[20], m[21] )

    EK  = np.sqrt( mkaon**2 + pK.Mag2() )
    m3pi1 = np.sqrt( E3pi1**2 - p3pi1.Mag2() )
    m3pi2 = np.sqrt( E3pi2**2 - p3pi2.Mag2() )
    
    # Initial values for xu
    if(init == 0): # MLP
        arr_m_1[0] = m[0]
        arr_m_2[0] = m[1]
        arr_m_3[0] = m[2]
        arr_m_4[0] = m[3]
        arr_m_5[0] = m[4]
        arr_m_6[0] = m[5]
        arr_m_7[0] = m[6]
        arr_m_8[0] = m[7]
        arr_m_9[0] = m[8]
        arr_m_10[0] = m[9]
        arr_m_11[0] = m[10]
        arr_m_12[0] = m[11]
        arr_m_13[0] = m[12]
        arr_m_14[0] = m[13]
        arr_m_15[0] = m[14]
        arr_m_16[0] = m[15]
        arr_m_17[0] = m[16]
        arr_m_18[0] = m[17]
        arr_m_19[0] = m[18]
        arr_m_20[0] = m[19]
        arr_m_21[0] = m[20]
        arr_m_22[0] = m[21]
        arr_Kp_RP_Z[0] = RPz
        arr_eventNumber[0] = (eventNumber%100)

        MLP_taup_PX = reader_taup_PX.EvaluateRegression("MLP")[0]
        MLP_taup_PY = reader_taup_PY.EvaluateRegression("MLP")[0]
        MLP_taup_PZ = reader_taup_PZ.EvaluateRegression("MLP")[0]
        MLP_taum_PX = reader_taum_PX.EvaluateRegression("MLP")[0]
        MLP_taum_PY = reader_taum_PY.EvaluateRegression("MLP")[0]
        MLP_taum_PZ = reader_taum_PZ.EvaluateRegression("MLP")[0]

        ptau1 = ROOT.Math.XYZVector( MLP_taup_PX, MLP_taup_PY, MLP_taup_PZ )
        ptau2 = ROOT.Math.XYZVector( MLP_taum_PX, MLP_taum_PY, MLP_taum_PZ )
        p3pi1 = ROOT.Math.XYZVector( m[6], m[7], m[8] )
        p3pi2 = ROOT.Math.XYZVector( m[13], m[14], m[15] )
        pnu1 = ptau1 - p3pi1
        pnu2 = ptau2 - p3pi2

        Etau1 = np.sqrt( mtau**2 + ptau1.Mag2() )
        Etau2 = np.sqrt( mtau**2 + ptau2.Mag2() )

        Enu1 = np.sqrt( pnu1.Mag2() )
        Enu2 = np.sqrt( mnu**2 + pnu2.Mag2() )

        pK = ROOT.Math.XYZVector( m[19], m[20], m[21] )
        EK = np.sqrt( mkaon**2 + pK.Mag2() )
        pB = ptau1 + ptau2 + pK
        EB = Etau1 + Etau2 + EK
        MB_squared = EB**2 - pB.Mag2()

        x0[dimM] = BV.x()
        x0[dimM+1] = BV.y()
        x0[dimM+2] = BV.z()
        x0[dimM+3] = pB.x()
        x0[dimM+4] = pB.y()
        x0[dimM+5] = pB.z()
        x0[dimM+6] = MB_squared
        x0[dimM+7] = ptau1.x()
        x0[dimM+8] = ptau1.y()
        x0[dimM+9] = ptau1.z()
        x0[dimM+10] = Etau1
        x0[dimM+11] = pnu1.x()
        x0[dimM+12] = pnu1.y()
        x0[dimM+13] = pnu1.z()
        x0[dimM+14] = Enu1
        x0[dimM+15] = ptau2.x()
        x0[dimM+16] = ptau2.y()
        x0[dimM+17] = ptau2.z()
        x0[dimM+18] = Etau2
        x0[dimM+19] = pnu2.x()
        x0[dimM+20] = pnu2.y()
        x0[dimM+21] = pnu2.z()
        x0[dimM+22] = Enu2

    elif(init == 1): # Anne
        PV_t = makeTransformation_point( pK, RP, PV, False )
        DV1_t = makeTransformation_point( pK, RP, DV1, False )
        p3pi1_t = makeTransformation_vec( pK, RP, p3pi1, False )
        DV2_t = makeTransformation_point( pK, RP, DV2, False )
        p3pi2_t = makeTransformation_vec( pK, RP, p3pi2, False )
        pK_t = makeTransformation_vec( pK, RP, pK, False )

        if(DV1_t.x() != 0):
            a1 = (DV1_t.y())/(DV1_t.x())
        else: 
            a1 = 1    
        if(DV2_t.x() != 0):
            a2 = (DV2_t.y())/(DV2_t.x())
        else:
            a2 = 1
        if(a2*PV_t.x() - PV_t.y() != 0):
            b = (PV_t.y() - a1*PV_t.x())/(a2*PV_t.x() - PV_t.y())
        else: 
            b = 1
        if(DV2_t.x() != 0):
            c = b*(DV1_t.x())/(DV2_t.x())
            d = b*((DV2_t.z() - DV1_t.z())/(DV2_t.x()))
        else: 
            c = 1
            d = 1
        if( (1+b)*DV1_t.x() - (1+c)*PV_t.x() != 0):
            e = ( (1+b)*(DV1_t.z() - PV_t.z()) + d*PV_t.x() )/( (1+b)*DV1_t.x() - (1+c)*PV_t.x() )
        else: 
            e = 1
        if((1+b)*DV1_t.x() - (1+c)*PV_t.x() != 0):
            f = ( PV_t.x()*np.sqrt( pK_t.Mag2() ) )/( (1+b)*DV1_t.x() - (1+c)*PV_t.x() )
        else:
            f = 1
        g = c*e + d
        h = f*c
        i = DV1_t.z() - e*DV1_t.x()
        j = f*DV1_t.x()

        x1 = p3pi1_t.x() + a1*p3pi1_t.y() + e*p3pi1_t.z()
        x2 = b*p3pi2_t.x() + a2*b*p3pi2_t.y() + g*p3pi2_t.z()

        p1 = 1 + a1**2 + e**2 - (x1/E3pi1)**2
        p2 = 2*e*f - ( mtau**2 + m3pi1**2 + 2*f*p3pi1_t.z() )*(x1/(E3pi1**2) )
        p3 = mtau**2 + f**2 - ( ( mtau**2 + m3pi1**2 + 2*f*p3pi1_t.z() )/(2*E3pi1) )**2 - mnu**2
        q1 = b**2 + (a2*b)**2 + g**2 - (x2/E3pi2)**2
        q2 = 2*g*h - ( mtau**2 + m3pi2**2 + 2*h*p3pi2_t.z() )*(x2/(E3pi2**2) )
        q3 =  mtau**2 + h**2 - ( ( mtau**2 + m3pi2**2 + 2*h*p3pi2_t.z() )/(2*E3pi2) )**2 - mnu**2

        if(p2*q1 - p1*q2 != 0):
            Ptau1x_t = (p1*q3 - p3*q1)/(p2*q1 - p1*q2)
        else:
            Ptau1x_t = 1
        Ptau1y_t = a1*Ptau1x_t
        Ptau1z_t = e*Ptau1x_t + f
        Ptau2x_t = b*Ptau1x_t
        Ptau2y_t = a2*b*Ptau1x_t
        Ptau2z_t = g*Ptau1x_t + h
        BVz_t = i - j*(1./Ptau1x_t)

        ptau1_t = ROOT.Math.XYZVector( Ptau1x_t, Ptau1y_t, Ptau1z_t )
        ptau2_t = ROOT.Math.XYZVector( Ptau2x_t, Ptau2y_t, Ptau2z_t )
        BV_t = ROOT.Math.XYZPoint(0, 0, BVz_t)

        ptau1 = makeTransformation_vec( pK, RP, ptau1_t, True )
        ptau2 = makeTransformation_vec( pK, RP, ptau2_t, True )
        BV_Anne = makeTransformation_point( pK, RP, BV_t, True )

        pnu1 = ptau1 - p3pi1
        pnu2 = ptau2 - p3pi2
        pB = ptau1 + ptau2 + pK

        Etau1 = np.sqrt( mtau**2 + ptau1.Mag2() )
        Etau2 = np.sqrt( mtau**2 + ptau2.Mag2() )
        Enu1 = np.sqrt( mnu**2 + pnu1.Mag2() )
        Enu2 = np.sqrt( mnu**2 + pnu2.Mag2() )

        EB = Etau1 + Etau2 + EK
        MB_squared = EB**2 - pB.Mag2()

        x0[dimM] = BV_Anne.x()
        x0[dimM+1] = BV_Anne.y()
        x0[dimM+2] = BV_Anne.z()
        x0[dimM+3] = pB.x()
        x0[dimM+4] = pB.y()
        x0[dimM+5] = pB.z()
        x0[dimM+6] = MB_squared
        x0[dimM+7] = ptau1.x()
        x0[dimM+8] = ptau1.y()
        x0[dimM+9] = ptau1.z()
        x0[dimM+10] = Etau1
        x0[dimM+11] = pnu1.x()
        x0[dimM+12] = pnu1.y()
        x0[dimM+13] = pnu1.z()
        x0[dimM+14] = Enu1
        x0[dimM+15] = ptau2.x()
        x0[dimM+16] = ptau2.y()
        x0[dimM+17] = ptau2.z()
        x0[dimM+18] = Etau2
        x0[dimM+19] = pnu2.x()
        x0[dimM+20] = pnu2.y()
        x0[dimM+21] = pnu2.z()
        x0[dimM+22] = Enu2

    elif(init == 2): # B->K* tautau initialisation taus direction based on vertices    
        # Magnitude of 3pi momenta
        p3pi1_mag = p3pi1.r()
        p3pi2_mag = p3pi2.r()

        # B+ flight direction
        bDir = (BV - PV).Unit()

        # Get K+ momentum perpendicular to the B+ flight direction
        pK_perp = pK - ROOT.Math.XYZVector((pK.Dot(bDir)) * bDir.x(), (pK.Dot(bDir)) * bDir.y(), (pK.Dot(bDir)) * bDir.z())

        # Get tau flight directions (from vertices)
        tau_dir1 = ROOT.Math.XYZVector(0,0,0)
        tau_dir2 = ROOT.Math.XYZVector(0,0,0)
        tau_dir1.SetXYZ( DV1.x() - BV.x(), DV1.y() - BV.y(), DV1.z() - BV.z() )
        tau_dir2.SetXYZ( DV2.x() - BV.x(), DV2.y() - BV.y(), DV2.z() - BV.z() )

        tau_dir1 = tau_dir1.Unit()
        tau_dir2 = tau_dir2.Unit()

        # Get tau direction unit vectors perpendicular to B+ flight direction
        tau_dir1_perp = (tau_dir1 - ROOT.Math.XYZVector((tau_dir1.Dot(bDir)) * bDir.x(), (tau_dir1.Dot(bDir)) * bDir.y(), (tau_dir1.Dot(bDir)) * bDir.z() )).Unit()
        tau_dir2_perp = (tau_dir2 - ROOT.Math.XYZVector((tau_dir2.Dot(bDir)) * bDir.x(), (tau_dir2.Dot(bDir)) * bDir.y(), (tau_dir2.Dot(bDir)) * bDir.z() )).Unit()

        # In plane perpendicular to B+ flight direction, get angles between tau momenta and K+ momentum
        cosphi1 = tau_dir1_perp.Dot(pK_perp.Unit())
        cosphi2 = tau_dir2_perp.Dot(pK_perp.Unit())

        phi1 = np.arccos(cosphi1)
        phi2 = np.arccos(cosphi2)

        # Calculate momentum component of taus in this plane
        pMag_tau1_perp = -1*pK_perp.R()/( cosphi1 + (cosphi2*(np.sin(phi1)/np.sin(phi2))) )
        pMag_tau2_perp = pMag_tau1_perp*(np.sin(phi1) / np.sin(phi2))

        p_tau1_perp = ROOT.Math.XYZVector(pMag_tau1_perp * tau_dir1_perp.x(), pMag_tau1_perp * tau_dir1_perp.y(), pMag_tau1_perp * tau_dir1_perp.z())
        p_tau2_perp = ROOT.Math.XYZVector(pMag_tau2_perp * tau_dir2_perp.x(), pMag_tau2_perp * tau_dir2_perp.y(), pMag_tau2_perp * tau_dir2_perp.z())

        # Get angles made by tau directions with B+ flight direction
        tau_B_cos1 = (tau_dir1.Unit()).Dot(bDir.Unit())
        tau_B_cos2 = (tau_dir2.Unit()).Dot(bDir.Unit())

        # Get tau momenta parallel to B+ flight direction
        pMag_tau1_long = np.abs(pMag_tau1_perp)*tau_B_cos1/np.sqrt( 1 - tau_B_cos1**2 )
        pMag_tau2_long = np.abs(pMag_tau2_perp)*tau_B_cos2/np.sqrt( 1 - tau_B_cos2**2 )

        # Total tau momentum vector
        ptau1 = p_tau1_perp + ROOT.Math.XYZVector(pMag_tau1_long * bDir.x(), pMag_tau1_long * bDir.y(), pMag_tau1_long * bDir.z())
        ptau2 = p_tau2_perp + ROOT.Math.XYZVector(pMag_tau2_long * bDir.x(), pMag_tau2_long * bDir.y(), pMag_tau2_long * bDir.z())

        # Get the rest of x
        pB = ptau1 + ptau2 + pK
        pnu1 = ptau1 - p3pi1
        pnu2 = ptau2 - p3pi2

        Etau1 = np.sqrt( mtau**2 + ptau1.Mag2() )
        Etau2 = np.sqrt( mtau**2 + ptau2.Mag2() )
        Enu1 = np.sqrt( mnu**2 + pnu1.Mag2() )
        Enu2 = np.sqrt( mnu**2 + pnu2.Mag2() )

        EB = Etau1 + Etau2 + EK
        MB_squared = EB**2 - pB.Mag2()

        x0[dimM] = BV.x()
        x0[dimM+1] = BV.y()
        x0[dimM+2] = BV.z()
        x0[dimM+3] = pB.x()
        x0[dimM+4] = pB.y()
        x0[dimM+5] = pB.z()
        x0[dimM+6] = MB_squared
        x0[dimM+7] = ptau1.x()
        x0[dimM+8] = ptau1.y()
        x0[dimM+9] = ptau1.z()
        x0[dimM+10] = Etau1
        x0[dimM+11] = pnu1.x()
        x0[dimM+12] = pnu1.y()
        x0[dimM+13] = pnu1.z()
        x0[dimM+14] = Enu1
        x0[dimM+15] = ptau2.x()
        x0[dimM+16] = ptau2.y()
        x0[dimM+17] = ptau2.z()
        x0[dimM+18] = Etau2
        x0[dimM+19] = pnu2.x()
        x0[dimM+20] = pnu2.y()
        x0[dimM+21] = pnu2.z()
        x0[dimM+22] = Enu2
    
    elif(init == 3): #  B->K* tautau initialisation taus direction based on visible 3pi momenta
        # Magnitude of 3pi momenta
        p3pi1_mag = p3pi1.r()
        p3pi2_mag = p3pi2.r()

        # B+ flight direction
        bDir = (BV - PV).Unit()

        # Get K+ momentum perpendicular to the B+ flight direction
        pK_perp = pK - ROOT.Math.XYZVector((pK.Dot(bDir)) * bDir.x(), (pK.Dot(bDir)) * bDir.y(), (pK.Dot(bDir)) * bDir.z())

        # Get tau flight directions (from pions visible momenta)
        tau_dir1 = ROOT.Math.XYZVector(0,0,0)
        tau_dir2 = ROOT.Math.XYZVector(0,0,0)
        tau_dir1.SetXYZ( p3pi1.x(), p3pi1.y(), p3pi1.z() )
        tau_dir2.SetXYZ( p3pi2.x(), p3pi2.y(), p3pi2.z() )

        tau_dir1 = tau_dir1.Unit()
        tau_dir2 = tau_dir2.Unit()

        # Get tau direction unit vectors perpendicular to B+ flight direction
        tau_dir1_perp = (tau_dir1 - ROOT.Math.XYZVector((tau_dir1.Dot(bDir)) * bDir.x(), (tau_dir1.Dot(bDir)) * bDir.y(), (tau_dir1.Dot(bDir)) * bDir.z() )).Unit()
        tau_dir2_perp = (tau_dir2 - ROOT.Math.XYZVector((tau_dir2.Dot(bDir)) * bDir.x(), (tau_dir2.Dot(bDir)) * bDir.y(), (tau_dir2.Dot(bDir)) * bDir.z() )).Unit()

        # In plane perpendicular to B+ flight direction, get angles between tau momenta and K+ momentum
        cosphi1 = tau_dir1_perp.Dot(pK_perp.Unit())
        cosphi2 = tau_dir2_perp.Dot(pK_perp.Unit())

        phi1 = np.arccos(cosphi1)
        phi2 = np.arccos(cosphi2)

        # Calculate momentum component of taus in this plane
        pMag_tau1_perp = -1*pK_perp.R()/( cosphi1 + (cosphi2*(np.sin(phi1)/np.sin(phi2))) )
        pMag_tau2_perp = pMag_tau1_perp*(np.sin(phi1) / np.sin(phi2))

        p_tau1_perp = ROOT.Math.XYZVector(pMag_tau1_perp * tau_dir1_perp.x(), pMag_tau1_perp * tau_dir1_perp.y(), pMag_tau1_perp * tau_dir1_perp.z())
        p_tau2_perp = ROOT.Math.XYZVector(pMag_tau2_perp * tau_dir2_perp.x(), pMag_tau2_perp * tau_dir2_perp.y(), pMag_tau2_perp * tau_dir2_perp.z())

        # Get angles made by tau directions with B+ flight direction
        tau_B_cos1 = (tau_dir1.Unit()).Dot(bDir.Unit())
        tau_B_cos2 = (tau_dir2.Unit()).Dot(bDir.Unit())

        # Get tau momenta parallel to B+ flight direction
        pMag_tau1_long = np.abs(pMag_tau1_perp)*tau_B_cos1/np.sqrt( 1 - tau_B_cos1**2 )
        pMag_tau2_long = np.abs(pMag_tau2_perp)*tau_B_cos2/np.sqrt( 1 - tau_B_cos2**2 )

        # Total tau momentum vector
        ptau1 = p_tau1_perp + ROOT.Math.XYZVector(pMag_tau1_long * bDir.x(), pMag_tau1_long * bDir.y(), pMag_tau1_long * bDir.z())
        ptau2 = p_tau2_perp + ROOT.Math.XYZVector(pMag_tau2_long * bDir.x(), pMag_tau2_long * bDir.y(), pMag_tau2_long * bDir.z())

        # Get the rest of x
        pB = ptau1 + ptau2 + pK
        pnu1 = ptau1 - p3pi1
        pnu2 = ptau2 - p3pi2

        Etau1 = np.sqrt( mtau**2 + ptau1.Mag2() )
        Etau2 = np.sqrt( mtau**2 + ptau2.Mag2() )
        Enu1 = np.sqrt( pnu1.Mag2() )
        Enu2 = np.sqrt( mnu**2 + pnu2.Mag2() )

        EB = Etau1 + Etau2 + EK
        MB_squared = EB**2 - pB.Mag2()

        x0[dimM] = BV.x()
        x0[dimM+1] = BV.y()
        x0[dimM+2] = BV.z()
        x0[dimM+3] = pB.x()
        x0[dimM+4] = pB.y()
        x0[dimM+5] = pB.z()
        x0[dimM+6] = MB_squared
        x0[dimM+7] = ptau1.x()
        x0[dimM+8] = ptau1.y()
        x0[dimM+9] = ptau1.z()
        x0[dimM+10] = Etau1
        x0[dimM+11] = pnu1.x()
        x0[dimM+12] = pnu1.y()
        x0[dimM+13] = pnu1.z()
        x0[dimM+14] = Enu1
        x0[dimM+15] = ptau2.x()
        x0[dimM+16] = ptau2.y()
        x0[dimM+17] = ptau2.z()
        x0[dimM+18] = Etau2
        x0[dimM+19] = pnu2.x()
        x0[dimM+20] = pnu2.y()
        x0[dimM+21] = pnu2.z()
        x0[dimM+22] = Enu2

    elif(init == 4): # RD initialisation
        u1 = (DV1 - BV).Unit()
        u2 = (DV2 - BV).Unit()

        # Use the maximum value of theta: neutrino takes the maximum portion of momentum from the tau
        theta1 = np.arcsin( ( np.sqrt( np.abs((mtau**2 + m3pi1**2 - mnu**2)**2 - 4*(mtau**2)*(m3pi1**2) )) )/( 2*mtau*np.sqrt( p3pi1.Mag2() ) ) )
        theta2 = np.arcsin( ( np.sqrt( np.abs((mtau**2 + m3pi2**2 - mnu**2)**2 - 4*(mtau**2)*(m3pi2**2) )) )/( 2*mtau*np.sqrt( p3pi2.Mag2() ) ) )

        ptau1_mag = ( (mtau**2 + m3pi1**2 - mnu**2)*np.sqrt(p3pi1.Mag2())*np.cos(theta1) )/( 2*( E3pi1**2 - p3pi1.Mag2()*np.cos(theta1)**2 ) )
        ptau2_mag = ( (mtau**2 + m3pi2**2 - mnu**2)*np.sqrt(p3pi2.Mag2())*np.cos(theta2) )/( 2*( E3pi2**2 - p3pi2.Mag2()*np.cos(theta2)**2 ) )

        ptau1 = ROOT.Math.XYZVector(ptau1_mag * u1.x(), ptau1_mag * u1.y(), ptau1_mag * u1.z())
        ptau2 = ROOT.Math.XYZVector(ptau2_mag * u2.x(), ptau2_mag * u2.y(), ptau2_mag * u2.z())
        pnu1 = ptau1 - p3pi1
        pnu2 = ptau2 - p3pi2

        Etau1 = np.sqrt( mtau**2 + ptau1.Mag2() )
        Etau2 = np.sqrt( mtau**2 + ptau2.Mag2() )
        Enu1 = np.sqrt( mnu**2 + pnu1.Mag2() )
        Enu2 = np.sqrt( mnu**2 + pnu2.Mag2() )

        pB = pK + ptau1 + ptau2
        EB = EK + Etau1 + Etau2
        MB_squared = EB**2 - pB.Mag2()

        x0[dimM] = BV.x()
        x0[dimM+1] = BV.y()
        x0[dimM+2] = BV.z()
        x0[dimM+3] = pB.x()
        x0[dimM+4] = pB.y()
        x0[dimM+5] = pB.z()
        x0[dimM+6] = MB_squared
        x0[dimM+7] = ptau1.x()
        x0[dimM+8] = ptau1.y()
        x0[dimM+9] = ptau1.z()
        x0[dimM+10] = Etau1
        x0[dimM+11] = pnu1.x()
        x0[dimM+12] = pnu1.y()
        x0[dimM+13] = pnu1.z()
        x0[dimM+14] = Enu1
        x0[dimM+15] = ptau2.x()
        x0[dimM+16] = ptau2.y()
        x0[dimM+17] = ptau2.z()
        x0[dimM+18] = Etau2
        x0[dimM+19] = pnu2.x()
        x0[dimM+20] = pnu2.y()
        x0[dimM+21] = pnu2.z()
        x0[dimM+22] = Enu2

    # elif(init == 5): # B->K* tautau initialisation tau+ direction from vertices  tau- direction from pions
    #     # Magnitude of 3pi momenta
    #     p3pi1_mag = p3pi1.r()
    #     p3pi2_mag = p3pi2.r()

    #     # B+ flight direction
    #     bDir = (BV - PV).Unit()

    #     # Get K+ momentum perpendicular to the B+ flight direction
    #     pK_perp = pK - ROOT.Math.XYZVector((pK.Dot(bDir)) * bDir.x(), (pK.Dot(bDir)) * bDir.y(), (pK.Dot(bDir)) * bDir.z())

    #     # Get tau flight directions (tau+ from vertices, tau- from pions)
    #     tau_dir1 = ROOT.Math.XYZVector(0,0,0)
    #     tau_dir2 = ROOT.Math.XYZVector(0,0,0)
    #     tau_dir1.SetXYZ( DV1.x() - BV.x(), DV1.y() - BV.y(), DV1.z() - BV.z() )
    #     tau_dir2.SetXYZ( p3pi2.x(), p3pi2.y(), p3pi2.z() )

    #     tau_dir1 = tau_dir1.Unit()
    #     tau_dir2 = tau_dir2.Unit()

    #     # Get tau direction unit vectors perpendicular to B+ flight direction
    #     tau_dir1_perp = (tau_dir1 - ROOT.Math.XYZVector((tau_dir1.Dot(bDir)) * bDir.x(), (tau_dir1.Dot(bDir)) * bDir.y(), (tau_dir1.Dot(bDir)) * bDir.z() )).Unit()
    #     tau_dir2_perp = (tau_dir2 - ROOT.Math.XYZVector((tau_dir2.Dot(bDir)) * bDir.x(), (tau_dir2.Dot(bDir)) * bDir.y(), (tau_dir2.Dot(bDir)) * bDir.z() )).Unit()

    #     # In plane perpendicular to B+ flight direction, get angles between tau momenta and K+ momentum
    #     cosphi1 = tau_dir1_perp.Dot(pK_perp.Unit())
    #     cosphi2 = tau_dir2_perp.Dot(pK_perp.Unit())

    #     phi1 = np.arccos(cosphi1)
    #     phi2 = np.arccos(cosphi2)

    #     # Calculate momentum component of taus in this plane
    #     pMag_tau1_perp = -1*pK_perp.R()/( cosphi1 + (cosphi2*(np.sin(phi1)/np.sin(phi2))) )
    #     pMag_tau2_perp = pMag_tau1_perp*(np.sin(phi1) / np.sin(phi2))

    #     p_tau1_perp = ROOT.Math.XYZVector(pMag_tau1_perp * tau_dir1_perp.x(), pMag_tau1_perp * tau_dir1_perp.y(), pMag_tau1_perp * tau_dir1_perp.z())
    #     p_tau2_perp = ROOT.Math.XYZVector(pMag_tau2_perp * tau_dir2_perp.x(), pMag_tau2_perp * tau_dir2_perp.y(), pMag_tau2_perp * tau_dir2_perp.z())

    #     # Get angles made by tau directions with B+ flight direction
    #     tau_B_cos1 = (tau_dir1.Unit()).Dot(bDir.Unit())
    #     tau_B_cos2 = (tau_dir2.Unit()).Dot(bDir.Unit())

    #     # Get tau momenta parallel to B+ flight direction
    #     pMag_tau1_long = np.abs(pMag_tau1_perp)*tau_B_cos1/np.sqrt( 1 - tau_B_cos1**2 )
    #     pMag_tau2_long = np.abs(pMag_tau2_perp)*tau_B_cos2/np.sqrt( 1 - tau_B_cos2**2 )

    #     # Total tau momentum vector
    #     ptau1 = p_tau1_perp + ROOT.Math.XYZVector(pMag_tau1_long * bDir.x(), pMag_tau1_long * bDir.y(), pMag_tau1_long * bDir.z())
    #     ptau2 = p_tau2_perp + ROOT.Math.XYZVector(pMag_tau2_long * bDir.x(), pMag_tau2_long * bDir.y(), pMag_tau2_long * bDir.z())

    #     # Get the rest of x
    #     pB = ptau1 + ptau2 + pK
    #     pnu1 = ptau1 - p3pi1
    #     pnu2 = ptau2 - p3pi2

    #     Etau1 = np.sqrt( mtau**2 + ptau1.Mag2() )
    #     Etau2 = np.sqrt( mtau**2 + ptau2.Mag2() )
    #     Enu1 = np.sqrt( mnu**2 + pnu1.Mag2() )
    #     Enu2 = np.sqrt( mnu**2 + pnu2.Mag2() )

    #     EB = Etau1 + Etau2 + EK
    #     MB_squared = EB**2 - pB.Mag2()

    #     x0[dimM] = BV.x()
    #     x0[dimM+1] = BV.y()
    #     x0[dimM+2] = BV.z()
    #     x0[dimM+3] = pB.x()
    #     x0[dimM+4] = pB.y()
    #     x0[dimM+5] = pB.z()
    #     x0[dimM+6] = MB_squared
    #     x0[dimM+7] = ptau1.x()
    #     x0[dimM+8] = ptau1.y()
    #     x0[dimM+9] = ptau1.z()
    #     x0[dimM+10] = Etau1
    #     x0[dimM+11] = pnu1.x()
    #     x0[dimM+12] = pnu1.y()
    #     x0[dimM+13] = pnu1.z()
    #     x0[dimM+14] = Enu1
    #     x0[dimM+15] = ptau2.x()
    #     x0[dimM+16] = ptau2.y()
    #     x0[dimM+17] = ptau2.z()
    #     x0[dimM+18] = Etau2
    #     x0[dimM+19] = pnu2.x()
    #     x0[dimM+20] = pnu2.y()
    #     x0[dimM+21] = pnu2.z()
    #     x0[dimM+22] = Enu2

    # elif(init == 6): # B->K* tautau initialisation tau- direction from vertices  tau+ direction from pions
    #     # Magnitude of 3pi momenta
    #     p3pi1_mag = p3pi1.r()
    #     p3pi2_mag = p3pi2.r()

    #     # B+ flight direction
    #     bDir = (BV - PV).Unit()

    #     # Get K+ momentum perpendicular to the B+ flight direction
    #     pK_perp = pK - ROOT.Math.XYZVector((pK.Dot(bDir)) * bDir.x(), (pK.Dot(bDir)) * bDir.y(), (pK.Dot(bDir)) * bDir.z())

    #     # Get tau flight directions (tau+ from vertices, tau- from pions)
    #     tau_dir1 = ROOT.Math.XYZVector(0,0,0)
    #     tau_dir2 = ROOT.Math.XYZVector(0,0,0)
    #     tau_dir1.SetXYZ( p3pi1.x(), p3pi1.y(), p3pi1.z() )
    #     tau_dir2.SetXYZ( DV2.x() - BV.x(), DV2.y() - BV.y(), DV2.z() - BV.z() )

    #     tau_dir1 = tau_dir1.Unit()
    #     tau_dir2 = tau_dir2.Unit()

    #     # Get tau direction unit vectors perpendicular to B+ flight direction
    #     tau_dir1_perp = (tau_dir1 - ROOT.Math.XYZVector((tau_dir1.Dot(bDir)) * bDir.x(), (tau_dir1.Dot(bDir)) * bDir.y(), (tau_dir1.Dot(bDir)) * bDir.z() )).Unit()
    #     tau_dir2_perp = (tau_dir2 - ROOT.Math.XYZVector((tau_dir2.Dot(bDir)) * bDir.x(), (tau_dir2.Dot(bDir)) * bDir.y(), (tau_dir2.Dot(bDir)) * bDir.z() )).Unit()

    #     # In plane perpendicular to B+ flight direction, get angles between tau momenta and K+ momentum
    #     cosphi1 = tau_dir1_perp.Dot(pK_perp.Unit())
    #     cosphi2 = tau_dir2_perp.Dot(pK_perp.Unit())

    #     phi1 = np.arccos(cosphi1)
    #     phi2 = np.arccos(cosphi2)

    #     # Calculate momentum component of taus in this plane
    #     pMag_tau1_perp = -1*pK_perp.R()/( cosphi1 + (cosphi2*(np.sin(phi1)/np.sin(phi2))) )
    #     pMag_tau2_perp = pMag_tau1_perp*(np.sin(phi1) / np.sin(phi2))

    #     p_tau1_perp = ROOT.Math.XYZVector(pMag_tau1_perp * tau_dir1_perp.x(), pMag_tau1_perp * tau_dir1_perp.y(), pMag_tau1_perp * tau_dir1_perp.z())
    #     p_tau2_perp = ROOT.Math.XYZVector(pMag_tau2_perp * tau_dir2_perp.x(), pMag_tau2_perp * tau_dir2_perp.y(), pMag_tau2_perp * tau_dir2_perp.z())

    #     # Get angles made by tau directions with B+ flight direction
    #     tau_B_cos1 = (tau_dir1.Unit()).Dot(bDir.Unit())
    #     tau_B_cos2 = (tau_dir2.Unit()).Dot(bDir.Unit())

    #     # Get tau momenta parallel to B+ flight direction
    #     pMag_tau1_long = np.abs(pMag_tau1_perp)*tau_B_cos1/np.sqrt( 1 - tau_B_cos1**2 )
    #     pMag_tau2_long = np.abs(pMag_tau2_perp)*tau_B_cos2/np.sqrt( 1 - tau_B_cos2**2 )

    #     # Total tau momentum vector
    #     ptau1 = p_tau1_perp + ROOT.Math.XYZVector(pMag_tau1_long * bDir.x(), pMag_tau1_long * bDir.y(), pMag_tau1_long * bDir.z())
    #     ptau2 = p_tau2_perp + ROOT.Math.XYZVector(pMag_tau2_long * bDir.x(), pMag_tau2_long * bDir.y(), pMag_tau2_long * bDir.z())

    #     # Get the rest of x
    #     pB = ptau1 + ptau2 + pK
    #     pnu1 = ptau1 - p3pi1
    #     pnu2 = ptau2 - p3pi2

    #     Etau1 = np.sqrt( mtau**2 + ptau1.Mag2() )
    #     Etau2 = np.sqrt( mtau**2 + ptau2.Mag2() )
    #     Enu1 = np.sqrt( mnu**2 + pnu1.Mag2() )
    #     Enu2 = np.sqrt( mnu**2 + pnu2.Mag2() )

    #     EB = Etau1 + Etau2 + EK
    #     MB_squared = EB**2 - pB.Mag2()

    #     x0[dimM] = BV.x()
    #     x0[dimM+1] = BV.y()
    #     x0[dimM+2] = BV.z()
    #     x0[dimM+3] = pB.x()
    #     x0[dimM+4] = pB.y()
    #     x0[dimM+5] = pB.z()
    #     x0[dimM+6] = MB_squared
    #     x0[dimM+7] = ptau1.x()
    #     x0[dimM+8] = ptau1.y()
    #     x0[dimM+9] = ptau1.z()
    #     x0[dimM+10] = Etau1
    #     x0[dimM+11] = pnu1.x()
    #     x0[dimM+12] = pnu1.y()
    #     x0[dimM+13] = pnu1.z()
    #     x0[dimM+14] = Enu1
    #     x0[dimM+15] = ptau2.x()
    #     x0[dimM+16] = ptau2.y()
    #     x0[dimM+17] = ptau2.z()
    #     x0[dimM+18] = Etau2
    #     x0[dimM+19] = pnu2.x()
    #     x0[dimM+20] = pnu2.y()
    #     x0[dimM+21] = pnu2.z()
    #     x0[dimM+22] = Enu2

    # elif(init == 7): # Mixes Marseille (tau+) and K* tau tau vertices (tau-)
    #     u1 = (DV1 - BV).Unit()
    #     u2 = (DV2 - BV).Unit()

    #     theta1 = np.arcsin( ( np.sqrt( (mtau**2 + m3pi1**2 - mnu**2)**2 - 4*(mtau**2)*(m3pi1**2) ) )/( 2*mtau*np.sqrt( p3pi1.Mag2() ) ) )

    #     ptau1_mag = ( (mtau**2 + m3pi1**2 - mnu**2)*np.sqrt(p3pi1.Mag2())*np.cos(theta1) )/( 2*( E3pi1**2 - p3pi1.Mag2()*np.cos(theta1)**2 ) )

    #     ptau1 = ROOT.Math.XYZVector(ptau1_mag * u1.x(), ptau1_mag * u1.y(), ptau1_mag * u1.z())

    #     bDir = (BV - PV).Unit()

    #     pK_perp = pK - ROOT.Math.XYZVector((pK.Dot(bDir)) * bDir.x(), (pK.Dot(bDir)) * bDir.y(), (pK.Dot(bDir)) * bDir.z())

    #     tau_dir1 = ROOT.Math.XYZVector(0,0,0)
    #     tau_dir2 = ROOT.Math.XYZVector(0,0,0)
    #     tau_dir1.SetXYZ( DV1.x() - BV.x(), DV1.y() - BV.y(), DV1.z() - BV.z() )
    #     tau_dir2.SetXYZ( DV2.x() - BV.x(), DV2.y() - BV.y(), DV2.z() - BV.z() )

    #     tau_dir1 = tau_dir1.Unit()
    #     tau_dir2 = tau_dir2.Unit()

    #     tau_dir1_perp = (tau_dir1 - ROOT.Math.XYZVector((tau_dir1.Dot(bDir)) * bDir.x(), (tau_dir1.Dot(bDir)) * bDir.y(), (tau_dir1.Dot(bDir)) * bDir.z() )).Unit()
    #     tau_dir2_perp = (tau_dir2 - ROOT.Math.XYZVector((tau_dir2.Dot(bDir)) * bDir.x(), (tau_dir2.Dot(bDir)) * bDir.y(), (tau_dir2.Dot(bDir)) * bDir.z() )).Unit()

    #     cosphi1 = tau_dir1_perp.Dot(pK_perp.Unit())
    #     cosphi2 = tau_dir2_perp.Dot(pK_perp.Unit())

    #     phi1 = np.arccos(cosphi1)
    #     phi2 = np.arccos(cosphi2)

    #     pMag_tau1_perp = -1*pK_perp.R()/( cosphi1 + (cosphi2*(np.sin(phi1)/np.sin(phi2))) )
    #     pMag_tau2_perp = pMag_tau1_perp*(np.sin(phi1) / np.sin(phi2))

    #     p_tau1_perp = ROOT.Math.XYZVector(pMag_tau1_perp * tau_dir1_perp.x(), pMag_tau1_perp * tau_dir1_perp.y(), pMag_tau1_perp * tau_dir1_perp.z())
    #     p_tau2_perp = ROOT.Math.XYZVector(pMag_tau2_perp * tau_dir2_perp.x(), pMag_tau2_perp * tau_dir2_perp.y(), pMag_tau2_perp * tau_dir2_perp.z())

    #     tau_B_cos1 = (tau_dir1.Unit()).Dot(bDir.Unit())
    #     tau_B_cos2 = (tau_dir2.Unit()).Dot(bDir.Unit())

    #     pMag_tau1_long = np.abs(pMag_tau1_perp)*tau_B_cos1/np.sqrt( 1 - tau_B_cos1**2 )
    #     pMag_tau2_long = np.abs(pMag_tau2_perp)*tau_B_cos2/np.sqrt( 1 - tau_B_cos2**2 )

    #     ptau2 = p_tau2_perp + ROOT.Math.XYZVector(pMag_tau2_long * bDir.x(), pMag_tau2_long * bDir.y(), pMag_tau2_long * bDir.z())

    #     pnu1 = ptau1 - p3pi1
    #     pnu2 = ptau2 - p3pi2

    #     Etau1 = np.sqrt( mtau**2 + ptau1.Mag2() )
    #     Etau2 = np.sqrt( mtau**2 + ptau2.Mag2() )
    #     Enu1 = np.sqrt( mnu**2 + pnu1.Mag2() )
    #     Enu2 = np.sqrt( mnu**2 + pnu2.Mag2() )

    #     pB = pK + ptau1 + ptau2
    #     EB = EK + Etau1 + Etau2
    #     MB_squared = EB**2 - pB.Mag2()

    #     x0[dimM] = BV.x()
    #     x0[dimM+1] = BV.y()
    #     x0[dimM+2] = BV.z()
    #     x0[dimM+3] = pB.x()
    #     x0[dimM+4] = pB.y()
    #     x0[dimM+5] = pB.z()
    #     x0[dimM+6] = MB_squared
    #     x0[dimM+7] = ptau1.x()
    #     x0[dimM+8] = ptau1.y()
    #     x0[dimM+9] = ptau1.z()
    #     x0[dimM+10] = Etau1
    #     x0[dimM+11] = pnu1.x()
    #     x0[dimM+12] = pnu1.y()
    #     x0[dimM+13] = pnu1.z()
    #     x0[dimM+14] = Enu1
    #     x0[dimM+15] = ptau2.x()
    #     x0[dimM+16] = ptau2.y()
    #     x0[dimM+17] = ptau2.z()
    #     x0[dimM+18] = Etau2
    #     x0[dimM+19] = pnu2.x()
    #     x0[dimM+20] = pnu2.y()
    #     x0[dimM+21] = pnu2.z()
    #     x0[dimM+22] = Enu2

    # elif(init == 8): # Mixes Marseille (tau-) and K* tau tau vertices (tau+)
    #     u1 = (DV1 - BV).Unit()
    #     u2 = (DV2 - BV).Unit()

    #     m3pi1 = np.sqrt( E3pi1**2 - p3pi1.Mag2() )
    #     m3pi2 = np.sqrt( E3pi2**2 - p3pi2.Mag2() )

    #     theta2 = np.arcsin( ( np.sqrt( (mtau**2 + m3pi2**2 - mnu**2)**2 - 4*(mtau**2)*(m3pi2**2) ) )/( 2*mtau*np.sqrt( p3pi2.Mag2() ) ) )

    #     ptau2_mag = ( (mtau**2 + m3pi2**2 - mnu**2)*np.sqrt(p3pi2.Mag2())*np.cos(theta2) )/( 2*( E3pi2**2 - p3pi2.Mag2()*np.cos(theta2)**2 ) )

    #     ptau2 = ROOT.Math.XYZVector(ptau2_mag * u2.x(), ptau2_mag * u2.y(), ptau2_mag * u2.z())

    #     bDir = (BV - PV).Unit()

    #     pK_perp = pK - ROOT.Math.XYZVector((pK.Dot(bDir)) * bDir.x(), (pK.Dot(bDir)) * bDir.y(), (pK.Dot(bDir)) * bDir.z())

    #     tau_dir1 = ROOT.Math.XYZVector(0,0,0)
    #     tau_dir2 = ROOT.Math.XYZVector(0,0,0)
    #     tau_dir1.SetXYZ( DV1.x() - BV.x(), DV1.y() - BV.y(), DV1.z() - BV.z() )
    #     tau_dir2.SetXYZ( DV2.x() - BV.x(), DV2.y() - BV.y(), DV2.z() - BV.z() )

    #     tau_dir1 = tau_dir1.Unit()
    #     tau_dir2 = tau_dir2.Unit()

    #     tau_dir1_perp = (tau_dir1 - ROOT.Math.XYZVector((tau_dir1.Dot(bDir)) * bDir.x(), (tau_dir1.Dot(bDir)) * bDir.y(), (tau_dir1.Dot(bDir)) * bDir.z() )).Unit()
    #     tau_dir2_perp = (tau_dir2 - ROOT.Math.XYZVector((tau_dir2.Dot(bDir)) * bDir.x(), (tau_dir2.Dot(bDir)) * bDir.y(), (tau_dir2.Dot(bDir)) * bDir.z() )).Unit()

    #     cosphi1 = tau_dir1_perp.Dot(pK_perp.Unit())
    #     cosphi2 = tau_dir2_perp.Dot(pK_perp.Unit())

    #     phi1 = np.arccos(cosphi1)
    #     phi2 = np.arccos(cosphi2)

    #     pMag_tau1_perp = -1*pK_perp.R()/( cosphi1 + (cosphi2*(np.sin(phi1)/np.sin(phi2))) )
    #     pMag_tau2_perp = pMag_tau1_perp*(np.sin(phi1) / np.sin(phi2))

    #     p_tau1_perp = ROOT.Math.XYZVector(pMag_tau1_perp * tau_dir1_perp.x(), pMag_tau1_perp * tau_dir1_perp.y(), pMag_tau1_perp * tau_dir1_perp.z())
    #     p_tau2_perp = ROOT.Math.XYZVector(pMag_tau2_perp * tau_dir2_perp.x(), pMag_tau2_perp * tau_dir2_perp.y(), pMag_tau2_perp * tau_dir2_perp.z())

    #     tau_B_cos1 = (tau_dir1.Unit()).Dot(bDir.Unit())
    #     tau_B_cos2 = (tau_dir2.Unit()).Dot(bDir.Unit())

    #     pMag_tau1_long = np.abs(pMag_tau1_perp)*tau_B_cos1/np.sqrt( 1 - tau_B_cos1**2 )
    #     pMag_tau2_long = np.abs(pMag_tau2_perp)*tau_B_cos2/np.sqrt( 1 - tau_B_cos2**2 )

    #     ptau1 = p_tau1_perp + ROOT.Math.XYZVector(pMag_tau1_long * bDir.x(), pMag_tau1_long * bDir.y(), pMag_tau1_long * bDir.z())

    #     pnu1 = ptau1 - p3pi1
    #     pnu2 = ptau2 - p3pi2

    #     Etau1 = np.sqrt( mtau**2 + ptau1.Mag2() )
    #     Etau2 = np.sqrt( mtau**2 + ptau2.Mag2() )
    #     Enu1 = np.sqrt( mnu**2 + pnu1.Mag2() )
    #     Enu2 = np.sqrt( mnu**2 + pnu2.Mag2() )

    #     pB = pK + ptau1 + ptau2
    #     EB = EK + Etau1 + Etau2
    #     MB_squared = EB**2 - pB.Mag2()

    #     x0[dimM] = BV.x()
    #     x0[dimM+1] = BV.y()
    #     x0[dimM+2] = BV.z()
    #     x0[dimM+3] = pB.x()
    #     x0[dimM+4] = pB.y()
    #     x0[dimM+5] = pB.z()
    #     x0[dimM+6] = MB_squared
    #     x0[dimM+7] = ptau1.x()
    #     x0[dimM+8] = ptau1.y()
    #     x0[dimM+9] = ptau1.z()
    #     x0[dimM+10] = Etau1
    #     x0[dimM+11] = pnu1.x()
    #     x0[dimM+12] = pnu1.y()
    #     x0[dimM+13] = pnu1.z()
    #     x0[dimM+14] = Enu1
    #     x0[dimM+15] = ptau2.x()
    #     x0[dimM+16] = ptau2.y()
    #     x0[dimM+17] = ptau2.z()
    #     x0[dimM+18] = Etau2
    #     x0[dimM+19] = pnu2.x()
    #     x0[dimM+20] = pnu2.y()
    #     x0[dimM+21] = pnu2.z()
    #     x0[dimM+22] = Enu2
    
    # elif(init == 9): # Mixes Marseille (tau+) and K* tau tau pions (tau-)
    #     u1 = (DV1 - BV).Unit()
    #     u2 = (DV2 - BV).Unit()

    #     theta1 = np.arcsin( ( np.sqrt( (mtau**2 + m3pi1**2 - mnu**2)**2 - 4*(mtau**2)*(m3pi1**2) ) )/( 2*mtau*np.sqrt( p3pi1.Mag2() ) ) )

    #     ptau1_mag = ( (mtau**2 + m3pi1**2 - mnu**2)*np.sqrt(p3pi1.Mag2())*np.cos(theta1) )/( 2*( E3pi1**2 - p3pi1.Mag2()*np.cos(theta1)**2 ) )

    #     ptau1 = ROOT.Math.XYZVector(ptau1_mag * u1.x(), ptau1_mag * u1.y(), ptau1_mag * u1.z())

    #     bDir = (BV - PV).Unit()

    #     pK_perp = pK - ROOT.Math.XYZVector((pK.Dot(bDir)) * bDir.x(), (pK.Dot(bDir)) * bDir.y(), (pK.Dot(bDir)) * bDir.z())

    #     tau_dir1 = ROOT.Math.XYZVector(0,0,0)
    #     tau_dir2 = ROOT.Math.XYZVector(0,0,0)
    #     tau_dir1.SetXYZ( p3pi1.x(), p3pi1.y(), p3pi1.z() )
    #     tau_dir2.SetXYZ( p3pi2.x(), p3pi2.y(), p3pi2.z() )

    #     tau_dir1 = tau_dir1.Unit()
    #     tau_dir2 = tau_dir2.Unit()

    #     tau_dir1_perp = (tau_dir1 - ROOT.Math.XYZVector((tau_dir1.Dot(bDir)) * bDir.x(), (tau_dir1.Dot(bDir)) * bDir.y(), (tau_dir1.Dot(bDir)) * bDir.z() )).Unit()
    #     tau_dir2_perp = (tau_dir2 - ROOT.Math.XYZVector((tau_dir2.Dot(bDir)) * bDir.x(), (tau_dir2.Dot(bDir)) * bDir.y(), (tau_dir2.Dot(bDir)) * bDir.z() )).Unit()

    #     cosphi1 = tau_dir1_perp.Dot(pK_perp.Unit())
    #     cosphi2 = tau_dir2_perp.Dot(pK_perp.Unit())

    #     phi1 = np.arccos(cosphi1)
    #     phi2 = np.arccos(cosphi2)

    #     pMag_tau1_perp = -1*pK_perp.R()/( cosphi1 + (cosphi2*(np.sin(phi1)/np.sin(phi2))) )
    #     pMag_tau2_perp = pMag_tau1_perp*(np.sin(phi1) / np.sin(phi2))

    #     p_tau1_perp = ROOT.Math.XYZVector(pMag_tau1_perp * tau_dir1_perp.x(), pMag_tau1_perp * tau_dir1_perp.y(), pMag_tau1_perp * tau_dir1_perp.z())
    #     p_tau2_perp = ROOT.Math.XYZVector(pMag_tau2_perp * tau_dir2_perp.x(), pMag_tau2_perp * tau_dir2_perp.y(), pMag_tau2_perp * tau_dir2_perp.z())

    #     tau_B_cos1 = (tau_dir1.Unit()).Dot(bDir.Unit())
    #     tau_B_cos2 = (tau_dir2.Unit()).Dot(bDir.Unit())

    #     pMag_tau1_long = np.abs(pMag_tau1_perp)*tau_B_cos1/np.sqrt( 1 - tau_B_cos1**2 )
    #     pMag_tau2_long = np.abs(pMag_tau2_perp)*tau_B_cos2/np.sqrt( 1 - tau_B_cos2**2 )

    #     ptau2 = p_tau2_perp + ROOT.Math.XYZVector(pMag_tau2_long * bDir.x(), pMag_tau2_long * bDir.y(), pMag_tau2_long * bDir.z())

    #     pnu1 = ptau1 - p3pi1
    #     pnu2 = ptau2 - p3pi2

    #     Etau1 = np.sqrt( mtau**2 + ptau1.Mag2() )
    #     Etau2 = np.sqrt( mtau**2 + ptau2.Mag2() )
    #     Enu1 = np.sqrt( mnu**2 + pnu1.Mag2() )
    #     Enu2 = np.sqrt( mnu**2 + pnu2.Mag2() )

    #     pB = pK + ptau1 + ptau2
    #     EB = EK + Etau1 + Etau2
    #     MB_squared = EB**2 - pB.Mag2()

    #     x0[dimM] = BV.x()
    #     x0[dimM+1] = BV.y()
    #     x0[dimM+2] = BV.z()
    #     x0[dimM+3] = pB.x()
    #     x0[dimM+4] = pB.y()
    #     x0[dimM+5] = pB.z()
    #     x0[dimM+6] = MB_squared
    #     x0[dimM+7] = ptau1.x()
    #     x0[dimM+8] = ptau1.y()
    #     x0[dimM+9] = ptau1.z()
    #     x0[dimM+10] = Etau1
    #     x0[dimM+11] = pnu1.x()
    #     x0[dimM+12] = pnu1.y()
    #     x0[dimM+13] = pnu1.z()
    #     x0[dimM+14] = Enu1
    #     x0[dimM+15] = ptau2.x()
    #     x0[dimM+16] = ptau2.y()
    #     x0[dimM+17] = ptau2.z()
    #     x0[dimM+18] = Etau2
    #     x0[dimM+19] = pnu2.x()
    #     x0[dimM+20] = pnu2.y()
    #     x0[dimM+21] = pnu2.z()
    #     x0[dimM+22] = Enu2

    # elif(init == 10): # Mixes Marseille (tau-) and K* tau tau pions (tau+)
    #     u1 = (DV1 - BV).Unit()
    #     u2 = (DV2 - BV).Unit()

    #     theta2 = np.arcsin( ( np.sqrt( (mtau**2 + m3pi2**2 - mnu**2)**2 - 4*(mtau**2)*(m3pi2**2) ) )/( 2*mtau*np.sqrt( p3pi2.Mag2() ) ) )

    #     ptau2_mag = ( (mtau**2 + m3pi2**2 - mnu**2)*np.sqrt(p3pi2.Mag2())*np.cos(theta2) )/( 2*( E3pi2*2 - p3pi2.Mag2()*np.cos(theta2)**2 ) )

    #     ptau2 = ROOT.Math.XYZVector(ptau2_mag * u2.x(), ptau2_mag * u2.y(), ptau2_mag * u2.z())

    #     bDir = (BV - PV).Unit()

    #     pK_perp = pK - ROOT.Math.XYZVector((pK.Dot(bDir)) * bDir.x(), (pK.Dot(bDir)) * bDir.y(), (pK.Dot(bDir)) * bDir.z())

    #     tau_dir1 = ROOT.Math.XYZVector(0,0,0)
    #     tau_dir2 = ROOT.Math.XYZVector(0,0,0)
    #     tau_dir1.SetXYZ( p3pi1.x(), p3pi1.y(), p3pi1.z() )
    #     tau_dir2.SetXYZ( p3pi2.x(), p3pi2.y(), p3pi2.z() )

    #     tau_dir1 = tau_dir1.Unit()
    #     tau_dir2 = tau_dir2.Unit()

    #     tau_dir1_perp = (tau_dir1 - ROOT.Math.XYZVector((tau_dir1.Dot(bDir)) * bDir.x(), (tau_dir1.Dot(bDir)) * bDir.y(), (tau_dir1.Dot(bDir)) * bDir.z() )).Unit()
    #     tau_dir2_perp = (tau_dir2 - ROOT.Math.XYZVector((tau_dir2.Dot(bDir)) * bDir.x(), (tau_dir2.Dot(bDir)) * bDir.y(), (tau_dir2.Dot(bDir)) * bDir.z() )).Unit()

    #     cosphi1 = tau_dir1_perp.Dot(pK_perp.Unit())
    #     cosphi2 = tau_dir2_perp.Dot(pK_perp.Unit())

    #     phi1 = np.arccos(cosphi1)
    #     phi2 = np.arccos(cosphi2)

    #     pMag_tau1_perp = -1*pK_perp.R()/( cosphi1 + (cosphi2*(np.sin(phi1)/np.sin(phi2))) )
    #     pMag_tau2_perp = pMag_tau1_perp*(np.sin(phi1) / np.sin(phi2))

    #     p_tau1_perp = ROOT.Math.XYZVector(pMag_tau1_perp * tau_dir1_perp.x(), pMag_tau1_perp * tau_dir1_perp.y(), pMag_tau1_perp * tau_dir1_perp.z())
    #     p_tau2_perp = ROOT.Math.XYZVector(pMag_tau2_perp * tau_dir2_perp.x(), pMag_tau2_perp * tau_dir2_perp.y(), pMag_tau2_perp * tau_dir2_perp.z())

    #     tau_B_cos1 = (tau_dir1.Unit()).Dot(bDir.Unit())
    #     tau_B_cos2 = (tau_dir2.Unit()).Dot(bDir.Unit())

    #     pMag_tau1_long = np.abs(pMag_tau1_perp)*tau_B_cos1/np.sqrt( 1 - tau_B_cos1**2 )
    #     pMag_tau2_long = np.abs(pMag_tau2_perp)*tau_B_cos2/np.sqrt( 1 - tau_B_cos2**2 )

    #     ptau1 = p_tau1_perp + ROOT.Math.XYZVector(pMag_tau1_long * bDir.x(), pMag_tau1_long * bDir.y(), pMag_tau1_long * bDir.z())

    #     pnu1 = ptau1 - p3pi1
    #     pnu2 = ptau2 - p3pi2

    #     Etau1 = np.sqrt( mtau**2 + ptau1.Mag2() )
    #     Etau2 = np.sqrt( mtau**2 + ptau2.Mag2() )
    #     Enu1 = np.sqrt( mnu**2 + pnu1.Mag2() )
    #     Enu2 = np.sqrt( mnu**2 + pnu2.Mag2() )

    #     pB = pK + ptau1 + ptau2
    #     EB = EK + Etau1 + Etau2
    #     MB_squared = EB**2 - pB.Mag2()

    #     x0[dimM] = BV.x()
    #     x0[dimM+1] = BV.y()
    #     x0[dimM+2] = BV.z()
    #     x0[dimM+3] = pB.x()
    #     x0[dimM+4] = pB.y()
    #     x0[dimM+5] = pB.z()
    #     x0[dimM+6] = MB_squared
    #     x0[dimM+7] = ptau1.x()
    #     x0[dimM+8] = ptau1.y()
    #     x0[dimM+9] = ptau1.z()
    #     x0[dimM+10] = Etau1
    #     x0[dimM+11] = pnu1.x()
    #     x0[dimM+12] = pnu1.y()
    #     x0[dimM+13] = pnu1.z()
    #     x0[dimM+14] = Enu1
    #     x0[dimM+15] = ptau2.x()
    #     x0[dimM+16] = ptau2.y()
    #     x0[dimM+17] = ptau2.z()
    #     x0[dimM+18] = Etau2
    #     x0[dimM+19] = pnu2.x()
    #     x0[dimM+20] = pnu2.y()
    #     x0[dimM+21] = pnu2.z()
    #     x0[dimM+22] = Enu2

    # elif(init == 11): # MLP, K*tautau pions
    #     ptau1 = ROOT.Math.XYZVector( MLP_taup_PX, MLP_taup_PY, MLP_taup_PZ )        

    #     # Magnitude of 3pi momenta
    #     p3pi2_mag = p3pi2.r()

    #     # B+ flight direction
    #     bDir = (BV - PV).Unit()

    #     # Get K+ momentum perpendicular to the B+ flight direction
    #     pK_perp = pK - ROOT.Math.XYZVector((pK.Dot(bDir)) * bDir.x(), (pK.Dot(bDir)) * bDir.y(), (pK.Dot(bDir)) * bDir.z())

    #     # Get tau flight directions (from pions visible momenta)
    #     tau_dir1 = ROOT.Math.XYZVector(0,0,0)
    #     tau_dir2 = ROOT.Math.XYZVector(0,0,0)
    #     tau_dir1.SetXYZ( p3pi1.x(), p3pi1.y(), p3pi1.z() )
    #     tau_dir2.SetXYZ( p3pi2.x(), p3pi2.y(), p3pi2.z() )

    #     tau_dir1 = tau_dir1.Unit()
    #     tau_dir2 = tau_dir2.Unit()

    #     # Get tau direction unit vectors perpendicular to B+ flight direction
    #     tau_dir1_perp = (tau_dir1 - ROOT.Math.XYZVector((tau_dir1.Dot(bDir)) * bDir.x(), (tau_dir1.Dot(bDir)) * bDir.y(), (tau_dir1.Dot(bDir)) * bDir.z() )).Unit()
    #     tau_dir2_perp = (tau_dir2 - ROOT.Math.XYZVector((tau_dir2.Dot(bDir)) * bDir.x(), (tau_dir2.Dot(bDir)) * bDir.y(), (tau_dir2.Dot(bDir)) * bDir.z() )).Unit()

    #     # In plane perpendicular to B+ flight direction, get angles between tau momenta and K+ momentum
    #     cosphi1 = tau_dir1_perp.Dot(pK_perp.Unit())
    #     cosphi2 = tau_dir2_perp.Dot(pK_perp.Unit())

    #     phi1 = np.arccos(cosphi1)
    #     phi2 = np.arccos(cosphi2)

    #     # Calculate momentum component of taus in this plane
    #     pMag_tau1_perp = -1*pK_perp.R()/( cosphi1 + (cosphi2*(np.sin(phi1)/np.sin(phi2))) )
    #     pMag_tau2_perp = pMag_tau1_perp*(np.sin(phi1) / np.sin(phi2))

    #     p_tau2_perp = ROOT.Math.XYZVector(pMag_tau2_perp * tau_dir2_perp.x(), pMag_tau2_perp * tau_dir2_perp.y(), pMag_tau2_perp * tau_dir2_perp.z())

    #     # Get angles made by tau directions with B+ flight direction
    #     tau_B_cos2 = (tau_dir2.Unit()).Dot(bDir.Unit())

    #     # Get tau momenta parallel to B+ flight direction
    #     pMag_tau2_long = np.abs(pMag_tau2_perp)*tau_B_cos2/np.sqrt( 1 - tau_B_cos2**2 )

    #     # Total tau momentum vector
    #     ptau2 = p_tau2_perp + ROOT.Math.XYZVector(pMag_tau2_long * bDir.x(), pMag_tau2_long * bDir.y(), pMag_tau2_long * bDir.z())

    #     pnu1 = ptau1 - p3pi1
    #     pnu2 = ptau2 - p3pi2

    #     Etau1 = np.sqrt( mtau**2 + ptau1.Mag2() )
    #     Etau2 = np.sqrt( mtau**2 + ptau2.Mag2() )

    #     Enu1 = np.sqrt( pnu1.Mag2() )
    #     Enu2 = np.sqrt( mnu**2 + pnu2.Mag2() )

    #     pK = ROOT.Math.XYZVector( m[19], m[20], m[21] )
    #     EK = np.sqrt( mkaon**2 + pK.Mag2() )
    #     pB = ptau1 + ptau2 + pK
    #     EB = Etau1 + Etau2 + EK
    #     MB_squared = EB**2 - pB.Mag2()

    #     x0[dimM] = BV.x()
    #     x0[dimM+1] = BV.y()
    #     x0[dimM+2] = BV.z()
    #     x0[dimM+3] = pB.x()
    #     x0[dimM+4] = pB.y()
    #     x0[dimM+5] = pB.z()
    #     x0[dimM+6] = MB_squared
    #     x0[dimM+7] = ptau1.x()
    #     x0[dimM+8] = ptau1.y()
    #     x0[dimM+9] = ptau1.z()
    #     x0[dimM+10] = Etau1
    #     x0[dimM+11] = pnu1.x()
    #     x0[dimM+12] = pnu1.y()
    #     x0[dimM+13] = pnu1.z()
    #     x0[dimM+14] = Enu1
    #     x0[dimM+15] = ptau2.x()
    #     x0[dimM+16] = ptau2.y()
    #     x0[dimM+17] = ptau2.z()
    #     x0[dimM+18] = Etau2
    #     x0[dimM+19] = pnu2.x()
    #     x0[dimM+20] = pnu2.y()
    #     x0[dimM+21] = pnu2.z()
    #     x0[dimM+22] = Enu2

    # elif(init == 12): # K*tautau (pions), MLP
    #     ptau2 = ROOT.Math.XYZVector( MLP_taum_PX, MLP_taum_PY, MLP_taum_PZ )

    #     # Magnitude of 3pi momenta
    #     p3pi1_mag = p3pi1.r()

    #     # B+ flight direction
    #     bDir = (BV - PV).Unit()

    #     # Get K+ momentum perpendicular to the B+ flight direction
    #     pK_perp = pK - ROOT.Math.XYZVector((pK.Dot(bDir)) * bDir.x(), (pK.Dot(bDir)) * bDir.y(), (pK.Dot(bDir)) * bDir.z())

    #     # Get tau flight directions (from pions visible momenta)
    #     tau_dir1 = ROOT.Math.XYZVector(0,0,0)
    #     tau_dir2 = ROOT.Math.XYZVector(0,0,0)
    #     tau_dir1.SetXYZ( p3pi1.x(), p3pi1.y(), p3pi1.z() )
    #     tau_dir2.SetXYZ( p3pi2.x(), p3pi2.y(), p3pi2.z() )

    #     tau_dir1 = tau_dir1.Unit()
    #     tau_dir2 = tau_dir2.Unit()

    #     # Get tau direction unit vectors perpendicular to B+ flight direction
    #     tau_dir1_perp = (tau_dir1 - ROOT.Math.XYZVector((tau_dir1.Dot(bDir)) * bDir.x(), (tau_dir1.Dot(bDir)) * bDir.y(), (tau_dir1.Dot(bDir)) * bDir.z() )).Unit()
    #     tau_dir2_perp = (tau_dir2 - ROOT.Math.XYZVector((tau_dir2.Dot(bDir)) * bDir.x(), (tau_dir2.Dot(bDir)) * bDir.y(), (tau_dir2.Dot(bDir)) * bDir.z() )).Unit()

    #     # In plane perpendicular to B+ flight direction, get angles between tau momenta and K+ momentum
    #     cosphi1 = tau_dir1_perp.Dot(pK_perp.Unit())
    #     cosphi2 = tau_dir2_perp.Dot(pK_perp.Unit())

    #     phi1 = np.arccos(cosphi1)
    #     phi2 = np.arccos(cosphi2)

    #     # Calculate momentum component of taus in this plane
    #     pMag_tau1_perp = -1*pK_perp.R()/( cosphi1 + (cosphi2*(np.sin(phi1)/np.sin(phi2))) )
    #     pMag_tau2_perp = pMag_tau1_perp*(np.sin(phi1) / np.sin(phi2))

    #     p_tau1_perp = ROOT.Math.XYZVector(pMag_tau1_perp * tau_dir1_perp.x(), pMag_tau1_perp * tau_dir1_perp.y(), pMag_tau1_perp * tau_dir1_perp.z())

    #     # Get angles made by tau directions with B+ flight direction
    #     tau_B_cos1 = (tau_dir1.Unit()).Dot(bDir.Unit())

    #     # Get tau momenta parallel to B+ flight direction
    #     pMag_tau1_long = np.abs(pMag_tau1_perp)*tau_B_cos1/np.sqrt( 1 - tau_B_cos1**2 )

    #     # Total tau momentum vector
    #     ptau1 = p_tau1_perp + ROOT.Math.XYZVector(pMag_tau1_long * bDir.x(), pMag_tau1_long * bDir.y(), pMag_tau1_long * bDir.z())

    #     pnu1 = ptau1 - p3pi1
    #     pnu2 = ptau2 - p3pi2

    #     Etau1 = np.sqrt( mtau**2 + ptau1.Mag2() )
    #     Etau2 = np.sqrt( mtau**2 + ptau2.Mag2() )

    #     Enu1 = np.sqrt( pnu1.Mag2() )
    #     Enu2 = np.sqrt( mnu**2 + pnu2.Mag2() )

    #     pK = ROOT.Math.XYZVector( m[19], m[20], m[21] )
    #     EK = np.sqrt( mkaon**2 + pK.Mag2() )
    #     pB = ptau1 + ptau2 + pK
    #     EB = Etau1 + Etau2 + EK
    #     MB_squared = EB**2 - pB.Mag2()

    #     x0[dimM] = BV.x()
    #     x0[dimM+1] = BV.y()
    #     x0[dimM+2] = BV.z()
    #     x0[dimM+3] = pB.x()
    #     x0[dimM+4] = pB.y()
    #     x0[dimM+5] = pB.z()
    #     x0[dimM+6] = MB_squared
    #     x0[dimM+7] = ptau1.x()
    #     x0[dimM+8] = ptau1.y()
    #     x0[dimM+9] = ptau1.z()
    #     x0[dimM+10] = Etau1
    #     x0[dimM+11] = pnu1.x()
    #     x0[dimM+12] = pnu1.y()
    #     x0[dimM+13] = pnu1.z()
    #     x0[dimM+14] = Enu1
    #     x0[dimM+15] = ptau2.x()
    #     x0[dimM+16] = ptau2.y()
    #     x0[dimM+17] = ptau2.z()
    #     x0[dimM+18] = Etau2
    #     x0[dimM+19] = pnu2.x()
    #     x0[dimM+20] = pnu2.y()
    #     x0[dimM+21] = pnu2.z()
    #     x0[dimM+22] = Enu2

    # elif(init == 13):
    #     ptau1 = ROOT.Math.XYZVector( MLP_taup_PX, MLP_taup_PY, MLP_taup_PZ )

    #     u2 = (DV2 - BV).Unit()
    #     theta2 = np.arcsin( ( np.sqrt( np.abs((mtau**2 + m3pi2**2 - mnu**2)**2 - 4*(mtau**2)*(m3pi2**2) )) )/( 2*mtau*np.sqrt( p3pi2.Mag2() ) ) )
    #     ptau2_mag = ( (mtau**2 + m3pi2**2 - mnu**2)*np.sqrt(p3pi2.Mag2())*np.cos(theta2) )/( 2*( E3pi2**2 - p3pi2.Mag2()*np.cos(theta2)**2 ) )
    #     ptau2 = ROOT.Math.XYZVector(ptau2_mag * u2.x(), ptau2_mag * u2.y(), ptau2_mag * u2.z())

    #     pnu1 = ptau1 - p3pi1
    #     pnu2 = ptau2 - p3pi2

    #     Etau1 = np.sqrt( mtau**2 + ptau1.Mag2() )
    #     Etau2 = np.sqrt( mtau**2 + ptau2.Mag2() )

    #     Enu1 = np.sqrt( pnu1.Mag2() )
    #     Enu2 = np.sqrt( mnu**2 + pnu2.Mag2() )

    #     pK = ROOT.Math.XYZVector( m[19], m[20], m[21] )
    #     EK = np.sqrt( mkaon**2 + pK.Mag2() )
    #     pB = ptau1 + ptau2 + pK
    #     EB = Etau1 + Etau2 + EK
    #     MB_squared = EB**2 - pB.Mag2()

    #     x0[dimM] = BV.x()
    #     x0[dimM+1] = BV.y()
    #     x0[dimM+2] = BV.z()
    #     x0[dimM+3] = pB.x()
    #     x0[dimM+4] = pB.y()
    #     x0[dimM+5] = pB.z()
    #     x0[dimM+6] = MB_squared
    #     x0[dimM+7] = ptau1.x()
    #     x0[dimM+8] = ptau1.y()
    #     x0[dimM+9] = ptau1.z()
    #     x0[dimM+10] = Etau1
    #     x0[dimM+11] = pnu1.x()
    #     x0[dimM+12] = pnu1.y()
    #     x0[dimM+13] = pnu1.z()
    #     x0[dimM+14] = Enu1
    #     x0[dimM+15] = ptau2.x()
    #     x0[dimM+16] = ptau2.y()
    #     x0[dimM+17] = ptau2.z()
    #     x0[dimM+18] = Etau2
    #     x0[dimM+19] = pnu2.x()
    #     x0[dimM+20] = pnu2.y()
    #     x0[dimM+21] = pnu2.z()
    #     x0[dimM+22] = Enu2

    # elif(init == 14):
    #     ptau2 = ROOT.Math.XYZVector( MLP_taum_PX, MLP_taum_PY, MLP_taum_PZ )

    #     u1 = (DV1 - BV).Unit()
    #     theta1 = np.arcsin( ( np.sqrt( np.abs((mtau**2 + m3pi1**2 - mnu**2)**2 - 4*(mtau**2)*(m3pi1**2) )) )/( 2*mtau*np.sqrt( p3pi1.Mag2() ) ) )
    #     ptau1_mag = ( (mtau**2 + m3pi1**2 - mnu**2)*np.sqrt(p3pi1.Mag2())*np.cos(theta1) )/( 2*( E3pi1**2 - p3pi1.Mag2()*np.cos(theta1)**2 ) )
    #     ptau1 = ROOT.Math.XYZVector(ptau1_mag * u1.x(), ptau1_mag * u1.y(), ptau1_mag * u1.z())

    #     pnu1 = ptau1 - p3pi1
    #     pnu2 = ptau2 - p3pi2

    #     Etau1 = np.sqrt( mtau**2 + ptau1.Mag2() )
    #     Etau2 = np.sqrt( mtau**2 + ptau2.Mag2() )
    #     Enu1 = np.sqrt( mnu**2 + pnu1.Mag2() )
    #     Enu2 = np.sqrt( mnu**2 + pnu2.Mag2() )

    #     pB = pK + ptau1 + ptau2
    #     EB = EK + Etau1 + Etau2
    #     MB_squared = EB**2 - pB.Mag2()

    #     x0[dimM] = BV.x()
    #     x0[dimM+1] = BV.y()
    #     x0[dimM+2] = BV.z()
    #     x0[dimM+3] = pB.x()
    #     x0[dimM+4] = pB.y()
    #     x0[dimM+5] = pB.z()
    #     x0[dimM+6] = MB_squared
    #     x0[dimM+7] = ptau1.x()
    #     x0[dimM+8] = ptau1.y()
    #     x0[dimM+9] = ptau1.z()
    #     x0[dimM+10] = Etau1
    #     x0[dimM+11] = pnu1.x()
    #     x0[dimM+12] = pnu1.y()
    #     x0[dimM+13] = pnu1.z()
    #     x0[dimM+14] = Enu1
    #     x0[dimM+15] = ptau2.x()
    #     x0[dimM+16] = ptau2.y()
    #     x0[dimM+17] = ptau2.z()
    #     x0[dimM+18] = Etau2
    #     x0[dimM+19] = pnu2.x()
    #     x0[dimM+20] = pnu2.y()
    #     x0[dimM+21] = pnu2.z()
    #     x0[dimM+22] = Enu2

    # elif(init == -2):
    #     # With the offline estimate for the BVz we can obtain all other unknown parameters
    #     BVz = BV.z()

    #     # BV must be in K+ trajectory
    #     BVx = RP.x() + (pK.x()/pK.z())*(BVz - RPz)
    #     BVy = RP.y() + (pK.y()/pK.z())*(BVz - RPz)

    #     # ptau1 must point back to the BV
    #     a1 = (DV1.x() - BVx)/(DV1.z() - BVz)
    #     a2 = (DV1.y() - BVy)/(DV1.z() - BVz)
    #     # ptau1x = a1*ptau1z
    #     # ptau1y = a2*ptau1z

    #     # ptau2 must point back to the BV
    #     b1 = (DV2.x() - BVx)/(DV2.z() - BVz)
    #     b2 = (DV2.y() - BVy)/(DV2.z() - BVz)
    #     # ptau2x = b1*ptau2z
    #     # ptau2y = b2*ptau2z

    #     # pB must point back to the PV
    #     c1 = (BVx - PV.x())/(BVz - PV.z())
    #     c2 = (BVy - PV.y())/(BVz - PV.z())

    #     d = (c1 - a1)/(b1 - c1)
    #     e = (c1*pK.z() - pK.x())/(b1 - c1)
    #     # ptau2z = d*ptau1z + e
        
    #     ptau1z = (c2*pK.z() - pK.y() + (c2 - b2)*e )/(a2 + b2*d - c2 - c2*d)
    #     ptau1x = a1*ptau1z
    #     ptau1y = a2*ptau1z
    #     ptau1 = ROOT.Math.XYZVector(ptau1x, ptau1y, ptau1z)

    #     ptau2z = d*ptau1z + e
    #     ptau2x = b1*ptau2z
    #     ptau2y = b2*ptau2z
    #     ptau2 = ROOT.Math.XYZVector(ptau2x, ptau2y, ptau2z)

    #     # Tau mass constraints
    #     # Etau1 = np.sqrt( mtau**2 + ptau1.Mag2() )
    #     # Etau2 = np.sqrt( mtau**2 + ptau2.Mag2() )

    #     # 3-momentum conservation in DV1 and DV2
    #     pnu1 = ptau1 - p3pi1
    #     pnu2 = ptau2 - p3pi2

    #     # Nu mass constraints
    #     # Enu1 = np.sqrt( mnu**2 + pnu1.Mag2() )
    #     # Enu2 = np.sqrt( mnu**2 + pnu2.Mag2() )

    #     m3pi1 = np.sqrt( E3pi1**2 - p3pi1.Mag2() )
    #     m3pi2 = np.sqrt( E3pi2**2 - p3pi2.Mag2() )

    #     Etau1 = ( mtau**2 + m3pi1**2 - mnu**2 + 2*( ptau1.Dot(p3pi1) ) )/(2*E3pi1)
    #     Etau2 = ( mtau**2 + m3pi2**2 - mnu**2 + 2*( ptau2.Dot(p3pi2) ) )/(2*E3pi2)

    #     Enu1 = ( mtau**2 + mnu**2 - m3pi1**2 + 2*( pnu1.Dot(p3pi1) ) )/(2*E3pi1)
    #     Enu2 = ( mtau**2 + mnu**2 - m3pi2**2 + 2*( pnu2.Dot(p3pi2) ) )/(2*E3pi2)

    #     # 4-momentum conservation in BV
    #     pB = pK + ptau1 + ptau2
    #     EB = EK + Etau1 + Etau2
    #     MB_squared = EB**2 - pB.Mag2()

    #     x0[dimM] = BVx
    #     x0[dimM+1] = BVy
    #     x0[dimM+2] = BVz
    #     x0[dimM+3] = pB.x()
    #     x0[dimM+4] = pB.y()
    #     x0[dimM+5] = pB.z()
    #     x0[dimM+6] = MB_squared
    #     x0[dimM+7] = ptau1.x()
    #     x0[dimM+8] = ptau1.y()
    #     x0[dimM+9] = ptau1.z()
    #     x0[dimM+10] = Etau1
    #     x0[dimM+11] = pnu1.x()
    #     x0[dimM+12] = pnu1.y()
    #     x0[dimM+13] = pnu1.z()
    #     x0[dimM+14] = Enu1
    #     x0[dimM+15] = ptau2.x()
    #     x0[dimM+16] = ptau2.y()
    #     x0[dimM+17] = ptau2.z()
    #     x0[dimM+18] = Etau2
    #     x0[dimM+19] = pnu2.x()
    #     x0[dimM+20] = pnu2.y()
    #     x0[dimM+21] = pnu2.z()
    #     x0[dimM+22] = Enu2

    # elif(init == -1):
    #     BV = ROOT.Math.XYZPoint( x_true[0], x_true[1], x_true[2] )
    #     pB = ROOT.Math.XYZVector( x_true[3], x_true[4], x_true[5] )
    #     Ptau1 = ROOT.Math.PxPyPzEVector( x_true[6], x_true[7], x_true[8], x_true[9] )
    #     P1 = ROOT.Math.PxPyPzEVector( x_true[10], x_true[11], x_true[12], x_true[13] )
    #     P2 = ROOT.Math.PxPyPzEVector( x_true[14], x_true[15], x_true[16], x_true[17] )
    #     P3 = ROOT.Math.PxPyPzEVector( x_true[18], x_true[19], x_true[20], x_true[21] )
    #     Ptau2 = ROOT.Math.PxPyPzEVector( x_true[22], x_true[23], x_true[24], x_true[25] )
    #     P4 = ROOT.Math.PxPyPzEVector( x_true[26], x_true[27], x_true[28], x_true[29] )
    #     P5 = ROOT.Math.PxPyPzEVector( x_true[30], x_true[31], x_true[32], x_true[33] )
    #     P6 = ROOT.Math.PxPyPzEVector( x_true[34], x_true[35], x_true[36], x_true[37] )

    #     x0[dimM] = BV.x()
    #     x0[dimM+1] = BV.y()
    #     x0[dimM+2] = BV.z()
    #     x0[dimM+3] = pB.x()
    #     x0[dimM+4] = pB.y()
    #     x0[dimM+5] = pB.z()
    #     x0[dimM+6] = (5279.41)**2
    #     x0[dimM+7] = Ptau1.x()
    #     x0[dimM+8] = Ptau1.y()
    #     x0[dimM+9] = Ptau1.z()
    #     x0[dimM+10] = Ptau1.t()
    #     x0[dimM+11] = Ptau1.x() - P1.x() - P2.x() - P3.x()
    #     x0[dimM+12] = Ptau1.y() - P1.y() - P2.y() - P3.y()
    #     x0[dimM+13] = Ptau1.z() - P1.z() - P2.z() - P3.z()
    #     x0[dimM+14] = Ptau1.t() - P1.t() - P2.t() - P3.t()
    #     x0[dimM+15] = Ptau2.x()
    #     x0[dimM+16] = Ptau2.y()
    #     x0[dimM+17] = Ptau2.z()
    #     x0[dimM+18] = Ptau2.t()
    #     x0[dimM+19] = Ptau2.x() - P4.x() - P5.x() - P6.x()
    #     x0[dimM+20] = Ptau2.y() - P4.y() - P5.y() - P6.y()
    #     x0[dimM+21] = Ptau2.z() - P4.z() - P5.z() - P6.z()
    #     x0[dimM+22] = Ptau2.t() - P4.t() - P5.t() - P6.t()

    # Initialise lambda (lb = 0 gives trivial solution)
    if(change_L == 0):
        for i in range(dimC):
            x0[dimM+dimX+i] = 0
    elif(change_L == 1):
        x0[dimM+dimX] = 0.001
        x0[dimM+dimX+1] = 0.001
        x0[dimM+dimX+2] = 0.001
        x0[dimM+dimX+3] = 0.001
        x0[dimM+dimX+4] = 0.001
        x0[dimM+dimX+5] = 0.001
        x0[dimM+dimX+6] = 0.0001
        x0[dimM+dimX+7] = 0.0001
        x0[dimM+dimX+8] = 0.0001
        x0[dimM+dimX+9] = 0.0001
        x0[dimM+dimX+10] = 0.01
        x0[dimM+dimX+11] = 0.01
        x0[dimM+dimX+12] = 0.01
        x0[dimM+dimX+13] = 0.01
        x0[dimM+dimX+14] = 0.0001
        x0[dimM+dimX+15] = 0.0001
        x0[dimM+dimX+16] = 0.01
        x0[dimM+dimX+17] = 0.01
        x0[dimM+dimX+18] = 0.01
        x0[dimM+dimX+19] = 0.01
        x0[dimM+dimX+20] = 0.001
        x0[dimM+dimX+21] = 0.001
        x0[dimM+dimX+22] = 0.0001
        x0[dimM+dimX+23] = 0
    elif(change_L == 2):
        x0[dimM+dimX] = 1./10000
        x0[dimM+dimX+1] = 1./10000
        x0[dimM+dimX+2] = 1./10000
        x0[dimM+dimX+3] = 1./10000
        x0[dimM+dimX+4] = 1./1000
        x0[dimM+dimX+5] = 1./1000
        x0[dimM+dimX+6] = 1./10000
        x0[dimM+dimX+7] = 1./10000
        x0[dimM+dimX+8] = 1./10000
        x0[dimM+dimX+9] = 1./10000
        x0[dimM+dimX+10] = 1./10000
        x0[dimM+dimX+11] = 1./10000
        x0[dimM+dimX+12] = 1./1000
        x0[dimM+dimX+13] = 1./1000
        x0[dimM+dimX+14] = 1./10000
        x0[dimM+dimX+15] = 1./10000
        x0[dimM+dimX+16] = 1./10000
        x0[dimM+dimX+17] = 1./10000
        x0[dimM+dimX+18] = 1./1000
        x0[dimM+dimX+19] = 1./1000
        x0[dimM+dimX+20] = 1./10000
        x0[dimM+dimX+21] = 1./10000
        x0[dimM+dimX+22] = 1./100000
        x0[dimM+dimX+23] = 1./100000

    return x0

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

def run_solver(init, year, species, line, BV_offline, params, use_generalised_region, max_iter, eps):
    global chi2_value, tol_value, status

    m = params[0]
    W = params[1]
    RPz = params[2]

    x0 =  x_initial_estimate(init, BV_offline, species, params)

    p3pi1z = x0[8]
    p3pi2z = x0[15]    
    ptau1z = x0[dimM+9]
    ptau2z = x0[dimM+17]

    if(ptau1z < p3pi1z):
        x0[dimM+9] = p3pi1z
    if(ptau2z < p3pi2z):
        x0[dimM+17] = p3pi2z

    try:
        run_fdfsolver(year, species, line, x0, params, use_generalised_region, max_iter, eps)
    except:
        print("ERROR: solver got stuck in iteration")
    
    chi2_value = chisquare_value(x,params)
    tol_value = absolute_sum(f)

    pass_second_test = second_derivative_test(x,params)

    if((status == 0) and (pass_second_test)):
        status = 0
    elif((status == 0) and (not pass_second_test)):
        status = 1
    else:
        status = 2

    # max_sum = np.max(np.abs(f))
    # if((max_sum < 0.1) and pass_second_test):
    #     status = 0
    # elif((max_sum < 0.1) and (not pass_second_test)):
    #     status = 1
    # else:
    #     status = 2

def lowest_chi2(year, species, line, BV_offline, params, init_list, max_iter, eps):
    global x, f, status, nIter, chi2_value, tol_value, init, N_local_minima, M_local_minima, Chi2_local_minima

    x_list = []
    f_list = []
    status_list = []
    nIter_list = []
    sum_list = []
    chi2_list = []
    mass_list = []

    N = len(init_list)

    for i in range(N):
        a = init_list[i]

        run_solver(a, year, species, line, BV_offline, params, True, max_iter, eps)
        if(status != 0):
            run_solver(a, year, species, line, BV_offline, params, False, max_iter, eps)
        
        x_list.append(x)
        f_list.append(f)
        status_list.append(status)
        nIter_list.append(nIter)
        sum_list.append(absolute_sum(f))
        chi2_list.append(chisquare_value(x,params))

        mass_squared = x[dimM+6]
        if(mass_squared > 0):
            mass = np.sqrt(mass_squared)
        else:
            mass = -np.sqrt( -mass_squared )

        mass_list.append(mass)

        # print("init = ", i, "status = ", status, "chi2 = ", chi2_list[i], "tol = ", sum_list[i], "mB = ", mass)

        nIter = 0

    # Return solution with lowest chi2
    all_fail = True
    for i in range(N):
        if( status_list[i] == 0 ):
            all_fail = False

    if(all_fail): # return solution with lowest value for the sum
        sum_index_sort = np.argsort(sum_list)
        i_min = sum_index_sort[0]
        N_local_minima[0] = 0
    else: # return solution with lowest chi2 value
        for i in range(N): # ingore failed fits
            if(status_list[i] != 0):
                chi2_list[i] = 1000000000000000000000000000000
                sum_list[i] = 1000000000000000000000000000000
    
        # Number of local minima 
        mass_pass = []
        for i in range(N):
            if(status_list[i] == 0):
                mass_pass.append(mass_list[i])
        mass_sort = np.sort(mass_pass)

        M_local_minima.append(mass_sort[0])
        N_local_minima[0] = 1
        for i in range(len(mass_sort)-1):
            if( int(mass_sort[i]) != int(mass_sort[i+1]) ):
                N_local_minima[0] += 1
                M_local_minima.append(mass_sort[i+1])

        for i in range(N_local_minima[0]):
            idx = np.where(mass_list == M_local_minima[i])[0][0]
            Chi2_local_minima.append(chi2_list[idx])

        chi2_index_sort = np.argsort(chi2_list)
        i_min = chi2_index_sort[0] 

    # chi2_min = 100000000000000000
    # f_min = 100
    # i_min = 0
    # if(all_fail): # If all fail, returns the one with the lowest value for the sum
    #     for i in range(N):
    #         if( sum_list[i] < f_min ):
    #             f_min = sum_list[i]
    #             i_min = i
    # else: # if one or more passes, from the ones that pass return the one that gives the lowest value for the chi^2
    #     for i in range(N):
    #         if( (status_list[i] == 0) and (chi2_list[i] < chi2_min) ):
    #             chi2_min = chi2_list[i]
    #             i_min = i

    x = x_list[i_min]
    f = f_list[i_min]
    status = status_list[i_min]
    nIter = nIter_list[i_min]
    tol_value = sum_list[i_min]
    chi2_value = chi2_list[i_min]
    init[0] = i_min

    x_list.clear()
    f_list.clear()
    status_list.clear()
    nIter_list.clear()
    sum_list.clear()
    chi2_list.clear()

def absolute_sum(F):
    abs_sum = 0
    for i in range(len(F)):
        abs_sum += np.abs(F[i])
    return abs_sum

def run_fdfsolver(year, species, line, x0, params, use_generalised_region, max_iter, eps):
    global x, f, status, nIter

    mysys = multiroots.gsl_multiroot_function_fdf(equations_f, equations_df, equations_fdf, params, dimM+dimX+dimC)
    if use_generalised_region:
        solver = multiroots.hybridsj(mysys, dimM+dimX+dimC)
    else:
        solver = multiroots.hybridj(mysys, dimM+dimX+dimC)

    tmp = numx.array(x0,)
    solver.set(tmp) 

    # print "# Testing solver ", solver.name() 
    for iter in range(max_iter):
        status = solver.iterate() 
        r = solver.root()
        x = solver.getx()
        f = solver.getf()
        status = multiroots.test_residual(f, eps)

        nIter += 1

        mass_squared = x[dimM+6]
        if(mass_squared > 0):
            mass = np.sqrt(mass_squared)
        else:
            mass = -np.sqrt( -mass_squared )

        chi2 = chisquare_value(x,params)
        tol  = absolute_sum(f)

        # print("iter = ", iter, "status = ", status, "chi2 = ", chi2, "tol = ", tol, "mB = ", mass, "max sum = ", max(f))

        if status == errno.GSL_SUCCESS:
            print("Converged")

        if status == errno.GSL_SUCCESS:
            break
    else:
        raise ValueError

def main(argv):
    global X, X_ERR, F, STATUS, NITER, MB, dMB, init, MB0, MB1, MB2, MB3, MB4
    global reader_taup_PX, reader_taup_PY, reader_taup_PZ, reader_taum_PX, reader_taum_PY, reader_taum_PZ
    global arr_m_1, arr_m_2, arr_m_3, arr_m_4, arr_m_5, arr_m_6, arr_m_7, arr_m_8, arr_m_9, arr_m_10, arr_m_11, arr_m_12, arr_m_13, arr_m_14, arr_m_15, arr_m_16, arr_m_17, arr_m_18, arr_m_19, arr_m_20, arr_m_21, arr_m_22, arr_Kp_RP_Z, arr_eventNumber 

    year = argv[1]
    species_str = argv[2]
    line = argv[3]
    # i_first = argv[4]
    # i_last = argv[5]

    year = int(year)
    species = int(species_str)
    line = int(line)
    # i_first = int(i_first)
    # i_last = int(i_last)

    RECO_files = "/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/201{0}/Species_{1}/pre_sel_tree.txt".format(year,species)

    isKtautau = False
    if((species == 10) or (species == 11) or (species == 12) or (species == 1) or (species == 2) or (species == 3)):
        isKtautau = True

    isD0D0K = False
    if((species == 9) or (species == 0) or (species == -1)):
        isD0D0K = True
    
    isDpDmK = False
    if((species == 4) or (species == 5) or (species == 6)):
        isDpDmK = True
    
    global mtau, mnu

    if isD0D0K: # B+ -> D0 D0 K+
        mtau = 1864.84
        mnu = 493.677
    if isDpDmK:
        mtau = 1869.66

    if(computeDerivatives):
        params = [m_symbols, W_symbols, RPz_symbol]

        print("Computing 1st derivatives")
        symbolic_f = equations_f_exp(params)
        print("Finished computing 1st derivatives")
        if isKtautau:
            with open("/panfs/felician/B2Ktautau/workflow/standalone_fitter/first_derivatives_ktautau.npy", "wb") as first_derivative_file:
                np.save(first_derivative_file, symbolic_f)
        if isDpDmK:
            with open("/panfs/felician/B2Ktautau/workflow/standalone_fitter/first_derivatives_dpdmk.npy", "wb") as first_derivative_file:
                np.save(first_derivative_file, symbolic_f)
        if isD0D0K:
            with open("/panfs/felician/B2Ktautau/workflow/standalone_fitter/first_derivatives_d0d0k.npy", "wb") as first_derivative_file:
                np.save(first_derivative_file, symbolic_f)

        print("Computing 2nd derivatives")
        symbolic_df = equations_df_exp(params)
        print("Finished computing 2nd derivatives")
        if isKtautau:
            with open("/panfs/felician/B2Ktautau/workflow/standalone_fitter/second_derivatives_ktautau.npy", "wb") as second_derivative_file:
                np.save(second_derivative_file, symbolic_df)
        if isDpDmK:
            with open("/panfs/felician/B2Ktautau/workflow/standalone_fitter/second_derivatives_dpdmk.npy", "wb") as second_derivative_file:
                np.save(second_derivative_file, symbolic_df)
        if isD0D0K:
            with open("/panfs/felician/B2Ktautau/workflow/standalone_fitter/second_derivatives_d0d0k.npy", "wb") as second_derivative_file:
                np.save(second_derivative_file, symbolic_df)

        print("Computing 2nd derivatives wrt m")
        symbolic_dfm = equations_dfm_exp(params)
        print("Finished computing 2nd derivatives wrt m")
        if isKtautau:
            with open("/panfs/felician/B2Ktautau/workflow/standalone_fitter/second_derivatives_m_ktautau.npy", "wb") as second_derivative_m_file:
                np.save(second_derivative_m_file, symbolic_dfm)
        if isDpDmK:
            with open("/panfs/felician/B2Ktautau/workflow/standalone_fitter/second_derivatives_m_dpdmk.npy", "wb") as second_derivative_m_file:
                np.save(second_derivative_m_file, symbolic_dfm)
        if isD0D0K:
            with open("/panfs/felician/B2Ktautau/workflow/standalone_fitter/second_derivatives_m_d0d0k.npy", "wb") as second_derivative_m_file:
                np.save(second_derivative_m_file, symbolic_dfm)

        print("Computing chi2 expression")
        symbolic_chi2 = chisquare(params)
        if isKtautau:
            with open("/panfs/felician/B2Ktautau/workflow/standalone_fitter/chisquare_ktautau.npy", "wb") as chi2_file:
                np.save(chi2_file, symbolic_chi2)
        if isDpDmK:
            with open("/panfs/felician/B2Ktautau/workflow/standalone_fitter/chisquare_dpdmk.npy", "wb") as chi2_file:
                np.save(chi2_file, symbolic_chi2)
        if isD0D0K:
            with open("/panfs/felician/B2Ktautau/workflow/standalone_fitter/chisquare_d0d0k.npy", "wb") as chi2_file:
                np.save(chi2_file, symbolic_chi2)

        print("Computing bordered Hessian")
        symbolic_bH = bordered_Hessian(params)
        print("Finished computing bordered Hessian")
        if isKtautau:
            with open("/panfs/felician/B2Ktautau/workflow/standalone_fitter/borderedHessian_ktautau.npy", "wb") as bordHess_file:
                np.save(bordHess_file, symbolic_bH)
        if isDpDmK:
            with open("/panfs/felician/B2Ktautau/workflow/standalone_fitter/borderedHessian_dpdmk.npy", "wb") as bordHess_file:
                np.save(bordHess_file, symbolic_bH)
        if isD0D0K:
            with open("/panfs/felician/B2Ktautau/workflow/standalone_fitter/borderedHessian_d0d0k.npy", "wb") as bordHess_file:
                np.save(bordHess_file, symbolic_bH)
            
        for i in range(dimC):
            print(symbolic_f[dimM+dimX+i])
                        
        quit()
    else:
        global lambdify_f, lambdify_df, lambdify_dfm, lambdify_chi2, lambdify_bH

        print("Loading 1st derivatives")
        if isKtautau:
            symbolic_f = np.load("/panfs/felician/B2Ktautau/workflow/standalone_fitter/first_derivatives_ktautau.npy")
        if isDpDmK:
            symbolic_f = np.load("/panfs/felician/B2Ktautau/workflow/standalone_fitter/first_derivatives_dpdmk.npy")
        if isD0D0K:
            symbolic_f = np.load("/panfs/felician/B2Ktautau/workflow/standalone_fitter/first_derivatives_d0d0k.npy")
        lambdify_f = sp.lambdify( (x_symbols, m_symbols, W_symbols, RPz_symbol), symbolic_f, "numpy")

        print("Loading 2nd derivatives")
        if isKtautau:
            symbolic_df = np.load("/panfs/felician/B2Ktautau/workflow/standalone_fitter/second_derivatives_ktautau.npy")
        if isDpDmK:
            symbolic_df = np.load("/panfs/felician/B2Ktautau/workflow/standalone_fitter/second_derivatives_dpdmk.npy") 
        if isD0D0K:
            symbolic_df = np.load("/panfs/felician/B2Ktautau/workflow/standalone_fitter/second_derivatives_d0d0k.npy") 
        lambdify_df = sp.lambdify( (x_symbols, m_symbols, W_symbols, RPz_symbol), symbolic_df, "numpy")

        print("Loading 2nd derivatives wrt m")
        if isKtautau:
            symbolic_dfm = np.load("/panfs/felician/B2Ktautau/workflow/standalone_fitter/second_derivatives_m_ktautau.npy")
        if isDpDmK:
            symbolic_dfm = np.load("/panfs/felician/B2Ktautau/workflow/standalone_fitter/second_derivatives_m_dpdmk.npy")
        if isD0D0K:
            symbolic_dfm = np.load("/panfs/felician/B2Ktautau/workflow/standalone_fitter/second_derivatives_m_d0d0k.npy")
        lambdify_dfm = sp.lambdify( (x_symbols, m_symbols, W_symbols, RPz_symbol), symbolic_dfm, "numpy")
        
        print("Loading chi2 expression")
        if isKtautau:
            symbolic_chi2 = np.load("/panfs/felician/B2Ktautau/workflow/standalone_fitter/chisquare_ktautau.npy")
        if isDpDmK:
            symbolic_chi2 = np.load("/panfs/felician/B2Ktautau/workflow/standalone_fitter/chisquare_dpdmk.npy")
        if isD0D0K:
            symbolic_chi2 = np.load("/panfs/felician/B2Ktautau/workflow/standalone_fitter/chisquare_d0d0k.npy")
        lambdify_chi2 = sp.lambdify( (x_symbols, m_symbols, W_symbols, RPz_symbol), symbolic_chi2, "numpy")

        print("Loading bordered Hessian expression")
        if isKtautau:
            symbolic_bH = np.load("/panfs/felician/B2Ktautau/workflow/standalone_fitter/borderedHessian_ktautau.npy")
        if isDpDmK:
            symbolic_bH = np.load("/panfs/felician/B2Ktautau/workflow/standalone_fitter/borderedHessian_dpdmk.npy")
        if isD0D0K:
            symbolic_bH = np.load("/panfs/felician/B2Ktautau/workflow/standalone_fitter/borderedHessian_d0d0k.npy")
        lambdify_bH = sp.lambdify( (x_symbols_2nd_test, m_symbols, W_symbols, RPz_symbol), symbolic_bH, "numpy")

    fc = ROOT.TFileCollection("fc", "fc", RECO_files, 1, line)
    t = ROOT.TChain("DecayTree")
    t.AddFileInfoList(fc.GetList())

    file = ROOT.TFile.Open("/panfs/felician/B2Ktautau/workflow/standalone_fitter/201{0}/Species_{1}/{2}.root".format(year,species,line), "UPDATE")
    tree = ROOT.TTree("DecayTree", "DecayTree")

    x_names = ["df_PVx", "df_PVy", "df_PVz", "df_DV1x", "df_DV1y", "df_DV1z", "df_3pi1_PX", "df_3pi1_PY", "df_3pi1_PZ", "df_3pi1_PE", "df_DV2x", "df_DV2y", "df_DV2z", "df_3pi2_PX", "df_3pi2_PY", "df_3pi2_PZ", "df_3pi2_PE", "df_RPx", "df_RPy", "df_Kp_PX", "df_Kp_PY", "df_Kp_PZ", "df_BVx", "df_BVy", "df_BVz", "df_Bp_PX", "df_Bp_PY", "df_Bp_PZ", "df_Bp_M2", "df_taup_PX", "df_taup_PY", "df_taup_PZ", "df_taup_PE", "df_antinutau_PX", "df_antinutau_PY", "df_antinutau_PZ", "df_antinutau_PE", "df_taum_PX", "df_taum_PY", "df_taum_PZ", "df_taum_PE", "df_nutau_PX", "df_nutau_PY", "df_nutau_PZ", "df_nutau_PE"]
    x_err_names = ["df_PVx_err", "df_PVy_err", "df_PVz_err", "df_DV1x_err", "df_DV1y_err", "df_DV1z_err", "df_3pi1_PX_err", "df_3pi1_PY_err", "df_3pi1_PZ_err", "df_3pi1_PE_err", "df_DV2x_err", "df_DV2y_err", "df_DV2z_err", "df_3pi2_PX_err", "df_3pi2_PY_err", "df_3pi2_PZ_err", "df_3pi2_PE_err", "df_RPx_err", "df_RPy_err", "df_Kp_PX_err", "df_Kp_PY_err", "df_Kp_PZ_err", "df_BVx_err", "df_BVy_err", "df_BVz_err", "df_Bp_PX_err", "df_Bp_PY_err", "df_Bp_PZ_err", "df_Bp_M2_err", "df_taup_PX_err", "df_taup_PY_err", "df_taup_PZ_err", "df_taup_PE_err", "df_antinutau_PX_err", "df_antinutau_PY_err", "df_antinutau_PZ_err", "df_antinutau_PE_err", "df_taum_PX_err", "df_taum_PY_err", "df_taum_PZ_err", "df_taum_PE_err", "df_nutau_PX_err", "df_nutau_PY_err", "df_nutau_PZ_err", "df_nutau_PE_err"]
    x_true_names = ["Bp_TRUEENDVERTEX_X", "Bp_TRUEENDVERTEX_Y", "Bp_TRUEENDVERTEX_Z", "Bp_TRUEP_X", "Bp_TRUEP_Y", "Bp_TRUEP_Z", "taup_TRUEP_X", "taup_TRUEP_Y", "taup_TRUEP_Z", "taup_TRUEP_E", "taup_pi1_TRUEP_X", "taup_pi1_TRUEP_Y", "taup_pi1_TRUEP_Z", "taup_pi1_TRUEP_E", "taup_pi2_TRUEP_X", "taup_pi2_TRUEP_Y", "taup_pi2_TRUEP_Z", "taup_pi2_TRUEP_E", "taup_pi3_TRUEP_X", "taup_pi3_TRUEP_Y", "taup_pi3_TRUEP_Z", "taup_pi3_TRUEP_E", "taum_TRUEP_X", "taum_TRUEP_Y", "taum_TRUEP_Z", "taum_TRUEP_E", "taum_pi1_TRUEP_X", "taum_pi1_TRUEP_Y", "taum_pi1_TRUEP_Z", "taum_pi1_TRUEP_E", "taum_pi2_TRUEP_X", "taum_pi2_TRUEP_Y", "taum_pi2_TRUEP_Z", "taum_pi2_TRUEP_E", "taum_pi3_TRUEP_X", "taum_pi3_TRUEP_Y", "taum_pi3_TRUEP_Z", "taum_pi3_TRUEP_E"]
    if isDpDmK:
        x_true_names = ["Bp_TRUEENDVERTEX_X", "Bp_TRUEENDVERTEX_Y", "Bp_TRUEENDVERTEX_Z", "Bp_TRUEP_X", "Bp_TRUEP_Y", "Bp_TRUEP_Z", "Dp_TRUEP_X", "Dp_TRUEP_Y", "Dp_TRUEP_Z", "Dp_TRUEP_E", "Dp_K_TRUEP_X", "Dp_K_TRUEP_Y", "Dp_K_TRUEP_Z", "Dp_K_TRUEP_E", "Dp_pi1_TRUEP_X", "Dp_pi1_TRUEP_Y", "Dp_pi1_TRUEP_Z", "Dp_pi1_TRUEP_E", "Dp_pi2_TRUEP_X", "Dp_pi2_TRUEP_Y", "Dp_pi2_TRUEP_Z", "Dp_pi2_TRUEP_E", "Dm_TRUEP_X", "Dm_TRUEP_Y", "Dm_TRUEP_Z", "Dm_TRUEP_E", "Dm_K_TRUEP_X", "Dm_K_TRUEP_Y", "Dm_K_TRUEP_Z", "Dm_K_TRUEP_E", "Dm_pi1_TRUEP_X", "Dm_pi1_TRUEP_Y", "Dm_pi1_TRUEP_Z", "Dm_pi1_TRUEP_E", "Dm_pi2_TRUEP_X", "Dm_pi2_TRUEP_Y", "Dm_pi2_TRUEP_Z", "Dm_pi2_TRUEP_E"]
    if isD0D0K:
        x_true_names = ["Bp_TRUEENDVERTEX_X", "Bp_TRUEENDVERTEX_Y", "Bp_TRUEENDVERTEX_Z", "Bp_TRUEP_X", "Bp_TRUEP_Y", "Bp_TRUEP_Z", "D0_TRUEP_X", "D0_TRUEP_Y", "D0_TRUEP_Z", "D0_TRUEP_E", "D0_pi1_TRUEP_X", "D0_pi1_TRUEP_Y", "D0_pi1_TRUEP_Z", "D0_pi1_TRUEP_E", "D0_pi2_TRUEP_X", "D0_pi2_TRUEP_Y", "D0_pi2_TRUEP_Z", "D0_pi2_TRUEP_E", "D0_pi3_TRUEP_X", "D0_pi3_TRUEP_Y", "D0_pi3_TRUEP_Z", "D0_pi3_TRUEP_E", "D0bar_TRUEP_X", "D0bar_TRUEP_Y", "D0bar_TRUEP_Z", "D0bar_TRUEP_E", "D0bar_pi1_TRUEP_X", "D0bar_pi1_TRUEP_Y", "D0bar_pi1_TRUEP_Z", "D0bar_pi1_TRUEP_E", "D0bar_pi2_TRUEP_X", "D0bar_pi2_TRUEP_Y", "D0bar_pi2_TRUEP_Z", "D0bar_pi2_TRUEP_E", "D0bar_pi3_TRUEP_X", "D0bar_pi3_TRUEP_Y", "D0bar_pi3_TRUEP_Z", "D0bar_pi3_TRUEP_E"]

    for i in range(dimM+dimX):
       tree.Branch(x_names[i], X[i], x_names[i]+"/D")
       tree.Branch(x_err_names[i], X_ERR[i], x_err_names[i]+"/D")
       # tree.Branch("df_det_{0}".format(i), D[i], "df_det_{0}/D".format(i))
    for i in range(dimM+dimX+dimC):
       tree.Branch("df_F_{0}".format(i), F[i], "df_F_{0}/D".format(i))
    tree.Branch("df_status", STATUS, "df_status/i")
    tree.Branch("df_nIter", NITER, "df_nIter/i")
    tree.Branch("df_Bp_M", MB, "df_Bp_M/D")
    tree.Branch("df_Bp_MERR", dMB, "df_Bp_MERR/D")
    tree.Branch("df_chi2", chi2, "df_chi2/D")
    tree.Branch("df_tolerance", tolerance, "df_tolerance/D")
    tree.Branch("df_init", init, "df_init/i")
    tree.Branch("df_N_local_min", N_local_minima, "df_N_local_min/i")
    tree.Branch("df_M_local_min", M_local, "df_M_local_min[df_N_local_min]/D")
    tree.Branch("df_chi2_local_min", Chi2_local, "df_chi2_local_min[df_N_local_min]/D")
    tree.Branch("df_time", theTime, "df_time/D")

    # tree.Branch("df_MB0", MB0, "df_MB0/D")
    # tree.Branch("df_MB1", MB1, "df_MB1/D")
    # tree.Branch("df_MB2", MB2, "df_MB2/D")
    # tree.Branch("df_MB3", MB3, "df_MB3/D")
    # tree.Branch("df_MB4", MB4, "df_MB4/D")

    # Load TMVA library (once)
    ROOT.TMVA.Tools.Instance()
    reader_taup_PX = ROOT.TMVA.Reader( "V:Color:Silent" )
    reader_taup_PY = ROOT.TMVA.Reader( "V:Color:Silent" ) 
    reader_taup_PZ = ROOT.TMVA.Reader( "V:Color:Silent" ) 
    reader_taum_PX = ROOT.TMVA.Reader( "V:Color:Silent" ) 
    reader_taum_PY = ROOT.TMVA.Reader( "V:Color:Silent" ) 
    reader_taum_PZ = ROOT.TMVA.Reader( "V:Color:Silent" ) 

    weightfile_taup_PX = ""
    weightfile_taup_PY = ""
    weightfile_taup_PZ = ""
    weightfile_taum_PX = ""
    weightfile_taum_PY = ""
    weightfile_taum_PZ = ""

    if((species == 1) or (species == 10) or (species == 11) or (species == 12) or (species == 2) or (species == 3)): # Ktautau
        weightfile_taup_PX = "/panfs/felician/MLP_weights/KTauTau_MLP_Train_taup_PX/dataset/weights/TMVARegression_taup_TRUEP_X_MLP.weights.xml"
        weightfile_taup_PY = "/panfs/felician/MLP_weights/KTauTau_MLP_Train_taup_PY/dataset/weights/TMVARegression_taup_TRUEP_Y_MLP.weights.xml"
        weightfile_taup_PZ = "/panfs/felician/MLP_weights/KTauTau_MLP_Train_taup_PZ/dataset/weights/TMVARegression_taup_TRUEP_Z_MLP.weights.xml"
        weightfile_taum_PX = "/panfs/felician/MLP_weights/KTauTau_MLP_Train_taum_PX/dataset/weights/TMVARegression_taum_TRUEP_X_MLP.weights.xml"
        weightfile_taum_PY = "/panfs/felician/MLP_weights/KTauTau_MLP_Train_taum_PY/dataset/weights/TMVARegression_taum_TRUEP_Y_MLP.weights.xml"
        weightfile_taum_PZ = "/panfs/felician/MLP_weights/KTauTau_MLP_Train_taum_PZ/dataset/weights/TMVARegression_taum_TRUEP_Z_MLP.weights.xml"
    elif( (species == 4) or (species == 5) or (species == 6)): # D+D-K+
        weightfile_taup_PX = "/panfs/felician/MLP_weights/DDK_MLP_Train_taup_PX/dataset/weights/TMVARegression_Dp_TRUEP_X_MLP.weights.xml"
        weightfile_taup_PY = "/panfs/felician/MLP_weights/DDK_MLP_Train_taup_PY/dataset/weights/TMVARegression_Dp_TRUEP_Y_MLP.weights.xml"
        weightfile_taup_PZ = "/panfs/felician/MLP_weights/DDK_MLP_Train_taup_PZ/dataset/weights/TMVARegression_Dp_TRUEP_Z_MLP.weights.xml"
        weightfile_taum_PX = "/panfs/felician/MLP_weights/DDK_MLP_Train_taum_PX/dataset/weights/TMVARegression_Dm_TRUEP_X_MLP.weights.xml"
        weightfile_taum_PY = "/panfs/felician/MLP_weights/DDK_MLP_Train_taum_PY/dataset/weights/TMVARegression_Dm_TRUEP_Y_MLP.weights.xml"
        weightfile_taum_PZ = "/panfs/felician/MLP_weights/DDK_MLP_Train_taum_PZ/dataset/weights/TMVARegression_Dm_TRUEP_Z_MLP.weights.xml"
    elif( (species == 9) or (species == 0) or (species == -1) ):
        weightfile_taup_PX = "/panfs/felician/MLP_weights/D0D0K_MLP_Train_D0bar_PX/dataset/weights/TMVARegression_D0bar_TRUEP_X_MLP.weights.xml"
        weightfile_taup_PY = "/panfs/felician/MLP_weights/D0D0K_MLP_Train_D0bar_PY/dataset/weights/TMVARegression_D0bar_TRUEP_Y_MLP.weights.xml"
        weightfile_taup_PZ = "/panfs/felician/MLP_weights/D0D0K_MLP_Train_D0bar_PZ/dataset/weights/TMVARegression_D0bar_TRUEP_Z_MLP.weights.xml"
        weightfile_taum_PX = "/panfs/felician/MLP_weights/D0D0K_MLP_Train_D0_PX/dataset/weights/TMVARegression_D0_TRUEP_X_MLP.weights.xml"
        weightfile_taum_PY = "/panfs/felician/MLP_weights/D0D0K_MLP_Train_D0_PY/dataset/weights/TMVARegression_D0_TRUEP_Y_MLP.weights.xml"
        weightfile_taum_PZ = "/panfs/felician/MLP_weights/D0D0K_MLP_Train_D0_PZ/dataset/weights/TMVARegression_D0_TRUEP_Z_MLP.weights.xml"

    arr_m_1 = array('f', [0])
    arr_m_2 = array('f', [0])
    arr_m_3 = array('f', [0])
    arr_m_4 = array('f', [0])
    arr_m_5 = array('f', [0])
    arr_m_6 = array('f', [0])
    arr_m_7 = array('f', [0])
    arr_m_8 = array('f', [0])
    arr_m_9 = array('f', [0])
    arr_m_10 = array('f', [0])
    arr_m_11 = array('f', [0])
    arr_m_12 = array('f', [0])
    arr_m_13 = array('f', [0])
    arr_m_14 = array('f', [0])
    arr_m_15 = array('f', [0])
    arr_m_16 = array('f', [0])
    arr_m_17 = array('f', [0])
    arr_m_18 = array('f', [0])
    arr_m_19 = array('f', [0])
    arr_m_20 = array('f', [0])
    arr_m_21 = array('f', [0])
    arr_m_22 = array('f', [0])
    arr_Kp_RP_Z = array('f', [0])
    arr_eventNumber = array('i', [0])

    if((species == 9) or (species == 0) or (species == -1)):
        reader_taup_PX.AddVariable("df_mprime_1", arr_m_1)
        reader_taup_PX.AddVariable("df_mprime_2", arr_m_2)
        reader_taup_PX.AddVariable("df_mprime_3", arr_m_3)
        reader_taup_PX.AddVariable("df_mprime_4", arr_m_4)
        reader_taup_PX.AddVariable("df_mprime_5", arr_m_5)
        reader_taup_PX.AddVariable("df_mprime_6", arr_m_6)
        reader_taup_PX.AddVariable("df_mprime_7", arr_m_7)
        reader_taup_PX.AddVariable("df_mprime_8", arr_m_8)
        reader_taup_PX.AddVariable("df_mprime_9", arr_m_9)
        reader_taup_PX.AddVariable("df_mprime_10", arr_m_10)
        reader_taup_PX.AddVariable("df_mprime_11", arr_m_11)
        reader_taup_PX.AddVariable("df_mprime_12", arr_m_12)
        reader_taup_PX.AddVariable("df_mprime_13", arr_m_13)
        reader_taup_PX.AddVariable("df_mprime_14", arr_m_14)
        reader_taup_PX.AddVariable("df_mprime_15", arr_m_15)
        reader_taup_PX.AddVariable("df_mprime_16", arr_m_16)
        reader_taup_PX.AddVariable("df_mprime_17", arr_m_17)
        reader_taup_PX.AddVariable("df_mprime_18", arr_m_18)
        reader_taup_PX.AddVariable("df_mprime_19", arr_m_19)
        reader_taup_PX.AddVariable("df_mprime_20", arr_m_20)
        reader_taup_PX.AddVariable("df_mprime_21", arr_m_21)
        reader_taup_PX.AddVariable("df_mprime_22", arr_m_22)
        reader_taup_PX.AddVariable("Kp_RP_Z", arr_Kp_RP_Z)
        reader_taup_PX.AddSpectator("eventNumber", arr_eventNumber)

        reader_taup_PY.AddVariable("df_mprime_1", arr_m_1)
        reader_taup_PY.AddVariable("df_mprime_2", arr_m_2)
        reader_taup_PY.AddVariable("df_mprime_3", arr_m_3)
        reader_taup_PY.AddVariable("df_mprime_4", arr_m_4)
        reader_taup_PY.AddVariable("df_mprime_5", arr_m_5)
        reader_taup_PY.AddVariable("df_mprime_6", arr_m_6)
        reader_taup_PY.AddVariable("df_mprime_7", arr_m_7)
        reader_taup_PY.AddVariable("df_mprime_8", arr_m_8)
        reader_taup_PY.AddVariable("df_mprime_9", arr_m_9)
        reader_taup_PY.AddVariable("df_mprime_10", arr_m_10)
        reader_taup_PY.AddVariable("df_mprime_11", arr_m_11)
        reader_taup_PY.AddVariable("df_mprime_12", arr_m_12)
        reader_taup_PY.AddVariable("df_mprime_13", arr_m_13)
        reader_taup_PY.AddVariable("df_mprime_14", arr_m_14)
        reader_taup_PY.AddVariable("df_mprime_15", arr_m_15)
        reader_taup_PY.AddVariable("df_mprime_16", arr_m_16)
        reader_taup_PY.AddVariable("df_mprime_17", arr_m_17)
        reader_taup_PY.AddVariable("df_mprime_18", arr_m_18)
        reader_taup_PY.AddVariable("df_mprime_19", arr_m_19)
        reader_taup_PY.AddVariable("df_mprime_20", arr_m_20)
        reader_taup_PY.AddVariable("df_mprime_21", arr_m_21)
        reader_taup_PY.AddVariable("df_mprime_22", arr_m_22)
        reader_taup_PY.AddVariable("Kp_RP_Z", arr_Kp_RP_Z)
        reader_taup_PY.AddSpectator("eventNumber", arr_eventNumber)

        reader_taup_PZ.AddVariable("df_mprime_1", arr_m_1)
        reader_taup_PZ.AddVariable("df_mprime_2", arr_m_2)
        reader_taup_PZ.AddVariable("df_mprime_3", arr_m_3)
        reader_taup_PZ.AddVariable("df_mprime_4", arr_m_4)
        reader_taup_PZ.AddVariable("df_mprime_5", arr_m_5)
        reader_taup_PZ.AddVariable("df_mprime_6", arr_m_6)
        reader_taup_PZ.AddVariable("df_mprime_7", arr_m_7)
        reader_taup_PZ.AddVariable("df_mprime_8", arr_m_8)
        reader_taup_PZ.AddVariable("df_mprime_9", arr_m_9)
        reader_taup_PZ.AddVariable("df_mprime_10", arr_m_10)
        reader_taup_PZ.AddVariable("df_mprime_11", arr_m_11)
        reader_taup_PZ.AddVariable("df_mprime_12", arr_m_12)
        reader_taup_PZ.AddVariable("df_mprime_13", arr_m_13)
        reader_taup_PZ.AddVariable("df_mprime_14", arr_m_14)
        reader_taup_PZ.AddVariable("df_mprime_15", arr_m_15)
        reader_taup_PZ.AddVariable("df_mprime_16", arr_m_16)
        reader_taup_PZ.AddVariable("df_mprime_17", arr_m_17)
        reader_taup_PZ.AddVariable("df_mprime_18", arr_m_18)
        reader_taup_PZ.AddVariable("df_mprime_19", arr_m_19)
        reader_taup_PZ.AddVariable("df_mprime_20", arr_m_20)
        reader_taup_PZ.AddVariable("df_mprime_21", arr_m_21)
        reader_taup_PZ.AddVariable("df_mprime_22", arr_m_22)
        reader_taup_PZ.AddVariable("Kp_RP_Z", arr_Kp_RP_Z)
        reader_taup_PZ.AddSpectator("eventNumber", arr_eventNumber)

        reader_taum_PX.AddVariable("df_mprime_1", arr_m_1)
        reader_taum_PX.AddVariable("df_mprime_2", arr_m_2)
        reader_taum_PX.AddVariable("df_mprime_3", arr_m_3)
        reader_taum_PX.AddVariable("df_mprime_4", arr_m_4)
        reader_taum_PX.AddVariable("df_mprime_5", arr_m_5)
        reader_taum_PX.AddVariable("df_mprime_6", arr_m_6)
        reader_taum_PX.AddVariable("df_mprime_7", arr_m_7)
        reader_taum_PX.AddVariable("df_mprime_8", arr_m_8)
        reader_taum_PX.AddVariable("df_mprime_9", arr_m_9)
        reader_taum_PX.AddVariable("df_mprime_10", arr_m_10)
        reader_taum_PX.AddVariable("df_mprime_11", arr_m_11)
        reader_taum_PX.AddVariable("df_mprime_12", arr_m_12)
        reader_taum_PX.AddVariable("df_mprime_13", arr_m_13)
        reader_taum_PX.AddVariable("df_mprime_14", arr_m_14)
        reader_taum_PX.AddVariable("df_mprime_15", arr_m_15)
        reader_taum_PX.AddVariable("df_mprime_16", arr_m_16)
        reader_taum_PX.AddVariable("df_mprime_17", arr_m_17)
        reader_taum_PX.AddVariable("df_mprime_18", arr_m_18)
        reader_taum_PX.AddVariable("df_mprime_19", arr_m_19)
        reader_taum_PX.AddVariable("df_mprime_20", arr_m_20)
        reader_taum_PX.AddVariable("df_mprime_21", arr_m_21)
        reader_taum_PX.AddVariable("df_mprime_22", arr_m_22)
        reader_taum_PX.AddVariable("Kp_RP_Z", arr_Kp_RP_Z)
        reader_taum_PX.AddSpectator("eventNumber", arr_eventNumber)

        reader_taum_PY.AddVariable("df_mprime_1", arr_m_1)
        reader_taum_PY.AddVariable("df_mprime_2", arr_m_2)
        reader_taum_PY.AddVariable("df_mprime_3", arr_m_3)
        reader_taum_PY.AddVariable("df_mprime_4", arr_m_4)
        reader_taum_PY.AddVariable("df_mprime_5", arr_m_5)
        reader_taum_PY.AddVariable("df_mprime_6", arr_m_6)
        reader_taum_PY.AddVariable("df_mprime_7", arr_m_7)
        reader_taum_PY.AddVariable("df_mprime_8", arr_m_8)
        reader_taum_PY.AddVariable("df_mprime_9", arr_m_9)
        reader_taum_PY.AddVariable("df_mprime_10", arr_m_10)
        reader_taum_PY.AddVariable("df_mprime_11", arr_m_11)
        reader_taum_PY.AddVariable("df_mprime_12", arr_m_12)
        reader_taum_PY.AddVariable("df_mprime_13", arr_m_13)
        reader_taum_PY.AddVariable("df_mprime_14", arr_m_14)
        reader_taum_PY.AddVariable("df_mprime_15", arr_m_15)
        reader_taum_PY.AddVariable("df_mprime_16", arr_m_16)
        reader_taum_PY.AddVariable("df_mprime_17", arr_m_17)
        reader_taum_PY.AddVariable("df_mprime_18", arr_m_18)
        reader_taum_PY.AddVariable("df_mprime_19", arr_m_19)
        reader_taum_PY.AddVariable("df_mprime_20", arr_m_20)
        reader_taum_PY.AddVariable("df_mprime_21", arr_m_21)
        reader_taum_PY.AddVariable("df_mprime_22", arr_m_22)
        reader_taum_PY.AddVariable("Kp_RP_Z", arr_Kp_RP_Z)
        reader_taum_PY.AddSpectator("eventNumber", arr_eventNumber)

        reader_taum_PZ.AddVariable("df_mprime_1", arr_m_1)
        reader_taum_PZ.AddVariable("df_mprime_2", arr_m_2)
        reader_taum_PZ.AddVariable("df_mprime_3", arr_m_3)
        reader_taum_PZ.AddVariable("df_mprime_4", arr_m_4)
        reader_taum_PZ.AddVariable("df_mprime_5", arr_m_5)
        reader_taum_PZ.AddVariable("df_mprime_6", arr_m_6)
        reader_taum_PZ.AddVariable("df_mprime_7", arr_m_7)
        reader_taum_PZ.AddVariable("df_mprime_8", arr_m_8)
        reader_taum_PZ.AddVariable("df_mprime_9", arr_m_9)
        reader_taum_PZ.AddVariable("df_mprime_10", arr_m_10)
        reader_taum_PZ.AddVariable("df_mprime_11", arr_m_11)
        reader_taum_PZ.AddVariable("df_mprime_12", arr_m_12)
        reader_taum_PZ.AddVariable("df_mprime_13", arr_m_13)
        reader_taum_PZ.AddVariable("df_mprime_14", arr_m_14)
        reader_taum_PZ.AddVariable("df_mprime_15", arr_m_15)
        reader_taum_PZ.AddVariable("df_mprime_16", arr_m_16)
        reader_taum_PZ.AddVariable("df_mprime_17", arr_m_17)
        reader_taum_PZ.AddVariable("df_mprime_18", arr_m_18)
        reader_taum_PZ.AddVariable("df_mprime_19", arr_m_19)
        reader_taum_PZ.AddVariable("df_mprime_20", arr_m_20)
        reader_taum_PZ.AddVariable("df_mprime_21", arr_m_21)
        reader_taum_PZ.AddVariable("df_mprime_22", arr_m_22)
        reader_taum_PZ.AddVariable("Kp_RP_Z", arr_Kp_RP_Z)
        reader_taum_PZ.AddSpectator("eventNumber", arr_eventNumber)
    else:
        reader_taup_PX.AddVariable("df_m_1", arr_m_1)
        reader_taup_PX.AddVariable("df_m_2", arr_m_2)
        reader_taup_PX.AddVariable("df_m_3", arr_m_3)
        reader_taup_PX.AddVariable("df_m_4", arr_m_4)
        reader_taup_PX.AddVariable("df_m_5", arr_m_5)
        reader_taup_PX.AddVariable("df_m_6", arr_m_6)
        reader_taup_PX.AddVariable("df_m_7", arr_m_7)
        reader_taup_PX.AddVariable("df_m_8", arr_m_8)
        reader_taup_PX.AddVariable("df_m_9", arr_m_9)
        reader_taup_PX.AddVariable("df_m_10", arr_m_10)
        reader_taup_PX.AddVariable("df_m_11", arr_m_11)
        reader_taup_PX.AddVariable("df_m_12", arr_m_12)
        reader_taup_PX.AddVariable("df_m_13", arr_m_13)
        reader_taup_PX.AddVariable("df_m_14", arr_m_14)
        reader_taup_PX.AddVariable("df_m_15", arr_m_15)
        reader_taup_PX.AddVariable("df_m_16", arr_m_16)
        reader_taup_PX.AddVariable("df_m_17", arr_m_17)
        reader_taup_PX.AddVariable("df_m_18", arr_m_18)
        reader_taup_PX.AddVariable("df_m_19", arr_m_19)
        reader_taup_PX.AddVariable("df_m_20", arr_m_20)
        reader_taup_PX.AddVariable("df_m_21", arr_m_21)
        reader_taup_PX.AddVariable("df_m_22", arr_m_22)
        reader_taup_PX.AddVariable("Kp_RP_Z", arr_Kp_RP_Z)
        reader_taup_PX.AddSpectator("eventNumber", arr_eventNumber)

        reader_taup_PY.AddVariable("df_m_1", arr_m_1)
        reader_taup_PY.AddVariable("df_m_2", arr_m_2)
        reader_taup_PY.AddVariable("df_m_3", arr_m_3)
        reader_taup_PY.AddVariable("df_m_4", arr_m_4)
        reader_taup_PY.AddVariable("df_m_5", arr_m_5)
        reader_taup_PY.AddVariable("df_m_6", arr_m_6)
        reader_taup_PY.AddVariable("df_m_7", arr_m_7)
        reader_taup_PY.AddVariable("df_m_8", arr_m_8)
        reader_taup_PY.AddVariable("df_m_9", arr_m_9)
        reader_taup_PY.AddVariable("df_m_10", arr_m_10)
        reader_taup_PY.AddVariable("df_m_11", arr_m_11)
        reader_taup_PY.AddVariable("df_m_12", arr_m_12)
        reader_taup_PY.AddVariable("df_m_13", arr_m_13)
        reader_taup_PY.AddVariable("df_m_14", arr_m_14)
        reader_taup_PY.AddVariable("df_m_15", arr_m_15)
        reader_taup_PY.AddVariable("df_m_16", arr_m_16)
        reader_taup_PY.AddVariable("df_m_17", arr_m_17)
        reader_taup_PY.AddVariable("df_m_18", arr_m_18)
        reader_taup_PY.AddVariable("df_m_19", arr_m_19)
        reader_taup_PY.AddVariable("df_m_20", arr_m_20)
        reader_taup_PY.AddVariable("df_m_21", arr_m_21)
        reader_taup_PY.AddVariable("df_m_22", arr_m_22)
        reader_taup_PY.AddVariable("Kp_RP_Z", arr_Kp_RP_Z)
        reader_taup_PY.AddSpectator("eventNumber", arr_eventNumber)

        reader_taup_PZ.AddVariable("df_m_1", arr_m_1)
        reader_taup_PZ.AddVariable("df_m_2", arr_m_2)
        reader_taup_PZ.AddVariable("df_m_3", arr_m_3)
        reader_taup_PZ.AddVariable("df_m_4", arr_m_4)
        reader_taup_PZ.AddVariable("df_m_5", arr_m_5)
        reader_taup_PZ.AddVariable("df_m_6", arr_m_6)
        reader_taup_PZ.AddVariable("df_m_7", arr_m_7)
        reader_taup_PZ.AddVariable("df_m_8", arr_m_8)
        reader_taup_PZ.AddVariable("df_m_9", arr_m_9)
        reader_taup_PZ.AddVariable("df_m_10", arr_m_10)
        reader_taup_PZ.AddVariable("df_m_11", arr_m_11)
        reader_taup_PZ.AddVariable("df_m_12", arr_m_12)
        reader_taup_PZ.AddVariable("df_m_13", arr_m_13)
        reader_taup_PZ.AddVariable("df_m_14", arr_m_14)
        reader_taup_PZ.AddVariable("df_m_15", arr_m_15)
        reader_taup_PZ.AddVariable("df_m_16", arr_m_16)
        reader_taup_PZ.AddVariable("df_m_17", arr_m_17)
        reader_taup_PZ.AddVariable("df_m_18", arr_m_18)
        reader_taup_PZ.AddVariable("df_m_19", arr_m_19)
        reader_taup_PZ.AddVariable("df_m_20", arr_m_20)
        reader_taup_PZ.AddVariable("df_m_21", arr_m_21)
        reader_taup_PZ.AddVariable("df_m_22", arr_m_22)
        reader_taup_PZ.AddVariable("Kp_RP_Z", arr_Kp_RP_Z)
        reader_taup_PZ.AddSpectator("eventNumber", arr_eventNumber)

        reader_taum_PX.AddVariable("df_m_1", arr_m_1)
        reader_taum_PX.AddVariable("df_m_2", arr_m_2)
        reader_taum_PX.AddVariable("df_m_3", arr_m_3)
        reader_taum_PX.AddVariable("df_m_4", arr_m_4)
        reader_taum_PX.AddVariable("df_m_5", arr_m_5)
        reader_taum_PX.AddVariable("df_m_6", arr_m_6)
        reader_taum_PX.AddVariable("df_m_7", arr_m_7)
        reader_taum_PX.AddVariable("df_m_8", arr_m_8)
        reader_taum_PX.AddVariable("df_m_9", arr_m_9)
        reader_taum_PX.AddVariable("df_m_10", arr_m_10)
        reader_taum_PX.AddVariable("df_m_11", arr_m_11)
        reader_taum_PX.AddVariable("df_m_12", arr_m_12)
        reader_taum_PX.AddVariable("df_m_13", arr_m_13)
        reader_taum_PX.AddVariable("df_m_14", arr_m_14)
        reader_taum_PX.AddVariable("df_m_15", arr_m_15)
        reader_taum_PX.AddVariable("df_m_16", arr_m_16)
        reader_taum_PX.AddVariable("df_m_17", arr_m_17)
        reader_taum_PX.AddVariable("df_m_18", arr_m_18)
        reader_taum_PX.AddVariable("df_m_19", arr_m_19)
        reader_taum_PX.AddVariable("df_m_20", arr_m_20)
        reader_taum_PX.AddVariable("df_m_21", arr_m_21)
        reader_taum_PX.AddVariable("df_m_22", arr_m_22)
        reader_taum_PX.AddVariable("Kp_RP_Z", arr_Kp_RP_Z)
        reader_taum_PX.AddSpectator("eventNumber", arr_eventNumber)

        reader_taum_PY.AddVariable("df_m_1", arr_m_1)
        reader_taum_PY.AddVariable("df_m_2", arr_m_2)
        reader_taum_PY.AddVariable("df_m_3", arr_m_3)
        reader_taum_PY.AddVariable("df_m_4", arr_m_4)
        reader_taum_PY.AddVariable("df_m_5", arr_m_5)
        reader_taum_PY.AddVariable("df_m_6", arr_m_6)
        reader_taum_PY.AddVariable("df_m_7", arr_m_7)
        reader_taum_PY.AddVariable("df_m_8", arr_m_8)
        reader_taum_PY.AddVariable("df_m_9", arr_m_9)
        reader_taum_PY.AddVariable("df_m_10", arr_m_10)
        reader_taum_PY.AddVariable("df_m_11", arr_m_11)
        reader_taum_PY.AddVariable("df_m_12", arr_m_12)
        reader_taum_PY.AddVariable("df_m_13", arr_m_13)
        reader_taum_PY.AddVariable("df_m_14", arr_m_14)
        reader_taum_PY.AddVariable("df_m_15", arr_m_15)
        reader_taum_PY.AddVariable("df_m_16", arr_m_16)
        reader_taum_PY.AddVariable("df_m_17", arr_m_17)
        reader_taum_PY.AddVariable("df_m_18", arr_m_18)
        reader_taum_PY.AddVariable("df_m_19", arr_m_19)
        reader_taum_PY.AddVariable("df_m_20", arr_m_20)
        reader_taum_PY.AddVariable("df_m_21", arr_m_21)
        reader_taum_PY.AddVariable("df_m_22", arr_m_22)
        reader_taum_PY.AddVariable("Kp_RP_Z", arr_Kp_RP_Z)
        reader_taum_PY.AddSpectator("eventNumber", arr_eventNumber)

        reader_taum_PZ.AddVariable("df_m_1", arr_m_1)
        reader_taum_PZ.AddVariable("df_m_2", arr_m_2)
        reader_taum_PZ.AddVariable("df_m_3", arr_m_3)
        reader_taum_PZ.AddVariable("df_m_4", arr_m_4)
        reader_taum_PZ.AddVariable("df_m_5", arr_m_5)
        reader_taum_PZ.AddVariable("df_m_6", arr_m_6)
        reader_taum_PZ.AddVariable("df_m_7", arr_m_7)
        reader_taum_PZ.AddVariable("df_m_8", arr_m_8)
        reader_taum_PZ.AddVariable("df_m_9", arr_m_9)
        reader_taum_PZ.AddVariable("df_m_10", arr_m_10)
        reader_taum_PZ.AddVariable("df_m_11", arr_m_11)
        reader_taum_PZ.AddVariable("df_m_12", arr_m_12)
        reader_taum_PZ.AddVariable("df_m_13", arr_m_13)
        reader_taum_PZ.AddVariable("df_m_14", arr_m_14)
        reader_taum_PZ.AddVariable("df_m_15", arr_m_15)
        reader_taum_PZ.AddVariable("df_m_16", arr_m_16)
        reader_taum_PZ.AddVariable("df_m_17", arr_m_17)
        reader_taum_PZ.AddVariable("df_m_18", arr_m_18)
        reader_taum_PZ.AddVariable("df_m_19", arr_m_19)
        reader_taum_PZ.AddVariable("df_m_20", arr_m_20)
        reader_taum_PZ.AddVariable("df_m_21", arr_m_21)
        reader_taum_PZ.AddVariable("df_m_22", arr_m_22)
        reader_taum_PZ.AddVariable("Kp_RP_Z", arr_Kp_RP_Z)
        reader_taum_PZ.AddSpectator("eventNumber", arr_eventNumber)

    reader_taup_PX.BookMVA("MLP", weightfile_taup_PX)
    reader_taup_PY.BookMVA("MLP", weightfile_taup_PY)
    reader_taup_PZ.BookMVA("MLP", weightfile_taup_PZ)
    reader_taum_PX.BookMVA("MLP", weightfile_taum_PX)
    reader_taum_PY.BookMVA("MLP", weightfile_taum_PY)
    reader_taum_PZ.BookMVA("MLP", weightfile_taum_PZ)

    m = np.zeros(dimM)
    V = np.zeros((dimM,dimM))

    Big_start = time.time()
    num_entries = t.GetEntries()

    for evt in range(num_entries):

        start = time.time()
        print("<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<", " evt = ", evt, "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<")
        global nIter, status

        entry = t.GetEntry(evt) 

        # Get vector m and covariance matrix V from TTree
        for i in range(dimM):

            if isD0D0K:
                m[i] = getattr(t, "df_mprime_{0}".format(i+1))
            else:
                m[i] = getattr(t, "df_m_{0}".format(i+1)) 

            for j in range(dimM):
                if isD0D0K:
                    if(j >= i):
                        V[i][j] = getattr(t, "df_Vprime_{0}_{1}".format(i+1,j+1)) 
                    else:
                        V[i][j] = getattr(t, "df_Vprime_{0}_{1}".format(j+1,i+1)) 
                else:
                    if(j >= i):
                        V[i][j] = getattr(t, "df_V_{0}_{1}".format(i+1,j+1))
                    else:
                        V[i][j] = getattr(t, "df_V_{0}_{1}".format(j+1,i+1))

        # true measured values
        # m[0] = getattr(t, "Bp_TRUEORIGINVERTEX_X")
        # m[1] = getattr(t, "Bp_TRUEORIGINVERTEX_Y")
        # m[2] = getattr(t, "Bp_TRUEORIGINVERTEX_Z")
        # m[3] = getattr(t, "taup_TRUEENDVERTEX_X")
        # m[4] = getattr(t, "taup_TRUEENDVERTEX_Y")
        # m[5] = getattr(t, "taup_TRUEENDVERTEX_Z")
        # pi1x = getattr(t, "taup_pi1_TRUEP_X")
        # pi1y = getattr(t, "taup_pi1_TRUEP_Y")
        # pi1z = getattr(t, "taup_pi1_TRUEP_Z")
        # E1 = getattr(t, "taup_pi1_TRUEP_E")
        # pi2x = getattr(t, "taup_pi2_TRUEP_X")
        # pi2y = getattr(t, "taup_pi2_TRUEP_Y")
        # pi2z = getattr(t, "taup_pi2_TRUEP_Z")
        # E2 = getattr(t, "taup_pi2_TRUEP_E")
        # pi3x = getattr(t, "taup_pi3_TRUEP_X")
        # pi3y = getattr(t, "taup_pi3_TRUEP_Y")
        # pi3z = getattr(t, "taup_pi3_TRUEP_Z")
        # E3 = getattr(t, "taup_pi3_TRUEP_E")
        # m[6] = pi1x + pi2x + pi3x
        # m[7] = pi1y + pi2y + pi3y
        # m[8] = pi1z + pi2z + pi3z
        # m[9] = E1 + E2 + E3
        # m[10] = getattr(t, "taum_TRUEENDVERTEX_X")
        # m[11] = getattr(t, "taum_TRUEENDVERTEX_Y")
        # m[12] = getattr(t, "taum_TRUEENDVERTEX_Z")
        # pi4x = getattr(t, "taum_pi1_TRUEP_X")
        # pi4y = getattr(t, "taum_pi1_TRUEP_Y")
        # pi4z = getattr(t, "taum_pi1_TRUEP_Z")
        # E4 = getattr(t, "taum_pi1_TRUEP_E")
        # pi5x = getattr(t, "taum_pi2_TRUEP_X")
        # pi5y = getattr(t, "taum_pi2_TRUEP_Y")
        # pi5z = getattr(t, "taum_pi2_TRUEP_Z")
        # E5 = getattr(t, "taum_pi2_TRUEP_E")
        # pi6x = getattr(t, "taum_pi3_TRUEP_X")
        # pi6y = getattr(t, "taum_pi3_TRUEP_Y")
        # pi6z = getattr(t, "taum_pi3_TRUEP_Z")
        # E6 = getattr(t, "taum_pi3_TRUEP_E")
        # m[13] = pi4x + pi5x + pi6x
        # m[14] = pi4y + pi5y + pi6y
        # m[15] = pi4z + pi5z + pi6z
        # m[16] = E4 + E5 + E6
        # m[17] = getattr(t, "Bp_TRUEENDVERTEX_X")
        # m[18] = getattr(t, "Bp_TRUEENDVERTEX_Y")
        # m[19] = getattr(t, "Kp_TRUEP_X")
        # m[20] = getattr(t, "Kp_TRUEP_Y")
        # m[21] = getattr(t, "Kp_TRUEP_Z")
        
        if((species == 10) or (species == 11) or (species == 12) or (species == 1) or (species == 4) or (species == 9)):
            global x_true
            for i in range(len(x_true_names)):
                x_true.append( getattr(t, x_true_names[i]) )
        
        RPz = getattr(t, "Kp_RP_Z")
        # RPz = getattr(t, "Bp_TRUEENDVERTEX_Z")
        eventNumber = getattr(t, "eventNumber")

        # Invert matrix V
        W = np.linalg.inv(V)
        # print(np.dot(W,V))
        # print(np.linalg.eig(V))

        BVx = getattr(t, "Bp_ENDVERTEX_X")
        BVy = getattr(t, "Bp_ENDVERTEX_Y")
        BVz = getattr(t, "Bp_ENDVERTEX_Z")
        # global BV_offline 
        BV_offline = ROOT.Math.XYZPoint(BVx, BVy, BVz)

        params = [m, W, RPz, eventNumber]

        global change_L
        change_L = 0

        # Single initialisation
        # init[0] = 14
        # run_solver(init[0], year, species, line, BV_offline, params, True, max_iter=10000, eps=0.000001)
        # if(status != 0):
        #     run_solver(init[0], year, species, line, BV_offline, params, False, max_iter=10000, eps=0.000001)
        # if(status != 0):
        #     change_L = 1
        #     run_solver(init[0], year, species, line, BV_offline, params, True, max_iter=10000, eps=0.000001)
        # if(status != 0):
        #     run_solver(init[0], year, species, line, BV_offline, params, False, max_iter=10000, eps=0.000001)
        # if(status != 0):
        #     change_L = 2
        #     run_solver(init[0], year, species, line, BV_offline, params, True, max_iter=10000, eps=0.000001)
        # if(status != 0):
        #     run_solver(init[0], year, species, line, BV_offline, params, False, max_iter=10000, eps=0.000001)

        # Lowest chi2 (DEFAULT)
        lowest_chi2(year, species, line, BV_offline, params, [0,1,2,3,4], max_iter=10000, eps=0.000001) 

        for i in range(N_local_minima[0]):
            M_local[i] =  M_local_minima[i]
            Chi2_local[i] = Chi2_local_minima[i]

        U = U_cov(x,params,V)

        for i in range(dimM+dimX):
            X[i][0] = x[i]
            X_ERR[i][0] = np.sqrt(U[i][i])
        for i in range(dimM+dimX+dimC):
            F[i][0] = f[i]

        chi2[0] =  chi2_value
        tolerance[0] = tol_value
        STATUS[0] = status
        NITER[0] = nIter

        mass_squared = x[dimM+6]
        if(mass_squared > 0):
            MB[0] = np.sqrt( mass_squared )
        else:
            MB[0] = - np.sqrt( -mass_squared )

        dMB_squared = np.sqrt(U[dimM+6][dimM+6])
        dMB[0] =  dMB_squared/(2*np.abs(MB[0]))

        print("init = ", init[0], "status = ", STATUS[0], "MB = ", MB[0], "dMB = ", dMB[0], "chi2 = ", chi2[0], "sum = ", tolerance[0], "N_local = ", N_local_minima[0])

        ######################################################################
        # MB0[0] = run_solver(0, year, species, line, BV_offline, params, True, max_iter=10000, eps=0.000001)
        # MB1[0] = run_solver(1, year, species, line, BV_offline, params, True, max_iter=10000, eps=0.000001)
        # MB2[0] = run_solver(2, year, species, line, BV_offline, params, True, max_iter=10000, eps=0.000001)
        # MB3[0] = run_solver(3, year, species, line, BV_offline, params, True, max_iter=10000, eps=0.000001)
        # MB4[0] = run_solver(4, year, species, line, BV_offline, params, True, max_iter=10000, eps=0.000001)

        end = time.time()
        theTime[0] = end-start
        print("time = ", theTime[0], " s")

        tree.Fill()
        x_true.clear()
        M_local_minima.clear()
        Chi2_local_minima.clear()
        nIter = 0

    Big_end = time.time()
    print("Run in {0} s".format(Big_end-Big_start))

    file.cd()
    tree.Write()
    file.Close()

if __name__ == "__main__":
    main(sys.argv)