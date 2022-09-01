import ROOT
import math
from math import pow
from array import array

useSmeared = True
def main():
    fileIn = ROOT.TFile.Open("fromMaria/B2Ktautau_tree.root")
    treeIn = fileIn.DecayTree

    nEntries = treeIn.GetEntries()

    print(f'treeIn has {nEntries} entries')

    #make the output file
    outFileName = "output1.root"
    if useSmeared:
        outFileName = "output_smeared.root"
    
    fileOut = ROOT.TFile(outFileName, 'RECREATE')
    treeOut = ROOT.TTree("DecayTree", "")
    
    B_M    = array('f', [0.])
    Tau1_M = array('f', [0.])
    Tau2_M = array('f', [0.])
    Tau1_P = array('f', [0.])
    Tau2_P = array('f', [0.])
    Tau1_P_orig = array('f', [0.])
    Tau2_P_orig = array('f', [0.])
    missing_mass_tau1 = array('f', [0.])
    missing_mass_tau2 = array('f', [0.])

    treeOut.Branch("B_M", B_M, "B_M/f")
    treeOut.Branch("Tau1_M", Tau1_M, "Tau1_M/f")
    treeOut.Branch("Tau2_M", Tau2_M, "Tau2_M/f")
    treeOut.Branch("Tau1_P", Tau1_P, "Tau1_P/f")
    treeOut.Branch("Tau2_P", Tau2_P, "Tau2_P/f")
    treeOut.Branch("Tau1_P_orig", Tau1_P_orig, "Tau1_P_orig/f")
    treeOut.Branch("Tau2_P_orig", Tau2_P_orig, "Tau2_P_orig/f")
    treeOut.Branch("Tau1_missingM", missing_mass_tau1, "Tau1_missingM/f")
    treeOut.Branch("Tau2_missingM", missing_mass_tau2, "Tau2_missingM/f")

    for iLoop in range(1000):
        if(iLoop%100 == 0): 
            print(f"iLoop = {iLoop}")
        
        treeIn.GetEntry(iLoop)

        PV  = ROOT.Math.XYZPoint(treeIn.Bp_0_origX_TRUE, treeIn.Bp_0_origY_TRUE, treeIn.Bp_0_origZ_TRUE)
        DV1 = ROOT.Math.XYZPoint(treeIn.taup_0_vtxX_TRUE, treeIn.taup_0_vtxY_TRUE, treeIn.taup_0_vtxZ_TRUE)
        DV2 = ROOT.Math.XYZPoint(treeIn.taum_0_vtxX_TRUE, treeIn.taum_0_vtxY_TRUE, treeIn.taum_0_vtxZ_TRUE)
        p_K = ROOT.Math.XYZVector(treeIn.Kp_0_PX_TRUE, treeIn.Kp_0_PY_TRUE, treeIn.Kp_0_PZ_TRUE)
        #for the moment, we cheat and use the true B decay vertex. 
        #later in the real data, have to find a way of extracting a track hit associated with the K+
        BV  = ROOT.Math.XYZPoint(treeIn.Bp_0_vtxX_TRUE, treeIn.Bp_0_vtxY_TRUE, treeIn.Bp_0_vtxZ_TRUE)
            
        #3pi variables
        p_tau1_pi1   = ROOT.Math.XYZVector(treeIn.pip_0_PX_TRUE, treeIn.pip_0_PY_TRUE, treeIn.pip_0_PZ_TRUE)
        p_tau1_pi2   = ROOT.Math.XYZVector(treeIn.pim_0_PX_TRUE, treeIn.pim_0_PY_TRUE, treeIn.pim_0_PZ_TRUE)
        p_tau1_pi3   = ROOT.Math.XYZVector(treeIn.pip_1_PX_TRUE, treeIn.pip_1_PY_TRUE, treeIn.pip_1_PZ_TRUE)

        p_tau2_pi1   = ROOT.Math.XYZVector(treeIn.pip_2_PX_TRUE, treeIn.pip_2_PY_TRUE, treeIn.pip_2_PZ_TRUE)
        p_tau2_pi2   = ROOT.Math.XYZVector(treeIn.pim_1_PX_TRUE, treeIn.pim_1_PY_TRUE, treeIn.pim_1_PZ_TRUE)
        p_tau2_pi3   = ROOT.Math.XYZVector(treeIn.pim_2_PX_TRUE, treeIn.pim_2_PY_TRUE, treeIn.pim_2_PZ_TRUE)

        p_3pi1_x = treeIn.pip_0_PX_TRUE + treeIn.pim_0_PX_TRUE + treeIn.pip_1_PX_TRUE
        p_3pi1_y = treeIn.pip_0_PY_TRUE + treeIn.pim_0_PY_TRUE + treeIn.pip_1_PY_TRUE
        p_3pi1_z = treeIn.pip_0_PZ_TRUE + treeIn.pim_0_PZ_TRUE + treeIn.pip_1_PZ_TRUE

        p_3pi2_x = treeIn.pip_2_PX_TRUE + treeIn.pim_1_PX_TRUE + treeIn.pim_2_PX_TRUE
        p_3pi2_y = treeIn.pip_2_PY_TRUE + treeIn.pim_1_PY_TRUE + treeIn.pim_2_PY_TRUE
        p_3pi2_z = treeIn.pip_2_PZ_TRUE + treeIn.pim_1_PZ_TRUE + treeIn.pim_2_PZ_TRUE

        E_3pi1   = treeIn.pip_0_E_TRUE + treeIn.pim_0_E_TRUE + treeIn.pip_1_E_TRUE
        E_3pi2   = treeIn.pip_2_E_TRUE + treeIn.pim_1_E_TRUE + treeIn.pim_2_E_TRUE

        if useSmeared:
            PV = ROOT.Math.XYZPoint(treeIn.Bp_0_origX, treeIn.Bp_0_origY, treeIn.Bp_0_origZ)
            DV1 = ROOT.Math.XYZPoint(treeIn.taup_0_vtxX, treeIn.taup_0_vtxY, treeIn.taup_0_vtxZ)
            DV2 = ROOT.Math.XYZPoint(treeIn.taum_0_vtxX, treeIn.taum_0_vtxY, treeIn.taum_0_vtxZ)
            BV  = ROOT.Math.XYZPoint(treeIn.Bp_0_vtxX, treeIn.Bp_0_vtxY, treeIn.Bp_0_vtxZ)

            p_tau1_pi1   = ROOT.Math.XYZVector(treeIn.pip_0_PX, treeIn.pip_0_PY, treeIn.pip_0_PZ)
            p_tau1_pi2   = ROOT.Math.XYZVector(treeIn.pim_0_PX, treeIn.pim_0_PY, treeIn.pim_0_PZ)
            p_tau1_pi3   = ROOT.Math.XYZVector(treeIn.pip_1_PX, treeIn.pip_1_PY, treeIn.pip_1_PZ)

            p_tau2_pi1   = ROOT.Math.XYZVector(treeIn.pip_2_PX, treeIn.pip_2_PY, treeIn.pip_2_PZ)
            p_tau2_pi2   = ROOT.Math.XYZVector(treeIn.pim_1_PX, treeIn.pim_1_PY, treeIn.pim_1_PZ)
            p_tau2_pi3   = ROOT.Math.XYZVector(treeIn.pim_2_PX, treeIn.pim_2_PY, treeIn.pim_2_PZ)

            p_3pi1_x = treeIn.pip_0_PX + treeIn.pim_0_PX + treeIn.pip_1_PX
            p_3pi1_y = treeIn.pip_0_PY + treeIn.pim_0_PY + treeIn.pip_1_PY
            p_3pi1_z = treeIn.pip_0_PZ + treeIn.pim_0_PZ + treeIn.pip_1_PZ

            p_3pi2_x = treeIn.pip_2_PX + treeIn.pim_1_PX + treeIn.pim_2_PX
            p_3pi2_y = treeIn.pip_2_PY + treeIn.pim_1_PY + treeIn.pim_2_PY
            p_3pi2_z = treeIn.pip_2_PZ + treeIn.pim_1_PZ + treeIn.pim_2_PZ

            E_3pi1   = treeIn.pip_0_E + treeIn.pim_0_E + treeIn.pip_1_E
            E_3pi2   = treeIn.pip_2_E + treeIn.pim_1_E + treeIn.pim_2_E

        m_3pi1   = math.sqrt(pow(E_3pi1,2) - pow(p_3pi1_x,2) - pow(p_3pi1_y,2) - pow(p_3pi1_z,2))
        m_3pi2   = math.sqrt(pow(E_3pi2,2) - pow(p_3pi2_x,2) - pow(p_3pi2_y,2) - pow(p_3pi2_z,2))

        #Apply the transformation. _t denotes the transformed versions of the vectors
        PV_t = makeTransformation(p_K, BV, PV, False)
        DV1_t = makeTransformation(p_K, BV, DV1, False)
        DV2_t = makeTransformation(p_K, BV, DV2, False)
        BV_t  = makeTransformation(p_K, BV, BV, False)

        p_K_t = makeTransformation(p_K, BV, p_K, False)

        p_tau1_pi1_t = makeTransformation(p_K, BV, p_tau1_pi1, False)
        p_tau1_pi2_t = makeTransformation(p_K, BV, p_tau1_pi2, False)
        p_tau1_pi3_t = makeTransformation(p_K, BV, p_tau1_pi3, False)

        p_tau2_pi1_t = makeTransformation(p_K, BV, p_tau2_pi1, False)
        p_tau2_pi2_t = makeTransformation(p_K, BV, p_tau2_pi2, False)
        p_tau2_pi3_t = makeTransformation(p_K, BV, p_tau2_pi3, False)

        p_3pi1_x_t = p_tau1_pi1_t.x() + p_tau1_pi2_t.x() + p_tau1_pi3_t.x()
        p_3pi1_y_t = p_tau1_pi1_t.y() + p_tau1_pi2_t.y() + p_tau1_pi3_t.y()
        p_3pi1_z_t = p_tau1_pi1_t.z() + p_tau1_pi2_t.z() + p_tau1_pi3_t.z()

        p_3pi2_x_t = p_tau2_pi1_t.x() + p_tau2_pi2_t.x() + p_tau2_pi3_t.x()
        p_3pi2_y_t = p_tau2_pi1_t.y() + p_tau2_pi2_t.y() + p_tau2_pi3_t.y()
        p_3pi2_z_t = p_tau2_pi1_t.z() + p_tau2_pi2_t.z() + p_tau2_pi3_t.z()

        m_3pi1_t   = math.sqrt(pow(E_3pi1,2) - pow(p_3pi1_x_t,2) - pow(p_3pi1_y_t,2) - pow(p_3pi1_z_t,2))
        m_3pi2_t   = math.sqrt(pow(E_3pi2,2) - pow(p_3pi2_x_t,2) - pow(p_3pi2_y_t,2) - pow(p_3pi2_z_t,2))

        # print(f"PV = ({PV.x():.2f},{PV.y():.2f},{PV.z():.2f}) --> ({PV_t.x():.2f},{PV_t.y():.2f},{PV_t.z():.2f})")
        # print(f"BV = ({BV.x():.2f},{BV.y():.2f},{BV.z():.2f}) --> ({BV_t.x():.2f},{BV_t.y():.2f},{BV_t.z():.2f})")
        # print(f"p_K = ({p_K.x():.2f},{p_K.y():.2f},{p_K.z():.2f}) --> ({p_K_t.x():.2f},{p_K_t.y():.2f},{p_K_t.z():.2f})")

        a1 = DV1_t.y()/DV1_t.x()
        a2 = DV2_t.y()/DV2_t.x()
        b  = (PV_t.y() - a1*PV_t.x())/(a2*PV_t.x() - PV_t.y())

        c1 = b * (DV2_t.z() - DV1_t.z())/DV2_t.x()
        c2 = b * DV1_t.x()/DV2_t.x()

        d1 = ((DV1_t.z() - PV_t.z())*(b + 1) + c1*PV_t.x())/(DV1_t.x()*(b + 1) - PV_t.x()*(1 + c2))
        #d2 = (p_K_t.r() *PV_t.x())/(DV1_t.x() * (b + 1) - PV_t.x()*(1 + c2))
        d2 = (p_K_t.z() *PV_t.x())/(DV1_t.x() * (b + 1) - PV_t.x()*(1 + c2))

        e1 = c1 + (c2 * d1)
        e2 = c2*d2

        #Now getting to the quadratic equations
        #Let X denote p_tau1_x
        #g*X^2 + h*X + i = 0
        #m*X^2 + n*X + o = 0

        mTau = 1.77686 #RapidSim uses GeV for the energy/momentum unit

        g = 1 + pow(a1,2) + pow(d1,2) - pow((p_3pi1_x_t + a1*p_3pi1_y_t + d1*p_3pi1_z_t)/E_3pi1, 2)
        h = 2*d1*d2 - (((pow(mTau,2) + pow(m_3pi1_t,2) + 2*d2*p_3pi1_z_t)*(p_3pi1_x_t + a1*p_3pi1_y_t + d1*p_3pi1_z_t))/(pow(E_3pi1,2)))
        i = pow(mTau, 2) + pow(d2, 2) - pow((pow(mTau, 2) + pow(m_3pi1_t, 2) + 2*d2*p_3pi1_z_t)/(2*E_3pi1), 2)

        m = (1 + pow(a2, 2))*pow(b, 2) + pow(e1, 2) - pow(((b*p_3pi2_x_t + a2*b*p_3pi2_y_t + e1*p_3pi2_z_t)/E_3pi2), 2)
        n = 2*e1*e2 - (((pow(mTau, 2) + pow(m_3pi2_t, 2) + 2*e2*p_3pi2_z_t)*(b*p_3pi2_x_t + a2*b*p_3pi2_y_t + e1*p_3pi2_z_t))/(pow(E_3pi2, 2)))
        o = pow(mTau, 2) + pow(e2, 2) - pow((pow(mTau, 2) + pow(m_3pi2_t, 2) + 2*e2*p_3pi2_z_t)/(2*E_3pi2), 2)

        disc1 = pow(h, 2) - 4*g*i
        disc2 = pow(n, 2) - 4*m*o

        sol1_eq1 = 0
        sol2_eq1 = 0
        sol1_eq2 = 0
        sol2_eq2 = 0

        if disc1 < 0:
            print(f"disc1 = {disc1:.2f}, -ve discriminant. Skipping to next iteration")
            print(f"pow(h, 2) = {pow(h, 2)}, 4gi = {4*g*i}")
            print(f"disc1/term1 = {disc1/pow(h, 2):.4f}, disc1/term2 = {disc1/(4*g*i):.4f}")
            print(f"iLoop = {iLoop}")
            print(f"disc1 = {disc1:.2f}")
            print(f"PV = ({PV.x():.2f},{PV.y():.2f},{PV.z():.2f}) --> ({PV_t.x():.2f},{PV_t.y():.2f},{PV_t.z():.2f})")
            print(f"BV = ({BV.x():.2f},{BV.y():.2f},{BV.z():.2f}) --> ({BV_t.x():.2f},{BV_t.y():.2f},{BV_t.z():.2f})")
            print(f"DV1 = ({DV1.x():.5f},{DV1.y():.5f},{DV1.z():.5f}) --> ({DV1_t.x():.5f},{DV1_t.y():.5f},{DV1_t.z():.5f})")
            print(f"DV2 = ({DV2.x():.5f},{DV2.y():.5f},{DV2.z():.5f}) --> ({DV2_t.x():.5f},{DV2_t.y():.5f},{DV2_t.z():.5f})")
            print(f"p_K = ({p_K.x():.2f},{p_K.y():.2f},{p_K.z():.2f}) --> ({p_K_t.x():.2f},{p_K_t.y():.2f},{p_K_t.z():.2f})")
            break
            #continue
        else:
            sol1_eq1 = (-h + math.sqrt(disc1))/(2*g)
            sol2_eq1 = (-h - math.sqrt(disc1))/(2*g)

        if disc2 < 0:
            print(f"disc2 = {disc2:.2f}, -ve discriminant. Skipping to next iteration")
            continue
        else:
            sol1_eq2 = (-1*n + math.sqrt(disc2))/(2*m)
            sol2_eq2 = (-1*n - math.sqrt(disc2))/(2*m)

        #print(f"disc1 = {disc1:.3f}, disc2 = {disc2:.3f}, sol1_eq1 = {sol1_eq1:.3f}, sol2_eq1 = {sol2_eq1:.3f}, sol1_eq2 = {sol1_eq2:.3f}, sol2_eq2 = {sol2_eq2:.3f}")

        #Match solutions for p_tau1_x from the two equations
        theSol = 0
        myTolerance = 0.001 #0.1%

        if math.fabs((sol1_eq2 - sol1_eq1)/sol1_eq1) < myTolerance:
            theSol = sol1_eq1
        elif math.fabs((sol2_eq2 - sol1_eq1)/sol1_eq1) < myTolerance:
            theSol = sol1_eq1
        elif math.fabs((sol1_eq2 - sol2_eq1)/sol2_eq1) < myTolerance:
            theSol = sol2_eq1
        elif math.fabs((sol2_eq2 - sol2_eq1)/sol2_eq1) < myTolerance:
            theSol = sol2_eq1
        else:
            print("None of the match pairs have a tolerance below 0.1%")
        
        #print(f"theSol = {theSol}")
        p_tau1_x = theSol
        p_tau2_x = b * p_tau1_x

        p_tau1_y = a1 * p_tau1_x
        p_tau2_y = a2 * p_tau2_x

        p_tau1_z = (d1 * p_tau1_x) + d2
        p_tau2_z = (e1 * p_tau1_x) + e2
        
        #New recalculated p_tau 3-vectors
        p_tau1 = ROOT.Math.XYZVector(p_tau1_x, p_tau1_y, p_tau1_z)
        p_tau2 = ROOT.Math.XYZVector(p_tau2_x, p_tau2_y, p_tau2_z)

        E_tau1 = math.sqrt(pow(mTau, 2) + p_tau1.Mag2())
        E_tau2 = math.sqrt(pow(mTau, 2) + p_tau2.Mag2())

        #make 4-vectors for K, tau1, tau2, and add them up to form the B
        P_K    = ROOT.Math.PxPyPzEVector(p_K_t.x(), p_K_t.y(), p_K_t.z(), treeIn.Kp_0_E_TRUE)
        P_tau1 = ROOT.Math.PxPyPzEVector(p_tau1_x, p_tau1_y, p_tau1_z, E_tau1)
        P_tau2 = ROOT.Math.PxPyPzEVector(p_tau2_x, p_tau2_y, p_tau2_z, E_tau2)

        P_B = P_K + P_tau1 + P_tau2

        B_M[0]    = P_B.M()
        Tau1_M[0] = P_tau1.M()
        Tau2_M[0] = P_tau2.M()

        #Now we want to apply the inverse transformation to the p_tau 3-vectors to take them back to the LHCb reference frame
        p_tau1_lhcb = makeTransformation(p_K, BV, p_tau1, True)
        p_tau2_lhcb = makeTransformation(p_K, BV, p_tau2, True)
        PV_lhcb     = makeTransformation(p_K, BV, PV_t, True)
        BV_lhcb     = makeTransformation(p_K, BV, BV_t, True)

        E_tau1_lhcb = math.sqrt(p_tau1_lhcb.Mag2() + pow(mTau, 2))
        E_tau2_lhcb = math.sqrt(p_tau2_lhcb.Mag2() + pow(mTau, 2))

        Tau1_P_orig[0] = math.sqrt(pow(p_3pi1_x ,2) + pow(p_3pi1_y ,2) + pow(p_3pi1_z ,2))
        Tau2_P_orig[0] = math.sqrt(pow(p_3pi2_x ,2) + pow(p_3pi2_y ,2) + pow(p_3pi2_z ,2))

        Tau1_P[0] = P_tau1.r()
        Tau2_P[0] = P_tau2.r()

        #missing mass
        missing_mass_sq_tau1 = pow(E_tau1_lhcb - E_3pi1, 2) - pow(p_tau1_lhcb.x() - p_3pi1_x, 2) - pow(p_tau1_lhcb.y() - p_3pi1_y, 2) - pow(p_tau1_lhcb.z() - p_3pi1_z, 2)
        missing_mass_sq_tau2 = pow(E_tau2_lhcb - E_3pi2, 2) - pow(p_tau2_lhcb.x() - p_3pi2_x, 2) - pow(p_tau2_lhcb.y() - p_3pi2_y, 2) - pow(p_tau2_lhcb.z() - p_3pi2_z, 2)

        #sometimes because of floating point precision, the missing mass squared becomes a very small negative number (needs some more thought)
        if missing_mass_sq_tau1 < 0:
            missing_mass_sq_tau1 *= -1
        if missing_mass_sq_tau2 < 0:
            missing_mass_sq_tau2 *= -1
        
        missing_mass_tau1[0] = math.sqrt(missing_mass_sq_tau1)
        missing_mass_tau2[0] = math.sqrt(missing_mass_sq_tau2)

        treeOut.Fill()
        
        #check to make sure inverse transform is working c
        # print(f"PV = ({PV.x():.2f},{PV.y():.2f},{PV.z():.2f}) --> ({PV_t.x():.2f},{PV_t.y():.2f},{PV_t.z():.2f}) --> ({PV_lhcb.x():.2f},{PV_lhcb.y():.2f},{PV_lhcb.z():.2f})")
        # print(f"BV = ({BV.x():.2f},{BV.y():.2f},{BV.z():.2f}) --> ({BV_t.x():.2f},{BV_t.y():.2f},{BV_t.z():.2f}) --> ({BV_lhcb.x():.2f},{BV_lhcb.y():.2f},{BV_lhcb.z():.2f})")

        # if abs(B_M[0]*1000-5279.2968) > 50:
        #     print(f"B_M = {B_M[0]*1000:.2f}, disc1 = {disc1:.3f}, disc2 = {disc2:.3f}, sol1_eq1 = {sol1_eq1:.3f}, sol2_eq1 = {sol2_eq1:.3f}, sol1_eq2 = {sol1_eq2:.3f}, sol2_eq2 = {sol2_eq2:.3f}, theSol = {theSol}")
        #     #print(f"B_M = {B_M[0]*1000:.2f}, Tau1_M = {Tau1_M[0]*1000:.2f}, Tau2_M = {Tau2_M[0]*1000:.2f}")
        # #print(f"m_B = {P_B.M()}")

    fileOut.cd()
    treeOut.Write("", ROOT.TObject.kOverwrite)
    fileOut.Close()
    fileIn.Close()

    
#define a function that does the transformation on a given vector
#original coordinate system has PV at origin

def makeTransformation(p_K, BV, theVector, invFlag): #theVec is the vector that needs to be transformed
    #we want to have a translation and two rotations
    #ROOT::Math::Translation3D only affects point vectors, not XYZ vectors.

    deltaX = BV.x()
    deltaY = BV.y()
    deltaZ = BV.z()

    myShift = ROOT.Math.Translation3D(-deltaX, -deltaY, -deltaZ) #ROOT::Math::Translation3D translates the point itself, not the coordinate system
    
    #Using the Euler angle formalism to define the rotations
    #See https://mathworld.wolfram.com/EulerAngles.html
    #https://root.cern.ch/doc/master/classROOT_1_1Math_1_1EulerAngles.html

    #Rotation about original Z axis to bring p_K into YZ plane. If Y component is +ve, rotation is clockwise, else anti-clockwise
    phi   = -1 * math.atan(p_K.x()/p_K.y())   
    #Clockwise rotation about new X axis to align Z axis with p_K
    theta = -1 * math.atan(p_K.Rho() * (p_K.y()/math.fabs(p_K.y())) /p_K.z()) 

    myRotation = ROOT.Math.EulerAngles(phi, theta, 0)

    #After this combined transformation, the BV will be at the origin. 
    #Will this be a problem?

    myTransformation = ROOT.Math.Transform3D(myRotation) * ROOT.Math.Transform3D(myShift)
    myTransformation_inverse = myTransformation.Inverse()

    if invFlag:
        theVector_t = myTransformation_inverse(theVector)
    else:
        theVector_t = myTransformation(theVector)

    return theVector_t
    #return myRotation(myShift(theVector))

#def makeTransformation_inverse(p_K, PV, theVector)

if __name__ == "__main__":
    main()