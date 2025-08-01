import sys
import numpy as np
import pandas as pd
from uncertainties import ufloat
import ROOT


def DD_cut(name, mother_ID):
    # Ds+ Ds-
    Ds_Ds = f"(abs(taup_TRUEID) == 431) && (abs(taum_TRUEID) == 431) && (abs(taup_MC_MOTHER_ID)=={mother_ID}) && (abs(taum_MC_MOTHER_ID)=={mother_ID})"
    Dsstar_Ds = f"((taup_TRUEID == 431) && (taup_MC_MOTHER_ID == 433) &&  (taup_MC_GD_MOTHER_ID == {mother_ID}) && (taum_TRUEID == -431) && (taum_MC_MOTHER_ID == {mother_ID})) || ((taum_TRUEID == 431) && (taum_MC_MOTHER_ID == 433) &&  (taum_MC_GD_MOTHER_ID == {mother_ID}) && (taup_TRUEID == -431) && (taup_MC_MOTHER_ID == {mother_ID})) || ((taup_TRUEID == -431) && (taup_MC_MOTHER_ID == -433) &&  (taup_MC_GD_MOTHER_ID == -{mother_ID}) && (taum_TRUEID == 431) && (taum_MC_MOTHER_ID == -{mother_ID})) || ((taum_TRUEID == -431) && (taum_MC_MOTHER_ID == -433) &&  (taum_MC_GD_MOTHER_ID == -{mother_ID}) && (taup_TRUEID == 431) && (taup_MC_MOTHER_ID == -{mother_ID}))"
    Ds_Dsstar = f"((taup_TRUEID == 431) && (taup_MC_MOTHER_ID == {mother_ID}) && (taum_TRUEID == -431) && (taum_MC_MOTHER_ID == -433) && (taum_MC_GD_MOTHER_ID == {mother_ID})) || ((taum_TRUEID == 431) && (taum_MC_MOTHER_ID == {mother_ID}) && (taup_TRUEID == -431) && (taup_MC_MOTHER_ID == -433) && (taup_MC_GD_MOTHER_ID == {mother_ID})) || ((taup_TRUEID == -431) && (taup_MC_MOTHER_ID == -{mother_ID}) && (taum_TRUEID == 431) && (taum_MC_MOTHER_ID == 433) && (taum_MC_GD_MOTHER_ID == -{mother_ID})) || ((taum_TRUEID == -431) && (taum_MC_MOTHER_ID == -{mother_ID}) && (taup_TRUEID == 431) && (taup_MC_MOTHER_ID == 433) && (taup_MC_GD_MOTHER_ID == -{mother_ID}))"
    Dsstar_Dsstar = f"(abs(taup_TRUEID) == 431) && (abs(taum_TRUEID) == 431) && (abs(taup_MC_MOTHER_ID)==433) && (abs(taup_MC_GD_MOTHER_ID)=={mother_ID}) && (abs(taum_MC_MOTHER_ID)==433) && (abs(taum_MC_GD_MOTHER_ID)=={mother_ID})"

    # D0bar D0
    D0_D0 = f"(abs(taup_TRUEID) == 421) && (abs(taum_TRUEID) == 421) && (abs(taup_MC_MOTHER_ID)=={mother_ID}) && (abs(taum_MC_MOTHER_ID)=={mother_ID})"
    D0star_D0 = f"((taup_TRUEID == -421) && (taup_MC_MOTHER_ID == -423) && (taup_MC_GD_MOTHER_ID == {mother_ID}) && (taum_TRUEID == 421) && (taum_MC_MOTHER_ID == {mother_ID})) || ((taum_TRUEID == -421) && (taum_MC_MOTHER_ID == -423) && (taum_MC_GD_MOTHER_ID == {mother_ID}) && (taup_TRUEID == 421) && (taup_MC_MOTHER_ID == {mother_ID})) || ((taup_TRUEID == 421) && (taup_MC_MOTHER_ID == 423) && (taup_MC_GD_MOTHER_ID == -{mother_ID}) && (taum_TRUEID == -421) && (taum_MC_MOTHER_ID == -{mother_ID})) || ((taum_TRUEID == 421) && (taum_MC_MOTHER_ID == 423) && (taum_MC_GD_MOTHER_ID == -{mother_ID}) && (taup_TRUEID == -421) && (taup_MC_MOTHER_ID == -{mother_ID}))"
    D0_D0star = f"((taup_TRUEID == -421) && (taup_MC_MOTHER_ID == {mother_ID}) && (taum_TRUEID == 421) && (taum_MC_MOTHER_ID == 423) && (taum_MC_GD_MOTHER_ID == {mother_ID})) || ((taum_TRUEID == -421) && (taum_MC_MOTHER_ID == {mother_ID}) && (taup_TRUEID == 421) && (taup_MC_MOTHER_ID == 423) && (taup_MC_GD_MOTHER_ID == {mother_ID})) || ((taup_TRUEID == 421) && (taup_MC_MOTHER_ID == -{mother_ID}) && (taum_TRUEID == -421) && (taum_MC_MOTHER_ID == -423) && (taum_MC_GD_MOTHER_ID == -{mother_ID})) || ((taum_TRUEID == 421) && (taum_MC_MOTHER_ID == -{mother_ID}) && (taup_TRUEID == -421) && (taup_MC_MOTHER_ID == -423) && (taup_MC_GD_MOTHER_ID == -{mother_ID}))"
    D0star_D0star = f"(abs(taup_TRUEID) == 421) && (abs(taum_TRUEID) == 421) && (abs(taup_MC_MOTHER_ID)==423) && (abs(taup_MC_GD_MOTHER_ID)=={mother_ID}) && (abs(taum_MC_MOTHER_ID)==423) && (abs(taum_MC_GD_MOTHER_ID)=={mother_ID})"

    # D+ D-
    Dp_Dm = f"(abs(taup_TRUEID) == 411) && (abs(taum_TRUEID) == 411) && (abs(taup_MC_MOTHER_ID)=={mother_ID}) && (abs(taum_MC_MOTHER_ID)=={mother_ID})"
    Dpstar_Dm = f"(((taup_TRUEID == 411) && (taum_TRUEID == -411) && (taup_MC_MOTHER_ID == 413) && (taup_MC_GD_MOTHER_ID == {mother_ID}) && (taum_MC_MOTHER_ID == {mother_ID})) || ((taup_TRUEID == 413) && (taum_TRUEID == -411) && (taup_MC_MOTHER_ID == {mother_ID}) && (taum_MC_MOTHER_ID == {mother_ID}))) || (((taum_TRUEID == 411) && (taup_TRUEID == -411) && (taum_MC_MOTHER_ID == 413) && (taum_MC_GD_MOTHER_ID == {mother_ID}) && (taup_MC_MOTHER_ID == {mother_ID})) || ((taum_TRUEID == 413) && (taup_TRUEID == -411) && (taum_MC_MOTHER_ID == {mother_ID}) && (taup_MC_MOTHER_ID == {mother_ID}))) || (((taup_TRUEID == -411) && (taum_TRUEID == 411) && (taup_MC_MOTHER_ID == -413) && (taup_MC_GD_MOTHER_ID == -{mother_ID}) && (taum_MC_MOTHER_ID == -{mother_ID})) || ((taup_TRUEID == -413) && (taum_TRUEID == 411) && (taup_MC_MOTHER_ID == -{mother_ID}) && (taum_MC_MOTHER_ID == -{mother_ID}))) || (((taum_TRUEID == -411) && (taup_TRUEID == 411) && (taum_MC_MOTHER_ID == -413) && (taum_MC_GD_MOTHER_ID == -{mother_ID}) && (taup_MC_MOTHER_ID == -{mother_ID})) || ((taum_TRUEID == -413) && (taup_TRUEID == 411) && (taum_MC_MOTHER_ID == -{mother_ID}) && (taup_MC_MOTHER_ID == -{mother_ID})))"
    Dp_Dmstar = f"(((taup_TRUEID == 411) && (taum_TRUEID == -411) && (taup_MC_MOTHER_ID == {mother_ID}) && (taum_MC_MOTHER_ID == -413) && (taum_MC_GD_MOTHER_ID == {mother_ID})) || ((taup_TRUEID == 411) && (taum_TRUEID == -413) && (taup_MC_MOTHER_ID == {mother_ID}) && (taum_MC_MOTHER_ID == {mother_ID}))) || (((taum_TRUEID == 411) && (taup_TRUEID == -411) && (taum_MC_MOTHER_ID == {mother_ID}) && (taup_MC_MOTHER_ID == -413) && (taup_MC_GD_MOTHER_ID == {mother_ID})) || ((taum_TRUEID == 411) && (taup_TRUEID == -413) && (taum_MC_MOTHER_ID == {mother_ID}) && (taup_MC_MOTHER_ID == {mother_ID}))) || (((taup_TRUEID == -411) && (taum_TRUEID == 411) && (taup_MC_MOTHER_ID == -{mother_ID}) && (taum_MC_MOTHER_ID == 413) && (taum_MC_GD_MOTHER_ID == -{mother_ID})) || ((taup_TRUEID == -411) && (taum_TRUEID == 413) && (taup_MC_MOTHER_ID == -{mother_ID}) && (taum_MC_MOTHER_ID == -{mother_ID}))) || (((taum_TRUEID == -411) && (taup_TRUEID == 411) && (taum_MC_MOTHER_ID == -{mother_ID}) && (taup_MC_MOTHER_ID == 413) && (taup_MC_GD_MOTHER_ID == -{mother_ID})) || ((taum_TRUEID == -411) && (taup_TRUEID == 413) && (taum_MC_MOTHER_ID == -{mother_ID}) && (taup_MC_MOTHER_ID == -{mother_ID})))"
    Dpstar_Dmstar = f"( (abs(taup_TRUEID) == 411) && (abs(taum_TRUEID) == 411) && (abs(taup_MC_MOTHER_ID)==413) && (abs(taup_MC_GD_MOTHER_ID)=={mother_ID}) && (abs(taum_MC_MOTHER_ID)==413) && (abs(taum_MC_GD_MOTHER_ID)=={mother_ID}) ) || ( (abs(taup_TRUEID) == 411) && (abs(taum_TRUEID) == 413) && (abs(taup_MC_MOTHER_ID)==413) && (abs(taup_MC_GD_MOTHER_ID)=={mother_ID}) && (abs(taum_MC_MOTHER_ID)=={mother_ID}) ) || ( (abs(taup_TRUEID) == 413) && (abs(taum_TRUEID) == 411) && (abs(taup_MC_MOTHER_ID)=={mother_ID}) && (abs(taum_MC_MOTHER_ID)==413) && (abs(taum_MC_GD_MOTHER_ID)=={mother_ID}) ) || ( (abs(taup_TRUEID) == 413) && (abs(taum_TRUEID) == 413) && (abs(taup_MC_MOTHER_ID)=={mother_ID}) && (abs(taum_MC_MOTHER_ID)=={mother_ID}) )"

    # D0 Ds
    D0_Ds = f"( (abs(taup_TRUEID) == 421) && (abs(taum_TRUEID) == 431) && (abs(taup_MC_MOTHER_ID)=={mother_ID}) && (abs(taum_MC_MOTHER_ID)=={mother_ID}) ) || ( (abs(taup_TRUEID) == 431) && (abs(taum_TRUEID) == 421) && (abs(taup_MC_MOTHER_ID)=={mother_ID}) && (abs(taum_MC_MOTHER_ID)=={mother_ID}) )"
    D0star_Ds = f"( (abs(taup_TRUEID) == 421) && (abs(taum_TRUEID) == 431) && (abs(taup_MC_MOTHER_ID)==423) && (abs(taup_MC_GD_MOTHER_ID)=={mother_ID}) && (abs(taum_MC_MOTHER_ID)=={mother_ID}) ) || ( (abs(taup_TRUEID) == 431) && (abs(taum_TRUEID) == 421) && (abs(taup_MC_MOTHER_ID)=={mother_ID}) && (abs(taum_MC_MOTHER_ID)==423) && (abs(taum_MC_GD_MOTHER_ID)=={mother_ID}) )"
    D0_Dsstar = f"( (abs(taup_TRUEID) == 421) && (abs(taum_TRUEID) == 431) && (abs(taup_MC_MOTHER_ID)=={mother_ID}) && (abs(taum_MC_MOTHER_ID)==433) && (abs(taum_MC_GD_MOTHER_ID)=={mother_ID}) ) || ( (abs(taup_TRUEID) == 431) && (abs(taum_TRUEID) == 421) && (abs(taup_MC_MOTHER_ID)==433) && (abs(taup_MC_GD_MOTHER_ID)=={mother_ID}) && (abs(taum_MC_MOTHER_ID)=={mother_ID}) )"
    D0star_Dsstar = f"( (abs(taup_TRUEID) == 421) && (abs(taum_TRUEID) == 431) && (abs(taup_MC_MOTHER_ID)==423) && (abs(taup_MC_GD_MOTHER_ID)=={mother_ID}) && (abs(taum_MC_MOTHER_ID)==433) && (abs(taum_MC_GD_MOTHER_ID)=={mother_ID}) ) || ( (abs(taup_TRUEID) == 431) && (abs(taum_TRUEID) == 421) && (abs(taup_MC_MOTHER_ID)==433) && (abs(taup_MC_GD_MOTHER_ID)=={mother_ID}) && (abs(taum_MC_MOTHER_ID)==423) && (abs(taum_MC_GD_MOTHER_ID)=={mother_ID}) )"

    # D+ D0
    Dp_D0 = f"( (abs(taup_TRUEID) == 411) && (abs(taum_TRUEID) == 421) && (abs(taup_MC_MOTHER_ID)=={mother_ID}) && (abs(taum_MC_MOTHER_ID)=={mother_ID}) ) || ( (abs(taup_TRUEID) == 421) && (abs(taum_TRUEID) == 411) && (abs(taup_MC_MOTHER_ID)=={mother_ID}) && (abs(taum_MC_MOTHER_ID)=={mother_ID}) )"
    Dpstar_D0 = f"( (abs(taup_TRUEID) == 411) && (abs(taum_TRUEID) == 421) && (abs(taup_MC_MOTHER_ID)==413) && (abs(taup_MC_GD_MOTHER_ID)=={mother_ID}) && (abs(taum_MC_MOTHER_ID)=={mother_ID}) ) || ( (abs(taup_TRUEID) == 421) && (abs(taum_TRUEID) == 411) && (abs(taup_MC_MOTHER_ID)=={mother_ID}) && (abs(taum_MC_MOTHER_ID)==413) && (abs(taum_MC_GD_MOTHER_ID)=={mother_ID}) ) || ( (abs(taup_TRUEID) == 413) && (abs(taum_TRUEID) == 421) && (abs(taup_MC_MOTHER_ID)=={mother_ID}) && (abs(taum_MC_MOTHER_ID)=={mother_ID}) ) || ( (abs(taup_TRUEID) == 421) && (abs(taum_TRUEID) == 413) && (abs(taup_MC_MOTHER_ID)=={mother_ID}) && (abs(taum_MC_MOTHER_ID)=={mother_ID}) )"
    Dp_D0star = f"( (abs(taup_TRUEID) == 411) && (abs(taum_TRUEID) == 421) && (abs(taup_MC_MOTHER_ID)=={mother_ID}) && (abs(taum_MC_MOTHER_ID)==423) && (abs(taum_MC_GD_MOTHER_ID)=={mother_ID}) ) || ( (abs(taup_TRUEID) == 421) && (abs(taum_TRUEID) == 411) && (abs(taup_MC_MOTHER_ID)==423) && (abs(taup_MC_GD_MOTHER_ID)=={mother_ID}) && (abs(taum_MC_MOTHER_ID)=={mother_ID}) )"
    Dpstar_D0star = f"( (abs(taup_TRUEID) == 411) && (abs(taum_TRUEID) == 421) && (abs(taup_MC_MOTHER_ID)==413) && (abs(taup_MC_GD_MOTHER_ID)=={mother_ID}) && (abs(taum_MC_MOTHER_ID)==423) && (abs(taum_MC_GD_MOTHER_ID)=={mother_ID}) ) || ( (abs(taup_TRUEID) == 421) && (abs(taum_TRUEID) == 411) && (abs(taup_MC_MOTHER_ID)==423) && (abs(taup_MC_GD_MOTHER_ID)=={mother_ID}) && (abs(taum_MC_MOTHER_ID)==413) && (abs(taum_MC_GD_MOTHER_ID)=={mother_ID}) ) || ( (abs(taup_TRUEID) == 413) && (abs(taum_TRUEID) == 421) && (abs(taup_MC_MOTHER_ID)=={mother_ID}) && (abs(taum_MC_MOTHER_ID)==423) && (abs(taum_MC_GD_MOTHER_ID)=={mother_ID}) ) || ( (abs(taup_TRUEID) == 421) && (abs(taum_TRUEID) == 413) && (abs(taup_MC_MOTHER_ID)==423) && (abs(taup_MC_GD_MOTHER_ID)=={mother_ID}) && (abs(taum_MC_MOTHER_ID)=={mother_ID}) )"

    # D+ Ds
    Dp_Ds = f"( (abs(taup_TRUEID) == 411) && (abs(taum_TRUEID) == 431) && (abs(taup_MC_MOTHER_ID)=={mother_ID}) && (abs(taum_MC_MOTHER_ID)=={mother_ID}) ) || ( (abs(taup_TRUEID) == 431) && (abs(taum_TRUEID) == 411) && (abs(taup_MC_MOTHER_ID)=={mother_ID}) && (abs(taum_MC_MOTHER_ID)=={mother_ID}) )"
    Dpstar_Ds = f"( (abs(taup_TRUEID) == 411) && (abs(taum_TRUEID) == 431) && (abs(taup_MC_MOTHER_ID)==413) && (abs(taup_MC_GD_MOTHER_ID)=={mother_ID}) && (abs(taum_MC_MOTHER_ID)=={mother_ID}) ) || ( (abs(taup_TRUEID) == 431) && (abs(taum_TRUEID) == 411) && (abs(taup_MC_MOTHER_ID)=={mother_ID}) && (abs(taum_MC_MOTHER_ID)==413) &&  (abs(taum_MC_GD_MOTHER_ID)=={mother_ID}) ) || ( (abs(taup_TRUEID) == 413) && (abs(taum_TRUEID) == 431) && (abs(taup_MC_MOTHER_ID)=={mother_ID}) && (abs(taum_MC_MOTHER_ID)=={mother_ID}) ) || ( (abs(taup_TRUEID) == 431) && (abs(taum_TRUEID) == 413) && (abs(taup_MC_MOTHER_ID)=={mother_ID}) && (abs(taum_MC_MOTHER_ID)=={mother_ID}) )"
    Dp_Dsstar = f"( (abs(taup_TRUEID) == 411) && (abs(taum_TRUEID) == 431) && (abs(taup_MC_MOTHER_ID)=={mother_ID}) && (abs(taum_MC_MOTHER_ID)==433) && (abs(taum_MC_GD_MOTHER_ID)=={mother_ID}) ) || ( (abs(taup_TRUEID) == 431) && (abs(taum_TRUEID) == 411) && (abs(taup_MC_MOTHER_ID)==433) && (abs(taup_MC_GD_MOTHER_ID)=={mother_ID}) && (abs(taum_MC_MOTHER_ID)=={mother_ID}) )"
    Dpstar_Dsstar = f"( (abs(taup_TRUEID) == 411) && (abs(taum_TRUEID) == 431) && (abs(taup_MC_MOTHER_ID)==413) && (abs(taup_MC_GD_MOTHER_ID)=={mother_ID}) && (abs(taum_MC_MOTHER_ID)==433) && (abs(taum_MC_GD_MOTHER_ID)=={mother_ID}) ) || ( (abs(taup_TRUEID) == 431) && (abs(taum_TRUEID) == 411) && (abs(taup_MC_MOTHER_ID)==433) && (abs(taup_MC_GD_MOTHER_ID)=={mother_ID}) && (abs(taum_MC_MOTHER_ID)==413) && (abs(taum_MC_GD_MOTHER_ID)=={mother_ID}) ) || ( (abs(taup_TRUEID) == 413) && (abs(taum_TRUEID) == 431) && (abs(taup_MC_MOTHER_ID)=={mother_ID}) && (abs(taum_MC_MOTHER_ID)==433) && (abs(taum_MC_GD_MOTHER_ID)=={mother_ID}) ) || ( (abs(taup_TRUEID) == 431) && (abs(taum_TRUEID) == 413) && (abs(taup_MC_MOTHER_ID)==433) && (abs(taup_MC_GD_MOTHER_ID)=={mother_ID}) &&  (abs(taum_MC_MOTHER_ID)=={mother_ID}) )"

    if(name == "Ds_Ds"):
        return Ds_Ds
    elif(name == "Dsstar_Ds"):
        return Dsstar_Ds
    elif(name == "Ds_Dsstar"):
        return Ds_Dsstar
    elif(name == "Dsstar_Dsstar"):
        return Dsstar_Dsstar
    elif(name == "D0_D0"):
        return D0_D0
    elif(name == "D0star_D0"):
        return D0star_D0
    elif(name == "D0_D0star"):
        return D0_D0star
    elif(name == "D0star_D0star"):
        return D0star_D0star
    elif(name == "Dp_Dm"): 
        return Dp_Dm
    elif(name == "Dpstar_Dm"):
        return Dpstar_Dm
    elif(name == "Dp_Dmstar"):
        return Dp_Dmstar
    elif(name == "Dp_Dmstar"):
        return Dp_Dmstar
    elif(name == "Dpstar_Dmstar"):
        return Dpstar_Dmstar
    elif(name == "D0_Ds"): 
        return D0_Ds
    elif(name == "D0star_Ds"):
        return D0star_Ds
    elif(name == "D0_Dsstar"):
        return D0_Dsstar
    elif(name == "D0star_Dsstar"):
        return D0star_Dsstar
    elif(name == "Dp_D0"):
        return Dp_D0
    elif(name == "Dpstar_D0"):
        return Dpstar_D0
    elif(name == "Dp_D0star"):
        return Dp_D0star
    elif(name == "Dpstar_D0star"):
        return Dpstar_D0star
    elif(name == "Dp_Ds"): 
        return Dp_Ds
    elif(name == "Dpstar_Ds"):
        return Dpstar_Ds
    elif(name == "Dp_Dsstar"):
        return Dp_Dsstar
    elif(name == "Dpstar_Dsstar"):
        return Dpstar_Dsstar


def truthMatch_cocktailMC(name):
    final_state = "(abs(taup_pi1_TRUEID) == 211) && (abs(taup_pi2_TRUEID) == 211) && (abs(taup_pi3_TRUEID) == 211) && (abs(taum_pi1_TRUEID) == 211) && (abs(taum_pi2_TRUEID) == 211) && (abs(taum_pi3_TRUEID) == 211) && (abs(Kp_TRUEID ) == 321)"

    # BuDDKp
    # B+ -> D0D0K+
    BuD0D0Kp = "(abs(Bp_TRUEID) == 521)"+" && "+final_state+" && "+DD_cut("D0_D0",521) 
    BuD0starD0Kp = "(abs(Bp_TRUEID) == 521)"+" && "+final_state+" && "+DD_cut("D0star_D0",521)
    BuD0D0starKp = "(abs(Bp_TRUEID) == 521)"+" && "+final_state+" && "+DD_cut("D0_D0star",521)
    BuD0starD0starKp = "(abs(Bp_TRUEID) == 521)"+" && "+final_state+" && "+DD_cut("D0star_D0star",521)

    theMap = {"BuD0D0Kp": BuD0D0Kp, "BuD0starD0Kp": BuD0starD0Kp, "BuD0D0starKp": BuD0D0starKp, "BuD0starD0starKp": BuD0starD0starKp}

    # B+ -> D+D-K+
    BuDpDmKp = "(abs(Bp_TRUEID) == 521)"+" && "+final_state+" && "+DD_cut("Dp_Dm",521)
    BuDpstarDmKp = "(abs(Bp_TRUEID) == 521)"+" && "+final_state+" && "+DD_cut("Dpstar_Dm",521)
    BuDpDmstarKp = "(abs(Bp_TRUEID) == 521)"+" && "+final_state+" && "+DD_cut("Dp_Dmstar",521)
    BuDpstarDmstarKp = "(abs(Bp_TRUEID) == 521)"+" && "+final_state+" && "+DD_cut("Dpstar_Dmstar",521)
    theMap["BuDpDmKp"] = BuDpDmKp
    theMap["BuDpstarDmKp"] = BuDpstarDmKp
    theMap["BuDpDmstarKp"] = BuDpDmstarKp
    theMap["BuDpstarDmstarKp"] = BuDpstarDmstarKp

    # B+ -> Ds+Ds-K+
    BuDsDsKp = "(abs(Bp_TRUEID) == 521)"+" && "+final_state+" && "+DD_cut("Ds_Ds",521)
    BuDsstarDsKp = "(abs(Bp_TRUEID) == 521)"+" && "+final_state+" && "+DD_cut("Dsstar_Ds",521)
    BuDsDsstarKp = "(abs(Bp_TRUEID) == 521)"+" && "+final_state+" && "+DD_cut("Ds_Dsstar",521)
    BuDsstarDsstarKp = "(abs(Bp_TRUEID) == 521)"+" && "+final_state+" && "+DD_cut("Dsstar_Dsstar",521)
    theMap["BuDsDsKp"] = BuDsDsKp
    theMap["BuDsstarDsKp"] = BuDsstarDsKp
    theMap["BuDsDsstarKp"] = BuDsDsstarKp
    theMap["BuDsstarDsstarKp"] = BuDsstarDsstarKp

    # BdDDKp
    # B0 -> D-D0K+
    BdDmD0Kp = "(abs(Bp_TRUEID) == 511)"+" && "+final_state+" && "+DD_cut("Dp_D0",511)
    BdDmstarD0Kp = "(abs(Bp_TRUEID) == 511)"+" && "+final_state+" && "+DD_cut("Dpstar_D0",511)
    BdDmD0starKp = "(abs(Bp_TRUEID) == 511)"+" && "+final_state+" && "+DD_cut("Dp_D0star",511)
    BdDmstarD0starKp = "(abs(Bp_TRUEID) == 511)"+" && "+final_state+" && "+DD_cut("Dpstar_D0star",511)
    theMap["BdDmD0Kp"] = BdDmD0Kp
    theMap["BdDmstarD0Kp"] = BdDmstarD0Kp
    theMap["BdDmD0starKp"] = BdDmD0starKp
    theMap["BdDmstarD0starKp"] = BdDmstarD0starKp

    # BsDDKp
    # Bs -> Ds-D0K+
    BsDsD0Kp = "(abs(Bp_TRUEID) == 531)"+" && "+final_state+" && "+DD_cut("D0_Ds",531)
    BsDsstarD0Kp = "(abs(Bp_TRUEID) == 531)"+" && "+final_state+" && "+DD_cut("D0_Dsstar",531)
    BsDsD0starKp = "(abs(Bp_TRUEID) == 531)"+" && "+final_state+" && "+DD_cut("D0star_Ds",531)
    BsDsstarD0starKp = "(abs(Bp_TRUEID) == 531)"+" && "+final_state+" && "+DD_cut("D0star_Dsstar",531)
    theMap["BsDsD0Kp"] = BsDsD0Kp
    theMap["BsDsstarD0Kp"] = BsDsstarD0Kp
    theMap["BsDsD0starKp"] = BsDsD0starKp
    theMap["BsDsstarD0starKp"] = BsDsstarD0starKp

    # BuDDK0
    # B+ -> D0D+K0
    BuD0DpK0 = final_state+" && "+DD_cut("Dp_D0",521)
    BuD0starDpK0 = final_state+" && "+DD_cut("Dp_D0star",521)
    BuD0DpstarK0 = final_state+" && "+DD_cut("Dpstar_D0",521)
    BuD0starDpstarK0 = final_state+" && "+DD_cut("Dpstar_D0star",521)
    theMap["BuD0DpK0"] = BuD0DpK0
    theMap["BuD0starDpK0"] = BuD0starDpK0
    theMap["BuD0DpstarK0"] = BuD0DpstarK0
    theMap["BuD0starDpstarK0"] = BuD0starDpstarK0

    # BdDDK0
    # B0 -> D+D-K0
    BdDpDmK0 = final_state+" && "+DD_cut("Dp_Dm",511)
    BdDpstarDmK0 = final_state+" && "+DD_cut("Dpstar_Dm",511)
    BdDpDmstarK0 = final_state+" && "+DD_cut("Dp_Dmstar",511)
    BdDpstarDmstarK0 = final_state+" && "+DD_cut("Dpstar_Dmstar",511)
    theMap["BdDpDmK0"] = BdDpDmK0
    theMap["BdDpstarDmK0"] = BdDpstarDmK0
    theMap["BdDpDmstarK0"] = BdDpDmstarK0
    theMap["BdDpstarDmstarK0"] = BdDpstarDmstarK0

    # B0 -> D0D0K0
    BdD0D0K0 = final_state+" && "+DD_cut("D0_D0",511)
    BdD0starD0K0 = final_state+" && "+DD_cut("D0star_D0",511)
    BdD0D0starK0 = final_state+" && "+DD_cut("D0_D0star",511)
    BdD0starD0starK0 = final_state+" && "+DD_cut("D0star_D0star",511)
    theMap["BdD0D0K0"] = BdD0D0K0
    theMap["BdD0starD0K0"] = BdD0starD0K0
    theMap["BdD0D0starK0"] = BdD0D0starK0
    theMap["BdD0starD0starK0"] = BdD0starD0starK0

    # BuDD
    # B+ -> D0Ds+
    BuD0Ds = final_state+" && "+DD_cut("D0_Ds",521)
    BuD0starDs = final_state+" && "+DD_cut("D0star_Ds",521)
    BuD0Dsstar = final_state+" && "+DD_cut("D0_Dsstar",521)
    BuD0starDsstar = final_state+" && "+DD_cut("D0star_Dsstar",521)
    theMap["BuD0Ds"] = BuD0Ds
    theMap["BuD0starDs"] = BuD0starDs
    theMap["BuD0Dsstar"] = BuD0Dsstar
    theMap["BuD0starDsstar"] = BuD0starDsstar

    # B+ -> D0D+
    BuD0Dp = final_state+" && "+DD_cut("Dp_D0",521)
    BuD0starDp = final_state+" && "+DD_cut("Dp_D0star",521)
    BuD0Dpstar = final_state+" && "+DD_cut("Dpstar_D0",521)
    BuD0starDpstar = final_state+" && "+DD_cut("Dpstar_D0star",521)
    theMap["BuD0Dp"] = BuD0Dp
    theMap["BuD0starDp"] = BuD0starDp
    theMap["BuD0Dpstar"] = BuD0Dpstar
    theMap["BuD0starDpstar"] = BuD0starDpstar


    # BdDD
    # B0 -> D0D0
    BdD0D0 = final_state+" && "+DD_cut("D0_D0",511)
    BdD0starD0 = final_state+" && "+DD_cut("D0star_D0",511)
    BdD0D0star = final_state+" && "+DD_cut("D0_D0star",511)
    BdD0starD0star = final_state+" && "+DD_cut("D0star_D0star",511)
    theMap["BdD0D0"] = BdD0D0
    theMap["BdD0starD0"] = BdD0starD0
    theMap["BdD0D0star"] = BdD0D0star
    theMap["BdD0starD0star"] = BdD0starD0star

    # B0 -> D+D-
    BdDpDm = final_state+" && "+DD_cut("Dp_Dm",511)
    BdDpstarDm = final_state+" && "+DD_cut("Dpstar_Dm",511)
    BdDpDmstar = final_state+" && "+DD_cut("Dp_Dmstar",511)
    BdDpstarDmstar = final_state+" && "+DD_cut("Dpstar_Dmstar",511)
    theMap["BdDpDm"] = BdDpDm
    theMap["BdDpstarDm"] = BdDpstarDm
    theMap["BdDpDmstar"] = BdDpDmstar
    theMap["BdDpstarDmstar"] = BdDpstarDmstar


    # B0 -> D-Ds+
    BdDpDs = final_state+" && "+DD_cut("Dp_Ds",511)
    BdDpstarDs = final_state+" && "+DD_cut("Dpstar_Ds",511)
    BdDpDsstar = final_state+" && "+DD_cut("Dp_Dsstar",511)
    BdDpstarDsstar = final_state+" && "+DD_cut("Dpstar_Dsstar",511)
    theMap["BdDpDs"] = BdDpDs
    theMap["BdDpstarDs"] = BdDpstarDs
    theMap["BdDpDsstar"] = BdDpDsstar
    theMap["BdDpstarDsstar"] = BdDpstarDsstar


    # B0 -> Ds+Ds-
    BdDsDs = final_state+" && "+DD_cut("Ds_Ds",511)
    BdDsstarDs = final_state+" && "+DD_cut("Dsstar_Ds",511)
    BdDsDsstar = final_state+" && "+DD_cut("Ds_Dsstar",511)
    BdDsstarDsstar = final_state+" && "+DD_cut("Dsstar_Dsstar",511)
    theMap["BdDsDs"] = BdDsDs
    theMap["BdDsstarDs"] = BdDsstarDs
    theMap["BdDsDsstar"] = BdDsDsstar
    theMap["BdDsstarDsstar"] = BdDsstarDsstar


    # Bs -> DD
    # Bs -> Ds+Ds-
    BsDsDs = final_state+" && "+DD_cut("Ds_Ds",531)
    BsDsstarDs = final_state+" && "+DD_cut("Dsstar_Ds",531)
    BsDsDsstar = final_state+" && "+DD_cut("Ds_Dsstar",531)
    BsDsstarDsstar = final_state+" && "+DD_cut("Dsstar_Dsstar",531)
    theMap["BsDsDs"] = BsDsDs
    theMap["BsDsstarDs"] = BsDsstarDs
    theMap["BsDsDsstar"] = BsDsDsstar
    theMap["BsDsstarDsstar"] = BsDsstarDsstar

    # Bs -> D-Ds+
    BsDpDs = final_state+" && "+DD_cut("Dp_Ds",531)
    BsDpstarDs = final_state+" && "+DD_cut("Dpstar_Ds",531)
    BsDpDsstar = final_state+" && "+DD_cut("Dp_Dsstar",531)
    BsDpstarDsstar = final_state+" && "+DD_cut("Dsstar_Dsstar",531)
    theMap["BsDpDs"] = BsDpDs
    theMap["BsDpstarDs"] = BsDpstarDs
    theMap["BsDpDsstar"] = BsDpDsstar
    theMap["BsDpstarDsstar"] = BsDpstarDsstar

    # Bs -> D+D-
    BsDpDm = final_state+" && "+DD_cut("Dp_Dm",531)
    BsDpstarDm = final_state+" && "+DD_cut("Dpstar_Dm",531)
    BsDpDmstar = final_state+" && "+DD_cut("Dp_Dmstar",531)
    BsDpstarDmstar = final_state+" && "+DD_cut("Dpstar_Dmstar",531)
    theMap["BsDpDm"] = BsDpDm
    theMap["BsDpstarDm"] = BsDpstarDm
    theMap["BsDpDmstar"] = BsDpDmstar
    theMap["BsDpstarDmstar"] = BsDpstarDmstar

    # Bs -> D0D0
    BsD0D0 = final_state+" && "+DD_cut("D0_D0",531)
    BsD0starD0 = final_state+" && "+DD_cut("D0star_D0",531)
    BsD0D0star = final_state+" && "+DD_cut("D0_D0star",531)
    BsD0starD0star = final_state+" && "+DD_cut("D0star_D0star",531)
    theMap["BsD0D0"] = BsD0D0
    theMap["BsD0starD0"] = BsD0starD0
    theMap["BsD0D0star"] = BsD0D0star
    theMap["BsD0starD0star"] = BsD0starD0star

    return theMap[name]


def create_tables(species, truthMatch, L0_trigger, HLT1_trigger, HLT2_trigger, trigger, rectangular_cuts_names, rectangular_cuts, pass_mass_fit, fit_range, mass_vetoes, mass_vetoes_cuts, best_cand):
    # Acceptance efficiency
    print("Acceptance efficiency")
    if(species == 100):
        acc_columns = ["All years (average MagUp MagDown) (\\%)"]
        acc_rows = ["$B^+ \\to D D K^+$", "$B^0 \\to D D K^+$", "$B^0_s \\to D D K^+$", "$B^+ \\to D D K^0$", "$B^+ \\to D D$"]
    else:
        acc_columns = ["Average of MagUp and MagDown (\\%)"]
        acc_rows = ["2016", "2017", "2018", "All"]

    if(species == 1):
        eps_acc_2016 = ufloat(5.318/100, 0.010/100)
        eps_acc_2017 = ufloat(5.321/100, 0.011/100)
        eps_acc_2018 = ufloat(5.330/100, 0.010/100)
        eps_acc = (28.5/165)*eps_acc_2016 + (60.0/165)*eps_acc_2017 + (76.5/165)*eps_acc_2018
    elif(species == 100):
        eps_acc_2016 = [ufloat(1.056/100, 0.002/100), ufloat(1.364/100, 0.002/100), ufloat(1.522/100, 0.003/100), ufloat(2.293/100, 0.004/100), ufloat(2.922/100, 0.005/100)]
        eps_acc_2017 = [ufloat(1.055/100, 0.002/100), ufloat(1.371/100, 0.003/100), ufloat(1.517/100, 0.003/100), ufloat(2.291/100, 0.004/100), ufloat(2.926/100, 0.005/100)]
        eps_acc_2018 = [ufloat(1.051/100, 0.002/100), ufloat(1.369/100, 0.003/100), ufloat(1.524/100, 0.003/100), ufloat(2.289/100, 0.004/100), ufloat(2.919/100, 0.005/100)]

        eps_acc = [(7.2/32.2)*eps_acc_2016[0] + (10.0/32.2)*eps_acc_2017[0] + (15.0/32.2)*eps_acc_2018[0], (7.2/32.2)*eps_acc_2016[1] + (10.0/32.2)*eps_acc_2017[1] + (15.0/32.2)*eps_acc_2018[1], (2.4/10.8)*eps_acc_2016[2] + (3.4/10.8)*eps_acc_2017[2] + (5.0/10.8)*eps_acc_2018[2], (14.3/64.3)*eps_acc_2016[3] + (20.0/64.3)*eps_acc_2017[3] + (30.0/64.3)*eps_acc_2018[3], (14.3/64.3)*eps_acc_2016[4] + (20.0/64.3)*eps_acc_2017[4] + (30.0/64.3)*eps_acc_2018[4]]
    else:
        eps_acc_2016 = ufloat(14.50/100, 0.10/100)
        eps_acc_2017 = ufloat(14.48/100, 0.10/100)
        eps_acc_2018 = ufloat(14.58/100, 0.10/100)
        eps_acc = (10.8/34.2)*eps_acc_2016 + (12./34.2)*eps_acc_2017 + (11.4/34.2)*eps_acc_2018

    acc_table = pd.DataFrame(index=acc_rows, columns=acc_columns)
    if(species == 100):
        acc_table.loc["$B^+ \\to D D K^+$", "All years (average MagUp MagDown) (\\%)"] = f"${eps_acc[0].nominal_value*100:.3f} \\pm {eps_acc[0].std_dev*100:.3f}$"
        acc_table.loc["$B^0 \\to D D K^+$", "All years (average MagUp MagDown) (\\%)"] = f"${eps_acc[1].nominal_value*100:.3f} \\pm {eps_acc[1].std_dev*100:.3f}$"
        acc_table.loc["$B^0_s \\to D D K^+$", "All years (average MagUp MagDown) (\\%)"] = f"${eps_acc[2].nominal_value*100:.3f} \\pm {eps_acc[2].std_dev*100:.3f}$"
        acc_table.loc["$B^+ \\to D D K^0$", "All years (average MagUp MagDown) (\\%)"] = f"${eps_acc[3].nominal_value*100:.3f} \\pm {eps_acc[3].std_dev*100:.3f}$"
        acc_table.loc["$B^+ \\to D D$", "All years (average MagUp MagDown) (\\%)"] = f"${eps_acc[4].nominal_value*100:.3f} \\pm {eps_acc[4].std_dev*100:.3f}$"
    else:
        acc_table.loc["2016", "Average of MagUp and MagDown (\\%)"] = f"${eps_acc_2016.nominal_value*100} \\pm {eps_acc_2016.std_dev*100}$"
        acc_table.loc["2017", "Average of MagUp and MagDown (\\%)"] = f"${eps_acc_2017.nominal_value*100} \\pm {eps_acc_2017.std_dev*100}$"
        acc_table.loc["2018", "Average of MagUp and MagDown (\\%)"] = f"${eps_acc_2018.nominal_value*100} \\pm {eps_acc_2018.std_dev*100}$"
        acc_table.loc["All", "Average of MagUp and MagDown (\\%)"] = f"${eps_acc.nominal_value*100:.3f} \\pm {eps_acc.std_dev*100:.3f}$"

    with open(f'/panfs/felician/B2Ktautau/workflow/selections_efficiency_tables/Species_{species}/acceptance_table.tex', 'w') as fout_acc:
        fout_acc.write(acc_table.to_latex())

    # Stripping filter efficiencies (for Ktautau only)
    print("Stripping efficiency")
    if(species == 100):
        strip_columns = ["All years (average MagUp MagDown) (\\%)"]
        strip_rows = ["$B^+ \\to D D K^+$", "$B^0 \\to D D K^+$", "$B^0_s \\to D D K^+$", "$B^+ \\to D D K^0$", "$B^+ \\to D D$"]

        eps_strip_2016 = [ufloat(1.719/100, 0.006/100), ufloat(1.591/100, 0.008/100), ufloat(2.75/100, 0.01/100), ufloat(0.785/100, 0.004/100), ufloat(1.102/100, 0.005/100)]
        eps_strip_2017 = [ufloat(1.760/100, 0.005/100), ufloat(1.643/100, 0.007/100), ufloat(2.81/100, 0.01/100), ufloat(0.806/100, 0.004/100), ufloat(1.140/100, 0.005/100)]
        eps_strip_2018 = [ufloat(1.763/100, 0.004/100), ufloat(1.625/100, 0.006/100), ufloat(2.799/100, 0.009/100), ufloat(0.805/100, 0.004/100), ufloat(1.127/100, 0.005/100)]
        eps_strip = [(7.2/32.2)*eps_strip_2016[0] + (10.0/32.2)*eps_strip_2017[0] + (15.0/32.2)*eps_strip_2018[0], (7.2/32.2)*eps_strip_2016[1] + (10.0/32.2)*eps_strip_2017[1] + (15.0/32.2)*eps_strip_2018[1], (2.4/10.8)*eps_strip_2016[2] + (3.4/10.8)*eps_strip_2017[2] + (5.0/10.8)*eps_strip_2018[2], (14.3/64.3)*eps_strip_2016[3] + (20.0/64.3)*eps_strip_2017[3] + (30.0/64.3)*eps_strip_2018[3], (14.3/64.3)*eps_strip_2016[4] + (20.0/64.3)*eps_strip_2017[4] + (30.0/64.3)*eps_strip_2018[4]]

        strip_table = pd.DataFrame(index=strip_rows, columns=strip_columns)
        strip_table.loc["$B^+ \\to D D K^+$", "All years (average MagUp MagDown) (\\%)"] = f"${eps_strip[0].nominal_value*100:.3f} \\pm {eps_strip[0].std_dev*100:.3f}$"
        strip_table.loc["$B^0 \\to D D K^+$", "All years (average MagUp MagDown) (\\%)"] = f"${eps_strip[1].nominal_value*100:.3f} \\pm {eps_strip[1].std_dev*100:.3f}$"
        strip_table.loc["$B^0_s \\to D D K^+$", "All years (average MagUp MagDown) (\\%)"] = f"${eps_strip[2].nominal_value*100:.3f} \\pm {eps_strip[2].std_dev*100:.3f}$"
        strip_table.loc["$B^+ \\to D D K^0$", "All years (average MagUp MagDown) (\\%)"] = f"${eps_strip[3].nominal_value*100:.3f} \\pm {eps_strip[3].std_dev*100:.3f}$"
        strip_table.loc["$B^+ \\to D D$", "All years (average MagUp MagDown) (\\%)"] = f"${eps_strip[4].nominal_value*100:.3f} \\pm {eps_strip[4].std_dev*100:.3f}$"

        with open(f'/panfs/felician/B2Ktautau/workflow/selections_efficiency_tables/Species_{species}/stripping_filter_table.tex', 'w') as fout_strip:
            fout_strip.write(strip_table.to_latex())

    elif(species == 1):
        strip_columns = ["Average of MagUp and MagDown (\\%)"]
        strip_rows = ["2016", "2017", "2018", "All"]

        eps_strip_2016 = ufloat(1.028/100, 0.003/100)
        eps_strip_2017 = ufloat(1.056/100, 0.003/100)
        eps_strip_2018 = ufloat(1.050/100, 0.002/100)
        eps_strip = (28.5/165)*eps_strip_2016 + (60.0/165)*eps_strip_2017 + (76.5/165)*eps_strip_2018

        strip_table = pd.DataFrame(index=strip_rows, columns=strip_columns)
        strip_table.loc["2016", "Average of MagUp and MagDown (\\%)"] = f"${eps_strip_2016.nominal_value*100} \\pm {eps_strip_2016.std_dev*100}$"
        strip_table.loc["2017", "Average of MagUp and MagDown (\\%)"] = f"${eps_strip_2017.nominal_value*100} \\pm {eps_strip_2017.std_dev*100}$"
        strip_table.loc["2018", "Average of MagUp and MagDown (\\%)"] = f"${eps_strip_2018.nominal_value*100} \\pm {eps_strip_2018.std_dev*100}$"
        strip_table.loc["All", "Average of MagUp and MagDown (\\%)"] = f"${eps_strip.nominal_value*100:.3f} \\pm {eps_strip.std_dev*100:.3f}$"

        with open(f'/panfs/felician/B2Ktautau/workflow/selections_efficiency_tables/Species_{species}/stripping_filter_table.tex', 'w') as fout_strip:
            fout_strip.write(strip_table.to_latex())

    np.save(f"/panfs/felician/B2Ktautau/workflow/selections_efficiency_tables/Species_{species}/acceptance_efficiency.npy", eps_acc)
    if(species == 1):
        np.save(f"/panfs/felician/B2Ktautau/workflow/selections_efficiency_tables/Species_{species}/stripping_efficiency.npy", eps_strip)


    if(species == 100):
        reco1_columns = ["Reconstruction efficiency I (\\%)"]
        reco1_rows = ["$B^+ \\to D D K^+$", "$B^0 \\to D D K^+$", "$B^0_s \\to D D K^+$", "$B^+ \\to D D K^0$", "$B^+ \\to D D$"]

        eps_reco1_2016 = [ufloat(99.8/100, 0.2/100), ufloat(99.8/100, 0.2/100), ufloat(99.7/100, 0.2/100), ufloat(99.8/100, 0.2/100), ufloat(99.9/100, 0.2/100)]
        eps_reco1_2017 = [ufloat(99.7/100, 0.1/100), ufloat(99.7/100, 0.2/100), ufloat(99.7/100, 0.1/100), ufloat(99.7/100, 0.2/100), ufloat(99.8/100, 0.2/100)]
        eps_reco1_2018 = [ufloat(99.7/100, 0.1/100), ufloat(99.7/100, 0.2/100), ufloat(99.7/100, 0.1/100), ufloat(99.8/100, 0.2/100), ufloat(99.8/100, 0.2/100)]
        eps_reco1 = [(7.2/32.2)*eps_reco1_2016[0] + (10.0/32.2)*eps_reco1_2017[0] + (15.0/32.2)*eps_reco1_2018[0], (7.2/32.2)*eps_reco1_2016[1] + (10.0/32.2)*eps_reco1_2017[1] + (15.0/32.2)*eps_reco1_2018[1], (2.4/10.8)*eps_reco1_2016[2] + (3.4/10.8)*eps_reco1_2017[2] + (5.0/10.8)*eps_reco1_2018[2], (14.3/64.3)*eps_reco1_2016[3] + (20.0/64.3)*eps_reco1_2017[3] + (30.0/64.3)*eps_reco1_2018[3], (14.3/64.3)*eps_reco1_2016[4] + (20.0/64.3)*eps_reco1_2017[4] + (30.0/64.3)*eps_reco1_2018[4]]

        reco1_table = pd.DataFrame(index=reco1_rows, columns=reco1_columns)
        reco1_table.loc["$B^+ \\to D D K^+$", "All years (average MagUp MagDown) (\\%)"] = f"${eps_reco1[0].nominal_value*100:.3f} \\pm {eps_reco1[0].std_dev*100:.3f}$"
        reco1_table.loc["$B^0 \\to D D K^+$", "All years (average MagUp MagDown) (\\%)"] = f"${eps_reco1[1].nominal_value*100:.3f} \\pm {eps_reco1[1].std_dev*100:.3f}$"
        reco1_table.loc["$B^0_s \\to D D K^+$", "All years (average MagUp MagDown) (\\%)"] = f"${eps_reco1[2].nominal_value*100:.3f} \\pm {eps_reco1[2].std_dev*100:.3f}$"
        reco1_table.loc["$B^+ \\to D D K^0$", "All years (average MagUp MagDown) (\\%)"] = f"${eps_reco1[3].nominal_value*100:.3f} \\pm {eps_reco1[3].std_dev*100:.3f}$"
        reco1_table.loc["$B^+ \\to D D$", "All years (average MagUp MagDown) (\\%)"] = f"${eps_reco1[4].nominal_value*100:.3f} \\pm {eps_reco1[4].std_dev*100:.3f}$"

        with open(f'/panfs/felician/B2Ktautau/workflow/selections_efficiency_tables/Species_{species}/reco_matching_eff_table.tex', 'w') as fout_reco1:
            fout_reco1.write(reco1_table.to_latex())

    # Reconstruction efficiency (truth-match reconstructed / accepted)
    print("Reconstruction efficiency")
    if(species == 1):
        fc_2016 = ROOT.TFileCollection("fc_2016", "fc_2016", "Files_on_grid/MC_2016.txt")
        fc_2017 = ROOT.TFileCollection("fc_2017", "fc_2017", "Files_on_grid/MC_2017.txt")
        fc_2018 = ROOT.TFileCollection("fc_2018", "fc_2018", "Files_on_grid/MC_2018.txt")

        # 2016
        t_gen_3pi3pi_2016 = ROOT.TChain("mc_ntuple_3pi_3pi/MCDecayTree")
        t_gen_3pi3pi_pi0_2016 = ROOT.TChain("mc_ntuple_3pi_3pipi0/MCDecayTree")
        t_gen_3pipi0_3pi_2016 = ROOT.TChain("mc_ntuple_3pipi0_3pi/MCDecayTree")
        t_gen_3pi3pi2pi0_2016 = ROOT.TChain("mc_ntuple_3pipi0_3pipi0/MCDecayTree")

        t_gen_3pi3pi_2016.AddFileInfoList(fc_2016.GetList())
        t_gen_3pi3pi_pi0_2016.AddFileInfoList(fc_2016.GetList())
        t_gen_3pipi0_3pi_2016.AddFileInfoList(fc_2016.GetList())
        t_gen_3pi3pi2pi0_2016.AddFileInfoList(fc_2016.GetList())

        N_gen_3pi3pi_2016 = t_gen_3pi3pi_2016.GetEntries()
        N_gen_3pi3pipi0_2016 = t_gen_3pi3pi_pi0_2016.GetEntries() + t_gen_3pipi0_3pi_2016.GetEntries()
        N_gen_3pi3pi2pi0_2016 = t_gen_3pi3pi2pi0_2016.GetEntries()
        N_gen_2016 = N_gen_3pi3pi_2016+N_gen_3pi3pipi0_2016+N_gen_3pi3pi2pi0_2016

        # 2017
        t_gen_3pi3pi_2017 = ROOT.TChain("mc_ntuple_3pi_3pi/MCDecayTree")
        t_gen_3pi3pi_pi0_2017 = ROOT.TChain("mc_ntuple_3pi_3pipi0/MCDecayTree")
        t_gen_3pipi0_3pi_2017 = ROOT.TChain("mc_ntuple_3pipi0_3pi/MCDecayTree")
        t_gen_3pi3pi2pi0_2017 = ROOT.TChain("mc_ntuple_3pipi0_3pipi0/MCDecayTree")

        t_gen_3pi3pi_2017.AddFileInfoList(fc_2017.GetList())
        t_gen_3pi3pi_pi0_2017.AddFileInfoList(fc_2017.GetList())
        t_gen_3pipi0_3pi_2017.AddFileInfoList(fc_2017.GetList())
        t_gen_3pi3pi2pi0_2017.AddFileInfoList(fc_2017.GetList())

        N_gen_3pi3pi_2017 = t_gen_3pi3pi_2017.GetEntries()
        N_gen_3pi3pipi0_2017 = t_gen_3pi3pi_pi0_2017.GetEntries() + t_gen_3pipi0_3pi_2017.GetEntries()
        N_gen_3pi3pi2pi0_2017 = t_gen_3pi3pi2pi0_2017.GetEntries()
        N_gen_2017 = N_gen_3pi3pi_2017+N_gen_3pi3pipi0_2017+N_gen_3pi3pi2pi0_2017

        # 2018
        t_gen_3pi3pi_2018 = ROOT.TChain("mc_ntuple_3pi_3pi/MCDecayTree")
        t_gen_3pi3pi_pi0_2018 = ROOT.TChain("mc_ntuple_3pi_3pipi0/MCDecayTree")
        t_gen_3pipi0_3pi_2018 = ROOT.TChain("mc_ntuple_3pipi0_3pi/MCDecayTree")
        t_gen_3pi3pi2pi0_2018 = ROOT.TChain("mc_ntuple_3pipi0_3pipi0/MCDecayTree")

        t_gen_3pi3pi_2018.AddFileInfoList(fc_2018.GetList())
        t_gen_3pi3pi_pi0_2018.AddFileInfoList(fc_2018.GetList())
        t_gen_3pipi0_3pi_2018.AddFileInfoList(fc_2018.GetList())
        t_gen_3pi3pi2pi0_2018.AddFileInfoList(fc_2018.GetList())

        N_gen_3pi3pi_2018 = t_gen_3pi3pi_2018.GetEntries()
        N_gen_3pi3pipi0_2018 = t_gen_3pi3pi_pi0_2018.GetEntries() + t_gen_3pipi0_3pi_2018.GetEntries()
        N_gen_3pi3pi2pi0_2018 = t_gen_3pi3pi2pi0_2018.GetEntries()
        N_gen_2018 = N_gen_3pi3pi_2018+N_gen_3pi3pipi0_2018+N_gen_3pi3pi2pi0_2018

        N_gen_3pi3pi = N_gen_3pi3pi_2016+N_gen_3pi3pi_2017+N_gen_3pi3pi_2018
        N_gen_3pi3pipi0 = N_gen_3pi3pipi0_2016+N_gen_3pi3pipi0_2017+N_gen_3pi3pipi0_2018
        N_gen_3pi3pi2pi0 = N_gen_3pi3pi2pi0_2016+N_gen_3pi3pi2pi0_2017+N_gen_3pi3pi2pi0_2018
        N_gen = N_gen_2016+N_gen_2017+N_gen_2018

        print("Ngen 3pi3pi = ", N_gen_3pi3pi)
        print("Ngen 3pi3pipi0 = ", N_gen_3pi3pipi0)
        print("Ngen 3pi3pi2pi0 = ", N_gen_3pi3pi2pi0)
        print("Ngen = ", N_gen)

    elif(species == 7):
        # fc_2016 = ROOT.TFileCollection("fc_2016", "fc_2016", "Files_on_grid/MC_D0Dps_2016.txt")
        # fc_2017 = ROOT.TFileCollection("fc_2017", "fc_2017", "Files_on_grid/MC_D0Dps_2017.txt")
        # fc_2018 = ROOT.TFileCollection("fc_2018", "fc_2018", "Files_on_grid/MC_D0Dps_2018.txt")

        # t_gen_2016 = ROOT.TChain("mc_ntuple/MCDecayTree")
        # t_gen_2017 = ROOT.TChain("mc_ntuple/MCDecayTree")
        # t_gen_2018 = ROOT.TChain("mc_ntuple/MCDecayTree")

        # t_gen_2016.AddFileInfoList(fc_2016.GetList())
        # t_gen_2017.AddFileInfoList(fc_2017.GetList())
        # t_gen_2018.AddFileInfoList(fc_2018.GetList())

        # N_gen = t_gen_2016.GetEntries() + t_gen_2017.GetEntries() + t_gen_2018.GetEntries()

        N_gen = 5350785+5782016+5763580

        print("N_gen = ", N_gen)
    elif(species == 100):
        N_gen_all = np.load('/panfs/felician/B2Ktautau/workflow/yield_estimates_inputs/eps_bkg_den.npy')

        n_BuDDKp = np.sum(N_gen_all[:3])
        n_BdDDKp = np.sum(N_gen_all[3])
        n_BsDDKp = np.sum(N_gen_all[4])
        n_BuDDK0 = np.sum(N_gen_all[5])
        n_BuDD = np.sum(N_gen_all[-2:])

        N_gen = [n_BuDDKp, n_BdDDKp, n_BsDDKp, n_BuDDK0, n_BuDD]

    np.save(f"/panfs/felician/B2Ktautau/workflow/selections_efficiency_tables/Species_{species}/N_gen.npy", N_gen)

    # Reconstructed w/ PIDGen vars tree
    fc_reco_pid_2016 = ROOT.TFileCollection("fc_reco_pid_2016", "fc_reco_pid_2016", f"/panfs/felician/B2Ktautau/workflow/PIDCalib/2016/Species_{species}/pid_corr.txt") 
    fc_reco_pid_2017 = ROOT.TFileCollection("fc_reco_pid_2017", "fc_reco_pid_2017", f"/panfs/felician/B2Ktautau/workflow/PIDCalib/2017/Species_{species}/pid_corr.txt") 
    fc_reco_pid_2018 = ROOT.TFileCollection("fc_reco_pid_2018", "fc_reco_pid_2018", f"/panfs/felician/B2Ktautau/workflow/PIDCalib/2018/Species_{species}/pid_corr.txt") 

    t_reco_pid_2016 = ROOT.TChain("DecayTree")
    t_reco_pid_2017 = ROOT.TChain("DecayTree")
    t_reco_pid_2018 = ROOT.TChain("DecayTree")

    t_reco_pid_2016.AddFileInfoList(fc_reco_pid_2016.GetList())
    t_reco_pid_2017.AddFileInfoList(fc_reco_pid_2017.GetList())
    t_reco_pid_2018.AddFileInfoList(fc_reco_pid_2018.GetList())

    if(species == 1):
        fc_comp_2016 = ROOT.TFileCollection("fc_comp_2016", "fc_comp_2016", "/panfs/felician/B2Ktautau/workflow/separate_reco_mc_components/2016/mc_components.txt")
        fc_comp_2017 = ROOT.TFileCollection("fc_comp_2017", "fc_comp_2017", "/panfs/felician/B2Ktautau/workflow/separate_reco_mc_components/2017/mc_components.txt")
        fc_comp_2018 = ROOT.TFileCollection("fc_comp_2018", "fc_comp_2018", "/panfs/felician/B2Ktautau/workflow/separate_reco_mc_components/2018/mc_components.txt")

        t_comp_2016 = ROOT.TChain("DecayTree")
        t_comp_2017 = ROOT.TChain("DecayTree")
        t_comp_2018 = ROOT.TChain("DecayTree")

        t_comp_2016.AddFileInfoList(fc_comp_2016.GetList())
        t_comp_2017.AddFileInfoList(fc_comp_2017.GetList())
        t_comp_2018.AddFileInfoList(fc_comp_2018.GetList())

        t_reco_pid_2016.AddFriend(t_comp_2016)
        t_reco_pid_2017.AddFriend(t_comp_2017)
        t_reco_pid_2018.AddFriend(t_comp_2018)

        N_reco_3pi3pi = t_reco_pid_2016.GetEntries(truthMatch+" && (component==0)") + t_reco_pid_2017.GetEntries(truthMatch+" && (component==0)") + t_reco_pid_2018.GetEntries(truthMatch+" && (component==0)")
        N_reco_3pi3pipi0 = t_reco_pid_2016.GetEntries(truthMatch+" && (component==1)") + t_reco_pid_2017.GetEntries(truthMatch+" && (component==1)") + t_reco_pid_2018.GetEntries(truthMatch+" && (component==1)")
        N_reco_3pi3pi2pi0 = t_reco_pid_2016.GetEntries(truthMatch+" && (component==2)") + t_reco_pid_2017.GetEntries(truthMatch+" && (component==2)") + t_reco_pid_2018.GetEntries(truthMatch+" && (component==2)")

        eps_reco_3pi3pi = N_reco_3pi3pi/N_gen
        up_reco_3pi3pi = ROOT.TEfficiency.Wilson(N_gen, N_reco_3pi3pi, 0.68, True)
        down_reco_3pi3pi = ROOT.TEfficiency.Wilson(N_gen, N_reco_3pi3pi, 0.68, False)
        eps_reco_3pi3pi_err = 0.5*(up_reco_3pi3pi - down_reco_3pi3pi)
        eps_reco_3pi3pi = ufloat(eps_reco_3pi3pi, eps_reco_3pi3pi_err)

        eps_reco_3pi3pipi0 = N_reco_3pi3pipi0/N_gen
        up_reco_3pi3pipi0 = ROOT.TEfficiency.Wilson(N_gen, N_reco_3pi3pipi0, 0.68, True)
        down_reco_3pi3pipi0 = ROOT.TEfficiency.Wilson(N_gen, N_reco_3pi3pipi0, 0.68, False)
        eps_reco_3pi3pipi0_err = 0.5*(up_reco_3pi3pipi0 - down_reco_3pi3pipi0)
        eps_reco_3pi3pipi0 = ufloat(eps_reco_3pi3pipi0, eps_reco_3pi3pipi0_err)

        eps_reco_3pi3pi2pi0 = N_reco_3pi3pi2pi0/N_gen
        up_reco_3pi3pi2pi0 = ROOT.TEfficiency.Wilson(N_gen, N_reco_3pi3pi2pi0, 0.68, True)
        down_reco_3pi3pi2pi0 = ROOT.TEfficiency.Wilson(N_gen, N_reco_3pi3pi2pi0, 0.68, False)
        eps_reco_3pi3pi2pi0_err = 0.5*(up_reco_3pi3pi2pi0 - down_reco_3pi3pi2pi0)
        eps_reco_3pi3pi2pi0 = ufloat(eps_reco_3pi3pi2pi0, eps_reco_3pi3pi2pi0_err)

        print("N_reco_3pi3pi = ", N_reco_3pi3pi)
        print("N_reco_3pi3pipi0 = ", N_reco_3pi3pipi0)
        print("N_reco_3pi3pi2pi0 = ", N_reco_3pi3pi2pi0)

    
    if(species == 100):
        fc_reco_pid_BuDDKp_2016 = ROOT.TFileCollection("fc_reco_pid_BuDDKp_2016", "fc_reco_pid_BuDDKp_2016", f"/panfs/felician/B2Ktautau/workflow/PIDCalib/2016/Species_100/pid_corr.txt") 
        fc_reco_pid_BuDDKp_2017 = ROOT.TFileCollection("fc_reco_pid_BuDDKp_2017", "fc_reco_pid_BuDDKp_2017", f"/panfs/felician/B2Ktautau/workflow/PIDCalib/2017/Species_100/pid_corr.txt") 
        fc_reco_pid_BuDDKp_2018 = ROOT.TFileCollection("fc_reco_pid_BuDDKp_2018", "fc_reco_pid_BuDDKp_2018", f"/panfs/felician/B2Ktautau/workflow/PIDCalib/2018/Species_100/pid_corr.txt") 

        t_reco_pid_BuDDKp_2016 = ROOT.TChain("DecayTree")
        t_reco_pid_BuDDKp_2017 = ROOT.TChain("DecayTree")
        t_reco_pid_BuDDKp_2018 = ROOT.TChain("DecayTree")

        t_reco_pid_BuDDKp_2016.AddFileInfoList(fc_reco_pid_BuDDKp_2016.GetList())
        t_reco_pid_BuDDKp_2017.AddFileInfoList(fc_reco_pid_BuDDKp_2017.GetList())
        t_reco_pid_BuDDKp_2018.AddFileInfoList(fc_reco_pid_BuDDKp_2018.GetList())     
    
        fc_reco_pid_BdDDKp_2016 = ROOT.TFileCollection("fc_reco_pid_BdDDKp_2016", "fc_reco_pid_BdDDKp_2016", f"/panfs/felician/B2Ktautau/workflow/PIDCalib/2016/Species_110/pid_corr.txt") 
        fc_reco_pid_BdDDKp_2017 = ROOT.TFileCollection("fc_reco_pid_BdDDKp_2017", "fc_reco_pid_BdDDKp_2017", f"/panfs/felician/B2Ktautau/workflow/PIDCalib/2017/Species_110/pid_corr.txt") 
        fc_reco_pid_BdDDKp_2018 = ROOT.TFileCollection("fc_reco_pid_BdDDKp_2018", "fc_reco_pid_BdDDKp_2018", f"/panfs/felician/B2Ktautau/workflow/PIDCalib/2018/Species_110/pid_corr.txt") 

        t_reco_pid_BdDDKp_2016 = ROOT.TChain("DecayTree")
        t_reco_pid_BdDDKp_2017 = ROOT.TChain("DecayTree")
        t_reco_pid_BdDDKp_2018 = ROOT.TChain("DecayTree")

        t_reco_pid_BdDDKp_2016.AddFileInfoList(fc_reco_pid_BdDDKp_2016.GetList())
        t_reco_pid_BdDDKp_2017.AddFileInfoList(fc_reco_pid_BdDDKp_2017.GetList())
        t_reco_pid_BdDDKp_2018.AddFileInfoList(fc_reco_pid_BdDDKp_2018.GetList())     

        fc_reco_pid_BsDDKp_2016 = ROOT.TFileCollection("fc_reco_pid_BsDDKp_2016", "fc_reco_pid_BsDDKp_2016", f"/panfs/felician/B2Ktautau/workflow/PIDCalib/2016/Species_120/pid_corr.txt") 
        fc_reco_pid_BsDDKp_2017 = ROOT.TFileCollection("fc_reco_pid_BsDDKp_2017", "fc_reco_pid_BsDDKp_2017", f"/panfs/felician/B2Ktautau/workflow/PIDCalib/2017/Species_120/pid_corr.txt") 
        fc_reco_pid_BsDDKp_2018 = ROOT.TFileCollection("fc_reco_pid_BsDDKp_2018", "fc_reco_pid_BsDDKp_2018", f"/panfs/felician/B2Ktautau/workflow/PIDCalib/2018/Species_120/pid_corr.txt") 

        t_reco_pid_BsDDKp_2016 = ROOT.TChain("DecayTree")
        t_reco_pid_BsDDKp_2017 = ROOT.TChain("DecayTree")
        t_reco_pid_BsDDKp_2018 = ROOT.TChain("DecayTree")

        t_reco_pid_BsDDKp_2016.AddFileInfoList(fc_reco_pid_BsDDKp_2016.GetList())
        t_reco_pid_BsDDKp_2017.AddFileInfoList(fc_reco_pid_BsDDKp_2017.GetList())
        t_reco_pid_BsDDKp_2018.AddFileInfoList(fc_reco_pid_BsDDKp_2018.GetList())     

        fc_reco_pid_BuDDK0_2016 = ROOT.TFileCollection("fc_reco_pid_BuDDK0_2016", "fc_reco_pid_BuDDK0_2016", f"/panfs/felician/B2Ktautau/workflow/PIDCalib/2016/Species_130/pid_corr.txt") 
        fc_reco_pid_BuDDK0_2017 = ROOT.TFileCollection("fc_reco_pid_BuDDK0_2017", "fc_reco_pid_BuDDK0_2017", f"/panfs/felician/B2Ktautau/workflow/PIDCalib/2017/Species_130/pid_corr.txt") 
        fc_reco_pid_BuDDK0_2018 = ROOT.TFileCollection("fc_reco_pid_BuDDK0_2018", "fc_reco_pid_BuDDK0_2018", f"/panfs/felician/B2Ktautau/workflow/PIDCalib/2018/Species_130/pid_corr.txt") 

        t_reco_pid_BuDDK0_2016 = ROOT.TChain("DecayTree")
        t_reco_pid_BuDDK0_2017 = ROOT.TChain("DecayTree")
        t_reco_pid_BuDDK0_2018 = ROOT.TChain("DecayTree")

        t_reco_pid_BuDDK0_2016.AddFileInfoList(fc_reco_pid_BuDDK0_2016.GetList())
        t_reco_pid_BuDDK0_2017.AddFileInfoList(fc_reco_pid_BuDDK0_2017.GetList())
        t_reco_pid_BuDDK0_2018.AddFileInfoList(fc_reco_pid_BuDDK0_2018.GetList())     

        fc_reco_pid_BuDD_2016 = ROOT.TFileCollection("fc_reco_pid_BuDD_2016", "fc_reco_pid_BuDD_2016", f"/panfs/felician/B2Ktautau/workflow/PIDCalib/2016/Species_150/pid_corr.txt") 
        fc_reco_pid_BuDD_2017 = ROOT.TFileCollection("fc_reco_pid_BuDD_2017", "fc_reco_pid_BuDD_2017", f"/panfs/felician/B2Ktautau/workflow/PIDCalib/2017/Species_150/pid_corr.txt") 
        fc_reco_pid_BuDD_2018 = ROOT.TFileCollection("fc_reco_pid_BuDD_2018", "fc_reco_pid_BuDD_2018", f"/panfs/felician/B2Ktautau/workflow/PIDCalib/2018/Species_150/pid_corr.txt") 

        t_reco_pid_BuDD_2016 = ROOT.TChain("DecayTree")
        t_reco_pid_BuDD_2017 = ROOT.TChain("DecayTree")
        t_reco_pid_BuDD_2018 = ROOT.TChain("DecayTree")

        t_reco_pid_BuDD_2016.AddFileInfoList(fc_reco_pid_BuDD_2016.GetList())
        t_reco_pid_BuDD_2017.AddFileInfoList(fc_reco_pid_BuDD_2017.GetList())
        t_reco_pid_BuDD_2018.AddFileInfoList(fc_reco_pid_BuDD_2018.GetList())     


    if(species == 100):
        truthMatch_BuDDKp_1 = truthMatch_cocktailMC("BuD0D0Kp")+" || "+truthMatch_cocktailMC("BuD0starD0Kp")+" || "+truthMatch_cocktailMC("BuD0D0starKp")+" || "+truthMatch_cocktailMC("BuD0starD0starKp")
        truthMatch_BuDDKp_2 = truthMatch_cocktailMC("BuDpDmKp")+" || "+truthMatch_cocktailMC("BuDpstarDmKp")+" || "+truthMatch_cocktailMC("BuDpDmstarKp")+" || "+truthMatch_cocktailMC("BuDpstarDmstarKp")
        truthMatch_BuDDKp_3 = truthMatch_cocktailMC("BuDsDsKp")+" || "+truthMatch_cocktailMC("BuDsstarDsKp")+" || "+truthMatch_cocktailMC("BuDsDsstarKp")+" || "+truthMatch_cocktailMC("BuDsstarDsstarKp")

        truthMatch_BdDDKp = truthMatch_cocktailMC("BdDmD0Kp")+" || "+truthMatch_cocktailMC("BdDmstarD0Kp")+" || "+truthMatch_cocktailMC("BdDmD0starKp")+" || "+truthMatch_cocktailMC("BdDmstarD0starKp")
        truthMatch_BsDDKp = truthMatch_cocktailMC("BsDsD0Kp")+" || "+truthMatch_cocktailMC("BsDsstarD0Kp")+" || "+truthMatch_cocktailMC("BsDsD0starKp")+" || "+truthMatch_cocktailMC("BsDsstarD0starKp")
        truthMatch_BuDDK0 = truthMatch_cocktailMC("BuD0DpK0")+" || "+truthMatch_cocktailMC("BuD0starDpK0")+" || "+truthMatch_cocktailMC("BuD0DpstarK0")+" || "+truthMatch_cocktailMC("BuD0starDpstarK0")
        truthMatch_BuDD = truthMatch_cocktailMC("BuD0Ds")+" || "+truthMatch_cocktailMC("BuD0starDs")+" || "+truthMatch_cocktailMC("BuD0Dsstar")+" || "+truthMatch_cocktailMC("BuD0starDsstar")+" || "+truthMatch_cocktailMC("BuD0Dp")+" || "+truthMatch_cocktailMC("BuD0starDp")+" || "+truthMatch_cocktailMC("BuD0Dpstar")+" || "+truthMatch_cocktailMC("BuD0starDpstar")

        N_reco_BuDDKp = t_reco_pid_BuDDKp_2016.GetEntries(truthMatch_BuDDKp_1) + t_reco_pid_BuDDKp_2017.GetEntries(truthMatch_BuDDKp_1) + t_reco_pid_BuDDKp_2018.GetEntries(truthMatch_BuDDKp_1) + t_reco_pid_BuDDKp_2016.GetEntries(truthMatch_BuDDKp_2) + t_reco_pid_BuDDKp_2017.GetEntries(truthMatch_BuDDKp_2) + t_reco_pid_BuDDKp_2018.GetEntries(truthMatch_BuDDKp_2) + t_reco_pid_BuDDKp_2016.GetEntries(truthMatch_BuDDKp_2) + t_reco_pid_BuDDKp_2017.GetEntries(truthMatch_BuDDKp_2) + t_reco_pid_BuDDKp_2018.GetEntries(truthMatch_BuDDKp_2)
        N_reco_BdDDKp = t_reco_pid_BdDDKp_2016.GetEntries(truthMatch_BdDDKp) + t_reco_pid_BdDDKp_2017.GetEntries(truthMatch_BdDDKp) + t_reco_pid_BdDDKp_2018.GetEntries(truthMatch_BdDDKp)
        N_reco_BsDDKp = t_reco_pid_BsDDKp_2016.GetEntries(truthMatch_BsDDKp) + t_reco_pid_BsDDKp_2017.GetEntries(truthMatch_BsDDKp) + t_reco_pid_BsDDKp_2018.GetEntries(truthMatch_BsDDKp)
        N_reco_BuDDK0 = t_reco_pid_BuDDK0_2016.GetEntries(truthMatch_BuDDK0) + t_reco_pid_BuDDK0_2017.GetEntries(truthMatch_BuDDK0) + t_reco_pid_BuDDK0_2018.GetEntries(truthMatch_BuDDK0)
        N_reco_BuDD = t_reco_pid_BuDD_2016.GetEntries(truthMatch_BuDD) + t_reco_pid_BuDD_2017.GetEntries(truthMatch_BuDD) + t_reco_pid_BuDD_2018.GetEntries(truthMatch_BuDD)

        up_reco_BuDDKp = ROOT.TEfficiency.Wilson(N_gen[0], N_reco_BuDDKp, 0.68, True)
        down_reco_BuDDKp = ROOT.TEfficiency.Wilson(N_gen[0], N_reco_BuDDKp, 0.68, False)
        eps_reco_err_BuDDKp = 0.5*(up_reco_BuDDKp - down_reco_BuDDKp)
        eps_reco_BuDDKp = ufloat(N_reco_BuDDKp/N_gen[0], eps_reco_err_BuDDKp)
                   
        up_reco_BdDDKp = ROOT.TEfficiency.Wilson(N_gen[1], N_reco_BdDDKp, 0.68, True)
        down_reco_BdDDKp = ROOT.TEfficiency.Wilson(N_gen[1], N_reco_BdDDKp, 0.68, False)
        eps_reco_err_BdDDKp = 0.5*(up_reco_BdDDKp - down_reco_BdDDKp)
        eps_reco_BdDDKp = ufloat(N_reco_BdDDKp/N_gen[1], eps_reco_err_BdDDKp)       
    
        up_reco_BsDDKp = ROOT.TEfficiency.Wilson(N_gen[2], N_reco_BsDDKp, 0.68, True)
        down_reco_BsDDKp = ROOT.TEfficiency.Wilson(N_gen[2], N_reco_BsDDKp, 0.68, False)
        eps_reco_err_BsDDKp = 0.5*(up_reco_BsDDKp- down_reco_BsDDKp)
        eps_reco_BsDDKp = ufloat(N_reco_BsDDKp/N_gen[2], eps_reco_err_BsDDKp)       
                                     
        up_reco_BuDDK0 = ROOT.TEfficiency.Wilson(N_gen[3], N_reco_BuDDK0, 0.68, True)
        down_reco_BuDDK0 = ROOT.TEfficiency.Wilson(N_gen[3], N_reco_BuDDK0, 0.68, False)
        eps_reco_err_BuDDK0 = 0.5*(up_reco_BuDDK0- down_reco_BuDDK0)
        eps_reco_BuDDK0 = ufloat(N_reco_BuDDK0/N_gen[3], eps_reco_err_BuDDK0)   

        up_reco_BuDD = ROOT.TEfficiency.Wilson(N_gen[4], N_reco_BuDD, 0.68, True)
        down_reco_BuDD = ROOT.TEfficiency.Wilson(N_gen[4], N_reco_BuDD, 0.68, False)
        eps_reco_err_BuDD = 0.5*(up_reco_BuDD- down_reco_BuDD)
        eps_reco_BuDD = ufloat(N_reco_BuDD/N_gen[4], eps_reco_err_BuDD)   

        eps_reco = [eps_reco_BuDDKp, eps_reco_BdDDKp, eps_reco_BsDDKp, eps_reco_BuDDK0, eps_reco_BuDD]
                   
    else:
        N_reco = t_reco_pid_2016.GetEntries(truthMatch) + t_reco_pid_2017.GetEntries(truthMatch) + t_reco_pid_2018.GetEntries(truthMatch)
        print("N_reco = ", N_reco)

        eps_reco = N_reco/N_gen
        up_reco = ROOT.TEfficiency.Wilson(N_gen, N_reco, 0.68, True)
        down_reco = ROOT.TEfficiency.Wilson(N_gen, N_reco, 0.68, False)
        eps_reco_err = 0.5*(up_reco - down_reco)
        eps_reco = ufloat(eps_reco, eps_reco_err)


    if(species == 100):
        reco_columns = ["Reconstruction efficiency II (\\%)"]
        reco_rows = ["$B^+ \\to D D K^+$", "$B^0 \\to D D K^+$", "$B^0_s \\to D D K^+$", "$B^+ \\to D D K^0$", "$B^+ \\to D D$"]

        reco_table = pd.DataFrame(index=reco_rows, columns=reco_columns)
        reco_table.loc["$B^+ \\to D D K^+$", "All years (\\%)"] = f"${eps_reco[0].nominal_value*100:.3f} \\pm {eps_reco[0].std_dev*100:.3f}$"
        reco_table.loc["$B^0 \\to D D K^+$", "All years (\\%)"] = f"${eps_reco[1].nominal_value*100:.3f} \\pm {eps_reco[1].std_dev*100:.3f}$"
        reco_table.loc["$B^0_s \\to D D K^+$", "All years (\\%)"] = f"${eps_reco[2].nominal_value*100:.3f} \\pm {eps_reco[2].std_dev*100:.3f}$"
        reco_table.loc["$B^+ \\to D D K^0$", "All years (\\%)"] = f"${eps_reco[3].nominal_value*100:.3f} \\pm {eps_reco[3].std_dev*100:.3f}$"
        reco_table.loc["$B^+ \\to D D$", "All years (\\%)"] = f"${eps_reco[4].nominal_value*100:.3f} \\pm {eps_reco[4].std_dev*100:.3f}$"
    
    else:
        reco_columns = ["Reconstruction efficiency (\\%)"]
        if(species == 1):
            reco_rows = ["$3\\pi 3\\pi$", "$3\\pi 3\\pi \\pi^0$", "$3\\pi 3\\pi 2\\pi^0$", "All"]
        else:
            reco_rows  = [""]

        reco_table = pd.DataFrame(index=reco_rows, columns=reco_columns)
        if(species == 1):
            reco_table.loc["$3\\pi 3\\pi$", "Reconstruction efficiency (\\%)"] = f"${eps_reco_3pi3pi.nominal_value*100:.2f} \\pm {eps_reco_3pi3pi.std_dev*100:.2f}$"
            reco_table.loc["$3\\pi 3\\pi \\pi^0$", "Reconstruction efficiency (\\%)"] = f"${eps_reco_3pi3pipi0.nominal_value*100:.2f} \\pm {eps_reco_3pi3pipi0.std_dev*100:.2f}$"
            reco_table.loc["$3\\pi 3\\pi 2\\pi^0$", "Reconstruction efficiency (\\%)"] = f"${eps_reco_3pi3pi2pi0.nominal_value*100:.2f} \\pm {eps_reco_3pi3pi2pi0.std_dev*100:.2f}$"
            reco_table.loc["All", "Reconstruction efficiency (\\%)"] = f"${eps_reco.nominal_value*100:.2f} \\pm {eps_reco.std_dev*100:.2f}$"
        else:
            reco_table.loc["", "Reconstruction efficiency (\\%)"] = f"${eps_reco.nominal_value*100:.3f} \\pm {eps_reco.std_dev*100:.3f}$"

    with open(f'/panfs/felician/B2Ktautau/workflow/selections_efficiency_tables/Species_{species}/reconstruction_efficiency_table.tex', 'w') as fout_reco:
        fout_reco.write(reco_table.to_latex())


    # Trigger efficiencies
    print("Trigger efficiency")
    if(species == 100):
        N_L0_BuDDKp = t_reco_pid_BuDDKp_2016.GetEntries(truthMatch_BuDDKp_1+" && "+L0_trigger) + t_reco_pid_BuDDKp_2017.GetEntries(truthMatch_BuDDKp_1+" && "+L0_trigger) + t_reco_pid_BuDDKp_2018.GetEntries(truthMatch_BuDDKp_1+" && "+L0_trigger) + t_reco_pid_BuDDKp_2016.GetEntries(truthMatch_BuDDKp_2+" && "+L0_trigger) + t_reco_pid_BuDDKp_2017.GetEntries(truthMatch_BuDDKp_2+" && "+L0_trigger) + t_reco_pid_BuDDKp_2018.GetEntries(truthMatch_BuDDKp_2+" && "+L0_trigger) + t_reco_pid_BuDDKp_2016.GetEntries(truthMatch_BuDDKp_3+" && "+L0_trigger) + t_reco_pid_BuDDKp_2017.GetEntries(truthMatch_BuDDKp_3+" && "+L0_trigger) + t_reco_pid_BuDDKp_2018.GetEntries(truthMatch_BuDDKp_3+" && "+L0_trigger)
        N_L0_BdDDKp = t_reco_pid_BdDDKp_2016.GetEntries(truthMatch_BdDDKp+" && "+L0_trigger) + t_reco_pid_BdDDKp_2017.GetEntries(truthMatch_BdDDKp+" && "+L0_trigger) + t_reco_pid_BdDDKp_2018.GetEntries(truthMatch_BdDDKp+" && "+L0_trigger)
        N_L0_BsDDKp = t_reco_pid_BsDDKp_2016.GetEntries(truthMatch_BsDDKp+" && "+L0_trigger) + t_reco_pid_BsDDKp_2017.GetEntries(truthMatch_BsDDKp+" && "+L0_trigger) + t_reco_pid_BsDDKp_2018.GetEntries(truthMatch_BsDDKp+" && "+L0_trigger)
        N_L0_BuDDK0 = t_reco_pid_BuDDK0_2016.GetEntries(truthMatch_BuDDK0+" && "+L0_trigger) + t_reco_pid_BuDDK0_2017.GetEntries(truthMatch_BuDDK0+" && "+L0_trigger) + t_reco_pid_BuDDK0_2018.GetEntries(truthMatch_BuDDK0+" && "+L0_trigger)
        N_L0_BuDD = t_reco_pid_BuDD_2016.GetEntries(truthMatch_BuDD+" && "+L0_trigger) + t_reco_pid_BuDD_2017.GetEntries(truthMatch_BuDD+" && "+L0_trigger) + t_reco_pid_BuDD_2018.GetEntries(truthMatch_BuDD+" && "+L0_trigger)

        up_L0_BuDDKp = ROOT.TEfficiency.Wilson(N_reco_BuDDKp, N_L0_BuDDKp, 0.68, True)
        down_L0_BuDDKp = ROOT.TEfficiency.Wilson(N_reco_BuDDKp, N_L0_BuDDKp, 0.68, False)
        eps_L0_BuDDKp_err = 0.5*(up_L0_BuDDKp - down_L0_BuDDKp)
        eps_L0_BuDDKp = ufloat(N_L0_BuDDKp/N_reco_BuDDKp, eps_L0_BuDDKp_err)

        up_L0_BdDDKp = ROOT.TEfficiency.Wilson(N_reco_BdDDKp, N_L0_BdDDKp, 0.68, True)
        down_L0_BdDDKp = ROOT.TEfficiency.Wilson(N_reco_BdDDKp, N_L0_BdDDKp, 0.68, False)
        eps_L0_BdDDKp_err = 0.5*(up_L0_BdDDKp - down_L0_BdDDKp)
        eps_L0_BdDDKp = ufloat(N_L0_BdDDKp/N_reco_BdDDKp, eps_L0_BdDDKp_err)

        up_L0_BsDDKp = ROOT.TEfficiency.Wilson(N_reco_BsDDKp, N_L0_BsDDKp, 0.68, True)
        down_L0_BsDDKp = ROOT.TEfficiency.Wilson(N_reco_BsDDKp, N_L0_BsDDKp, 0.68, False)
        eps_L0_BsDDKp_err = 0.5*(up_L0_BsDDKp - down_L0_BsDDKp)
        eps_L0_BsDDKp = ufloat(N_L0_BsDDKp/N_reco_BsDDKp, eps_L0_BsDDKp_err)

        up_L0_BuDDK0 = ROOT.TEfficiency.Wilson(N_reco_BuDDK0, N_L0_BuDDK0, 0.68, True)
        down_L0_BuDDK0 = ROOT.TEfficiency.Wilson(N_reco_BuDDK0, N_L0_BuDDK0, 0.68, False)
        eps_L0_BuDDK0_err = 0.5*(up_L0_BuDDK0 - down_L0_BuDDK0)
        eps_L0_BuDDK0 = ufloat(N_L0_BuDDK0/N_reco_BuDDK0, eps_L0_BuDDK0_err)

        up_L0_BuDD = ROOT.TEfficiency.Wilson(N_reco_BuDD, N_L0_BuDD, 0.68, True)
        down_L0_BuDD = ROOT.TEfficiency.Wilson(N_reco_BuDD, N_L0_BuDD, 0.68, False)
        eps_L0_BuDD_err = 0.5*(up_L0_BuDD - down_L0_BuDD)
        eps_L0_BuDD = ufloat(N_L0_BuDD/N_reco_BuDD, eps_L0_BuDD_err)

        eps_L0 = [eps_L0_BuDDKp, eps_L0_BdDDKp, eps_L0_BsDDKp, eps_L0_BuDDK0, eps_L0_BuDD]

        N_HLT1_BuDDKp = t_reco_pid_BuDDKp_2016.GetEntries(truthMatch_BuDDKp_1+" && "+L0_trigger+" && "+HLT1_trigger) + t_reco_pid_BuDDKp_2017.GetEntries(truthMatch_BuDDKp_1+" && "+L0_trigger+" && "+HLT1_trigger) + t_reco_pid_BuDDKp_2018.GetEntries(truthMatch_BuDDKp_1+" && "+L0_trigger+" && "+HLT1_trigger) + t_reco_pid_BuDDKp_2016.GetEntries(truthMatch_BuDDKp_2+" && "+L0_trigger+" && "+HLT1_trigger) + t_reco_pid_BuDDKp_2017.GetEntries(truthMatch_BuDDKp_2+" && "+L0_trigger+" && "+HLT1_trigger) + t_reco_pid_BuDDKp_2018.GetEntries(truthMatch_BuDDKp_2+" && "+L0_trigger+" && "+HLT1_trigger) + t_reco_pid_BuDDKp_2016.GetEntries(truthMatch_BuDDKp_3+" && "+L0_trigger+" && "+HLT1_trigger) + t_reco_pid_BuDDKp_2017.GetEntries(truthMatch_BuDDKp_3+" && "+L0_trigger+" && "+HLT1_trigger) + t_reco_pid_BuDDKp_2018.GetEntries(truthMatch_BuDDKp_3+" && "+L0_trigger+" && "+HLT1_trigger)
        N_HLT1_BdDDKp = t_reco_pid_BdDDKp_2016.GetEntries(truthMatch_BdDDKp+" && "+L0_trigger+" && "+HLT1_trigger) + t_reco_pid_BdDDKp_2017.GetEntries(truthMatch_BdDDKp+" && "+L0_trigger+" && "+HLT1_trigger) + t_reco_pid_BdDDKp_2018.GetEntries(truthMatch_BdDDKp+" && "+L0_trigger+" && "+HLT1_trigger)
        N_HLT1_BsDDKp = t_reco_pid_BsDDKp_2016.GetEntries(truthMatch_BsDDKp+" && "+L0_trigger+" && "+HLT1_trigger) + t_reco_pid_BsDDKp_2017.GetEntries(truthMatch_BsDDKp+" && "+L0_trigger+" && "+HLT1_trigger) + t_reco_pid_BsDDKp_2018.GetEntries(truthMatch_BsDDKp+" && "+L0_trigger+" && "+HLT1_trigger)
        N_HLT1_BuDDK0 = t_reco_pid_BuDDK0_2016.GetEntries(truthMatch_BuDDK0+" && "+L0_trigger+" && "+HLT1_trigger) + t_reco_pid_BuDDK0_2017.GetEntries(truthMatch_BuDDK0+" && "+L0_trigger+" && "+HLT1_trigger) + t_reco_pid_BuDDK0_2018.GetEntries(truthMatch_BuDDK0+" && "+L0_trigger+" && "+HLT1_trigger)
        N_HLT1_BuDD = t_reco_pid_BuDD_2016.GetEntries(truthMatch_BuDD+" && "+L0_trigger+" && "+HLT1_trigger) + t_reco_pid_BuDD_2017.GetEntries(truthMatch_BuDD+" && "+L0_trigger+" && "+HLT1_trigger) + t_reco_pid_BuDD_2018.GetEntries(truthMatch_BuDD+" && "+L0_trigger+" && "+HLT1_trigger)

        up_HLT1_BuDDKp = ROOT.TEfficiency.Wilson(N_reco_BuDDKp, N_HLT1_BuDDKp, 0.68, True)
        down_HLT1_BuDDKp = ROOT.TEfficiency.Wilson(N_reco_BuDDKp, N_HLT1_BuDDKp, 0.68, False)
        eps_HLT1_BuDDKp_err = 0.5*(up_HLT1_BuDDKp - down_HLT1_BuDDKp)
        eps_HLT1_BuDDKp = ufloat(N_HLT1_BuDDKp/N_reco_BuDDKp, eps_HLT1_BuDDKp_err)

        up_HLT1_BdDDKp = ROOT.TEfficiency.Wilson(N_reco_BdDDKp, N_HLT1_BdDDKp, 0.68, True)
        down_HLT1_BdDDKp = ROOT.TEfficiency.Wilson(N_reco_BdDDKp, N_HLT1_BdDDKp, 0.68, False)
        eps_HLT1_BdDDKp_err = 0.5*(up_HLT1_BdDDKp - down_HLT1_BdDDKp)
        eps_HLT1_BdDDKp = ufloat(N_HLT1_BdDDKp/N_reco_BdDDKp, eps_HLT1_BdDDKp_err)

        up_HLT1_BsDDKp = ROOT.TEfficiency.Wilson(N_reco_BsDDKp, N_HLT1_BsDDKp, 0.68, True)
        down_HLT1_BsDDKp = ROOT.TEfficiency.Wilson(N_reco_BsDDKp, N_HLT1_BsDDKp, 0.68, False)
        eps_HLT1_BsDDKp_err = 0.5*(up_HLT1_BsDDKp - down_HLT1_BsDDKp)
        eps_HLT1_BsDDKp = ufloat(N_HLT1_BsDDKp/N_reco_BsDDKp, eps_HLT1_BsDDKp_err)

        up_HLT1_BuDDK0 = ROOT.TEfficiency.Wilson(N_reco_BuDDK0, N_HLT1_BuDDK0, 0.68, True)
        down_HLT1_BuDDK0 = ROOT.TEfficiency.Wilson(N_reco_BuDDK0, N_HLT1_BuDDK0, 0.68, False)
        eps_HLT1_BuDDK0_err = 0.5*(up_HLT1_BuDDK0 - down_HLT1_BuDDK0)
        eps_HLT1_BuDDK0 = ufloat(N_HLT1_BuDDK0/N_reco_BuDDK0, eps_HLT1_BuDDK0_err)

        up_HLT1_BuDD = ROOT.TEfficiency.Wilson(N_reco_BuDD, N_HLT1_BuDD, 0.68, True)
        down_HLT1_BuDD = ROOT.TEfficiency.Wilson(N_reco_BuDD, N_HLT1_BuDD, 0.68, False)
        eps_HLT1_BuDD_err = 0.5*(up_HLT1_BuDD - down_HLT1_BuDD)
        eps_HLT1_BuDD = ufloat(N_HLT1_BuDD/N_reco_BuDD, eps_HLT1_BuDD_err)

        eps_HLT1 = [eps_HLT1_BuDDKp, eps_HLT1_BdDDKp, eps_HLT1_BsDDKp, eps_HLT1_BuDDK0, eps_HLT1_BuDD]

        N_trigger_BuDDKp = t_reco_pid_BuDDKp_2016.GetEntries(truthMatch_BuDDKp_1+" && "+trigger) + t_reco_pid_BuDDKp_2017.GetEntries(truthMatch_BuDDKp_1+" && "+trigger) + t_reco_pid_BuDDKp_2018.GetEntries(truthMatch_BuDDKp_1+" && "+trigger) + t_reco_pid_BuDDKp_2016.GetEntries(truthMatch_BuDDKp_2+" && "+trigger) + t_reco_pid_BuDDKp_2017.GetEntries(truthMatch_BuDDKp_2+" && "+trigger) + t_reco_pid_BuDDKp_2018.GetEntries(truthMatch_BuDDKp_2+" && "+trigger) + t_reco_pid_BuDDKp_2016.GetEntries(truthMatch_BuDDKp_3+" && "+trigger) + t_reco_pid_BuDDKp_2017.GetEntries(truthMatch_BuDDKp_3+" && "+trigger) + t_reco_pid_BuDDKp_2018.GetEntries(truthMatch_BuDDKp_3+" && "+trigger)
        N_trigger_BdDDKp = t_reco_pid_BdDDKp_2016.GetEntries(truthMatch_BdDDKp+" && "+trigger) + t_reco_pid_BdDDKp_2017.GetEntries(truthMatch_BdDDKp+" && "+trigger) + t_reco_pid_BdDDKp_2018.GetEntries(truthMatch_BdDDKp+" && "+trigger)
        N_trigger_BsDDKp = t_reco_pid_BsDDKp_2016.GetEntries(truthMatch_BsDDKp+" && "+trigger) + t_reco_pid_BsDDKp_2017.GetEntries(truthMatch_BsDDKp+" && "+trigger) + t_reco_pid_BsDDKp_2018.GetEntries(truthMatch_BsDDKp+" && "+trigger)
        N_trigger_BuDDK0 = t_reco_pid_BuDDK0_2016.GetEntries(truthMatch_BuDDK0+" && "+trigger) + t_reco_pid_BuDDK0_2017.GetEntries(truthMatch_BuDDK0+" && "+trigger) + t_reco_pid_BuDDK0_2018.GetEntries(truthMatch_BuDDK0+" && "+trigger)
        N_trigger_BuDD = t_reco_pid_BuDD_2016.GetEntries(truthMatch_BuDD+" && "+trigger) + t_reco_pid_BuDD_2017.GetEntries(truthMatch_BuDD+" && "+trigger) + t_reco_pid_BuDD_2018.GetEntries(truthMatch_BuDD+" && "+trigger)

        up_trigger_BuDDKp = ROOT.TEfficiency.Wilson(N_reco_BuDDKp, N_trigger_BuDDKp, 0.68, True)
        down_trigger_BuDDKp = ROOT.TEfficiency.Wilson(N_reco_BuDDKp, N_trigger_BuDDKp, 0.68, False)
        eps_trigger_BuDDKp_err = 0.5*(up_trigger_BuDDKp - down_trigger_BuDDKp)
        eps_trigger_BuDDKp = ufloat(N_trigger_BuDDKp/N_reco_BuDDKp, eps_trigger_BuDDKp_err)

        up_trigger_BdDDKp = ROOT.TEfficiency.Wilson(N_reco_BdDDKp, N_trigger_BdDDKp, 0.68, True)
        down_trigger_BdDDKp = ROOT.TEfficiency.Wilson(N_reco_BdDDKp, N_trigger_BdDDKp, 0.68, False)
        eps_trigger_BdDDKp_err = 0.5*(up_trigger_BdDDKp - down_trigger_BdDDKp)
        eps_trigger_BdDDKp = ufloat(N_trigger_BdDDKp/N_reco_BdDDKp, eps_trigger_BdDDKp_err)

        up_trigger_BsDDKp = ROOT.TEfficiency.Wilson(N_reco_BsDDKp, N_trigger_BsDDKp, 0.68, True)
        down_trigger_BsDDKp = ROOT.TEfficiency.Wilson(N_reco_BsDDKp, N_trigger_BsDDKp, 0.68, False)
        eps_trigger_BsDDKp_err = 0.5*(up_trigger_BsDDKp - down_trigger_BsDDKp)
        eps_trigger_BsDDKp = ufloat(N_trigger_BsDDKp/N_reco_BsDDKp, eps_trigger_BsDDKp_err)

        up_trigger_BuDDK0 = ROOT.TEfficiency.Wilson(N_reco_BuDDK0, N_trigger_BuDDK0, 0.68, True)
        down_trigger_BuDDK0 = ROOT.TEfficiency.Wilson(N_reco_BuDDK0, N_trigger_BuDDK0, 0.68, False)
        eps_trigger_BuDDK0_err = 0.5*(up_trigger_BuDDK0 - down_trigger_BuDDK0)
        eps_trigger_BuDDK0 = ufloat(N_trigger_BuDDK0/N_reco_BuDDK0, eps_trigger_BuDDK0_err)

        up_trigger_BuDD = ROOT.TEfficiency.Wilson(N_reco_BuDD, N_trigger_BuDD, 0.68, True)
        down_trigger_BuDD = ROOT.TEfficiency.Wilson(N_reco_BuDD, N_trigger_BuDD, 0.68, False)
        eps_trigger_BuDD_err = 0.5*(up_trigger_BuDD - down_trigger_BuDD)
        eps_trigger_BuDD = ufloat(N_trigger_BuDD/N_reco_BuDD, eps_trigger_BuDD_err)

        eps_trigger = [eps_trigger_BuDDKp, eps_trigger_BdDDKp, eps_trigger_BsDDKp, eps_trigger_BuDDK0, eps_trigger_BuDD]

        trigger_columns = ["L0 (\\%)", "L0+HLT1 (\\%)", "Trigger (\\%)"]
        trigger_rows = ["$B^+ \\to D D K^+$", "$B^0 \\to D D K^+$", "$B^0_s \\to D D K^+$", "$B^+ \\to D D K^0$", "$B^+ \\to D D$"]

        trigger_table = pd.DataFrame(index=trigger_rows, columns=trigger_columns)
        trigger_table.loc["$B^+ \\to D D K^+$", "L0 (\\%)"] = f"${eps_L0[0].nominal_value*100:.2f} \\pm {eps_L0[0].std_dev*100:.2f}$"
        trigger_table.loc["$B^0 \\to D D K^+$", "L0 (\\%)"] = f"${eps_L0[1].nominal_value*100:.2f} \\pm {eps_L0[1].std_dev*100:.2f}$"
        trigger_table.loc["$B^0_s \\to D D K^+$", "L0 (\\%)"] = f"${eps_L0[2].nominal_value*100:.2f} \\pm {eps_L0[2].std_dev*100:.2f}$"
        trigger_table.loc["$B^+ \\to D D K^0$", "L0 (\\%)"] = f"${eps_L0[3].nominal_value*100:.2f} \\pm {eps_L0[3].std_dev*100:.2f}$"
        trigger_table.loc["$B^+ \\to D D$", "L0 (\\%)"] = f"${eps_L0[4].nominal_value*100:.2f} \\pm {eps_L0[4].std_dev*100:.2f}$"

        trigger_table.loc["$B^+ \\to D D K^+$", "L0+HLT1 (\\%)"] = f"${eps_HLT1[0].nominal_value*100:.2f} \\pm {eps_HLT1[0].std_dev*100:.2f}$"
        trigger_table.loc["$B^0 \\to D D K^+$", "L0+HLT1 (\\%)"] = f"${eps_HLT1[1].nominal_value*100:.2f} \\pm {eps_HLT1[1].std_dev*100:.2f}$"
        trigger_table.loc["$B^0_s \\to D D K^+$", "L0+HLT1 (\\%)"] = f"${eps_HLT1[2].nominal_value*100:.2f} \\pm {eps_HLT1[2].std_dev*100:.2f}$"
        trigger_table.loc["$B^+ \\to D D K^0$", "L0+HLT1 (\\%)"] = f"${eps_HLT1[3].nominal_value*100:.2f} \\pm {eps_HLT1[3].std_dev*100:.2f}$"
        trigger_table.loc["$B^+ \\to D D$", "L0+HLT1 (\\%)"] = f"${eps_HLT1[4].nominal_value*100:.2f} \\pm {eps_HLT1[4].std_dev*100:.2f}$"

        trigger_table.loc["$B^+ \\to D D K^+$", "Trigger (\\%)"] = f"${eps_trigger[0].nominal_value*100:.2f} \\pm {eps_trigger[0].std_dev*100:.2f}$"
        trigger_table.loc["$B^0 \\to D D K^+$", "Trigger (\\%)"] = f"${eps_trigger[1].nominal_value*100:.2f} \\pm {eps_trigger[1].std_dev*100:.2f}$"
        trigger_table.loc["$B^0_s \\to D D K^+$", "Trigger (\\%)"] = f"${eps_trigger[2].nominal_value*100:.2f} \\pm {eps_trigger[2].std_dev*100:.2f}$"
        trigger_table.loc["$B^+ \\to D D K^0$", "Trigger (\\%)"] = f"${eps_trigger[3].nominal_value*100:.2f} \\pm {eps_trigger[3].std_dev*100:.2f}$"
        trigger_table.loc["$B^+ \\to D D$", "Trigger (\\%)"] = f"${eps_trigger[4].nominal_value*100:.2f} \\pm {eps_trigger[4].std_dev*100:.2f}$"

    else:
        if(species == 1):
            N_L0_trigger_3pi3pi = t_reco_pid_2016.GetEntries(truthMatch+" && "+L0_trigger+" && (component == 0)") + t_reco_pid_2017.GetEntries(truthMatch+" && "+L0_trigger+" && (component == 0)") + t_reco_pid_2018.GetEntries(truthMatch+" && "+L0_trigger+" && (component == 0)")
            N_L0_trigger_3pi3pipi0 = t_reco_pid_2016.GetEntries(truthMatch+" && "+L0_trigger+" && (component == 1)") + t_reco_pid_2017.GetEntries(truthMatch+" && "+L0_trigger+" && (component == 1)") + t_reco_pid_2018.GetEntries(truthMatch+" && "+L0_trigger+" && (component == 1)")
            N_L0_trigger_3pi3pi2pi0 = t_reco_pid_2016.GetEntries(truthMatch+" && "+L0_trigger+" && (component == 2)") + t_reco_pid_2017.GetEntries(truthMatch+" && "+L0_trigger+" && (component == 2)") + t_reco_pid_2018.GetEntries(truthMatch+" && "+L0_trigger+" && (component == 2)")

            eps_L0_3pi3pi = N_L0_trigger_3pi3pi/N_reco_3pi3pi
            up_L0_3pi3pi = ROOT.TEfficiency.Wilson(N_reco_3pi3pi, N_L0_trigger_3pi3pi, 0.68, True)
            down_L0_3pi3pi = ROOT.TEfficiency.Wilson(N_reco_3pi3pi, N_L0_trigger_3pi3pi, 0.68, False)
            eps_L0_3pi3pi_err = 0.5*(up_L0_3pi3pi - down_L0_3pi3pi)
            eps_L0_3pi3pi = ufloat(eps_L0_3pi3pi, eps_L0_3pi3pi_err)

            eps_L0_3pi3pipi0 = N_L0_trigger_3pi3pipi0/N_reco_3pi3pipi0
            up_L0_3pi3pipi0 = ROOT.TEfficiency.Wilson(N_reco_3pi3pipi0, N_L0_trigger_3pi3pipi0, 0.68, True)
            down_L0_3pi3pipi0 = ROOT.TEfficiency.Wilson(N_reco_3pi3pipi0, N_L0_trigger_3pi3pipi0, 0.68, False)
            eps_L0_3pi3pipi0_err = 0.5*(up_L0_3pi3pipi0 - down_L0_3pi3pipi0)
            eps_L0_3pi3pipi0 = ufloat(eps_L0_3pi3pipi0, eps_L0_3pi3pipi0_err)

            eps_L0_3pi3pi2pi0 = N_L0_trigger_3pi3pi2pi0/N_reco_3pi3pi2pi0
            up_L0_3pi3pi2pi0 = ROOT.TEfficiency.Wilson(N_reco_3pi3pi2pi0, N_L0_trigger_3pi3pi2pi0, 0.68, True)
            down_L0_3pi3pi2pi0 = ROOT.TEfficiency.Wilson(N_reco_3pi3pi2pi0, N_L0_trigger_3pi3pi2pi0, 0.68, False)
            eps_L0_3pi3pi2pi0_err = 0.5*(up_L0_3pi3pi2pi0 - down_L0_3pi3pi2pi0)
            eps_L0_3pi3pi2pi0 = ufloat(eps_L0_3pi3pi2pi0, eps_L0_3pi3pi2pi0_err)

            N_HLT1_trigger_3pi3pi = t_reco_pid_2016.GetEntries(truthMatch+" && "+L0_trigger+" && "+HLT1_trigger+" && (component == 0)") + t_reco_pid_2017.GetEntries(truthMatch+" && "+L0_trigger+" && "+HLT1_trigger+" && (component == 0)") + t_reco_pid_2018.GetEntries(truthMatch+" && "+L0_trigger+" && "+HLT1_trigger+" && (component == 0)")
            N_HLT1_trigger_3pi3pipi0 = t_reco_pid_2016.GetEntries(truthMatch+" && "+L0_trigger+" && "+HLT1_trigger+" && (component == 1)") + t_reco_pid_2017.GetEntries(truthMatch+" && "+L0_trigger+" && "+HLT1_trigger+" && (component == 1)") + t_reco_pid_2018.GetEntries(truthMatch+" && "+L0_trigger+" && "+HLT1_trigger+" && (component == 1)")
            N_HLT1_trigger_3pi3pi2pi0 = t_reco_pid_2016.GetEntries(truthMatch+" && "+L0_trigger+" && "+HLT1_trigger+" && (component == 2)") + t_reco_pid_2017.GetEntries(truthMatch+" && "+L0_trigger+" && "+HLT1_trigger+" && (component == 2)") + t_reco_pid_2018.GetEntries(truthMatch+" && "+L0_trigger+" && "+HLT1_trigger+" && (component == 2)")

            eps_HLT1_3pi3pi = N_HLT1_trigger_3pi3pi/N_reco_3pi3pi
            up_HLT1_3pi3pi = ROOT.TEfficiency.Wilson(N_reco_3pi3pi, N_HLT1_trigger_3pi3pi, 0.68, True)
            down_HLT1_3pi3pi = ROOT.TEfficiency.Wilson(N_reco_3pi3pi, N_HLT1_trigger_3pi3pi, 0.68, False)
            eps_HLT1_3pi3pi_err = 0.5*(up_HLT1_3pi3pi - down_HLT1_3pi3pi)
            eps_HLT1_3pi3pi = ufloat(eps_HLT1_3pi3pi, eps_HLT1_3pi3pi_err)

            eps_HLT1_3pi3pipi0 = N_HLT1_trigger_3pi3pipi0/N_reco_3pi3pipi0
            up_HLT1_3pi3pipi0 = ROOT.TEfficiency.Wilson(N_reco_3pi3pipi0, N_HLT1_trigger_3pi3pipi0, 0.68, True)
            down_HLT1_3pi3pipi0 = ROOT.TEfficiency.Wilson(N_reco_3pi3pipi0, N_HLT1_trigger_3pi3pipi0, 0.68, False)
            eps_HLT1_3pi3pipi0_err = 0.5*(up_HLT1_3pi3pipi0 - down_HLT1_3pi3pipi0)
            eps_HLT1_3pi3pipi0 = ufloat(eps_HLT1_3pi3pipi0, eps_HLT1_3pi3pipi0_err)

            eps_HLT1_3pi3pi2pi0 = N_HLT1_trigger_3pi3pi2pi0/N_reco_3pi3pi2pi0
            up_HLT1_3pi3pi2pi0 = ROOT.TEfficiency.Wilson(N_reco_3pi3pi2pi0, N_HLT1_trigger_3pi3pi2pi0, 0.68, True)
            down_HLT1_3pi3pi2pi0 = ROOT.TEfficiency.Wilson(N_reco_3pi3pi2pi0, N_HLT1_trigger_3pi3pi2pi0, 0.68, False)
            eps_HLT1_3pi3pi2pi0_err = 0.5*(up_HLT1_3pi3pi2pi0 - down_HLT1_3pi3pi2pi0)
            eps_HLT1_3pi3pi2pi0 = ufloat(eps_HLT1_3pi3pi2pi0, eps_HLT1_3pi3pi2pi0_err)

            N_trigger_3pi3pi = t_reco_pid_2016.GetEntries(truthMatch+" && "+trigger+" && (component==0)") + t_reco_pid_2017.GetEntries(truthMatch+" && "+trigger+" && (component==0)") + t_reco_pid_2018.GetEntries(truthMatch+" && "+trigger+" && (component==0)")
            N_trigger_3pi3pipi0 = t_reco_pid_2016.GetEntries(truthMatch+" && "+trigger+" && (component==1)") + t_reco_pid_2017.GetEntries(truthMatch+" && "+trigger+" && (component==1)") + t_reco_pid_2018.GetEntries(truthMatch+" && "+trigger+" && (component==1)")
            N_trigger_3pi3pi2pi0 = t_reco_pid_2016.GetEntries(truthMatch+" && "+trigger+" && (component==2)") + t_reco_pid_2017.GetEntries(truthMatch+" && "+trigger+" && (component==2)") + t_reco_pid_2018.GetEntries(truthMatch+" && "+trigger+" && (component==2)")

            eps_trigger_3pi3pi = N_trigger_3pi3pi/N_reco_3pi3pi
            up_trigger_3pi3pi = ROOT.TEfficiency.Wilson(N_reco_3pi3pi, N_trigger_3pi3pi, 0.68, True)
            down_trigger_3pi3pi = ROOT.TEfficiency.Wilson(N_reco_3pi3pi, N_trigger_3pi3pi, 0.68, False)
            eps_trigger_3pi3pi_err = 0.5*(up_trigger_3pi3pi - down_trigger_3pi3pi)
            eps_trigger_3pi3pi = ufloat(eps_trigger_3pi3pi, eps_trigger_3pi3pi_err)

            eps_trigger_3pi3pipi0 = N_trigger_3pi3pipi0/N_reco_3pi3pipi0
            up_trigger_3pi3pipi0 = ROOT.TEfficiency.Wilson(N_reco_3pi3pipi0, N_trigger_3pi3pipi0, 0.68, True)
            down_trigger_3pi3pipi0 = ROOT.TEfficiency.Wilson(N_reco_3pi3pipi0, N_trigger_3pi3pipi0, 0.68, False)
            eps_trigger_3pi3pipi0_err = 0.5*(up_trigger_3pi3pipi0 - down_trigger_3pi3pipi0)
            eps_trigger_3pi3pipi0 = ufloat(eps_trigger_3pi3pipi0, eps_trigger_3pi3pipi0_err)

            eps_trigger_3pi3pi2pi0 = N_trigger_3pi3pi2pi0/N_reco_3pi3pi2pi0
            up_trigger_3pi3pi2pi0 = ROOT.TEfficiency.Wilson(N_reco_3pi3pi2pi0, N_trigger_3pi3pi2pi0, 0.68, True)
            down_trigger_3pi3pi2pi0 = ROOT.TEfficiency.Wilson(N_reco_3pi3pi2pi0, N_trigger_3pi3pi2pi0, 0.68, False)
            eps_trigger_3pi3pi2pi0_err = 0.5*(up_trigger_3pi3pi2pi0 - down_trigger_3pi3pi2pi0)
            eps_trigger_3pi3pi2pi0 = ufloat(eps_trigger_3pi3pi2pi0, eps_trigger_3pi3pi2pi0_err)

        N_L0_trigger = t_reco_pid_2016.GetEntries(truthMatch+" && "+L0_trigger) + t_reco_pid_2017.GetEntries(truthMatch+" && "+L0_trigger) + t_reco_pid_2018.GetEntries(truthMatch+" && "+L0_trigger)
        N_HLT1_trigger = t_reco_pid_2016.GetEntries(truthMatch+" && "+L0_trigger+" && "+HLT1_trigger) + t_reco_pid_2017.GetEntries(truthMatch+" && "+L0_trigger+" && "+HLT1_trigger) + t_reco_pid_2018.GetEntries(truthMatch+" && "+L0_trigger+" && "+HLT1_trigger)
        N_trigger = t_reco_pid_2016.GetEntries(truthMatch+" && "+trigger) + t_reco_pid_2017.GetEntries(truthMatch+" && "+trigger) + t_reco_pid_2018.GetEntries(truthMatch+" && "+trigger)

        print("L0 = ", N_L0_trigger)
        print("L0+HLT1 = ", N_HLT1_trigger)
        print("Trigger = ", N_trigger)

        eps_L0 = N_L0_trigger/N_reco
        up_L0 = ROOT.TEfficiency.Wilson(N_reco, N_L0_trigger, 0.68, True)
        down_L0 = ROOT.TEfficiency.Wilson(N_reco, N_L0_trigger, 0.68, False)
        eps_L0_err = 0.5*(up_L0 - down_L0)
        eps_L0 = ufloat(eps_L0, eps_L0_err)

        eps_HLT1 = N_HLT1_trigger/N_reco
        up_HLT1 = ROOT.TEfficiency.Wilson(N_reco, N_HLT1_trigger, 0.68, True)
        down_HLT1 = ROOT.TEfficiency.Wilson(N_reco, N_HLT1_trigger, 0.68, False)
        eps_HLT1_err = 0.5*(up_HLT1 - down_HLT1)
        eps_HLT1 = ufloat(eps_HLT1, eps_HLT1_err)

        eps_trigger = N_trigger/N_reco
        up_trigger = ROOT.TEfficiency.Wilson(N_reco, N_trigger, 0.68, True)
        down_trigger = ROOT.TEfficiency.Wilson(N_reco, N_trigger, 0.68, False)
        eps_trigger_err = 0.5*(up_trigger - down_trigger)
        eps_trigger = ufloat(eps_trigger, eps_trigger_err)

        trigger_columns = ["L0 (\\%)", "L0+HLT1 (\\%)", "Trigger (\\%)"]
        if(species == 1):
            trigger_rows = ["$3\\pi 3\\pi$", "$3\\pi 3\\pi \\pi^0$", "$3\\pi 3\\pi 2\\pi^0$", "All"]
        else:
            trigger_rows  = [""]

        trigger_table = pd.DataFrame(index=trigger_rows, columns=trigger_columns)
        if(species == 1):
            trigger_table.loc["$3\\pi 3\\pi$", "L0 (\\%)"] = f"${eps_L0_3pi3pi.nominal_value*100:.2f} \\pm {eps_L0_3pi3pi.std_dev*100:.2f}$"
            trigger_table.loc["$3\\pi 3\\pi \\pi^0$", "L0 (\\%)"] = f"${eps_L0_3pi3pipi0.nominal_value*100:.2f} \\pm {eps_L0_3pi3pipi0.std_dev*100:.2f}$"
            trigger_table.loc["$3\\pi 3\\pi 2\\pi^0$", "L0 (\\%)"] = f"${eps_L0_3pi3pi2pi0.nominal_value*100:.2f} \\pm {eps_L0_3pi3pi2pi0.std_dev*100:.2f}$"
            trigger_table.loc["All", "L0 (\\%)"] = f"${eps_L0.nominal_value*100:.2f} \\pm {eps_L0.std_dev*100:.2f}$"

            trigger_table.loc["$3\\pi 3\\pi$", "L0+HLT1 (\\%)"] = f"${eps_HLT1_3pi3pi.nominal_value*100:.2f} \\pm {eps_HLT1_3pi3pi.std_dev*100:.2f}$"
            trigger_table.loc["$3\\pi 3\\pi \\pi^0$", "L0+HLT1 (\\%)"] = f"${eps_HLT1_3pi3pipi0.nominal_value*100:.2f} \\pm {eps_HLT1_3pi3pipi0.std_dev*100:.2f}$"
            trigger_table.loc["$3\\pi 3\\pi 2\\pi^0$", "L0+HLT1 (\\%)"] = f"${eps_HLT1_3pi3pi2pi0.nominal_value*100:.2f} \\pm {eps_HLT1_3pi3pi2pi0.std_dev*100:.2f}$"
            trigger_table.loc["All", "L0+HLT1 (\\%)"] = f"${eps_HLT1.nominal_value*100:.2f} \\pm {eps_HLT1.std_dev*100:.2f}$"

            trigger_table.loc["$3\\pi 3\\pi$", "Trigger (\\%)"] = f"${eps_trigger_3pi3pi.nominal_value*100:.2f} \\pm {eps_trigger_3pi3pi.std_dev*100:.2f}$"
            trigger_table.loc["$3\\pi 3\\pi \\pi^0$", "Trigger (\\%)"] = f"${eps_trigger_3pi3pipi0.nominal_value*100:.2f} \\pm {eps_trigger_3pi3pipi0.std_dev*100:.2f}$"
            trigger_table.loc["$3\\pi 3\\pi 2\\pi^0$", "Trigger (\\%)"] = f"${eps_trigger_3pi3pi2pi0.nominal_value*100:.2f} \\pm {eps_trigger_3pi3pi2pi0.std_dev*100:.2f}$"
            trigger_table.loc["All", "Trigger (\\%)"] = f"${eps_trigger.nominal_value*100:.2f} \\pm {eps_trigger.std_dev*100:.2f}$"

        else:
            trigger_table.loc["", "L0 (\\%)"] = f"${eps_L0.nominal_value*100:.2f} \\pm {eps_L0.std_dev*100:.2f}$"
            trigger_table.loc["", "L0+HLT1 (\\%)"] = f"${eps_HLT1.nominal_value*100:.2f} \\pm {eps_HLT1.std_dev*100:.2f}$"
            trigger_table.loc["", "Trigger (\\%)"] = f"${eps_trigger.nominal_value*100:.2f} \\pm {eps_trigger.std_dev*100:.2f}$"

    with open(f'/panfs/felician/B2Ktautau/workflow/selections_efficiency_tables/Species_{species}/trigger_efficiency_table.tex', 'w') as fout_trigger:
        fout_trigger.write(trigger_table.to_latex())


    # Rectangular cuts
    print("Rectangular cuts efficiency")
    eps_rect_cuts_3pi3pi = []
    eps_rect_cuts_3pi3pipi0 = []
    eps_rect_cuts_3pi3pi2pi0 = []
    eps_rect_cuts_MC = []
    eps_rect_cuts_RS_data = []
    eps_rect_cuts_WS_data = []

    if(species == 100):
        all_rectangular_cuts = ""

        for i in range(len(rectangular_cuts_names)):
            print(rectangular_cuts_names[i])
            if(i == len(rectangular_cuts_names)-1):
                all_rectangular_cuts += rectangular_cuts[i]
            else:
                all_rectangular_cuts += rectangular_cuts[i]+" && "

        fc_reco_BuDDKp_2016 = ROOT.TFileCollection("fc_reco_BuDDKp_2016", "fc_reco_BuDDKp_2016", "/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2016/Species_1000/pre_sel_tree.txt")
        fc_reco_BuDDKp_2017 = ROOT.TFileCollection("fc_reco_BuDDKp_2017", "fc_reco_BuDDKp_2017", "/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2017/Species_1000/pre_sel_tree.txt")
        fc_reco_BuDDKp_2018 = ROOT.TFileCollection("fc_reco_BuDDKp_2018", "fc_reco_BuDDKp_2018", "/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2018/Species_1000/pre_sel_tree.txt")

        t_reco_BuDDKp_2016 = ROOT.TChain("DecayTree")
        t_reco_BuDDKp_2017 = ROOT.TChain("DecayTree")
        t_reco_BuDDKp_2018 = ROOT.TChain("DecayTree")

        t_reco_BuDDKp_2016.AddFileInfoList(fc_reco_BuDDKp_2016.GetList())
        t_reco_BuDDKp_2017.AddFileInfoList(fc_reco_BuDDKp_2017.GetList())
        t_reco_BuDDKp_2018.AddFileInfoList(fc_reco_BuDDKp_2018.GetList())

        fc_reco_BdDDKp_2016 = ROOT.TFileCollection("fc_reco_BdDDKp_2016", "fc_reco_BdDDKp_2016", "/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2016/Species_1100/pre_sel_tree.txt")
        fc_reco_BdDDKp_2017 = ROOT.TFileCollection("fc_reco_BdDDKp_2017", "fc_reco_BdDDKp_2017", "/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2017/Species_1100/pre_sel_tree.txt")
        fc_reco_BdDDKp_2018 = ROOT.TFileCollection("fc_reco_BdDDKp_2018", "fc_reco_BdDDKp_2018", "/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2018/Species_1100/pre_sel_tree.txt")

        t_reco_BdDDKp_2016 = ROOT.TChain("DecayTree")
        t_reco_BdDDKp_2017 = ROOT.TChain("DecayTree")
        t_reco_BdDDKp_2018 = ROOT.TChain("DecayTree")

        t_reco_BdDDKp_2016.AddFileInfoList(fc_reco_BdDDKp_2016.GetList())
        t_reco_BdDDKp_2017.AddFileInfoList(fc_reco_BdDDKp_2017.GetList())
        t_reco_BdDDKp_2018.AddFileInfoList(fc_reco_BdDDKp_2018.GetList())

        fc_reco_BsDDKp_2016 = ROOT.TFileCollection("fc_reco_BsDDKp_2016", "fc_reco_BsDDKp_2016", "/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2016/Species_1200/pre_sel_tree.txt")
        fc_reco_BsDDKp_2017 = ROOT.TFileCollection("fc_reco_BsDDKp_2017", "fc_reco_BsDDKp_2017", "/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2017/Species_1200/pre_sel_tree.txt")
        fc_reco_BsDDKp_2018 = ROOT.TFileCollection("fc_reco_BsDDKp_2018", "fc_reco_BsDDKp_2018", "/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2018/Species_1200/pre_sel_tree.txt")

        t_reco_BsDDKp_2016 = ROOT.TChain("DecayTree")
        t_reco_BsDDKp_2017 = ROOT.TChain("DecayTree")
        t_reco_BsDDKp_2018 = ROOT.TChain("DecayTree")

        t_reco_BsDDKp_2016.AddFileInfoList(fc_reco_BsDDKp_2016.GetList())
        t_reco_BsDDKp_2017.AddFileInfoList(fc_reco_BsDDKp_2017.GetList())
        t_reco_BsDDKp_2018.AddFileInfoList(fc_reco_BsDDKp_2018.GetList())

        fc_reco_BuDDK0_2016 = ROOT.TFileCollection("fc_reco_BuDDK0_2016", "fc_reco_BuDDK0_2016", "/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2016/Species_1300/pre_sel_tree.txt")
        fc_reco_BuDDK0_2017 = ROOT.TFileCollection("fc_reco_BuDDK0_2017", "fc_reco_BuDDK0_2017", "/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2017/Species_1300/pre_sel_tree.txt")
        fc_reco_BuDDK0_2018 = ROOT.TFileCollection("fc_reco_BuDDK0_2018", "fc_reco_BuDDK0_2018", "/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2018/Species_1300/pre_sel_tree.txt")

        t_reco_BuDDK0_2016 = ROOT.TChain("DecayTree")
        t_reco_BuDDK0_2017 = ROOT.TChain("DecayTree")
        t_reco_BuDDK0_2018 = ROOT.TChain("DecayTree")

        t_reco_BuDDK0_2016.AddFileInfoList(fc_reco_BuDDK0_2016.GetList())
        t_reco_BuDDK0_2017.AddFileInfoList(fc_reco_BuDDK0_2017.GetList())
        t_reco_BuDDK0_2018.AddFileInfoList(fc_reco_BuDDK0_2018.GetList())

        fc_reco_BuDD_2016 = ROOT.TFileCollection("fc_reco_BuDD_2016", "fc_reco_BuDD_2016", "/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2016/Species_1500/pre_sel_tree.txt")
        fc_reco_BuDD_2017 = ROOT.TFileCollection("fc_reco_BuDD_2017", "fc_reco_BuDD_2017", "/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2017/Species_1500/pre_sel_tree.txt")
        fc_reco_BuDD_2018 = ROOT.TFileCollection("fc_reco_BuDD_2018", "fc_reco_BuDD_2018", "/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2018/Species_1500/pre_sel_tree.txt")

        t_reco_BuDD_2016 = ROOT.TChain("DecayTree")
        t_reco_BuDD_2017 = ROOT.TChain("DecayTree")
        t_reco_BuDD_2018 = ROOT.TChain("DecayTree")

        t_reco_BuDD_2016.AddFileInfoList(fc_reco_BuDD_2016.GetList())
        t_reco_BuDD_2017.AddFileInfoList(fc_reco_BuDD_2017.GetList())
        t_reco_BuDD_2018.AddFileInfoList(fc_reco_BuDD_2018.GetList())

        fc_gsl_BuDDKp_2016 = ROOT.TFileCollection("fc_gsl_BuDDKp_2016", "fc_gsl_BuDDKp_2016", "/panfs/felician/B2Ktautau/workflow/standalone_fitter/2016/Species_1000/fit_results.txt")
        fc_gsl_BuDDKp_2017 = ROOT.TFileCollection("fc_gsl_BuDDKp_2017", "fc_gsl_BuDDKp_2017", "/panfs/felician/B2Ktautau/workflow/standalone_fitter/2017/Species_1000/fit_results.txt")
        fc_gsl_BuDDKp_2018 = ROOT.TFileCollection("fc_gsl_BuDDKp_2018", "fc_gsl_BuDDKp_2018", "/panfs/felician/B2Ktautau/workflow/standalone_fitter/2018/Species_1000/fit_results.txt")

        t_gsl_BuDDKp_2016 = ROOT.TChain("DecayTree")
        t_gsl_BuDDKp_2017 = ROOT.TChain("DecayTree")
        t_gsl_BuDDKp_2018 = ROOT.TChain("DecayTree")

        t_gsl_BuDDKp_2016.AddFileInfoList(fc_gsl_BuDDKp_2016.GetList())
        t_gsl_BuDDKp_2017.AddFileInfoList(fc_gsl_BuDDKp_2017.GetList())
        t_gsl_BuDDKp_2018.AddFileInfoList(fc_gsl_BuDDKp_2018.GetList())

        fc_gsl_BdDDKp_2016 = ROOT.TFileCollection("fc_gsl_BdDDKp_2016", "fc_gsl_BdDDKp_2016", "/panfs/felician/B2Ktautau/workflow/standalone_fitter/2016/Species_1100/fit_results.txt")
        fc_gsl_BdDDKp_2017 = ROOT.TFileCollection("fc_gsl_BdDDKp_2017", "fc_gsl_BdDDKp_2017", "/panfs/felician/B2Ktautau/workflow/standalone_fitter/2017/Species_1100/fit_results.txt")
        fc_gsl_BdDDKp_2018 = ROOT.TFileCollection("fc_gsl_BdDDKp_2018", "fc_gsl_BdDDKp_2018", "/panfs/felician/B2Ktautau/workflow/standalone_fitter/2018/Species_1100/fit_results.txt")

        t_gsl_BdDDKp_2016 = ROOT.TChain("DecayTree")
        t_gsl_BdDDKp_2017 = ROOT.TChain("DecayTree")
        t_gsl_BdDDKp_2018 = ROOT.TChain("DecayTree")

        t_gsl_BdDDKp_2016.AddFileInfoList(fc_gsl_BdDDKp_2016.GetList())
        t_gsl_BdDDKp_2017.AddFileInfoList(fc_gsl_BdDDKp_2017.GetList())
        t_gsl_BdDDKp_2018.AddFileInfoList(fc_gsl_BdDDKp_2018.GetList())

        fc_gsl_BsDDKp_2016 = ROOT.TFileCollection("fc_gsl_BsDDKp_2016", "fc_gsl_BsDDKp_2016", "/panfs/felician/B2Ktautau/workflow/standalone_fitter/2016/Species_1200/fit_results.txt")
        fc_gsl_BsDDKp_2017 = ROOT.TFileCollection("fc_gsl_BsDDKp_2017", "fc_gsl_BsDDKp_2017", "/panfs/felician/B2Ktautau/workflow/standalone_fitter/2017/Species_1200/fit_results.txt")
        fc_gsl_BsDDKp_2018 = ROOT.TFileCollection("fc_gsl_BsDDKp_2018", "fc_gsl_BsDDKp_2018", "/panfs/felician/B2Ktautau/workflow/standalone_fitter/2018/Species_1200/fit_results.txt")

        t_gsl_BsDDKp_2016 = ROOT.TChain("DecayTree")
        t_gsl_BsDDKp_2017 = ROOT.TChain("DecayTree")
        t_gsl_BsDDKp_2018 = ROOT.TChain("DecayTree")

        t_gsl_BsDDKp_2016.AddFileInfoList(fc_gsl_BsDDKp_2016.GetList())
        t_gsl_BsDDKp_2017.AddFileInfoList(fc_gsl_BsDDKp_2017.GetList())
        t_gsl_BsDDKp_2018.AddFileInfoList(fc_gsl_BsDDKp_2018.GetList())

        fc_gsl_BuDDK0_2016 = ROOT.TFileCollection("fc_gsl_BuDDK0_2016", "fc_gsl_BuDDK0_2016", "/panfs/felician/B2Ktautau/workflow/standalone_fitter/2016/Species_1300/fit_results.txt")
        fc_gsl_BuDDK0_2017 = ROOT.TFileCollection("fc_gsl_BuDDK0_2017", "fc_gsl_BuDDK0_2017", "/panfs/felician/B2Ktautau/workflow/standalone_fitter/2017/Species_1300/fit_results.txt")
        fc_gsl_BuDDK0_2018 = ROOT.TFileCollection("fc_gsl_BuDDK0_2018", "fc_gsl_BuDDK0_2018", "/panfs/felician/B2Ktautau/workflow/standalone_fitter/2018/Species_1300/fit_results.txt")

        t_gsl_BuDDK0_2016 = ROOT.TChain("DecayTree")
        t_gsl_BuDDK0_2017 = ROOT.TChain("DecayTree")
        t_gsl_BuDDK0_2018 = ROOT.TChain("DecayTree")

        t_gsl_BuDDK0_2016.AddFileInfoList(fc_gsl_BuDDK0_2016.GetList())
        t_gsl_BuDDK0_2017.AddFileInfoList(fc_gsl_BuDDK0_2017.GetList())
        t_gsl_BuDDK0_2018.AddFileInfoList(fc_gsl_BuDDK0_2018.GetList())

        fc_gsl_BuDD_2016 = ROOT.TFileCollection("fc_gsl_BuDD_2016", "fc_gsl_BuDD_2016", "/panfs/felician/B2Ktautau/workflow/standalone_fitter/2016/Species_1500/fit_results.txt")
        fc_gsl_BuDD_2017 = ROOT.TFileCollection("fc_gsl_BuDD_2017", "fc_gsl_BuDD_2017", "/panfs/felician/B2Ktautau/workflow/standalone_fitter/2017/Species_1500/fit_results.txt")
        fc_gsl_BuDD_2018 = ROOT.TFileCollection("fc_gsl_BuDD_2018", "fc_gsl_BuDD_2018", "/panfs/felician/B2Ktautau/workflow/standalone_fitter/2018/Species_1500/fit_results.txt")

        t_gsl_BuDD_2016 = ROOT.TChain("DecayTree")
        t_gsl_BuDD_2017 = ROOT.TChain("DecayTree")
        t_gsl_BuDD_2018 = ROOT.TChain("DecayTree")

        t_gsl_BuDD_2016.AddFileInfoList(fc_gsl_BuDD_2016.GetList())
        t_gsl_BuDD_2017.AddFileInfoList(fc_gsl_BuDD_2017.GetList())
        t_gsl_BuDD_2018.AddFileInfoList(fc_gsl_BuDD_2018.GetList())

        fc_mass_BuDDKp_2016 =  ROOT.TFileCollection("fc_mass_BuDDKp_2016", "fc_mass_BuDDKp_2016", "/panfs/felician/B2Ktautau/workflow/create_invariant_mass_tree/2016/Species_1000/invariant_mass_tree.txt")
        fc_mass_BuDDKp_2017 =  ROOT.TFileCollection("fc_mass_BuDDKp_2017", "fc_mass_BuDDKp_2017", "/panfs/felician/B2Ktautau/workflow/create_invariant_mass_tree/2017/Species_1000/invariant_mass_tree.txt")
        fc_mass_BuDDKp_2018 =  ROOT.TFileCollection("fc_mass_BuDDKp_2018", "fc_mass_BuDDKp_2018", "/panfs/felician/B2Ktautau/workflow/create_invariant_mass_tree/2018/Species_1000/invariant_mass_tree.txt")

        t_mass_BuDDKp_2016 = ROOT.TChain("DecayTree")
        t_mass_BuDDKp_2017 = ROOT.TChain("DecayTree")
        t_mass_BuDDKp_2018 = ROOT.TChain("DecayTree")

        t_mass_BuDDKp_2016.AddFileInfoList(fc_mass_BuDDKp_2016.GetList())
        t_mass_BuDDKp_2017.AddFileInfoList(fc_mass_BuDDKp_2017.GetList())
        t_mass_BuDDKp_2018.AddFileInfoList(fc_mass_BuDDKp_2018.GetList())

        fc_mass_BdDDKp_2016 =  ROOT.TFileCollection("fc_mass_BdDDKp_2016", "fc_mass_BdDDKp_2016", "/panfs/felician/B2Ktautau/workflow/create_invariant_mass_tree/2016/Species_1100/invariant_mass_tree.txt")
        fc_mass_BdDDKp_2017 =  ROOT.TFileCollection("fc_mass_BdDDKp_2017", "fc_mass_BdDDKp_2017", "/panfs/felician/B2Ktautau/workflow/create_invariant_mass_tree/2017/Species_1100/invariant_mass_tree.txt")
        fc_mass_BdDDKp_2018 =  ROOT.TFileCollection("fc_mass_BdDDKp_2018", "fc_mass_BdDDKp_2018", "/panfs/felician/B2Ktautau/workflow/create_invariant_mass_tree/2018/Species_1100/invariant_mass_tree.txt")

        t_mass_BdDDKp_2016 = ROOT.TChain("DecayTree")
        t_mass_BdDDKp_2017 = ROOT.TChain("DecayTree")
        t_mass_BdDDKp_2018 = ROOT.TChain("DecayTree")

        t_mass_BdDDKp_2016.AddFileInfoList(fc_mass_BdDDKp_2016.GetList())
        t_mass_BdDDKp_2017.AddFileInfoList(fc_mass_BdDDKp_2017.GetList())
        t_mass_BdDDKp_2018.AddFileInfoList(fc_mass_BdDDKp_2018.GetList())

        fc_mass_BsDDKp_2016 =  ROOT.TFileCollection("fc_mass_BsDDKp_2016", "fc_mass_BsDDKp_2016", "/panfs/felician/B2Ktautau/workflow/create_invariant_mass_tree/2016/Species_1200/invariant_mass_tree.txt")
        fc_mass_BsDDKp_2017 =  ROOT.TFileCollection("fc_mass_BsDDKp_2017", "fc_mass_BsDDKp_2017", "/panfs/felician/B2Ktautau/workflow/create_invariant_mass_tree/2017/Species_1200/invariant_mass_tree.txt")
        fc_mass_BsDDKp_2018 =  ROOT.TFileCollection("fc_mass_BsDDKp_2018", "fc_mass_BsDDKp_2018", "/panfs/felician/B2Ktautau/workflow/create_invariant_mass_tree/2018/Species_1200/invariant_mass_tree.txt")

        t_mass_BsDDKp_2016 = ROOT.TChain("DecayTree")
        t_mass_BsDDKp_2017 = ROOT.TChain("DecayTree")
        t_mass_BsDDKp_2018 = ROOT.TChain("DecayTree")

        t_mass_BsDDKp_2016.AddFileInfoList(fc_mass_BsDDKp_2016.GetList())
        t_mass_BsDDKp_2017.AddFileInfoList(fc_mass_BsDDKp_2017.GetList())
        t_mass_BsDDKp_2018.AddFileInfoList(fc_mass_BsDDKp_2018.GetList())

        fc_mass_BuDDK0_2016 =  ROOT.TFileCollection("fc_mass_BuDDK0_2016", "fc_mass_BuDDK0_2016", "/panfs/felician/B2Ktautau/workflow/create_invariant_mass_tree/2016/Species_1300/invariant_mass_tree.txt")
        fc_mass_BuDDK0_2017 =  ROOT.TFileCollection("fc_mass_BuDDK0_2017", "fc_mass_BuDDK0_2017", "/panfs/felician/B2Ktautau/workflow/create_invariant_mass_tree/2017/Species_1300/invariant_mass_tree.txt")
        fc_mass_BuDDK0_2018 =  ROOT.TFileCollection("fc_mass_BuDDK0_2018", "fc_mass_BuDDK0_2018", "/panfs/felician/B2Ktautau/workflow/create_invariant_mass_tree/2018/Species_1300/invariant_mass_tree.txt")

        t_mass_BuDDK0_2016 = ROOT.TChain("DecayTree")
        t_mass_BuDDK0_2017 = ROOT.TChain("DecayTree")
        t_mass_BuDDK0_2018 = ROOT.TChain("DecayTree")

        t_mass_BuDDK0_2016.AddFileInfoList(fc_mass_BuDDK0_2016.GetList())
        t_mass_BuDDK0_2017.AddFileInfoList(fc_mass_BuDDK0_2017.GetList())
        t_mass_BuDDK0_2018.AddFileInfoList(fc_mass_BuDDK0_2018.GetList())

        fc_mass_BuDD_2016 =  ROOT.TFileCollection("fc_mass_BuDD_2016", "fc_mass_BuDD_2016", "/panfs/felician/B2Ktautau/workflow/create_invariant_mass_tree/2016/Species_1500/invariant_mass_tree.txt")
        fc_mass_BuDD_2017 =  ROOT.TFileCollection("fc_mass_BuDD_2017", "fc_mass_BuDD_2017", "/panfs/felician/B2Ktautau/workflow/create_invariant_mass_tree/2017/Species_1500/invariant_mass_tree.txt")
        fc_mass_BuDD_2018 =  ROOT.TFileCollection("fc_mass_BuDD_2018", "fc_mass_BuDD_2018", "/panfs/felician/B2Ktautau/workflow/create_invariant_mass_tree/2018/Species_1500/invariant_mass_tree.txt")

        t_mass_BuDD_2016 = ROOT.TChain("DecayTree")
        t_mass_BuDD_2017 = ROOT.TChain("DecayTree")
        t_mass_BuDD_2018 = ROOT.TChain("DecayTree")

        t_mass_BuDD_2016.AddFileInfoList(fc_mass_BuDD_2016.GetList())
        t_mass_BuDD_2017.AddFileInfoList(fc_mass_BuDD_2017.GetList())
        t_mass_BuDD_2018.AddFileInfoList(fc_mass_BuDD_2018.GetList())

        t_reco_BuDDKp_2016.AddFriend(t_gsl_BuDDKp_2016)
        t_reco_BuDDKp_2016.AddFriend(t_mass_BuDDKp_2016)

        t_reco_BuDDKp_2017.AddFriend(t_gsl_BuDDKp_2017)
        t_reco_BuDDKp_2017.AddFriend(t_mass_BuDDKp_2017)

        t_reco_BuDDKp_2018.AddFriend(t_gsl_BuDDKp_2018)
        t_reco_BuDDKp_2018.AddFriend(t_mass_BuDDKp_2018)

        t_reco_BdDDKp_2016.AddFriend(t_gsl_BdDDKp_2016)
        t_reco_BdDDKp_2016.AddFriend(t_mass_BdDDKp_2016)

        t_reco_BdDDKp_2017.AddFriend(t_gsl_BdDDKp_2017)
        t_reco_BdDDKp_2017.AddFriend(t_mass_BdDDKp_2017)

        t_reco_BdDDKp_2018.AddFriend(t_gsl_BdDDKp_2018)
        t_reco_BdDDKp_2018.AddFriend(t_mass_BdDDKp_2018)

        t_reco_BsDDKp_2016.AddFriend(t_gsl_BsDDKp_2016)
        t_reco_BsDDKp_2016.AddFriend(t_mass_BsDDKp_2016)

        t_reco_BsDDKp_2017.AddFriend(t_gsl_BsDDKp_2017)
        t_reco_BsDDKp_2017.AddFriend(t_mass_BsDDKp_2017)

        t_reco_BsDDKp_2018.AddFriend(t_gsl_BsDDKp_2018)
        t_reco_BsDDKp_2018.AddFriend(t_mass_BsDDKp_2018)

        t_reco_BuDDK0_2016.AddFriend(t_gsl_BuDDK0_2016)
        t_reco_BuDDK0_2016.AddFriend(t_mass_BuDDK0_2016)

        t_reco_BuDDK0_2017.AddFriend(t_gsl_BuDDK0_2017)
        t_reco_BuDDK0_2017.AddFriend(t_mass_BuDDK0_2017)

        t_reco_BuDDK0_2018.AddFriend(t_gsl_BuDDK0_2018)
        t_reco_BuDDK0_2018.AddFriend(t_mass_BuDDK0_2018)

        t_reco_BuDD_2016.AddFriend(t_gsl_BuDD_2016)
        t_reco_BuDD_2016.AddFriend(t_mass_BuDD_2016)

        t_reco_BuDD_2017.AddFriend(t_gsl_BuDD_2017)
        t_reco_BuDD_2017.AddFriend(t_mass_BuDD_2017)

        t_reco_BuDD_2018.AddFriend(t_gsl_BuDD_2018)
        t_reco_BuDD_2018.AddFriend(t_mass_BuDD_2018)

        N_rect_cuts_BuDDKp = t_reco_BuDDKp_2016.GetEntries(trigger+" && "+all_rectangular_cuts) + t_reco_BuDDKp_2017.GetEntries(trigger+" && "+all_rectangular_cuts) + t_reco_BuDDKp_2018.GetEntries(trigger+" && "+all_rectangular_cuts)
        N_rect_cuts_BdDDKp = t_reco_BdDDKp_2016.GetEntries(trigger+" && "+all_rectangular_cuts) + t_reco_BdDDKp_2017.GetEntries(trigger+" && "+all_rectangular_cuts) + t_reco_BdDDKp_2018.GetEntries(trigger+" && "+all_rectangular_cuts)
        N_rect_cuts_BsDDKp = t_reco_BsDDKp_2016.GetEntries(trigger+" && "+all_rectangular_cuts) + t_reco_BsDDKp_2017.GetEntries(trigger+" && "+all_rectangular_cuts) + t_reco_BsDDKp_2018.GetEntries(trigger+" && "+all_rectangular_cuts)
        N_rect_cuts_BuDDK0 = t_reco_BuDDK0_2016.GetEntries(trigger+" && "+all_rectangular_cuts) + t_reco_BuDDK0_2017.GetEntries(trigger+" && "+all_rectangular_cuts) + t_reco_BuDDK0_2018.GetEntries(trigger+" && "+all_rectangular_cuts)
        N_rect_cuts_BuDD = t_reco_BuDD_2016.GetEntries(trigger+" && "+all_rectangular_cuts) + t_reco_BuDD_2017.GetEntries(trigger+" && "+all_rectangular_cuts) + t_reco_BuDD_2018.GetEntries(trigger+" && "+all_rectangular_cuts)

        up_rect_cuts_BuDDKp = ROOT.TEfficiency.Wilson(N_trigger_BuDDKp, N_rect_cuts_BuDDKp, 0.68, True)
        down_rect_cuts_BuDDKp = ROOT.TEfficiency.Wilson(N_trigger_BuDDKp, N_rect_cuts_BuDDKp, 0.68, False)
        eps_rect_cuts_BuDDKp_err = 0.5*(up_rect_cuts_BuDDKp - down_rect_cuts_BuDDKp)
        eps_rect_cuts_BuDDKp = ufloat(N_rect_cuts_BuDDKp/N_trigger_BuDDKp, eps_rect_cuts_BuDDKp_err)

        up_rect_cuts_BdDDKp = ROOT.TEfficiency.Wilson(N_trigger_BdDDKp, N_rect_cuts_BdDDKp, 0.68, True)
        down_rect_cuts_BdDDKp = ROOT.TEfficiency.Wilson(N_trigger_BdDDKp, N_rect_cuts_BdDDKp, 0.68, False)
        eps_rect_cuts_BdDDKp_err = 0.5*(up_rect_cuts_BdDDKp - down_rect_cuts_BdDDKp)
        eps_rect_cuts_BdDDKp = ufloat(N_rect_cuts_BdDDKp/N_trigger_BdDDKp, eps_rect_cuts_BdDDKp_err)

        up_rect_cuts_BsDDKp = ROOT.TEfficiency.Wilson(N_trigger_BsDDKp, N_rect_cuts_BsDDKp, 0.68, True)
        down_rect_cuts_BsDDKp = ROOT.TEfficiency.Wilson(N_trigger_BsDDKp, N_rect_cuts_BsDDKp, 0.68, False)
        eps_rect_cuts_BsDDKp_err = 0.5*(up_rect_cuts_BsDDKp - down_rect_cuts_BsDDKp)
        eps_rect_cuts_BsDDKp = ufloat(N_rect_cuts_BsDDKp/N_trigger_BsDDKp, eps_rect_cuts_BsDDKp_err)

        up_rect_cuts_BuDDK0 = ROOT.TEfficiency.Wilson(N_trigger_BuDDK0, N_rect_cuts_BuDDK0, 0.68, True)
        down_rect_cuts_BuDDK0 = ROOT.TEfficiency.Wilson(N_trigger_BuDDK0, N_rect_cuts_BuDDK0, 0.68, False)
        eps_rect_cuts_BuDDK0_err = 0.5*(up_rect_cuts_BuDDK0 - down_rect_cuts_BuDDK0)
        eps_rect_cuts_BuDDK0 = ufloat(N_rect_cuts_BuDDK0/N_trigger_BuDDK0, eps_rect_cuts_BuDDK0_err)

        up_rect_cuts_BuDD = ROOT.TEfficiency.Wilson(N_trigger_BuDD, N_rect_cuts_BuDD, 0.68, True)
        down_rect_cuts_BuDD = ROOT.TEfficiency.Wilson(N_trigger_BuDD, N_rect_cuts_BuDD, 0.68, False)
        eps_rect_cuts_BuDD_err = 0.5*(up_rect_cuts_BuDD - down_rect_cuts_BuDD)
        eps_rect_cuts_BuDD = ufloat(N_rect_cuts_BuDD/N_trigger_BuDD, eps_rect_cuts_BuDD_err)

        eps_rect_cuts = [eps_rect_cuts_BuDDKp, eps_rect_cuts_BdDDKp, eps_rect_cuts_BsDDKp, eps_rect_cuts_BuDDK0, eps_rect_cuts_BuDD]

        rect_cuts_columns = ["Rectangular cuts (\\%)"]
        rect_cuts_rows = ["$B^+ \\to D D K^+$", "$B^0 \\to D D K^+$", "$B^0_s \\to D D K^+$", "$B^+ \\to D D K^0$", "$B^+ \\to D D$"]

        rect_cut_table_1 = pd.DataFrame(index=rect_cuts_rows, columns=rect_cuts_columns)
        rect_cut_table_1.loc["$B^+ \\to D D K^+$", "Rectangular cuts (\\%)"] = f"${eps_rect_cuts[0].nominal_value*100:.2f} \\pm {eps_rect_cuts[0].std_dev*100:.2f}$"
        rect_cut_table_1.loc["$B^0 \\to D D K^+$", "Rectangular cuts (\\%)"] = f"${eps_rect_cuts[1].nominal_value*100:.2f} \\pm {eps_rect_cuts[1].std_dev*100:.2f}$"
        rect_cut_table_1.loc["$B^0_s \\to D D K^+$", "Rectangular cuts (\\%)"] = f"${eps_rect_cuts[2].nominal_value*100:.2f} \\pm {eps_rect_cuts[2].std_dev*100:.2f}$"
        rect_cut_table_1.loc["$B^+ \\to D D K^0$", "Rectangular cuts (\\%)"] = f"${eps_rect_cuts[3].nominal_value*100:.2f} \\pm {eps_rect_cuts[3].std_dev*100:.2f}$"
        rect_cut_table_1.loc["$B^+ \\to D D$", "Rectangular cuts (\\%)"] = f"${eps_rect_cuts[4].nominal_value*100:.2f} \\pm {eps_rect_cuts[4].std_dev*100:.2f}$"

    
    else:
        if(species == 1):
            fc_data_2016_noRectCuts = ROOT.TFileCollection("fc_data_2016_noRectCuts", "fc_data_2016_noRectCuts", "Files_on_grid/data_2016_noSel_DTF.txt")
            fc_data_2017_noRectCuts = ROOT.TFileCollection("fc_data_2017_noRectCuts", "fc_data_2017_noRectCuts", "Files_on_grid/data_2017_noSel_DTF.txt")
            fc_data_2018_noRectCuts = ROOT.TFileCollection("fc_data_2018_noRectCuts", "fc_data_2018_noRectCuts", "Files_on_grid/data_2018_noSel_DTF.txt")

            t_rs_data_2016_noRectCuts = ROOT.TChain("ntuple/DecayTree")
            t_rs_data_2017_noRectCuts = ROOT.TChain("ntuple/DecayTree")
            t_rs_data_2018_noRectCuts = ROOT.TChain("ntuple/DecayTree")

            t_rs_data_2016_noRectCuts.AddFileInfoList(fc_data_2016_noRectCuts.GetList())
            t_rs_data_2017_noRectCuts.AddFileInfoList(fc_data_2017_noRectCuts.GetList())
            t_rs_data_2018_noRectCuts.AddFileInfoList(fc_data_2018_noRectCuts.GetList())

            t_ws_data_2016_noRectCuts = ROOT.TChain("ntuple_SS/DecayTree")
            t_ws_data_2017_noRectCuts = ROOT.TChain("ntuple_SS/DecayTree")
            t_ws_data_2018_noRectCuts = ROOT.TChain("ntuple_SS/DecayTree")
            
            t_ws_data_2016_noRectCuts.AddFileInfoList(fc_data_2016_noRectCuts.GetList())
            t_ws_data_2017_noRectCuts.AddFileInfoList(fc_data_2017_noRectCuts.GetList())
            t_ws_data_2018_noRectCuts.AddFileInfoList(fc_data_2018_noRectCuts.GetList())
        else:
            fc_data_2016_noRectCuts = ROOT.TFileCollection("fc_data_2016_noRectCuts", "fc_data_2016_noRectCuts", "Files_on_grid/data_D0Dps_2016.txt", 50)
            fc_data_2017_noRectCuts = ROOT.TFileCollection("fc_data_2017_noRectCuts", "fc_data_2017_noRectCuts", "Files_on_grid/data_D0Dps_2017.txt", 50)
            fc_data_2018_noRectCuts = ROOT.TFileCollection("fc_data_2018_noRectCuts", "fc_data_2018_noRectCuts", "Files_on_grid/data_D0Dps_2018.txt", 50)

            t_rs_data_2016_noRectCuts = ROOT.TChain("ntuple/DecayTree")
            t_rs_data_2017_noRectCuts = ROOT.TChain("ntuple/DecayTree")
            t_rs_data_2018_noRectCuts = ROOT.TChain("ntuple/DecayTree")

            t_rs_data_2016_noRectCuts.AddFileInfoList(fc_data_2016_noRectCuts.GetList())
            t_rs_data_2017_noRectCuts.AddFileInfoList(fc_data_2017_noRectCuts.GetList())
            t_rs_data_2018_noRectCuts.AddFileInfoList(fc_data_2018_noRectCuts.GetList())

        all_rectangular_cuts_MC = ""
        all_rectangular_cuts_data = ""

        for i in range(len(rectangular_cuts_names)):
            print(rectangular_cuts_names[i])
            if(species == 1):
                N_rect_cuts_3pi3pi = t_reco_pid_2016.GetEntries(truthMatch+" && "+trigger+" && "+rectangular_cuts[i]+" && (component == 0)") + t_reco_pid_2017.GetEntries(truthMatch+" && "+trigger+" && "+rectangular_cuts[i]+" && (component == 0)") + t_reco_pid_2018.GetEntries(truthMatch+" && "+trigger+" && "+rectangular_cuts[i]+" && (component == 0)") 
                N_rect_cuts_3pi3pipi0 = t_reco_pid_2016.GetEntries(truthMatch+" && "+trigger+" && "+rectangular_cuts[i]+" && (component == 1)") + t_reco_pid_2017.GetEntries(truthMatch+" && "+trigger+" && "+rectangular_cuts[i]+" && (component == 1)") + t_reco_pid_2018.GetEntries(truthMatch+" && "+trigger+" && "+rectangular_cuts[i]+" && (component == 1)") 
                N_rect_cuts_3pi3pi2pi0 = t_reco_pid_2016.GetEntries(truthMatch+" && "+trigger+" && "+rectangular_cuts[i]+" && (component == 2)") + t_reco_pid_2017.GetEntries(truthMatch+" && "+trigger+" && "+rectangular_cuts[i]+" && (component == 2)") + t_reco_pid_2018.GetEntries(truthMatch+" && "+trigger+" && "+rectangular_cuts[i]+" && (component == 2)") 
                    
                up_rect_cuts_3pi3pi = ROOT.TEfficiency.Wilson(N_trigger_3pi3pi, N_rect_cuts_3pi3pi, 0.68, True)
                down_rect_cuts_3pi3pi = ROOT.TEfficiency.Wilson(N_trigger_3pi3pi, N_rect_cuts_3pi3pi, 0.68, False)
                eps_rect_cuts_3pi3pi_err = 0.5*(up_rect_cuts_3pi3pi - down_rect_cuts_3pi3pi)
                eps_rect_cuts_3pi3pi.append( ufloat(N_rect_cuts_3pi3pi/N_trigger_3pi3pi, eps_rect_cuts_3pi3pi_err) )

                up_rect_cuts_3pi3pipi0 = ROOT.TEfficiency.Wilson(N_trigger_3pi3pipi0, N_rect_cuts_3pi3pipi0, 0.68, True)
                down_rect_cuts_3pi3pipi0 = ROOT.TEfficiency.Wilson(N_trigger_3pi3pipi0, N_rect_cuts_3pi3pipi0, 0.68, False)
                eps_rect_cuts_3pi3pipi0_err = 0.5*(up_rect_cuts_3pi3pipi0 - down_rect_cuts_3pi3pipi0)
                eps_rect_cuts_3pi3pipi0.append( ufloat(N_rect_cuts_3pi3pipi0/N_trigger_3pi3pipi0, eps_rect_cuts_3pi3pipi0_err) )

                up_rect_cuts_3pi3pi2pi0 = ROOT.TEfficiency.Wilson(N_trigger_3pi3pi2pi0, N_rect_cuts_3pi3pi2pi0, 0.68, True)
                down_rect_cuts_3pi3pi2pi0 = ROOT.TEfficiency.Wilson(N_trigger_3pi3pi2pi0, N_rect_cuts_3pi3pi2pi0, 0.68, False)
                eps_rect_cuts_3pi3pi2pi0_err = 0.5*(up_rect_cuts_3pi3pi2pi0 - down_rect_cuts_3pi3pi2pi0)
                eps_rect_cuts_3pi3pi2pi0.append( ufloat(N_rect_cuts_3pi3pi2pi0/N_trigger_3pi3pi2pi0, eps_rect_cuts_3pi3pi2pi0_err) )

                cut_name = rectangular_cuts[i]
                if(cut_name == "(Kp_ProbNNk_pidgen_default > 0.2)"):
                    cut_name = "(Kp_ProbNNk > 0.2)"
                if(cut_name == "(taup_pi1_ProbNNpi_pidgen_default > 0.55) && (taup_pi2_ProbNNpi_pidgen_default > 0.55) && (taup_pi3_ProbNNpi_pidgen_default > 0.55) && (taum_pi1_ProbNNpi_pidgen_default > 0.55) && (taum_pi2_ProbNNpi_pidgen_default > 0.55) && (taum_pi3_ProbNNpi_pidgen_default > 0.55)"):
                    cut_name = "(taup_pi1_ProbNNpi > 0.55) && (taup_pi2_ProbNNpi > 0.55) && (taup_pi3_ProbNNpi > 0.55) && (taum_pi1_ProbNNpi > 0.55) && (taum_pi2_ProbNNpi > 0.55) && (taum_pi3_ProbNNpi > 0.55)"

                N_rect_cuts_RS_data = t_rs_data_2016_noRectCuts.GetEntries(trigger+" && "+cut_name) + t_rs_data_2017_noRectCuts.GetEntries(trigger+" && "+cut_name) + t_rs_data_2018_noRectCuts.GetEntries(trigger+" && "+cut_name)
                N_rect_cuts_WS_data = t_ws_data_2016_noRectCuts.GetEntries(trigger+" && "+cut_name) + t_ws_data_2017_noRectCuts.GetEntries(trigger+" && "+cut_name) + t_ws_data_2018_noRectCuts.GetEntries(trigger+" && "+cut_name)

                N_trigger_RS_data = t_rs_data_2016_noRectCuts.GetEntries(trigger) + t_rs_data_2017_noRectCuts.GetEntries(trigger) + t_rs_data_2018_noRectCuts.GetEntries(trigger)
                N_trigger_WS_data = t_ws_data_2016_noRectCuts.GetEntries(trigger) + t_ws_data_2017_noRectCuts.GetEntries(trigger) + t_ws_data_2018_noRectCuts.GetEntries(trigger)
                    
                if(i == len(rectangular_cuts_names)-1):
                    all_rectangular_cuts_MC += rectangular_cuts[i]
                    all_rectangular_cuts_data += cut_name
                else:
                    all_rectangular_cuts_MC += rectangular_cuts[i]+" && "
                    all_rectangular_cuts_data += cut_name+" && "

            else:
                cut_name = rectangular_cuts[i]
                if(cut_name == "(D0bar_K_ProbNNk_pidgen_default > 0.55)"):
                    cut_name = "(D0bar_K_ProbNNk > 0.55)"
                if(cut_name == "(Dsp_K1_ProbNNk_pidgen_default > 0.55)"):
                    cut_name = "(Dsp_K1_ProbNNk > 0.55)"
                if(cut_name == "(Dsp_K2_ProbNNk_pidgen_default > 0.55)"):
                    cut_name = "(Dsp_K2_ProbNNk > 0.55)"
                if(cut_name == "(D0bar_pi_ProbNNpi_pidgen_default > 0.55)"):
                    cut_name = "(D0bar_pi_ProbNNpi > 0.55)"
                if(cut_name == "(Dsp_pi_ProbNNpi_pidgen_default > 0.55)"):
                    cut_name = "(Dsp_pi_ProbNNpi > 0.55)"

                N_rect_cuts_RS_data = t_rs_data_2016_noRectCuts.GetEntries(trigger+" && "+cut_name) + t_rs_data_2017_noRectCuts.GetEntries(trigger+" && "+cut_name) + t_rs_data_2018_noRectCuts.GetEntries(trigger+" && "+cut_name)
                N_rect_cuts_WS_data = t_rs_data_2016_noRectCuts.GetEntries(trigger+" && "+cut_name+" && (Bp_dtf_M[0] > 5320)") + t_rs_data_2017_noRectCuts.GetEntries(trigger+" && "+cut_name+" && (Bp_dtf_M[0] > 5320)") + t_rs_data_2018_noRectCuts.GetEntries(trigger+" && "+cut_name+" && (Bp_dtf_M[0] > 5320)")

                N_trigger_RS_data = t_rs_data_2016_noRectCuts.GetEntries(trigger) + t_rs_data_2017_noRectCuts.GetEntries(trigger) + t_rs_data_2018_noRectCuts.GetEntries(trigger)
                N_trigger_WS_data = t_rs_data_2016_noRectCuts.GetEntries(trigger+" && (Bp_dtf_M[0] > 5320)") + t_rs_data_2017_noRectCuts.GetEntries(trigger+" && (Bp_dtf_M[0] > 5320)") + t_rs_data_2018_noRectCuts.GetEntries(trigger+" && (Bp_dtf_M[0] > 5320)")

                if(i == len(rectangular_cuts_names)-1):
                    all_rectangular_cuts_MC += rectangular_cuts[i]
                    all_rectangular_cuts_data += cut_name
                else:
                    all_rectangular_cuts_MC += rectangular_cuts[i]+" && "
                    all_rectangular_cuts_data += cut_name+" && "

            N_rect_cuts_MC = t_reco_pid_2016.GetEntries(truthMatch+" && "+trigger+" && "+rectangular_cuts[i]) + t_reco_pid_2017.GetEntries(truthMatch+" && "+trigger+" && "+rectangular_cuts[i]) + t_reco_pid_2018.GetEntries(truthMatch+" && "+trigger+" && "+rectangular_cuts[i]) 
            N_rect_cuts_MC_all = t_reco_pid_2016.GetEntries(truthMatch+" && "+trigger+" && "+all_rectangular_cuts_MC) + t_reco_pid_2017.GetEntries(truthMatch+" && "+trigger+" && "+all_rectangular_cuts_MC) + t_reco_pid_2018.GetEntries(truthMatch+" && "+trigger+" && "+all_rectangular_cuts_MC) 
            N_rect_cuts_RS_data_all = t_rs_data_2016_noRectCuts.GetEntries(trigger+" && "+all_rectangular_cuts_data) + t_rs_data_2017_noRectCuts.GetEntries(trigger+" && "+all_rectangular_cuts_data) + t_rs_data_2018_noRectCuts.GetEntries(trigger+" && "+all_rectangular_cuts_data)

            if(species == 1):
                N_rect_cuts_3pi3pi_all = t_reco_pid_2016.GetEntries(truthMatch+" && "+trigger+" && "+all_rectangular_cuts_MC+" && (component == 0)") + t_reco_pid_2017.GetEntries(truthMatch+" && "+trigger+" && "+all_rectangular_cuts_MC+" && (component == 0)") + t_reco_pid_2018.GetEntries(truthMatch+" && "+trigger+" && "+all_rectangular_cuts_MC+" && (component == 0)")
                N_rect_cuts_3pi3pipi0_all = t_reco_pid_2016.GetEntries(truthMatch+" && "+trigger+" && "+all_rectangular_cuts_MC+" && (component == 1)") + t_reco_pid_2017.GetEntries(truthMatch+" && "+trigger+" && "+all_rectangular_cuts_MC+" && (component == 1)") + t_reco_pid_2018.GetEntries(truthMatch+" && "+trigger+" && "+all_rectangular_cuts_MC+" && (component == 1)")
                N_rect_cuts_3pi3pi2pi0_all = t_reco_pid_2016.GetEntries(truthMatch+" && "+trigger+" && "+all_rectangular_cuts_MC+" && (component == 2)") + t_reco_pid_2017.GetEntries(truthMatch+" && "+trigger+" && "+all_rectangular_cuts_MC+" && (component == 2)") + t_reco_pid_2018.GetEntries(truthMatch+" && "+trigger+" && "+all_rectangular_cuts_MC+" && (component == 2)")
            
                N_rect_cuts_WS_data_all = t_ws_data_2016_noRectCuts.GetEntries(trigger+" && "+all_rectangular_cuts_data) + t_ws_data_2017_noRectCuts.GetEntries(trigger+" && "+all_rectangular_cuts_data) + t_ws_data_2018_noRectCuts.GetEntries(trigger+" && "+all_rectangular_cuts_data)  
            else:
                N_rect_cuts_WS_data_all = t_rs_data_2016_noRectCuts.GetEntries(trigger+" && "+all_rectangular_cuts_data+" && (Bp_dtf_M[0] > 5320)") + t_rs_data_2017_noRectCuts.GetEntries(trigger+" && "+all_rectangular_cuts_data+" && (Bp_dtf_M[0] > 5320)") + t_rs_data_2018_noRectCuts.GetEntries(trigger+" && "+all_rectangular_cuts_data+" && (Bp_dtf_M[0] > 5320)")

            up_rect_cuts_MC = ROOT.TEfficiency.Wilson(N_trigger, N_rect_cuts_MC, 0.68, True)
            down_rect_cuts_MC = ROOT.TEfficiency.Wilson(N_trigger, N_rect_cuts_MC, 0.68, False)
            eps_rect_cuts_MC_err = 0.5*(up_rect_cuts_MC - down_rect_cuts_MC)
            eps_rect_cuts_MC.append( ufloat(N_rect_cuts_MC/N_trigger, eps_rect_cuts_MC_err) )

            up_rect_cuts_RS_data = ROOT.TEfficiency.Wilson(N_trigger_RS_data, N_rect_cuts_RS_data, 0.68, True)
            down_rect_cuts_RS_data = ROOT.TEfficiency.Wilson(N_trigger_RS_data, N_rect_cuts_RS_data, 0.68, False)
            eps_rect_cuts_RS_data_err = 0.5*(up_rect_cuts_RS_data - down_rect_cuts_RS_data)
            eps_rect_cuts_RS_data.append( ufloat(N_rect_cuts_RS_data/N_trigger_RS_data, eps_rect_cuts_RS_data_err) )
        
            up_rect_cuts_WS_data = ROOT.TEfficiency.Wilson(N_trigger_WS_data, N_rect_cuts_WS_data, 0.68, True)
            down_rect_cuts_WS_data = ROOT.TEfficiency.Wilson(N_trigger_WS_data, N_rect_cuts_WS_data, 0.68, False)
            eps_rect_cuts_WS_data_err = 0.5*(up_rect_cuts_WS_data - down_rect_cuts_WS_data)
            eps_rect_cuts_WS_data.append( ufloat(N_rect_cuts_WS_data/N_trigger_WS_data, eps_rect_cuts_WS_data_err) )

        if(species == 1):
            up_rect_cuts_3pi3pi_all = ROOT.TEfficiency.Wilson(N_trigger_3pi3pi, N_rect_cuts_3pi3pi_all, 0.68, True)
            down_rect_cuts_3pi3pi_all = ROOT.TEfficiency.Wilson(N_trigger_3pi3pi, N_rect_cuts_3pi3pi_all, 0.68, False)
            eps_rect_cuts_3pi3pi_all_err = 0.5*(up_rect_cuts_3pi3pi_all - down_rect_cuts_3pi3pi_all)
            eps_rect_cuts_3pi3pi_all =  ufloat(N_rect_cuts_3pi3pi_all/N_trigger_3pi3pi, eps_rect_cuts_3pi3pi_all_err) 

            up_rect_cuts_3pi3pipi0_all = ROOT.TEfficiency.Wilson(N_trigger_3pi3pipi0, N_rect_cuts_3pi3pipi0_all, 0.68, True)
            down_rect_cuts_3pi3pipi0_all = ROOT.TEfficiency.Wilson(N_trigger_3pi3pipi0, N_rect_cuts_3pi3pipi0_all, 0.68, False)
            eps_rect_cuts_3pi3pipi0_all_err = 0.5*(up_rect_cuts_3pi3pipi0_all - down_rect_cuts_3pi3pipi0_all)
            eps_rect_cuts_3pi3pipi0_all = ufloat(N_rect_cuts_3pi3pipi0_all/N_trigger_3pi3pipi0, eps_rect_cuts_3pi3pipi0_all_err) 

            up_rect_cuts_3pi3pi2pi0_all = ROOT.TEfficiency.Wilson(N_trigger_3pi3pi2pi0, N_rect_cuts_3pi3pi2pi0_all, 0.68, True)
            down_rect_cuts_3pi3pi2pi0_all = ROOT.TEfficiency.Wilson(N_trigger_3pi3pi2pi0, N_rect_cuts_3pi3pi2pi0_all, 0.68, False)
            eps_rect_cuts_3pi3pi2pi0_all_err = 0.5*(up_rect_cuts_3pi3pi2pi0_all - down_rect_cuts_3pi3pi2pi0_all)
            eps_rect_cuts_3pi3pi2pi0_all = ufloat(N_rect_cuts_3pi3pi2pi0_all/N_trigger_3pi3pi2pi0, eps_rect_cuts_3pi3pi2pi0_all_err) 

        print("N_rect_cuts MC: ", N_rect_cuts_MC_all)
        print("N_rect_cuts RS: ", N_rect_cuts_RS_data_all)
        print("N_rect_cuts WS: ", N_rect_cuts_WS_data_all)

        up_rect_cuts_MC_all = ROOT.TEfficiency.Wilson(N_trigger, N_rect_cuts_MC_all, 0.68, True)
        down_rect_cuts_MC_all = ROOT.TEfficiency.Wilson(N_trigger, N_rect_cuts_MC_all, 0.68, False)
        eps_rect_cuts_MC_all_err = 0.5*(up_rect_cuts_MC_all - down_rect_cuts_MC_all)
        eps_rect_cuts_MC_all = ufloat(N_rect_cuts_MC_all/N_trigger, eps_rect_cuts_MC_all_err) 

        up_rect_cuts_RS_data_all = ROOT.TEfficiency.Wilson(N_trigger_RS_data, N_rect_cuts_RS_data_all, 0.68, True)
        down_rect_cuts_RS_data_all = ROOT.TEfficiency.Wilson(N_trigger_RS_data, N_rect_cuts_RS_data_all, 0.68, False)
        eps_rect_cuts_RS_data_all_err = 0.5*(up_rect_cuts_RS_data_all - down_rect_cuts_RS_data_all)
        eps_rect_cuts_RS_data_all = ufloat(N_rect_cuts_RS_data_all/N_trigger_RS_data, eps_rect_cuts_RS_data_all_err) 

        up_rect_cuts_WS_data_all = ROOT.TEfficiency.Wilson(N_trigger_WS_data, N_rect_cuts_WS_data_all, 0.68, True)
        down_rect_cuts_WS_data_all = ROOT.TEfficiency.Wilson(N_trigger_WS_data, N_rect_cuts_WS_data_all, 0.68, False)
        eps_rect_cuts_WS_data_all_err = 0.5*(up_rect_cuts_WS_data_all - down_rect_cuts_WS_data_all)
        eps_rect_cuts_WS_data_all = ufloat(N_rect_cuts_WS_data_all/N_trigger_WS_data, eps_rect_cuts_WS_data_all_err) 


        if(species == 1):
            rect_cuts_columns_1 = ["$3\\pi 3\\pi$", "$3\\pi 3\\pi \\pi^0$", "$3\\pi 3\\pi 2\\pi^0$", "All"]
            rect_cuts_columns_2 = ["MC", "RS data", "WS data"]
        else:
            rect_cuts_columns_1 = ["MC", "Data", "Data upper sideband"]

        rect_cut_table_1 = pd.DataFrame(index=rectangular_cuts_names+["All"], columns=rect_cuts_columns_1)
        if(species == 1):
            rect_cut_table_2 = pd.DataFrame(index=rectangular_cuts_names, columns=rect_cuts_columns_2)


        if(species == 1):
            for i in range(len(rectangular_cuts_names)):
                rect_cut_table_1.loc[rectangular_cuts_names[i], "$3\\pi 3\\pi$"] = f"${eps_rect_cuts_3pi3pi[i].nominal_value*100:.2f} \\pm {eps_rect_cuts_3pi3pi[i].std_dev*100:.2f}$"
                rect_cut_table_1.loc[rectangular_cuts_names[i], "$3\\pi 3\\pi \\pi^0$"] = f"${eps_rect_cuts_3pi3pipi0[i].nominal_value*100:.2f} \\pm {eps_rect_cuts_3pi3pipi0[i].std_dev*100:.2f}$"
                rect_cut_table_1.loc[rectangular_cuts_names[i], "$3\\pi 3\\pi 2\\pi^0$"] = f"${eps_rect_cuts_3pi3pi2pi0[i].nominal_value*100:.2f} \\pm {eps_rect_cuts_3pi3pi2pi0[i].std_dev*100:.2f}$"
                rect_cut_table_1.loc[rectangular_cuts_names[i], "All"] = f"${eps_rect_cuts_MC[i].nominal_value*100:.2f} \\pm {eps_rect_cuts_MC[i].std_dev*100:.2f}$"

                rect_cut_table_2.loc[rectangular_cuts_names[i], "MC"] = f"${eps_rect_cuts_MC[i].nominal_value*100:.2f} \\pm {eps_rect_cuts_MC[i].std_dev*100:.2f}$"
                rect_cut_table_2.loc[rectangular_cuts_names[i], "RS data"] = f"${eps_rect_cuts_RS_data[i].nominal_value*100:.2f} \\pm {eps_rect_cuts_RS_data[i].std_dev*100:.2f}$"
                rect_cut_table_2.loc[rectangular_cuts_names[i], "WS data"] = f"${eps_rect_cuts_WS_data[i].nominal_value*100:.2f} \\pm {eps_rect_cuts_WS_data[i].std_dev*100:.2f}$"
        
            rect_cut_table_1.loc["All", "$3\\pi 3\\pi$"] = f"${eps_rect_cuts_3pi3pi_all.nominal_value*100:.2f} \\pm {eps_rect_cuts_3pi3pi_all.std_dev*100:.2f}$"
            rect_cut_table_1.loc["All", "$3\\pi 3\\pi \\pi^0$"] = f"${eps_rect_cuts_3pi3pipi0_all.nominal_value*100:.2f} \\pm {eps_rect_cuts_3pi3pipi0_all.std_dev*100:.2f}$"
            rect_cut_table_1.loc["All", "$3\\pi 3\\pi 2\\pi^0$"] = f"${eps_rect_cuts_3pi3pi2pi0_all.nominal_value*100:.2f} \\pm {eps_rect_cuts_3pi3pi2pi0_all.std_dev*100:.2f}$"
            rect_cut_table_1.loc["All", "All"] = f"${eps_rect_cuts_MC_all.nominal_value*100:.2f} \\pm {eps_rect_cuts_MC_all.std_dev*100:.2f}$"

            rect_cut_table_2.loc["All", "MC"] = f"${eps_rect_cuts_MC_all.nominal_value*100:.2f} \\pm {eps_rect_cuts_MC_all.std_dev*100:.2f}$"
            rect_cut_table_2.loc["All", "RS data"] = f"${eps_rect_cuts_RS_data_all.nominal_value*100:.2f} \\pm {eps_rect_cuts_RS_data_all.std_dev*100:.2f}$"
            rect_cut_table_2.loc["All", "WS data"] = f"${eps_rect_cuts_WS_data_all.nominal_value*100:.2f} \\pm {eps_rect_cuts_WS_data_all.std_dev*100:.2f}$"

        else:
            for i in range(len(rectangular_cuts_names)):
                rect_cut_table_1.loc[rectangular_cuts_names[i], "MC"] = f"${eps_rect_cuts_MC[i].nominal_value*100:.2f} \\pm {eps_rect_cuts_MC[i].std_dev*100:.2f}$"
                rect_cut_table_1.loc[rectangular_cuts_names[i], "Data"] = f"${eps_rect_cuts_RS_data[i].nominal_value*100:.2f} \\pm {eps_rect_cuts_RS_data[i].std_dev*100:.2f}$"
                rect_cut_table_1.loc[rectangular_cuts_names[i], "Data upper sideband"] = f"${eps_rect_cuts_WS_data[i].nominal_value*100:.2f} \\pm {eps_rect_cuts_WS_data[i].std_dev*100:.2f}$"

            rect_cut_table_1.loc["All", "MC"] = f"${eps_rect_cuts_MC_all.nominal_value*100:.2f} \\pm {eps_rect_cuts_MC_all.std_dev*100:.2f}$"
            rect_cut_table_1.loc["All", "Data"] = f"${eps_rect_cuts_RS_data_all.nominal_value*100:.2f} \\pm {eps_rect_cuts_RS_data_all.std_dev*100:.2f}$"
            rect_cut_table_1.loc["All", "Data upper sideband"] = f"${eps_rect_cuts_WS_data_all.nominal_value*100:.2f} \\pm {eps_rect_cuts_WS_data_all.std_dev*100:.2f}$"


    with open(f'/panfs/felician/B2Ktautau/workflow/selections_efficiency_tables/Species_{species}/rectangular_cuts_efficiency_table.tex', 'w') as fout_rect_cuts:
        fout_rect_cuts.write(rect_cut_table_1.to_latex())
    
    if(species == 1):
        with open(f'/panfs/felician/B2Ktautau/workflow/selections_efficiency_tables/Species_{species}/rectangular_cuts_efficiency_table_2.tex', 'w') as fout_rect_cuts_2:
            fout_rect_cuts_2.write(rect_cut_table_2.to_latex())

    # Pass fitter
    print("Pass mass fit efficiency")
    if(species == 100):
        N_pass_fit_BuDDKp = t_reco_BuDDKp_2016.GetEntries(trigger+" && "+all_rectangular_cuts+" && "+pass_mass_fit) + t_reco_BuDDKp_2017.GetEntries(trigger+" && "+all_rectangular_cuts+" && "+pass_mass_fit) + t_reco_BuDDKp_2018.GetEntries(trigger+" && "+all_rectangular_cuts+" && "+pass_mass_fit) 
        N_pass_fit_BdDDKp = t_reco_BdDDKp_2016.GetEntries(trigger+" && "+all_rectangular_cuts+" && "+pass_mass_fit) + t_reco_BdDDKp_2017.GetEntries(trigger+" && "+all_rectangular_cuts+" && "+pass_mass_fit) + t_reco_BdDDKp_2018.GetEntries(trigger+" && "+all_rectangular_cuts+" && "+pass_mass_fit)
        N_pass_fit_BsDDKp = t_reco_BsDDKp_2016.GetEntries(trigger+" && "+all_rectangular_cuts+" && "+pass_mass_fit) + t_reco_BsDDKp_2017.GetEntries(trigger+" && "+all_rectangular_cuts+" && "+pass_mass_fit) + t_reco_BsDDKp_2018.GetEntries(trigger+" && "+all_rectangular_cuts+" && "+pass_mass_fit)
        N_pass_fit_BuDDK0 = t_reco_BuDDK0_2016.GetEntries(trigger+" && "+all_rectangular_cuts+" && "+pass_mass_fit) + t_reco_BuDDK0_2017.GetEntries(trigger+" && "+all_rectangular_cuts+" && "+pass_mass_fit) + t_reco_BuDDK0_2018.GetEntries(trigger+" && "+all_rectangular_cuts+" && "+pass_mass_fit)
        N_pass_fit_BuDD = t_reco_BuDD_2016.GetEntries(trigger+" && "+all_rectangular_cuts+" && "+pass_mass_fit) + t_reco_BuDD_2017.GetEntries(trigger+" && "+all_rectangular_cuts+" && "+pass_mass_fit) + t_reco_BuDD_2018.GetEntries(trigger+" && "+all_rectangular_cuts+" && "+pass_mass_fit)

        up_pass_fit_BuDDKp = ROOT.TEfficiency.Wilson(N_rect_cuts_BuDDKp, N_pass_fit_BuDDKp, 0.68, True)
        down_pass_fit_BuDDKp = ROOT.TEfficiency.Wilson(N_rect_cuts_BuDDKp, N_pass_fit_BuDDKp, 0.68, False)
        eps_pass_fit_BuDDKp_err = 0.5*(up_pass_fit_BuDDKp - down_pass_fit_BuDDKp)
        eps_pass_fit_BuDDKp = ufloat(N_pass_fit_BuDDKp/N_rect_cuts_BuDDKp, eps_pass_fit_BuDDKp_err)

        up_pass_fit_BdDDKp = ROOT.TEfficiency.Wilson(N_rect_cuts_BdDDKp, N_pass_fit_BdDDKp, 0.68, True)
        down_pass_fit_BdDDKp = ROOT.TEfficiency.Wilson(N_rect_cuts_BdDDKp, N_pass_fit_BdDDKp, 0.68, False)
        eps_pass_fit_BdDDKp_err = 0.5*(up_pass_fit_BdDDKp - down_pass_fit_BdDDKp)
        eps_pass_fit_BdDDKp = ufloat(N_pass_fit_BdDDKp/N_rect_cuts_BdDDKp, eps_pass_fit_BdDDKp_err)

        up_pass_fit_BsDDKp = ROOT.TEfficiency.Wilson(N_rect_cuts_BsDDKp, N_pass_fit_BsDDKp, 0.68, True)
        down_pass_fit_BsDDKp = ROOT.TEfficiency.Wilson(N_rect_cuts_BsDDKp, N_pass_fit_BsDDKp, 0.68, False)
        eps_pass_fit_BsDDKp_err = 0.5*(up_pass_fit_BsDDKp - down_pass_fit_BsDDKp)
        eps_pass_fit_BsDDKp = ufloat(N_pass_fit_BsDDKp/N_rect_cuts_BsDDKp, eps_pass_fit_BsDDKp_err)

        up_pass_fit_BuDDK0 = ROOT.TEfficiency.Wilson(N_rect_cuts_BuDDK0, N_pass_fit_BuDDK0, 0.68, True)
        down_pass_fit_BuDDK0 = ROOT.TEfficiency.Wilson(N_rect_cuts_BuDDK0, N_pass_fit_BuDDK0, 0.68, False)
        eps_pass_fit_BuDDK0_err = 0.5*(up_pass_fit_BuDDK0 - down_pass_fit_BuDDK0)
        eps_pass_fit_BuDDK0 = ufloat(N_pass_fit_BuDDK0/N_rect_cuts_BuDDK0, eps_pass_fit_BuDDK0_err)

        up_pass_fit_BuDD = ROOT.TEfficiency.Wilson(N_rect_cuts_BuDD, N_pass_fit_BuDD, 0.68, True)
        down_pass_fit_BuDD = ROOT.TEfficiency.Wilson(N_rect_cuts_BuDD, N_pass_fit_BuDD, 0.68, False)
        eps_pass_fit_BuDD_err = 0.5*(up_pass_fit_BuDD - down_pass_fit_BuDD)
        eps_pass_fit_BuDD = ufloat(N_pass_fit_BuDD/N_rect_cuts_BuDD, eps_pass_fit_BuDD_err)

        eps_pass_fit = [eps_pass_fit_BuDDKp, eps_pass_fit_BdDDKp, eps_pass_fit_BsDDKp, eps_pass_fit_BuDDK0, eps_pass_fit_BuDD]

        pass_fit_columns = ["Pass mass fit (\\%)"]
        pass_fit_rows = ["$B^+ \\to D D K^+$", "$B^0 \\to D D K^+$", "$B^0_s \\to D D K^+$", "$B^+ \\to D D K^0$", "$B^+ \\to D D$"]

        pass_fit_table = pd.DataFrame(index=pass_fit_rows, columns=pass_fit_columns)
        pass_fit_table.loc["$B^+ \\to D D K^+$", "Pass mass fit (\\%)"] = f"${eps_pass_fit[0].nominal_value*100:.2f} \\pm {eps_pass_fit[0].std_dev*100:.2f}$"
        pass_fit_table.loc["$B^0 \\to D D K^+$", "Pass mass fit (\\%)"] = f"${eps_pass_fit[1].nominal_value*100:.2f} \\pm {eps_pass_fit[1].std_dev*100:.2f}$"
        pass_fit_table.loc["$B^0_s \\to D D K^+$", "Pass mass fit (\\%)"] = f"${eps_pass_fit[2].nominal_value*100:.2f} \\pm {eps_pass_fit[2].std_dev*100:.2f}$"
        pass_fit_table.loc["$B^+ \\to D D K^0$", "Pass mass fit (\\%)"] = f"${eps_pass_fit[3].nominal_value*100:.2f} \\pm {eps_pass_fit[3].std_dev*100:.2f}$"
        pass_fit_table.loc["$B^+ \\to D D$", "Pass mass fit (\\%)"] = f"${eps_pass_fit[4].nominal_value*100:.2f} \\pm {eps_pass_fit[4].std_dev*100:.2f}$"


    else:
        if(species == 1):
            fc_reco_mc_2016 = ROOT.TFileCollection("fc_reco_mc_2016", "fc_reco_mc_2016", "/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2016/Species_10/pre_sel_tree.txt")
            fc_reco_mc_2017 = ROOT.TFileCollection("fc_reco_mc_2017", "fc_reco_mc_2017", "/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2017/Species_10/pre_sel_tree.txt")
            fc_reco_mc_2018 = ROOT.TFileCollection("fc_reco_mc_2018", "fc_reco_mc_2018", "/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2018/Species_10/pre_sel_tree.txt")

            t_reco_mc_2016 = ROOT.TChain("DecayTree")
            t_reco_mc_2017 = ROOT.TChain("DecayTree")
            t_reco_mc_2018 = ROOT.TChain("DecayTree")

            t_reco_mc_2016.AddFileInfoList(fc_reco_mc_2016.GetList())
            t_reco_mc_2017.AddFileInfoList(fc_reco_mc_2017.GetList())
            t_reco_mc_2018.AddFileInfoList(fc_reco_mc_2018.GetList())

            fc_gsl_mc_2016 = ROOT.TFileCollection("fc_gsl_mc_2016", "fc_gsl_mc_2016", "/panfs/felician/B2Ktautau/workflow/standalone_fitter/2016/Species_10/fit_results.txt")
            fc_gsl_mc_2017 = ROOT.TFileCollection("fc_gsl_mc_2017", "fc_gsl_mc_2017", "/panfs/felician/B2Ktautau/workflow/standalone_fitter/2017/Species_10/fit_results.txt")
            fc_gsl_mc_2018 = ROOT.TFileCollection("fc_gsl_mc_2018", "fc_gsl_mc_2018", "/panfs/felician/B2Ktautau/workflow/standalone_fitter/2018/Species_10/fit_results.txt")

            t_gsl_mc_2016 = ROOT.TChain("DecayTree")
            t_gsl_mc_2017 = ROOT.TChain("DecayTree")
            t_gsl_mc_2018 = ROOT.TChain("DecayTree")

            t_gsl_mc_2016.AddFileInfoList(fc_gsl_mc_2016.GetList())
            t_gsl_mc_2017.AddFileInfoList(fc_gsl_mc_2017.GetList())
            t_gsl_mc_2018.AddFileInfoList(fc_gsl_mc_2018.GetList())

            fc_mass_mc_2016 =  ROOT.TFileCollection("fc_mass_mc_2016", "fc_mass_mc_2016", "/panfs/felician/B2Ktautau/workflow/create_invariant_mass_tree/2016/Species_10/invariant_mass_tree.txt")
            fc_mass_mc_2017 =  ROOT.TFileCollection("fc_mass_mc_2017", "fc_mass_mc_2017", "/panfs/felician/B2Ktautau/workflow/create_invariant_mass_tree/2017/Species_10/invariant_mass_tree.txt")
            fc_mass_mc_2018 =  ROOT.TFileCollection("fc_mass_mc_2018", "fc_mass_mc_2018", "/panfs/felician/B2Ktautau/workflow/create_invariant_mass_tree/2018/Species_10/invariant_mass_tree.txt")

            t_mass_mc_2016 = ROOT.TChain("DecayTree")
            t_mass_mc_2017 = ROOT.TChain("DecayTree")
            t_mass_mc_2018 = ROOT.TChain("DecayTree")

            t_mass_mc_2016.AddFileInfoList(fc_mass_mc_2016.GetList())
            t_mass_mc_2017.AddFileInfoList(fc_mass_mc_2017.GetList())
            t_mass_mc_2018.AddFileInfoList(fc_mass_mc_2018.GetList())

            t_reco_mc_2016.AddFriend(t_gsl_mc_2016)
            t_reco_mc_2016.AddFriend(t_mass_mc_2016)

            t_reco_mc_2017.AddFriend(t_gsl_mc_2017)
            t_reco_mc_2017.AddFriend(t_mass_mc_2017)

            t_reco_mc_2018.AddFriend(t_gsl_mc_2018)
            t_reco_mc_2018.AddFriend(t_mass_mc_2018)

            N_pass_fitter_3pi3pi = t_reco_mc_2016.GetEntries(truthMatch+" && "+trigger+" && "+all_rectangular_cuts_MC+" && "+pass_mass_fit+" && (component == 0)") + t_reco_mc_2017.GetEntries(truthMatch+" && "+trigger+" && "+all_rectangular_cuts_MC+" && "+pass_mass_fit+" && (component == 0)") + t_reco_mc_2018.GetEntries(truthMatch+" && "+trigger+" && "+all_rectangular_cuts_MC+" && "+pass_mass_fit+" && (component == 0)")
            N_pass_fitter_3pi3pipi0 = t_reco_mc_2016.GetEntries(truthMatch+" && "+trigger+" && "+all_rectangular_cuts_MC+" && "+pass_mass_fit+" && (component == 1)") + t_reco_mc_2017.GetEntries(truthMatch+" && "+trigger+" && "+all_rectangular_cuts_MC+" && "+pass_mass_fit+" && (component == 1)") + t_reco_mc_2018.GetEntries(truthMatch+" && "+trigger+" && "+all_rectangular_cuts_MC+" && "+pass_mass_fit+" && (component == 1)")
            N_pass_fitter_3pi3pi2pi0 = t_reco_mc_2016.GetEntries(truthMatch+" && "+trigger+" && "+all_rectangular_cuts_MC+" && "+pass_mass_fit+" && (component == 2)") + t_reco_mc_2017.GetEntries(truthMatch+" && "+trigger+" && "+all_rectangular_cuts_MC+" && "+pass_mass_fit+" && (component == 2)") + t_reco_mc_2018.GetEntries(truthMatch+" && "+trigger+" && "+all_rectangular_cuts_MC+" && "+pass_mass_fit+" && (component == 2)")

            N_pass_cuts_3pi3pi = t_reco_mc_2016.GetEntries(truthMatch+" && "+trigger+" && "+all_rectangular_cuts_MC+" && (component == 0)") + t_reco_mc_2017.GetEntries(truthMatch+" && "+trigger+" && "+all_rectangular_cuts_MC+" && (component == 0)") + t_reco_mc_2018.GetEntries(truthMatch+" && "+trigger+" && "+all_rectangular_cuts_MC+" && (component == 0)")
            N_pass_cuts_3pi3pipi0 = t_reco_mc_2016.GetEntries(truthMatch+" && "+trigger+" && "+all_rectangular_cuts_MC+" && (component == 1)") + t_reco_mc_2017.GetEntries(truthMatch+" && "+trigger+" && "+all_rectangular_cuts_MC+" && (component == 1)") + t_reco_mc_2018.GetEntries(truthMatch+" && "+trigger+" && "+all_rectangular_cuts_MC+" && (component == 1)")
            N_pass_cuts_3pi3pi2pi0 = t_reco_mc_2016.GetEntries(truthMatch+" && "+trigger+" && "+all_rectangular_cuts_MC+" && (component == 2)") + t_reco_mc_2017.GetEntries(truthMatch+" && "+trigger+" && "+all_rectangular_cuts_MC+" && (component == 2)") + t_reco_mc_2018.GetEntries(truthMatch+" && "+trigger+" && "+all_rectangular_cuts_MC+" && (component == 2)")

            up_pass_fit_3pi3pi = ROOT.TEfficiency.Wilson(N_pass_cuts_3pi3pi, N_pass_fitter_3pi3pi, 0.68, True)
            down_pass_fit_3pi3pi = ROOT.TEfficiency.Wilson(N_pass_cuts_3pi3pi, N_pass_fitter_3pi3pi, 0.68, False)
            eps_pass_fit_3pi3pi_err = 0.5*(up_pass_fit_3pi3pi - down_pass_fit_3pi3pi)
            eps_pass_fit_3pi3pi = ufloat(N_pass_fitter_3pi3pi/N_pass_cuts_3pi3pi, eps_pass_fit_3pi3pi_err)

            up_pass_fit_3pi3pipi0 = ROOT.TEfficiency.Wilson(N_pass_cuts_3pi3pipi0, N_pass_fitter_3pi3pipi0, 0.68, True)
            down_pass_fit_3pi3pipi0 = ROOT.TEfficiency.Wilson(N_pass_cuts_3pi3pipi0, N_pass_fitter_3pi3pipi0, 0.68, False)
            eps_pass_fit_3pi3pipi0_err = 0.5*(up_pass_fit_3pi3pipi0 - down_pass_fit_3pi3pipi0)
            eps_pass_fit_3pi3pipi0 = ufloat(N_pass_fitter_3pi3pipi0/N_pass_cuts_3pi3pipi0, eps_pass_fit_3pi3pipi0_err)

            up_pass_fit_3pi3pi2pi0 = ROOT.TEfficiency.Wilson(N_pass_cuts_3pi3pi2pi0, N_pass_fitter_3pi3pi2pi0, 0.68, True)
            down_pass_fit_3pi3pi2pi0 = ROOT.TEfficiency.Wilson(N_pass_cuts_3pi3pi2pi0, N_pass_fitter_3pi3pi2pi0, 0.68, False)
            eps_pass_fit_3pi3pi2pi0_err = 0.5*(up_pass_fit_3pi3pi2pi0 - down_pass_fit_3pi3pi2pi0)
            eps_pass_fit_3pi3pi2pi0 = ufloat(N_pass_fitter_3pi3pi2pi0/N_pass_cuts_3pi3pi2pi0, eps_pass_fit_3pi3pi2pi0_err)

            fc_reco_rs_2016 = ROOT.TFileCollection("fc_reco_rs_2016", "fc_reco_rs_2016", "/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2016/Species_2/pre_sel_tree.txt", 10)
            fc_reco_rs_2017 = ROOT.TFileCollection("fc_reco_rs_2017", "fc_reco_rs_2017", "/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2017/Species_2/pre_sel_tree.txt", 10)
            fc_reco_rs_2018 = ROOT.TFileCollection("fc_reco_rs_2018", "fc_reco_rs_2018", "/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2018/Species_2/pre_sel_tree.txt", 10)

            t_reco_rs_2016 = ROOT.TChain("DecayTree")
            t_reco_rs_2017 = ROOT.TChain("DecayTree")
            t_reco_rs_2018 = ROOT.TChain("DecayTree")

            t_reco_rs_2016.AddFileInfoList(fc_reco_rs_2016.GetList())
            t_reco_rs_2017.AddFileInfoList(fc_reco_rs_2017.GetList())
            t_reco_rs_2018.AddFileInfoList(fc_reco_rs_2018.GetList())

            fc_gsl_rs_2016 = ROOT.TFileCollection("fc_gsl_rs_2016", "fc_gsl_rs_2016", "/panfs/felician/B2Ktautau/workflow/standalone_fitter/2016/Species_2/fit_results.txt", 10)
            fc_gsl_rs_2017 = ROOT.TFileCollection("fc_gsl_rs_2017", "fc_gsl_rs_2017", "/panfs/felician/B2Ktautau/workflow/standalone_fitter/2017/Species_2/fit_results.txt", 10)
            fc_gsl_rs_2018 = ROOT.TFileCollection("fc_gsl_rs_2018", "fc_gsl_rs_2018", "/panfs/felician/B2Ktautau/workflow/standalone_fitter/2018/Species_2/fit_results.txt", 10)

            t_gsl_rs_2016 = ROOT.TChain("DecayTree")
            t_gsl_rs_2017 = ROOT.TChain("DecayTree")
            t_gsl_rs_2018 = ROOT.TChain("DecayTree")

            t_gsl_rs_2016.AddFileInfoList(fc_gsl_rs_2016.GetList())
            t_gsl_rs_2017.AddFileInfoList(fc_gsl_rs_2017.GetList())
            t_gsl_rs_2018.AddFileInfoList(fc_gsl_rs_2018.GetList())

            fc_mass_rs_2016 =  ROOT.TFileCollection("fc_mass_rs_2016", "fc_mass_rs_2016", "/panfs/felician/B2Ktautau/workflow/create_invariant_mass_tree/2016/Species_2/invariant_mass_tree.txt", 10)
            fc_mass_rs_2017 =  ROOT.TFileCollection("fc_mass_rs_2017", "fc_mass_rs_2017", "/panfs/felician/B2Ktautau/workflow/create_invariant_mass_tree/2017/Species_2/invariant_mass_tree.txt", 10)
            fc_mass_rs_2018 =  ROOT.TFileCollection("fc_mass_rs_2018", "fc_mass_rs_2018", "/panfs/felician/B2Ktautau/workflow/create_invariant_mass_tree/2018/Species_2/invariant_mass_tree.txt", 10)

            t_mass_rs_2016 = ROOT.TChain("DecayTree")
            t_mass_rs_2017 = ROOT.TChain("DecayTree")
            t_mass_rs_2018 = ROOT.TChain("DecayTree")

            t_mass_rs_2016.AddFileInfoList(fc_mass_rs_2016.GetList())
            t_mass_rs_2017.AddFileInfoList(fc_mass_rs_2017.GetList())
            t_mass_rs_2018.AddFileInfoList(fc_mass_rs_2018.GetList())

            t_reco_rs_2016.AddFriend(t_gsl_rs_2016)
            t_reco_rs_2016.AddFriend(t_mass_rs_2016)

            t_reco_rs_2017.AddFriend(t_gsl_rs_2017)
            t_reco_rs_2017.AddFriend(t_mass_rs_2017)

            t_reco_rs_2018.AddFriend(t_gsl_rs_2018)
            t_reco_rs_2018.AddFriend(t_mass_rs_2018)

            fc_reco_ws_2016 = ROOT.TFileCollection("fc_reco_ws_2016", "fc_reco_ws_2016", "/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2016/Species_3/pre_sel_tree.txt", 10)
            fc_reco_ws_2017 = ROOT.TFileCollection("fc_reco_ws_2017", "fc_reco_ws_2017", "/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2017/Species_3/pre_sel_tree.txt", 10)
            fc_reco_ws_2018 = ROOT.TFileCollection("fc_reco_ws_2018", "fc_reco_ws_2018", "/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2018/Species_3/pre_sel_tree.txt", 10)

            t_reco_ws_2016 = ROOT.TChain("DecayTree")
            t_reco_ws_2017 = ROOT.TChain("DecayTree")
            t_reco_ws_2018 = ROOT.TChain("DecayTree")

            t_reco_ws_2016.AddFileInfoList(fc_reco_ws_2016.GetList())
            t_reco_ws_2017.AddFileInfoList(fc_reco_ws_2017.GetList())
            t_reco_ws_2018.AddFileInfoList(fc_reco_ws_2018.GetList())

            fc_gsl_ws_2016 = ROOT.TFileCollection("fc_gsl_ws_2016", "fc_gsl_ws_2016", "/panfs/felician/B2Ktautau/workflow/standalone_fitter/2016/Species_3/fit_results.txt", 10)
            fc_gsl_ws_2017 = ROOT.TFileCollection("fc_gsl_ws_2017", "fc_gsl_ws_2017", "/panfs/felician/B2Ktautau/workflow/standalone_fitter/2017/Species_3/fit_results.txt", 10)
            fc_gsl_ws_2018 = ROOT.TFileCollection("fc_gsl_ws_2018", "fc_gsl_ws_2018", "/panfs/felician/B2Ktautau/workflow/standalone_fitter/2018/Species_3/fit_results.txt", 10)

            t_gsl_ws_2016 = ROOT.TChain("DecayTree")
            t_gsl_ws_2017 = ROOT.TChain("DecayTree")
            t_gsl_ws_2018 = ROOT.TChain("DecayTree")

            t_gsl_ws_2016.AddFileInfoList(fc_gsl_ws_2016.GetList())
            t_gsl_ws_2017.AddFileInfoList(fc_gsl_ws_2017.GetList())
            t_gsl_ws_2018.AddFileInfoList(fc_gsl_ws_2018.GetList())

            fc_mass_ws_2016 =  ROOT.TFileCollection("fc_mass_ws_2016", "fc_mass_ws_2016", "/panfs/felician/B2Ktautau/workflow/create_invariant_mass_tree/2016/Species_3/invariant_mass_tree.txt", 10)
            fc_mass_ws_2017 =  ROOT.TFileCollection("fc_mass_ws_2017", "fc_mass_ws_2017", "/panfs/felician/B2Ktautau/workflow/create_invariant_mass_tree/2017/Species_3/invariant_mass_tree.txt", 10)
            fc_mass_ws_2018 =  ROOT.TFileCollection("fc_mass_ws_2018", "fc_mass_ws_2018", "/panfs/felician/B2Ktautau/workflow/create_invariant_mass_tree/2018/Species_3/invariant_mass_tree.txt", 10)

            t_mass_ws_2016 = ROOT.TChain("DecayTree")
            t_mass_ws_2017 = ROOT.TChain("DecayTree")
            t_mass_ws_2018 = ROOT.TChain("DecayTree")

            t_mass_ws_2016.AddFileInfoList(fc_mass_ws_2016.GetList())
            t_mass_ws_2017.AddFileInfoList(fc_mass_ws_2017.GetList())
            t_mass_ws_2018.AddFileInfoList(fc_mass_ws_2018.GetList())

            t_reco_ws_2016.AddFriend(t_gsl_ws_2016)
            t_reco_ws_2016.AddFriend(t_mass_ws_2016)

            t_reco_ws_2017.AddFriend(t_gsl_ws_2017)
            t_reco_ws_2017.AddFriend(t_mass_ws_2017)

            t_reco_ws_2018.AddFriend(t_gsl_ws_2018)
            t_reco_ws_2018.AddFriend(t_mass_ws_2018)

            N_pass_fitter_WS = t_reco_ws_2016.GetEntries(trigger+" && "+all_rectangular_cuts_data+" && "+pass_mass_fit) + t_reco_ws_2017.GetEntries(trigger+" && "+all_rectangular_cuts_data+" && "+pass_mass_fit) + t_reco_ws_2018.GetEntries(trigger+" && "+all_rectangular_cuts_data+" && "+pass_mass_fit)
            N_cuts_WS = t_reco_ws_2016.GetEntries(trigger+" && "+all_rectangular_cuts_data) + t_reco_ws_2017.GetEntries(trigger+" && "+all_rectangular_cuts_data) + t_reco_ws_2018.GetEntries(trigger+" && "+all_rectangular_cuts_data)

            up_pass_fit_WS = ROOT.TEfficiency.Wilson(N_cuts_WS, N_pass_fitter_WS, 0.68, True)
            down_pass_fit_WS = ROOT.TEfficiency.Wilson(N_cuts_WS, N_pass_fitter_WS, 0.68, False)
            eps_pass_fit_WS_err = 0.5*(up_pass_fit_WS - down_pass_fit_WS)
            eps_pass_fit_WS = ufloat(N_pass_fitter_WS/N_cuts_WS, eps_pass_fit_WS_err)
        
        else:
            fc_reco_mc_2016 = ROOT.TFileCollection("fc_reco_mc_2016", "fc_reco_mc_2016", "/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2016/Species_72/pre_sel_tree.txt")
            fc_reco_mc_2017 = ROOT.TFileCollection("fc_reco_mc_2017", "fc_reco_mc_2017", "/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2017/Species_72/pre_sel_tree.txt")
            fc_reco_mc_2018 = ROOT.TFileCollection("fc_reco_mc_2018", "fc_reco_mc_2018", "/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2018/Species_72/pre_sel_tree.txt")

            t_reco_mc_2016 = ROOT.TChain("DecayTree")
            t_reco_mc_2017 = ROOT.TChain("DecayTree")
            t_reco_mc_2018 = ROOT.TChain("DecayTree")

            t_reco_mc_2016.AddFileInfoList(fc_reco_mc_2016.GetList())
            t_reco_mc_2017.AddFileInfoList(fc_reco_mc_2017.GetList())
            t_reco_mc_2018.AddFileInfoList(fc_reco_mc_2018.GetList())

            fc_best_cand_mc_2016 = ROOT.TFileCollection("fc_best_cand_mc_2016", "fc_best_cand_mc_2016", "/panfs/felician/B2Ktautau/workflow/multiple_events/2016/Species_72/multiple_events.txt")
            fc_best_cand_mc_2017 = ROOT.TFileCollection("fc_best_cand_mc_2017", "fc_best_cand_mc_2017", "/panfs/felician/B2Ktautau/workflow/multiple_events/2017/Species_72/multiple_events.txt")
            fc_best_cand_mc_2018 = ROOT.TFileCollection("fc_best_cand_mc_2018", "fc_best_cand_mc_2018", "/panfs/felician/B2Ktautau/workflow/multiple_events/2018/Species_72/multiple_events.txt")

            t_best_cand_mc_2016 = ROOT.TChain("DecayTree")
            t_best_cand_mc_2017 = ROOT.TChain("DecayTree")
            t_best_cand_mc_2018 = ROOT.TChain("DecayTree")

            t_best_cand_mc_2016.AddFileInfoList(fc_best_cand_mc_2016.GetList())
            t_best_cand_mc_2017.AddFileInfoList(fc_best_cand_mc_2017.GetList())
            t_best_cand_mc_2018.AddFileInfoList(fc_best_cand_mc_2018.GetList())

            t_reco_mc_2016.AddFriend(t_best_cand_mc_2016)
            t_reco_mc_2017.AddFriend(t_best_cand_mc_2017)
            t_reco_mc_2018.AddFriend(t_best_cand_mc_2018)

            fc_reco_rs_2016 = ROOT.TFileCollection("fc_reco_rs_2016", "fc_reco_rs_2016", "/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2016/Species_8/pre_sel_tree.txt", 100)
            fc_reco_rs_2017 = ROOT.TFileCollection("fc_reco_rs_2017", "fc_reco_rs_2017", "/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2017/Species_8/pre_sel_tree.txt", 100)
            fc_reco_rs_2018 = ROOT.TFileCollection("fc_reco_rs_2018", "fc_reco_rs_2018", "/panfs/felician/B2Ktautau/workflow/create_pre_selection_tree/2018/Species_8/pre_sel_tree.txt", 100)

            t_reco_rs_2016 = ROOT.TChain("DecayTree")
            t_reco_rs_2017 = ROOT.TChain("DecayTree")
            t_reco_rs_2018 = ROOT.TChain("DecayTree")

            t_reco_rs_2016.AddFileInfoList(fc_reco_rs_2016.GetList())
            t_reco_rs_2017.AddFileInfoList(fc_reco_rs_2017.GetList())
            t_reco_rs_2018.AddFileInfoList(fc_reco_rs_2018.GetList())

            fc_best_cand_rs_2016 = ROOT.TFileCollection("fc_best_cand_rs_2016", "fc_best_cand_rs_2016", "/panfs/felician/B2Ktautau/workflow/multiple_events/2016/Species_8/multiple_events.txt", 100)
            fc_best_cand_rs_2017 = ROOT.TFileCollection("fc_best_cand_rs_2017", "fc_best_cand_rs_2017", "/panfs/felician/B2Ktautau/workflow/multiple_events/2017/Species_8/multiple_events.txt", 100)
            fc_best_cand_rs_2018 = ROOT.TFileCollection("fc_best_cand_rs_2018", "fc_best_cand_rs_2018", "/panfs/felician/B2Ktautau/workflow/multiple_events/2018/Species_8/multiple_events.txt", 100)

            t_best_cand_rs_2016 = ROOT.TChain("DecayTree")
            t_best_cand_rs_2017 = ROOT.TChain("DecayTree")
            t_best_cand_rs_2018 = ROOT.TChain("DecayTree")

            t_best_cand_rs_2016.AddFileInfoList(fc_best_cand_rs_2016.GetList())
            t_best_cand_rs_2017.AddFileInfoList(fc_best_cand_rs_2017.GetList())
            t_best_cand_rs_2018.AddFileInfoList(fc_best_cand_rs_2018.GetList())

            t_reco_rs_2016.AddFriend(t_best_cand_rs_2016)
            t_reco_rs_2017.AddFriend(t_best_cand_rs_2017)
            t_reco_rs_2018.AddFriend(t_best_cand_rs_2018)


        if(species == 1):
            N_pass_fitter_MC = t_reco_mc_2016.GetEntries(truthMatch+" && "+trigger+" && "+all_rectangular_cuts_MC+" && "+pass_mass_fit) + t_reco_mc_2017.GetEntries(truthMatch+" && "+trigger+" && "+all_rectangular_cuts_MC+" && "+pass_mass_fit) + t_reco_mc_2018.GetEntries(truthMatch+" && "+trigger+" && "+all_rectangular_cuts_MC+" && "+pass_mass_fit)
            N_cuts_MC = t_reco_mc_2016.GetEntries(truthMatch+" && "+trigger+" && "+all_rectangular_cuts_MC) + t_reco_mc_2017.GetEntries(truthMatch+" && "+trigger+" && "+all_rectangular_cuts_MC) + t_reco_mc_2018.GetEntries(truthMatch+" && "+trigger+" && "+all_rectangular_cuts_MC)
        
            N_pass_fitter_RS = t_reco_rs_2016.GetEntries(trigger+" && "+all_rectangular_cuts_data+" && "+pass_mass_fit) + t_reco_rs_2017.GetEntries(trigger+" && "+all_rectangular_cuts_data+" && "+pass_mass_fit) + t_reco_rs_2018.GetEntries(trigger+" && "+all_rectangular_cuts_data+" && "+pass_mass_fit)
            N_cuts_RS = t_reco_rs_2016.GetEntries(trigger+" && "+all_rectangular_cuts_data) + t_reco_rs_2017.GetEntries(trigger+" && "+all_rectangular_cuts_data) + t_reco_rs_2018.GetEntries(trigger+" && "+all_rectangular_cuts_data)

            up_pass_fit_MC = ROOT.TEfficiency.Wilson(N_cuts_MC, N_pass_fitter_MC, 0.68, True)
            down_pass_fit_MC = ROOT.TEfficiency.Wilson(N_cuts_MC, N_pass_fitter_MC, 0.68, False)
            eps_pass_fit_MC_err = 0.5*(up_pass_fit_MC - down_pass_fit_MC)
            eps_pass_fit_MC = ufloat(N_pass_fitter_MC/N_cuts_MC, eps_pass_fit_MC_err)

            up_pass_fit_RS = ROOT.TEfficiency.Wilson(N_cuts_RS, N_pass_fitter_RS, 0.68, True)
            down_pass_fit_RS = ROOT.TEfficiency.Wilson(N_cuts_RS, N_pass_fitter_RS, 0.68, False)
            eps_pass_fit_RS_err = 0.5*(up_pass_fit_RS - down_pass_fit_RS)
            eps_pass_fit_RS = ufloat(N_pass_fitter_RS/N_cuts_RS, eps_pass_fit_RS_err)
        else:
            N_pass_fitter_MC = t_reco_pid_2016.GetEntries(truthMatch+" && "+trigger+" && "+all_rectangular_cuts_MC+" && "+pass_mass_fit) + t_reco_pid_2017.GetEntries(truthMatch+" && "+trigger+" && "+all_rectangular_cuts_MC+" && "+pass_mass_fit) + t_reco_pid_2018.GetEntries(truthMatch+" && "+trigger+" && "+all_rectangular_cuts_MC+" && "+pass_mass_fit)
            N_pass_fitter_RS = t_rs_data_2016_noRectCuts.GetEntries(trigger+" && "+all_rectangular_cuts_data+" && "+pass_mass_fit) + t_rs_data_2017_noRectCuts.GetEntries(trigger+" && "+all_rectangular_cuts_data+" && "+pass_mass_fit) + t_rs_data_2018_noRectCuts.GetEntries(trigger+" && "+all_rectangular_cuts_data+" && "+pass_mass_fit)

            up_pass_fit_MC = ROOT.TEfficiency.Wilson(N_rect_cuts_MC_all, N_pass_fitter_MC, 0.68, True)
            down_pass_fit_MC = ROOT.TEfficiency.Wilson(N_rect_cuts_MC_all, N_pass_fitter_MC, 0.68, False)
            eps_pass_fit_MC_err = 0.5*(up_pass_fit_MC - down_pass_fit_MC)
            eps_pass_fit_MC = ufloat(N_pass_fitter_MC/N_rect_cuts_MC_all, eps_pass_fit_MC_err)

            up_pass_fit_RS = ROOT.TEfficiency.Wilson(N_rect_cuts_RS_data_all, N_pass_fitter_RS, 0.68, True)
            down_pass_fit_RS = ROOT.TEfficiency.Wilson(N_rect_cuts_RS_data_all, N_pass_fitter_RS, 0.68, False)
            eps_pass_fit_RS_err = 0.5*(up_pass_fit_RS - down_pass_fit_RS)
            eps_pass_fit_RS = ufloat(N_pass_fitter_RS/N_rect_cuts_RS_data_all, eps_pass_fit_RS_err)

        print("N_pass_fit MC:", N_pass_fitter_MC)
        print("N_pass_fit RS:", N_pass_fitter_RS)

        pass_fit_columns = ["Efficiency (\\%)"]
        if(species == 1):
            pass_fit_rows = ["$3\\pi 3\\pi$", "$3\\pi 3\\pi \\pi^0$", "$3\\pi 3\\pi 2\\pi^0$", "All MC", "RS data", "WS data"]
        else:
            pass_fit_rows  = ["MC", "Data"]

        pass_fit_table = pd.DataFrame(index=pass_fit_rows, columns=pass_fit_columns)
        if(species == 1):
            pass_fit_table.loc["$3\\pi 3\\pi$", "Efficiency (\\%)"] = f"${eps_pass_fit_3pi3pi.nominal_value*100:.2f} \\pm {eps_pass_fit_3pi3pi.std_dev*100:.2f}$"
            pass_fit_table.loc["$3\\pi 3\\pi \\pi^0$"] = f"${eps_pass_fit_3pi3pipi0.nominal_value*100:.2f} \\pm {eps_pass_fit_3pi3pipi0.std_dev*100:.2f}$"
            pass_fit_table.loc["$3\\pi 3\\pi 2\\pi^0$"] = f"${eps_pass_fit_3pi3pi2pi0.nominal_value*100:.2f} \\pm {eps_pass_fit_3pi3pi2pi0.std_dev*100:.2f}$"
            pass_fit_table.loc["All MC"] = f"${eps_pass_fit_MC.nominal_value*100:.2f} \\pm {eps_pass_fit_MC.std_dev*100:.2f}$"
            pass_fit_table.loc["RS data"] = f"${eps_pass_fit_RS.nominal_value*100:.2f} \\pm {eps_pass_fit_RS.std_dev*100:.2f}$"
            pass_fit_table.loc["WS data"] = f"${eps_pass_fit_WS.nominal_value*100:.2f} \\pm {eps_pass_fit_WS.std_dev*100:.2f}$"
        else:
            pass_fit_table.loc["MC", "Efficiency (\\%)"] = f"${eps_pass_fit_MC.nominal_value*100:.2f} \\pm {eps_pass_fit_MC.std_dev*100:.2f}$"
            pass_fit_table.loc["Data", "Efficiency (\\%)"] = f"${eps_pass_fit_RS.nominal_value*100:.2f} \\pm {eps_pass_fit_RS.std_dev*100:.2f}$"

    with open(f'/panfs/felician/B2Ktautau/workflow/selections_efficiency_tables/Species_{species}/pass_mass_fit_efficiency_table.tex', 'w') as fout_pass_fit:
        fout_pass_fit.write(pass_fit_table.to_latex())


    # Fit region
    print("Fit region efficiency")
    if(species == 100):
        N_fit_region_BuDDKp = t_reco_BuDDKp_2016.GetEntries(trigger+" && "+all_rectangular_cuts+" && "+pass_mass_fit+" && "+fit_range) + t_reco_BuDDKp_2017.GetEntries(trigger+" && "+all_rectangular_cuts+" && "+pass_mass_fit+" && "+fit_range) + t_reco_BuDDKp_2018.GetEntries(trigger+" && "+all_rectangular_cuts+" && "+pass_mass_fit+" && "+fit_range) 
        N_fit_region_BdDDKp = t_reco_BdDDKp_2016.GetEntries(trigger+" && "+all_rectangular_cuts+" && "+pass_mass_fit+" && "+fit_range) + t_reco_BdDDKp_2017.GetEntries(trigger+" && "+all_rectangular_cuts+" && "+pass_mass_fit+" && "+fit_range) + t_reco_BdDDKp_2018.GetEntries(trigger+" && "+all_rectangular_cuts+" && "+pass_mass_fit+" && "+fit_range)
        N_fit_region_BsDDKp = t_reco_BsDDKp_2016.GetEntries(trigger+" && "+all_rectangular_cuts+" && "+pass_mass_fit+" && "+fit_range) + t_reco_BsDDKp_2017.GetEntries(trigger+" && "+all_rectangular_cuts+" && "+pass_mass_fit+" && "+fit_range) + t_reco_BsDDKp_2018.GetEntries(trigger+" && "+all_rectangular_cuts+" && "+pass_mass_fit+" && "+fit_range)
        N_fit_region_BuDDK0 = t_reco_BuDDK0_2016.GetEntries(trigger+" && "+all_rectangular_cuts+" && "+pass_mass_fit+" && "+fit_range) + t_reco_BuDDK0_2017.GetEntries(trigger+" && "+all_rectangular_cuts+" && "+pass_mass_fit+" && "+fit_range) + t_reco_BuDDK0_2018.GetEntries(trigger+" && "+all_rectangular_cuts+" && "+pass_mass_fit+" && "+fit_range)
        N_fit_region_BuDD = t_reco_BuDD_2016.GetEntries(trigger+" && "+all_rectangular_cuts+" && "+pass_mass_fit+" && "+fit_range) + t_reco_BuDD_2017.GetEntries(trigger+" && "+all_rectangular_cuts+" && "+pass_mass_fit+" && "+fit_range) + t_reco_BuDD_2018.GetEntries(trigger+" && "+all_rectangular_cuts+" && "+pass_mass_fit+" && "+fit_range)

        up_fit_region_BuDDKp = ROOT.TEfficiency.Wilson(N_pass_fit_BuDDKp, N_fit_region_BuDDKp, 0.68, True)
        down_fit_region_BuDDKp = ROOT.TEfficiency.Wilson(N_pass_fit_BuDDKp, N_fit_region_BuDDKp, 0.68, False)
        eps_fit_region_BuDDKp_err = 0.5*(up_fit_region_BuDDKp - down_fit_region_BuDDKp)
        eps_fit_region_BuDDKp = ufloat(N_fit_region_BuDDKp/N_pass_fit_BuDDKp, eps_fit_region_BuDDKp_err)

        up_fit_region_BdDDKp = ROOT.TEfficiency.Wilson(N_pass_fit_BdDDKp, N_fit_region_BdDDKp, 0.68, True)
        down_fit_region_BdDDKp = ROOT.TEfficiency.Wilson(N_pass_fit_BdDDKp, N_fit_region_BdDDKp, 0.68, False)
        eps_fit_region_BdDDKp_err = 0.5*(up_fit_region_BdDDKp - down_fit_region_BdDDKp)
        eps_fit_region_BdDDKp = ufloat(N_fit_region_BdDDKp/N_pass_fit_BdDDKp, eps_fit_region_BdDDKp_err)

        up_fit_region_BsDDKp = ROOT.TEfficiency.Wilson(N_pass_fit_BsDDKp, N_fit_region_BsDDKp, 0.68, True)
        down_fit_region_BsDDKp = ROOT.TEfficiency.Wilson(N_pass_fit_BsDDKp, N_fit_region_BsDDKp, 0.68, False)
        eps_fit_region_BsDDKp_err = 0.5*(up_fit_region_BsDDKp - down_fit_region_BsDDKp)
        eps_fit_region_BsDDKp = ufloat(N_fit_region_BsDDKp/N_pass_fit_BsDDKp, eps_fit_region_BsDDKp_err)

        up_fit_region_BuDDK0 = ROOT.TEfficiency.Wilson(N_pass_fit_BuDDK0, N_fit_region_BuDDK0, 0.68, True)
        down_fit_region_BuDDK0 = ROOT.TEfficiency.Wilson(N_pass_fit_BuDDK0, N_fit_region_BuDDK0, 0.68, False)
        eps_fit_region_BuDDK0_err = 0.5*(up_fit_region_BuDDK0 - down_fit_region_BuDDK0)
        eps_fit_region_BuDDK0 = ufloat(N_fit_region_BuDDK0/N_pass_fit_BuDDK0, eps_fit_region_BuDDK0_err)

        up_fit_region_BuDD = ROOT.TEfficiency.Wilson(N_pass_fit_BuDD, N_fit_region_BuDD, 0.68, True)
        down_fit_region_BuDD = ROOT.TEfficiency.Wilson(N_pass_fit_BuDD, N_fit_region_BuDD, 0.68, False)
        eps_fit_region_BuDD_err = 0.5*(up_fit_region_BuDD - down_fit_region_BuDD)
        eps_fit_region_BuDD = ufloat(N_fit_region_BuDD/N_pass_fit_BuDD, eps_fit_region_BuDD_err)

        eps_fit_region = [eps_fit_region_BuDDKp, eps_fit_region_BdDDKp, eps_fit_region_BsDDKp, eps_fit_region_BuDDK0, eps_fit_region_BuDD]

        fit_region_columns = ["Fit region efficiency (\\%)"]
        fit_region_rows = ["$B^+ \\to D D K^+$", "$B^0 \\to D D K^+$", "$B^0_s \\to D D K^+$", "$B^+ \\to D D K^0$", "$B^+ \\to D D$"]

        fit_region_table = pd.DataFrame(index=fit_region_rows, columns=fit_region_columns)
        fit_region_table.loc["$B^+ \\to D D K^+$", "Fit region efficiency (\\%)"] = f"${eps_fit_region[0].nominal_value*100:.2f} \\pm {eps_fit_region[0].std_dev*100:.2f}$"
        fit_region_table.loc["$B^0 \\to D D K^+$", "Fit region efficiency (\\%)"] = f"${eps_fit_region[1].nominal_value*100:.2f} \\pm {eps_fit_region[1].std_dev*100:.2f}$"
        fit_region_table.loc["$B^0_s \\to D D K^+$", "Fit region efficiency (\\%)"] = f"${eps_fit_region[2].nominal_value*100:.2f} \\pm {eps_fit_region[2].std_dev*100:.2f}$"
        fit_region_table.loc["$B^+ \\to D D K^0$", "Fit region efficiency (\\%)"] = f"${eps_fit_region[3].nominal_value*100:.2f} \\pm {eps_fit_region[3].std_dev*100:.2f}$"
        fit_region_table.loc["$B^+ \\to D D$", "Fit region efficiency (\\%)"] = f"${eps_fit_region[4].nominal_value*100:.2f} \\pm {eps_fit_region[4].std_dev*100:.2f}$"

    else:
        if(species == 1):
            N_fit_region_3pi3pi = t_reco_mc_2016.GetEntries(truthMatch+" && "+trigger+" && "+all_rectangular_cuts_MC+" && "+pass_mass_fit+" && "+fit_range+" && (component == 0)") + t_reco_mc_2017.GetEntries(truthMatch+" && "+trigger+" && "+all_rectangular_cuts_MC+" && "+pass_mass_fit+" && "+fit_range+" && (component == 0)") + t_reco_mc_2018.GetEntries(truthMatch+" && "+trigger+" && "+all_rectangular_cuts_MC+" && "+pass_mass_fit+" && "+fit_range+" && (component == 0)")
            N_fit_region_3pi3pipi0 = t_reco_mc_2016.GetEntries(truthMatch+" && "+trigger+" && "+all_rectangular_cuts_MC+" && "+pass_mass_fit+" && "+fit_range+" && (component == 1)") + t_reco_mc_2017.GetEntries(truthMatch+" && "+trigger+" && "+all_rectangular_cuts_MC+" && "+pass_mass_fit+" && "+fit_range+" && (component == 1)") + t_reco_mc_2018.GetEntries(truthMatch+" && "+trigger+" && "+all_rectangular_cuts_MC+" && "+pass_mass_fit+" && "+fit_range+" && (component == 1)")
            N_fit_region_3pi3pi2pi0 = t_reco_mc_2016.GetEntries(truthMatch+" && "+trigger+" && "+all_rectangular_cuts_MC+" && "+pass_mass_fit+" && "+fit_range+" && (component == 2)") + t_reco_mc_2017.GetEntries(truthMatch+" && "+trigger+" && "+all_rectangular_cuts_MC+" && "+pass_mass_fit+" && "+fit_range+" && (component == 2)") + t_reco_mc_2018.GetEntries(truthMatch+" && "+trigger+" && "+all_rectangular_cuts_MC+" && "+pass_mass_fit+" && "+fit_range+" && (component == 2)")

            up_fit_region_3pi3pi = ROOT.TEfficiency.Wilson(N_pass_fitter_3pi3pi, N_fit_region_3pi3pi, 0.68, True)
            down_fit_region_3pi3pi = ROOT.TEfficiency.Wilson(N_pass_fitter_3pi3pi, N_fit_region_3pi3pi, 0.68, False)
            eps_fit_region_3pi3pi_err = 0.5*(up_fit_region_3pi3pi - down_fit_region_3pi3pi)
            eps_fit_region_3pi3pi = ufloat(N_fit_region_3pi3pi/N_pass_fitter_3pi3pi, eps_fit_region_3pi3pi_err)

            up_fit_region_3pi3pipi0 = ROOT.TEfficiency.Wilson(N_pass_fitter_3pi3pipi0, N_fit_region_3pi3pipi0, 0.68, True)
            down_fit_region_3pi3pipi0 = ROOT.TEfficiency.Wilson(N_pass_fitter_3pi3pipi0, N_fit_region_3pi3pipi0, 0.68, False)
            eps_fit_region_3pi3pipi0_err = 0.5*(up_fit_region_3pi3pipi0 - down_fit_region_3pi3pipi0)
            eps_fit_region_3pi3pipi0 = ufloat(N_fit_region_3pi3pipi0/N_pass_fitter_3pi3pipi0, eps_fit_region_3pi3pipi0_err)

            up_fit_region_3pi3pi2pi0 = ROOT.TEfficiency.Wilson(N_pass_fitter_3pi3pi2pi0, N_fit_region_3pi3pi2pi0, 0.68, True)
            down_fit_region_3pi3pi2pi0 = ROOT.TEfficiency.Wilson(N_pass_fitter_3pi3pi2pi0, N_fit_region_3pi3pi2pi0, 0.68, False)
            eps_fit_region_3pi3pi2pi0_err = 0.5*(up_fit_region_3pi3pi2pi0 - down_fit_region_3pi3pi2pi0)
            eps_fit_region_3pi3pi2pi0 = ufloat(N_fit_region_3pi3pi2pi0/N_pass_fitter_3pi3pi2pi0, eps_fit_region_3pi3pi2pi0_err)

            N_fit_region_WS = t_reco_ws_2016.GetEntries(trigger+" && "+all_rectangular_cuts_data+" && "+pass_mass_fit+" && "+fit_range) + t_reco_ws_2017.GetEntries(trigger+" && "+all_rectangular_cuts_data+" && "+pass_mass_fit+" && "+fit_range) + t_reco_ws_2018.GetEntries(trigger+" && "+all_rectangular_cuts_data+" && "+pass_mass_fit+" && "+fit_range)

            up_fit_region_WS = ROOT.TEfficiency.Wilson(N_pass_fitter_WS, N_fit_region_WS, 0.68, True)
            down_fit_region_WS = ROOT.TEfficiency.Wilson(N_pass_fitter_WS, N_fit_region_WS, 0.68, False)
            eps_fit_region_WS_err = 0.5*(up_fit_region_WS - down_fit_region_WS)
            eps_fit_region_WS = ufloat(N_fit_region_WS/N_pass_fitter_WS, eps_fit_region_WS_err)
        

        if(species == 1):
            N_fit_region_MC = t_reco_mc_2016.GetEntries(truthMatch+" && "+trigger+" && "+all_rectangular_cuts_MC+" && "+pass_mass_fit+" && "+fit_range) + t_reco_mc_2017.GetEntries(truthMatch+" && "+trigger+" && "+all_rectangular_cuts_MC+" && "+pass_mass_fit+" && "+fit_range) + t_reco_mc_2018.GetEntries(truthMatch+" && "+trigger+" && "+all_rectangular_cuts_MC+" && "+pass_mass_fit+" && "+fit_range)
            N_fit_region_RS = t_reco_rs_2016.GetEntries(trigger+" && "+all_rectangular_cuts_data+" && "+pass_mass_fit+" && "+fit_range) + t_reco_rs_2017.GetEntries(trigger+" && "+all_rectangular_cuts_data+" && "+pass_mass_fit+" && "+fit_range) + t_reco_rs_2018.GetEntries(trigger+" && "+all_rectangular_cuts_data+" && "+pass_mass_fit+" && "+fit_range)
        else:
            N_fit_region_MC = t_reco_pid_2016.GetEntries(truthMatch+" && "+trigger+" && "+all_rectangular_cuts_MC+" && "+pass_mass_fit+" && "+fit_range) + t_reco_pid_2017.GetEntries(truthMatch+" && "+trigger+" && "+all_rectangular_cuts_MC+" && "+pass_mass_fit+" && "+fit_range) + t_reco_pid_2018.GetEntries(truthMatch+" && "+trigger+" && "+all_rectangular_cuts_MC+" && "+pass_mass_fit+" && "+fit_range)
            N_fit_region_RS = t_rs_data_2016_noRectCuts.GetEntries(trigger+" && "+all_rectangular_cuts_data+" && "+pass_mass_fit+" && "+fit_range) + t_rs_data_2017_noRectCuts.GetEntries(trigger+" && "+all_rectangular_cuts_data+" && "+pass_mass_fit+" && "+fit_range) + t_rs_data_2018_noRectCuts.GetEntries(trigger+" && "+all_rectangular_cuts_data+" && "+pass_mass_fit+" && "+fit_range)

        up_fit_region_MC = ROOT.TEfficiency.Wilson(N_pass_fitter_MC, N_fit_region_MC, 0.68, True)
        down_fit_region_MC = ROOT.TEfficiency.Wilson(N_pass_fitter_MC, N_fit_region_MC, 0.68, False)
        eps_fit_region_MC_err = 0.5*(up_fit_region_MC - down_fit_region_MC)
        eps_fit_region_MC = ufloat(N_fit_region_MC/N_pass_fitter_MC, eps_fit_region_MC_err)

        up_fit_region_RS = ROOT.TEfficiency.Wilson(N_pass_fitter_RS, N_fit_region_RS, 0.68, True)
        down_fit_region_RS = ROOT.TEfficiency.Wilson(N_pass_fitter_RS, N_fit_region_RS, 0.68, False)
        eps_fit_region_RS_err = 0.5*(up_fit_region_RS - down_fit_region_RS)
        eps_fit_region_RS = ufloat(N_fit_region_RS/N_pass_fitter_RS, eps_fit_region_RS_err)

        fit_region_columns = ["Efficiency (\\%)"]
        if(species == 1):
            fit_region_rows = ["$3\\pi 3\\pi$", "$3\\pi 3\\pi \\pi^0$", "$3\\pi 3\\pi 2\\pi^0$", "All MC", "RS data", "WS data"]
        else:
            fit_region_rows  = ["MC", "Data"]

        fit_region_table = pd.DataFrame(index=fit_region_rows, columns=fit_region_columns)
        if(species == 1):
            fit_region_table.loc["$3\\pi 3\\pi$", "Efficiency (\\%)"] = f"${eps_fit_region_3pi3pi.nominal_value*100:.2f} \\pm {eps_fit_region_3pi3pi.std_dev*100:.2f}$"
            fit_region_table.loc["$3\\pi 3\\pi \\pi^0$"] = f"${eps_fit_region_3pi3pipi0.nominal_value*100:.2f} \\pm {eps_fit_region_3pi3pipi0.std_dev*100:.2f}$"
            fit_region_table.loc["$3\\pi 3\\pi 2\\pi^0$"] = f"${eps_fit_region_3pi3pi2pi0.nominal_value*100:.2f} \\pm {eps_fit_region_3pi3pi2pi0.std_dev*100:.2f}$"
            fit_region_table.loc["All MC"] = f"${eps_fit_region_MC.nominal_value*100:.2f} \\pm {eps_fit_region_MC.std_dev*100:.2f}$"
            fit_region_table.loc["RS data"] = f"${eps_fit_region_RS.nominal_value*100:.2f} \\pm {eps_fit_region_RS.std_dev*100:.2f}$"
            fit_region_table.loc["WS data"] = f"${eps_fit_region_WS.nominal_value*100:.2f} \\pm {eps_fit_region_WS.std_dev*100:.2f}$"
        else:
            fit_region_table.loc["MC", "Efficiency (\\%)"] = f"${eps_fit_region_MC.nominal_value*100:.2f} \\pm {eps_fit_region_MC.std_dev*100:.2f}$"
            fit_region_table.loc["Data", "Efficiency (\\%)"] = f"${eps_fit_region_RS.nominal_value*100:.2f} \\pm {eps_fit_region_RS.std_dev*100:.2f}$"

    with open(f'/panfs/felician/B2Ktautau/workflow/selections_efficiency_tables/Species_{species}/fit_region_efficiency_table.tex', 'w') as fout_fit_region:
        fout_fit_region.write(fit_region_table.to_latex())


    # Mass vetoes
    print("Mass vetoes efficiency")
    eps_mass_vetoes_3pi3pi = []
    eps_mass_vetoes_3pi3pipi0 = []
    eps_mass_vetoes_3pi3pi2pi0 = []
    eps_mass_vetoes_MC = []
    eps_mass_vetoes_RS = []
    eps_mass_vetoes_WS = [] 

    if(species == 100):
        N_mass_vetoes_BuDDKp = t_reco_BuDDKp_2016.GetEntries(trigger+" && "+all_rectangular_cuts+" && "+pass_mass_fit+" && "+fit_range+" && "+mass_vetoes) + t_reco_BuDDKp_2017.GetEntries(trigger+" && "+all_rectangular_cuts+" && "+pass_mass_fit+" && "+fit_range+" && "+mass_vetoes) + t_reco_BuDDKp_2018.GetEntries(trigger+" && "+all_rectangular_cuts+" && "+pass_mass_fit+" && "+fit_range+" && "+mass_vetoes) 
        N_mass_vetoes_BdDDKp = t_reco_BdDDKp_2016.GetEntries(trigger+" && "+all_rectangular_cuts+" && "+pass_mass_fit+" && "+fit_range+" && "+mass_vetoes) + t_reco_BdDDKp_2017.GetEntries(trigger+" && "+all_rectangular_cuts+" && "+pass_mass_fit+" && "+fit_range+" && "+mass_vetoes) + t_reco_BdDDKp_2018.GetEntries(trigger+" && "+all_rectangular_cuts+" && "+pass_mass_fit+" && "+fit_range+" && "+mass_vetoes)
        N_mass_vetoes_BsDDKp = t_reco_BsDDKp_2016.GetEntries(trigger+" && "+all_rectangular_cuts+" && "+pass_mass_fit+" && "+fit_range+" && "+mass_vetoes) + t_reco_BsDDKp_2017.GetEntries(trigger+" && "+all_rectangular_cuts+" && "+pass_mass_fit+" && "+fit_range+" && "+mass_vetoes) + t_reco_BsDDKp_2018.GetEntries(trigger+" && "+all_rectangular_cuts+" && "+pass_mass_fit+" && "+fit_range+" && "+mass_vetoes)
        N_mass_vetoes_BuDDK0 = t_reco_BuDDK0_2016.GetEntries(trigger+" && "+all_rectangular_cuts+" && "+pass_mass_fit+" && "+fit_range+" && "+mass_vetoes) + t_reco_BuDDK0_2017.GetEntries(trigger+" && "+all_rectangular_cuts+" && "+pass_mass_fit+" && "+fit_range+" && "+mass_vetoes) + t_reco_BuDDK0_2018.GetEntries(trigger+" && "+all_rectangular_cuts+" && "+pass_mass_fit+" && "+fit_range+" && "+mass_vetoes)
        N_mass_vetoes_BuDD = t_reco_BuDD_2016.GetEntries(trigger+" && "+all_rectangular_cuts+" && "+pass_mass_fit+" && "+fit_range+" && "+mass_vetoes) + t_reco_BuDD_2017.GetEntries(trigger+" && "+all_rectangular_cuts+" && "+pass_mass_fit+" && "+fit_range+" && "+mass_vetoes) + t_reco_BuDD_2018.GetEntries(trigger+" && "+all_rectangular_cuts+" && "+pass_mass_fit+" && "+fit_range+" && "+mass_vetoes)

        up_mass_vetoes_BuDDKp = ROOT.TEfficiency.Wilson(N_fit_region_BuDDKp, N_mass_vetoes_BuDDKp, 0.68, True)
        down_mass_vetoes_BuDDKp = ROOT.TEfficiency.Wilson(N_fit_region_BuDDKp, N_mass_vetoes_BuDDKp, 0.68, False)
        eps_mass_vetoes_BuDDKp_err = 0.5*(up_mass_vetoes_BuDDKp - down_mass_vetoes_BuDDKp)
        eps_mass_vetoes_BuDDKp = ufloat(N_mass_vetoes_BuDDKp/N_fit_region_BuDDKp, eps_mass_vetoes_BuDDKp_err)

        up_mass_vetoes_BdDDKp = ROOT.TEfficiency.Wilson(N_fit_region_BdDDKp, N_mass_vetoes_BdDDKp, 0.68, True)
        down_mass_vetoes_BdDDKp = ROOT.TEfficiency.Wilson(N_fit_region_BdDDKp, N_mass_vetoes_BdDDKp, 0.68, False)
        eps_mass_vetoes_BdDDKp_err = 0.5*(up_mass_vetoes_BdDDKp - down_mass_vetoes_BdDDKp)
        eps_mass_vetoes_BdDDKp = ufloat(N_mass_vetoes_BdDDKp/N_fit_region_BdDDKp, eps_mass_vetoes_BdDDKp_err)

        up_mass_vetoes_BsDDKp = ROOT.TEfficiency.Wilson(N_fit_region_BsDDKp, N_mass_vetoes_BsDDKp, 0.68, True)
        down_mass_vetoes_BsDDKp = ROOT.TEfficiency.Wilson(N_fit_region_BsDDKp, N_mass_vetoes_BsDDKp, 0.68, False)
        eps_mass_vetoes_BsDDKp_err = 0.5*(up_mass_vetoes_BsDDKp - down_mass_vetoes_BsDDKp)
        eps_mass_vetoes_BsDDKp = ufloat(N_mass_vetoes_BsDDKp/N_fit_region_BsDDKp, eps_mass_vetoes_BsDDKp_err)

        up_mass_vetoes_BuDDK0 = ROOT.TEfficiency.Wilson(N_fit_region_BuDDK0, N_mass_vetoes_BuDDK0, 0.68, True)
        down_mass_vetoes_BuDDK0 = ROOT.TEfficiency.Wilson(N_fit_region_BuDDK0, N_mass_vetoes_BuDDK0, 0.68, False)
        eps_mass_vetoes_BuDDK0_err = 0.5*(up_mass_vetoes_BuDDK0 - down_mass_vetoes_BuDDK0)
        eps_mass_vetoes_BuDDK0 = ufloat(N_mass_vetoes_BuDDK0/N_fit_region_BuDDK0, eps_mass_vetoes_BuDDK0_err)

        up_mass_vetoes_BuDD = ROOT.TEfficiency.Wilson(N_fit_region_BuDD, N_mass_vetoes_BuDD, 0.68, True)
        down_mass_vetoes_BuDD = ROOT.TEfficiency.Wilson(N_fit_region_BuDD, N_mass_vetoes_BuDD, 0.68, False)
        eps_mass_vetoes_BuDD_err = 0.5*(up_mass_vetoes_BuDD - down_mass_vetoes_BuDD)
        eps_mass_vetoes_BuDD = ufloat(N_mass_vetoes_BuDD/N_fit_region_BuDD, eps_mass_vetoes_BuDD_err)

        eps_mass_vetoes = [eps_mass_vetoes_BuDDKp, eps_mass_vetoes_BdDDKp, eps_mass_vetoes_BsDDKp, eps_mass_vetoes_BuDDK0, eps_mass_vetoes_BuDD]

        mass_vetoes_columns = ["Mass vetoes efficiency (\\%)"]
        mass_vetoes_rows = ["$B^+ \\to D D K^+$", "$B^0 \\to D D K^+$", "$B^0_s \\to D D K^+$", "$B^+ \\to D D K^0$", "$B^+ \\to D D$"]

        mass_vetoes_table = pd.DataFrame(index=mass_vetoes_rows, columns=mass_vetoes_columns)
        mass_vetoes_table.loc["$B^+ \\to D D K^+$", "Mass vetoes efficiency (\\%)"] = f"${eps_mass_vetoes[0].nominal_value*100:.2f} \\pm {eps_mass_vetoes[0].std_dev*100:.2f}$"
        mass_vetoes_table.loc["$B^0 \\to D D K^+$", "Mass vetoes efficiency (\\%)"] = f"${eps_mass_vetoes[1].nominal_value*100:.2f} \\pm {eps_mass_vetoes[1].std_dev*100:.2f}$"
        mass_vetoes_table.loc["$B^0_s \\to D D K^+$", "Mass vetoes efficiency (\\%)"] = f"${eps_mass_vetoes[2].nominal_value*100:.2f} \\pm {eps_mass_vetoes[2].std_dev*100:.2f}$"
        mass_vetoes_table.loc["$B^+ \\to D D K^0$", "Mass vetoes efficiency (\\%)"] = f"${eps_mass_vetoes[3].nominal_value*100:.2f} \\pm {eps_mass_vetoes[3].std_dev*100:.2f}$"
        mass_vetoes_table.loc["$B^+ \\to D D$", "Mass vetoes efficiency (\\%)"] = f"${eps_mass_vetoes[4].nominal_value*100:.2f} \\pm {eps_mass_vetoes[4].std_dev*100:.2f}$"

    elif(species == 1):
        for i in range(len(mass_vetoes_cuts)):
            N_mass_vetoes_3pi3pi = t_reco_mc_2016.GetEntries(truthMatch+" && "+trigger+" && "+all_rectangular_cuts_MC+" && "+pass_mass_fit+" && "+fit_range+" && (component == 0)"+" && "+mass_vetoes_cuts[i]) + t_reco_mc_2017.GetEntries(truthMatch+" && "+trigger+" && "+all_rectangular_cuts_MC+" && "+pass_mass_fit+" && "+fit_range+" && (component == 0)"+" && "+mass_vetoes_cuts[i]) + t_reco_mc_2018.GetEntries(truthMatch+" && "+trigger+" && "+all_rectangular_cuts_MC+" && "+pass_mass_fit+" && "+fit_range+" && (component == 0)"+" && "+mass_vetoes_cuts[i])
            N_mass_vetoes_3pi3pipi0 = t_reco_mc_2016.GetEntries(truthMatch+" && "+trigger+" && "+all_rectangular_cuts_MC+" && "+pass_mass_fit+" && "+fit_range+" && (component == 1)"+" && "+mass_vetoes_cuts[i]) + t_reco_mc_2017.GetEntries(truthMatch+" && "+trigger+" && "+all_rectangular_cuts_MC+" && "+pass_mass_fit+" && "+fit_range+" && (component == 1)"+" && "+mass_vetoes_cuts[i]) + t_reco_mc_2018.GetEntries(truthMatch+" && "+trigger+" && "+all_rectangular_cuts_MC+" && "+pass_mass_fit+" && "+fit_range+" && (component == 1)"+" && "+mass_vetoes_cuts[i])
            N_mass_vetoes_3pi3pi2pi0 = t_reco_mc_2016.GetEntries(truthMatch+" && "+trigger+" && "+all_rectangular_cuts_MC+" && "+pass_mass_fit+" && "+fit_range+" && (component == 2)"+" && "+mass_vetoes_cuts[i]) + t_reco_mc_2017.GetEntries(truthMatch+" && "+trigger+" && "+all_rectangular_cuts_MC+" && "+pass_mass_fit+" && "+fit_range+" && (component == 2)"+" && "+mass_vetoes_cuts[i]) + t_reco_mc_2018.GetEntries(truthMatch+" && "+trigger+" && "+all_rectangular_cuts_MC+" && "+pass_mass_fit+" && "+fit_range+" && (component == 2)"+" && "+mass_vetoes_cuts[i])
            N_mass_vetoes_MC = t_reco_mc_2016.GetEntries(truthMatch+" && "+trigger+" && "+all_rectangular_cuts_MC+" && "+pass_mass_fit+" && "+fit_range+" && "+mass_vetoes_cuts[i]) + t_reco_mc_2017.GetEntries(truthMatch+" && "+trigger+" && "+all_rectangular_cuts_MC+" && "+pass_mass_fit+" && "+fit_range+" && "+mass_vetoes_cuts[i]) + t_reco_mc_2018.GetEntries(truthMatch+" && "+trigger+" && "+all_rectangular_cuts_MC+" && "+pass_mass_fit+" && "+fit_range+" && "+mass_vetoes_cuts[i])
            N_mass_vetoes_RS = t_reco_rs_2016.GetEntries(trigger+" && "+all_rectangular_cuts_data+" && "+pass_mass_fit+" && "+fit_range+" && "+mass_vetoes_cuts[i]) + t_reco_rs_2017.GetEntries(trigger+" && "+all_rectangular_cuts_data+" && "+pass_mass_fit+" && "+fit_range+" && "+mass_vetoes_cuts[i]) + t_reco_rs_2018.GetEntries(trigger+" && "+all_rectangular_cuts_data+" && "+pass_mass_fit+" && "+fit_range+" && "+mass_vetoes_cuts[i])
            N_mass_vetoes_WS = t_reco_ws_2016.GetEntries(trigger+" && "+all_rectangular_cuts_data+" && "+pass_mass_fit+" && "+fit_range+" && "+mass_vetoes_cuts[i]) + t_reco_ws_2017.GetEntries(trigger+" && "+all_rectangular_cuts_data+" && "+pass_mass_fit+" && "+fit_range+" && "+mass_vetoes_cuts[i]) + t_reco_ws_2018.GetEntries(trigger+" && "+all_rectangular_cuts_data+" && "+pass_mass_fit+" && "+fit_range+" && "+mass_vetoes_cuts[i])

            up_mass_vetoes_3pi3pi = ROOT.TEfficiency.Wilson(N_fit_region_3pi3pi, N_mass_vetoes_3pi3pi, 0.68, True)
            down_mass_vetoes_3pi3pi = ROOT.TEfficiency.Wilson(N_fit_region_3pi3pi, N_mass_vetoes_3pi3pi, 0.68, False)
            eps_mass_vetoes_3pi3pi_err = 0.5*(up_mass_vetoes_3pi3pi - down_mass_vetoes_3pi3pi)
            eps_mass_vetoes_3pi3pi.append( ufloat(N_mass_vetoes_3pi3pi/N_fit_region_3pi3pi, eps_mass_vetoes_3pi3pi_err) )

            up_mass_vetoes_3pi3pipi0 = ROOT.TEfficiency.Wilson(N_fit_region_3pi3pipi0, N_mass_vetoes_3pi3pipi0, 0.68, True)
            down_mass_vetoes_3pi3pipi0 = ROOT.TEfficiency.Wilson(N_fit_region_3pi3pipi0, N_mass_vetoes_3pi3pipi0, 0.68, False)
            eps_mass_vetoes_3pi3pipi0_err = 0.5*(up_mass_vetoes_3pi3pipi0 - down_mass_vetoes_3pi3pipi0)
            eps_mass_vetoes_3pi3pipi0.append( ufloat(N_mass_vetoes_3pi3pipi0/N_fit_region_3pi3pipi0, eps_mass_vetoes_3pi3pipi0_err) )

            up_mass_vetoes_3pi3pi2pi0 = ROOT.TEfficiency.Wilson(N_fit_region_3pi3pi2pi0, N_mass_vetoes_3pi3pi2pi0, 0.68, True)
            down_mass_vetoes_3pi3pi2pi0 = ROOT.TEfficiency.Wilson(N_fit_region_3pi3pi2pi0, N_mass_vetoes_3pi3pi2pi0, 0.68, False)
            eps_mass_vetoes_3pi3pi2pi0_err = 0.5*(up_mass_vetoes_3pi3pi2pi0 - down_mass_vetoes_3pi3pi2pi0)
            eps_mass_vetoes_3pi3pi2pi0.append( ufloat(N_mass_vetoes_3pi3pi2pi0/N_fit_region_3pi3pi2pi0, eps_mass_vetoes_3pi3pi2pi0_err) )

            up_mass_vetoes_MC = ROOT.TEfficiency.Wilson(N_fit_region_MC, N_mass_vetoes_MC, 0.68, True)
            down_mass_vetoes_MC = ROOT.TEfficiency.Wilson(N_fit_region_MC, N_mass_vetoes_MC, 0.68, False)
            eps_mass_vetoes_MC_err = 0.5*(up_mass_vetoes_MC - down_mass_vetoes_MC)
            eps_mass_vetoes_MC.append( ufloat(N_mass_vetoes_MC/N_fit_region_MC, eps_mass_vetoes_MC_err) )

            up_mass_vetoes_RS = ROOT.TEfficiency.Wilson(N_fit_region_RS, N_mass_vetoes_RS, 0.68, True)
            down_mass_vetoes_RS = ROOT.TEfficiency.Wilson(N_fit_region_RS, N_mass_vetoes_RS, 0.68, False)
            eps_mass_vetoes_RS_err = 0.5*(up_mass_vetoes_RS - down_mass_vetoes_RS)
            eps_mass_vetoes_RS.append( ufloat(N_mass_vetoes_RS/N_fit_region_RS, eps_mass_vetoes_RS_err) )

            up_mass_vetoes_WS = ROOT.TEfficiency.Wilson(N_fit_region_WS, N_mass_vetoes_WS, 0.68, True)
            down_mass_vetoes_WS = ROOT.TEfficiency.Wilson(N_fit_region_WS, N_mass_vetoes_WS, 0.68, False)
            eps_mass_vetoes_WS_err = 0.5*(up_mass_vetoes_WS - down_mass_vetoes_WS)
            eps_mass_vetoes_WS.append( ufloat(N_mass_vetoes_WS/N_fit_region_WS, eps_mass_vetoes_WS_err) )

        N_all_mass_vetoes_3pi3pi = t_reco_mc_2016.GetEntries(truthMatch+" && "+trigger+" && "+all_rectangular_cuts_MC+" && "+pass_mass_fit+" && "+fit_range+" && (component == 0)"+" && "+mass_vetoes) + t_reco_mc_2017.GetEntries(truthMatch+" && "+trigger+" && "+all_rectangular_cuts_MC+" && "+pass_mass_fit+" && "+fit_range+" && (component == 0)"+" && "+mass_vetoes) + t_reco_mc_2018.GetEntries(truthMatch+" && "+trigger+" && "+all_rectangular_cuts_MC+" && "+pass_mass_fit+" && "+fit_range+" && (component == 0)"+" && "+mass_vetoes)
        N_all_mass_vetoes_3pi3pipi0 = t_reco_mc_2016.GetEntries(truthMatch+" && "+trigger+" && "+all_rectangular_cuts_MC+" && "+pass_mass_fit+" && "+fit_range+" && (component == 1)"+" && "+mass_vetoes) + t_reco_mc_2017.GetEntries(truthMatch+" && "+trigger+" && "+all_rectangular_cuts_MC+" && "+pass_mass_fit+" && "+fit_range+" && (component == 1)"+" && "+mass_vetoes) + t_reco_mc_2018.GetEntries(truthMatch+" && "+trigger+" && "+all_rectangular_cuts_MC+" && "+pass_mass_fit+" && "+fit_range+" && (component == 1)"+" && "+mass_vetoes)
        N_all_mass_vetoes_3pi3pi2pi0 = t_reco_mc_2016.GetEntries(truthMatch+" && "+trigger+" && "+all_rectangular_cuts_MC+" && "+pass_mass_fit+" && "+fit_range+" && (component == 2)"+" && "+mass_vetoes) + t_reco_mc_2017.GetEntries(truthMatch+" && "+trigger+" && "+all_rectangular_cuts_MC+" && "+pass_mass_fit+" && "+fit_range+" && (component == 2)"+" && "+mass_vetoes) + t_reco_mc_2018.GetEntries(truthMatch+" && "+trigger+" && "+all_rectangular_cuts_MC+" && "+pass_mass_fit+" && "+fit_range+" && (component == 2)"+" && "+mass_vetoes)
        N_all_mass_vetoes_MC = t_reco_mc_2016.GetEntries(truthMatch+" && "+trigger+" && "+all_rectangular_cuts_MC+" && "+pass_mass_fit+" && "+fit_range+" && "+mass_vetoes) + t_reco_mc_2017.GetEntries(truthMatch+" && "+trigger+" && "+all_rectangular_cuts_MC+" && "+pass_mass_fit+" && "+fit_range+" && "+mass_vetoes) + t_reco_mc_2018.GetEntries(truthMatch+" && "+trigger+" && "+all_rectangular_cuts_MC+" && "+pass_mass_fit+" && "+fit_range+" && "+mass_vetoes)
        N_all_mass_vetoes_RS = t_reco_rs_2016.GetEntries(trigger+" && "+all_rectangular_cuts_data+" && "+pass_mass_fit+" && "+fit_range+" && "+mass_vetoes) + t_reco_rs_2017.GetEntries(trigger+" && "+all_rectangular_cuts_data+" && "+pass_mass_fit+" && "+fit_range+" && "+mass_vetoes) + t_reco_rs_2018.GetEntries(trigger+" && "+all_rectangular_cuts_data+" && "+pass_mass_fit+" && "+fit_range+" && "+mass_vetoes)
        N_all_mass_vetoes_WS = t_reco_ws_2016.GetEntries(trigger+" && "+all_rectangular_cuts_data+" && "+pass_mass_fit+" && "+fit_range+" && "+mass_vetoes) + t_reco_ws_2017.GetEntries(trigger+" && "+all_rectangular_cuts_data+" && "+pass_mass_fit+" && "+fit_range+" && "+mass_vetoes) + t_reco_ws_2018.GetEntries(trigger+" && "+all_rectangular_cuts_data+" && "+pass_mass_fit+" && "+fit_range+" && "+mass_vetoes)

        up_all_mass_vetoes_3pi3pi = ROOT.TEfficiency.Wilson(N_fit_region_3pi3pi, N_all_mass_vetoes_3pi3pi, 0.68, True)
        down_all_mass_vetoes_3pi3pi = ROOT.TEfficiency.Wilson(N_fit_region_3pi3pi, N_all_mass_vetoes_3pi3pi, 0.68, False)
        eps_all_mass_vetoes_3pi3pi_err = 0.5*(up_all_mass_vetoes_3pi3pi - down_all_mass_vetoes_3pi3pi)
        eps_all_mass_vetoes_3pi3pi =  ufloat(N_all_mass_vetoes_3pi3pi/N_fit_region_3pi3pi, eps_all_mass_vetoes_3pi3pi_err) 

        up_all_mass_vetoes_3pi3pipi0 = ROOT.TEfficiency.Wilson(N_fit_region_3pi3pipi0, N_all_mass_vetoes_3pi3pipi0, 0.68, True)
        down_all_mass_vetoes_3pi3pipi0 = ROOT.TEfficiency.Wilson(N_fit_region_3pi3pipi0, N_all_mass_vetoes_3pi3pipi0, 0.68, False)
        eps_all_mass_vetoes_3pi3pipi0_err = 0.5*(up_all_mass_vetoes_3pi3pipi0 - down_all_mass_vetoes_3pi3pipi0)
        eps_all_mass_vetoes_3pi3pipi0 = ufloat(N_all_mass_vetoes_3pi3pipi0/N_fit_region_3pi3pipi0, eps_all_mass_vetoes_3pi3pipi0_err) 

        up_all_mass_vetoes_3pi3pi2pi0 = ROOT.TEfficiency.Wilson(N_fit_region_3pi3pi2pi0, N_all_mass_vetoes_3pi3pi2pi0, 0.68, True)
        down_all_mass_vetoes_3pi3pi2pi0 = ROOT.TEfficiency.Wilson(N_fit_region_3pi3pi2pi0, N_all_mass_vetoes_3pi3pi2pi0, 0.68, False)
        eps_all_mass_vetoes_3pi3pi2pi0_err = 0.5*(up_all_mass_vetoes_3pi3pi2pi0 - down_all_mass_vetoes_3pi3pi2pi0)
        eps_all_mass_vetoes_3pi3pi2pi0 =  ufloat(N_all_mass_vetoes_3pi3pi2pi0/N_fit_region_3pi3pi2pi0, eps_all_mass_vetoes_3pi3pi2pi0_err) 

        up_all_mass_vetoes_MC = ROOT.TEfficiency.Wilson(N_fit_region_MC, N_all_mass_vetoes_MC, 0.68, True)
        down_all_mass_vetoes_MC = ROOT.TEfficiency.Wilson(N_fit_region_MC, N_all_mass_vetoes_MC, 0.68, False)
        eps_all_mass_vetoes_MC_err = 0.5*(up_all_mass_vetoes_MC - down_all_mass_vetoes_MC)
        eps_all_mass_vetoes_MC =  ufloat(N_all_mass_vetoes_MC/N_fit_region_MC, eps_all_mass_vetoes_MC_err) 

        up_all_mass_vetoes_RS = ROOT.TEfficiency.Wilson(N_fit_region_RS, N_all_mass_vetoes_RS, 0.68, True)
        down_all_mass_vetoes_RS = ROOT.TEfficiency.Wilson(N_fit_region_RS, N_all_mass_vetoes_RS, 0.68, False)
        eps_all_mass_vetoes_RS_err = 0.5*(up_all_mass_vetoes_RS - down_all_mass_vetoes_RS)
        eps_all_mass_vetoes_RS =  ufloat(N_all_mass_vetoes_RS/N_fit_region_RS, eps_all_mass_vetoes_RS_err) 

        up_all_mass_vetoes_WS = ROOT.TEfficiency.Wilson(N_fit_region_WS, N_all_mass_vetoes_WS, 0.68, True)
        down_all_mass_vetoes_WS = ROOT.TEfficiency.Wilson(N_fit_region_WS, N_all_mass_vetoes_WS, 0.68, False)
        eps_all_mass_vetoes_WS_err = 0.5*(up_all_mass_vetoes_WS - down_all_mass_vetoes_WS)
        eps_all_mass_vetoes_WS =  ufloat(N_all_mass_vetoes_WS/N_fit_region_WS, eps_all_mass_vetoes_WS_err) 

        mass_vetoes_columns = ["$3\\pi 3\\pi$", "$3\\pi 3\\pi \\pi^0$", "$3\\pi 3\\pi 2\\pi^0$", "All MC", "RS data", "WS data"]
        mass_vetoes_rows = ["$D^0$ (2 particles)", "$D^-$ (3 particles)", "$D^0$ (4 particles)", "$D^{*-}$ (5 particles)", "$IP > 0.25$\\,mm ($K^{*0}$)", "$IP > 0.15$\\,mm ($D^{0}$)",  "All"]

        mass_vetoes_table = pd.DataFrame(index=mass_vetoes_rows, columns=mass_vetoes_columns)
        mass_vetoes_table.loc["$D^0$ (2 particles)", "$3\\pi 3\\pi$"] = f"${eps_mass_vetoes_3pi3pi[0].nominal_value*100:.2f} \\pm {eps_mass_vetoes_3pi3pi[0].std_dev*100:.2f}$"
        mass_vetoes_table.loc["$D^-$ (3 particles)", "$3\\pi 3\\pi$"] = f"${eps_mass_vetoes_3pi3pi[1].nominal_value*100:.2f} \\pm {eps_mass_vetoes_3pi3pi[1].std_dev*100:.2f}$"
        mass_vetoes_table.loc["$D^0$ (4 particles)", "$3\\pi 3\\pi$"] = f"${eps_mass_vetoes_3pi3pi[2].nominal_value*100:.2f} \\pm {eps_mass_vetoes_3pi3pi[2].std_dev*100:.2f}$"
        mass_vetoes_table.loc["$D^{*-}$ (5 particles)", "$3\\pi 3\\pi$"] = f"${eps_mass_vetoes_3pi3pi[3].nominal_value*100:.2f} \\pm {eps_mass_vetoes_3pi3pi[3].std_dev*100:.2f}$"
        mass_vetoes_table.loc["$IP > 0.25$\\,mm ($K^{*0}$)", "$3\\pi 3\\pi$"] = f"${eps_mass_vetoes_3pi3pi[4].nominal_value*100:.2f} \\pm {eps_mass_vetoes_3pi3pi[4].std_dev*100:.2f}$"
        mass_vetoes_table.loc["$IP > 0.15$\\,mm ($D^{0}$)", "$3\\pi 3\\pi$"] = f"${eps_mass_vetoes_3pi3pi[5].nominal_value*100:.2f} \\pm {eps_mass_vetoes_3pi3pi[5].std_dev*100:.2f}$"

        mass_vetoes_table.loc["$D^0$ (2 particles)", "$3\\pi 3\\pi \\pi^0$"] = f"${eps_mass_vetoes_3pi3pipi0[0].nominal_value*100:.2f} \\pm {eps_mass_vetoes_3pi3pipi0[0].std_dev*100:.2f}$"
        mass_vetoes_table.loc["$D^-$ (3 particles)", "$3\\pi 3\\pi \\pi^0$"] = f"${eps_mass_vetoes_3pi3pipi0[1].nominal_value*100:.2f} \\pm {eps_mass_vetoes_3pi3pipi0[1].std_dev*100:.2f}$"
        mass_vetoes_table.loc["$D^0$ (4 particles)", "$3\\pi 3\\pi \\pi^0$"] = f"${eps_mass_vetoes_3pi3pipi0[2].nominal_value*100:.2f} \\pm {eps_mass_vetoes_3pi3pipi0[2].std_dev*100:.2f}$"
        mass_vetoes_table.loc["$D^{*-}$ (5 particles)", "$3\\pi 3\\pi \\pi^0$"] = f"${eps_mass_vetoes_3pi3pipi0[3].nominal_value*100:.2f} \\pm {eps_mass_vetoes_3pi3pipi0[3].std_dev*100:.2f}$"
        mass_vetoes_table.loc["$IP > 0.25$\\,mm ($K^{*0}$)", "$3\\pi 3\\pi \\pi^0$"] = f"${eps_mass_vetoes_3pi3pipi0[4].nominal_value*100:.2f} \\pm {eps_mass_vetoes_3pi3pipi0[4].std_dev*100:.2f}$"
        mass_vetoes_table.loc["$IP > 0.15$\\,mm ($D^{0}$)", "$3\\pi 3\\pi \\pi^0$"] = f"${eps_mass_vetoes_3pi3pipi0[5].nominal_value*100:.2f} \\pm {eps_mass_vetoes_3pi3pipi0[5].std_dev*100:.2f}$"

        mass_vetoes_table.loc["$D^0$ (2 particles)", "$3\\pi 3\\pi 2\\pi^0$"] = f"${eps_mass_vetoes_3pi3pi2pi0[0].nominal_value*100:.2f} \\pm {eps_mass_vetoes_3pi3pi2pi0[0].std_dev*100:.2f}$"
        mass_vetoes_table.loc["$D^-$ (3 particles)", "$3\\pi 3\\pi 2\\pi^0$"] = f"${eps_mass_vetoes_3pi3pi2pi0[1].nominal_value*100:.2f} \\pm {eps_mass_vetoes_3pi3pi2pi0[1].std_dev*100:.2f}$"
        mass_vetoes_table.loc["$D^0$ (4 particles)", "$3\\pi 3\\pi 2\\pi^0$"] = f"${eps_mass_vetoes_3pi3pi2pi0[2].nominal_value*100:.2f} \\pm {eps_mass_vetoes_3pi3pi2pi0[2].std_dev*100:.2f}$"
        mass_vetoes_table.loc["$D^{*-}$ (5 particles)", "$3\\pi 3\\pi 2\\pi^0$"] = f"${eps_mass_vetoes_3pi3pi2pi0[3].nominal_value*100:.2f} \\pm {eps_mass_vetoes_3pi3pi2pi0[3].std_dev*100:.2f}$"
        mass_vetoes_table.loc["$IP > 0.25$\\,mm ($K^{*0}$)", "$3\\pi 3\\pi 2\\pi^0$"] = f"${eps_mass_vetoes_3pi3pipi0[4].nominal_value*100:.2f} \\pm {eps_mass_vetoes_3pi3pipi0[4].std_dev*100:.2f}$"
        mass_vetoes_table.loc["$IP > 0.15$\\,mm ($D^{0}$)", "$3\\pi 3\\pi 2\\pi^0$"] = f"${eps_mass_vetoes_3pi3pi2pi0[5].nominal_value*100:.2f} \\pm {eps_mass_vetoes_3pi3pi2pi0[5].std_dev*100:.2f}$"

        mass_vetoes_table.loc["$D^0$ (2 particles)", "All MC"] = f"${eps_mass_vetoes_MC[0].nominal_value*100:.2f} \\pm {eps_mass_vetoes_MC[0].std_dev*100:.2f}$"
        mass_vetoes_table.loc["$D^-$ (3 particles)", "All MC"] = f"${eps_mass_vetoes_MC[1].nominal_value*100:.2f} \\pm {eps_mass_vetoes_MC[1].std_dev*100:.2f}$"
        mass_vetoes_table.loc["$D^0$ (4 particles)", "All MC"] = f"${eps_mass_vetoes_MC[2].nominal_value*100:.2f} \\pm {eps_mass_vetoes_MC[2].std_dev*100:.2f}$"
        mass_vetoes_table.loc["$D^{*-}$ (5 particles)", "All MC"] = f"${eps_mass_vetoes_MC[3].nominal_value*100:.2f} \\pm {eps_mass_vetoes_MC[3].std_dev*100:.2f}$"
        mass_vetoes_table.loc["$IP > 0.25$\\,mm ($K^{*0}$)", "All MC"] = f"${eps_mass_vetoes_MC[4].nominal_value*100:.2f} \\pm {eps_mass_vetoes_MC[4].std_dev*100:.2f}$"
        mass_vetoes_table.loc["$IP > 0.15$\\,mm ($D^{0}$)", "All MC"] = f"${eps_mass_vetoes_MC[5].nominal_value*100:.2f} \\pm {eps_mass_vetoes_MC[5].std_dev*100:.2f}$"

        mass_vetoes_table.loc["$D^0$ (2 particles)", "RS data"] = f"${eps_mass_vetoes_RS[0].nominal_value*100:.2f} \\pm {eps_mass_vetoes_RS[0].std_dev*100:.2f}$"
        mass_vetoes_table.loc["$D^-$ (3 particles)", "RS data"] = f"${eps_mass_vetoes_RS[1].nominal_value*100:.2f} \\pm {eps_mass_vetoes_RS[1].std_dev*100:.2f}$"
        mass_vetoes_table.loc["$D^0$ (4 particles)", "RS data"] = f"${eps_mass_vetoes_RS[2].nominal_value*100:.2f} \\pm {eps_mass_vetoes_RS[2].std_dev*100:.2f}$"
        mass_vetoes_table.loc["$D^{*-}$ (5 particles)", "RS data"] = f"${eps_mass_vetoes_RS[3].nominal_value*100:.2f} \\pm {eps_mass_vetoes_RS[3].std_dev*100:.2f}$"
        mass_vetoes_table.loc["$IP > 0.25$\\,mm ($K^{*0}$)", "RS data"] = f"${eps_mass_vetoes_RS[4].nominal_value*100:.2f} \\pm {eps_mass_vetoes_RS[4].std_dev*100:.2f}$"
        mass_vetoes_table.loc["$IP > 0.15$\\,mm ($D^{0}$)", "RS data"] = f"${eps_mass_vetoes_RS[5].nominal_value*100:.2f} \\pm {eps_mass_vetoes_RS[5].std_dev*100:.2f}$"

        mass_vetoes_table.loc["$D^0$ (2 particles)", "WS data"] = f"${eps_mass_vetoes_WS[0].nominal_value*100:.2f} \\pm {eps_mass_vetoes_WS[0].std_dev*100:.2f}$"
        mass_vetoes_table.loc["$D^-$ (3 particles)", "WS data"] = f"${eps_mass_vetoes_WS[1].nominal_value*100:.2f} \\pm {eps_mass_vetoes_WS[1].std_dev*100:.2f}$"
        mass_vetoes_table.loc["$D^0$ (4 particles)", "WS data"] = f"${eps_mass_vetoes_WS[2].nominal_value*100:.2f} \\pm {eps_mass_vetoes_WS[2].std_dev*100:.2f}$"
        mass_vetoes_table.loc["$D^{*-}$ (5 particles)", "WS data"] = f"${eps_mass_vetoes_WS[3].nominal_value*100:.2f} \\pm {eps_mass_vetoes_WS[3].std_dev*100:.2f}$"
        mass_vetoes_table.loc["$IP > 0.25$\\,mm ($K^{*0}$)", "WS data"] = f"${eps_mass_vetoes_WS[4].nominal_value*100:.2f} \\pm {eps_mass_vetoes_WS[4].std_dev*100:.2f}$"
        mass_vetoes_table.loc["$IP > 0.15$\\,mm ($D^{0}$)", "WS data"] = f"${eps_mass_vetoes_WS[5].nominal_value*100:.2f} \\pm {eps_mass_vetoes_WS[5].std_dev*100:.2f}$"

        mass_vetoes_table.loc["All", "$3\\pi 3\\pi$"] = f"${eps_all_mass_vetoes_3pi3pi.nominal_value*100:.2f} \\pm {eps_all_mass_vetoes_3pi3pi.std_dev*100:.2f}$"
        mass_vetoes_table.loc["All", "$3\\pi 3\\pi \\pi^0$"] = f"${eps_all_mass_vetoes_3pi3pipi0.nominal_value*100:.2f} \\pm {eps_all_mass_vetoes_3pi3pipi0.std_dev*100:.2f}$"
        mass_vetoes_table.loc["All", "$3\\pi 3\\pi 2\\pi^0$"] = f"${eps_all_mass_vetoes_3pi3pi2pi0.nominal_value*100:.2f} \\pm {eps_all_mass_vetoes_3pi3pi2pi0.std_dev*100:.2f}$"
        mass_vetoes_table.loc["All", "All MC"] = f"${eps_all_mass_vetoes_MC.nominal_value*100:.2f} \\pm {eps_all_mass_vetoes_MC.std_dev*100:.2f}$"
        mass_vetoes_table.loc["All", "RS data"] = f"${eps_all_mass_vetoes_RS.nominal_value*100:.2f} \\pm {eps_all_mass_vetoes_RS.std_dev*100:.2f}$"
        mass_vetoes_table.loc["All", "WS data"] = f"${eps_all_mass_vetoes_WS.nominal_value*100:.2f} \\pm {eps_all_mass_vetoes_WS.std_dev*100:.2f}$"


    if((species == 100) or (species == 1)):
        with open(f'/panfs/felician/B2Ktautau/workflow/selections_efficiency_tables/Species_{species}/mass_vetoes_efficiency_table.tex', 'w') as fout_mass_vetoes:
            fout_mass_vetoes.write(mass_vetoes_table.to_latex())

    # BDT efficiency (BDT > 0.999)
    if(species == 1):
        f_mc_post_sel = ROOT.TFile("/panfs/felician/B2Ktautau/workflow/create_post_selection_tree/Species_10/post_sel_tree_bdt_0.0.root")
        t_mc_post_sel = f_mc_post_sel.Get("DecayTree")
        f_mc_best_cand = ROOT.TFile("/panfs/felician/B2Ktautau/workflow/candidate_selection/Species_10/candidate_bdt_0.999.root")
        t_mc_best_cand = f_mc_best_cand.Get("DecayTree")

        f_rs_post_sel = ROOT.TFile("/panfs/felician/B2Ktautau/workflow/create_post_selection_tree/Species_2/post_sel_tree_bdt_0.0.root")
        t_rs_post_sel = f_rs_post_sel.Get("DecayTree")
        f_rs_best_cand = ROOT.TFile("/panfs/felician/B2Ktautau/workflow/candidate_selection/Species_2/candidate_bdt_0.999.root")
        t_rs_best_cand = f_rs_best_cand.Get("DecayTree")

        f_ws_post_sel = ROOT.TFile("/panfs/felician/B2Ktautau/workflow/create_post_selection_tree/Species_3/post_sel_tree_bdt_0.0.root")
        t_ws_post_sel = f_ws_post_sel.Get("DecayTree")
        f_ws_best_cand = ROOT.TFile("/panfs/felician/B2Ktautau/workflow/candidate_selection/Species_3/candidate_bdt_0.999.root")
        t_ws_best_cand = f_ws_best_cand.Get("DecayTree")

        N_bdt_3pi3pi = t_mc_post_sel.GetEntries(truthMatch+" && "+trigger+" && "+all_rectangular_cuts_MC+" && "+pass_mass_fit+" && "+fit_range+" && (component == 0)"+" && "+mass_vetoes+" && (BDT > 0.999)")
        N_bdt_3pi3pipi0 = t_mc_post_sel.GetEntries(truthMatch+" && "+trigger+" && "+all_rectangular_cuts_MC+" && "+pass_mass_fit+" && "+fit_range+" && (component == 1)"+" && "+mass_vetoes+" && (BDT > 0.999)")
        N_bdt_3pi3pi2pi0 = t_mc_post_sel.GetEntries(truthMatch+" && "+trigger+" && "+all_rectangular_cuts_MC+" && "+pass_mass_fit+" && "+fit_range+" && (component == 2)"+" && "+mass_vetoes+" && (BDT > 0.999)")
        N_bdt_MC = t_mc_post_sel.GetEntries(trigger+" && "+all_rectangular_cuts_MC+" && "+pass_mass_fit+" && "+fit_range+" && "+mass_vetoes+" && (BDT > 0.999)")
        N_bdt_RS = t_rs_post_sel.GetEntries(trigger+" && "+all_rectangular_cuts_data+" && "+pass_mass_fit+" && "+fit_range+" && "+mass_vetoes+" && (BDT > 0.999)")
        N_bdt_WS = t_ws_post_sel.GetEntries(trigger+" && "+all_rectangular_cuts_data+" && "+pass_mass_fit+" && "+fit_range+" && "+mass_vetoes+" && (BDT > 0.999)")

        up_bdt_3pi3pi = ROOT.TEfficiency.Wilson(N_all_mass_vetoes_3pi3pi, N_bdt_3pi3pi, 0.68, True)
        down_bdt_3pi3pi = ROOT.TEfficiency.Wilson(N_all_mass_vetoes_3pi3pi, N_bdt_3pi3pi, 0.68, False)
        eps_bdt_3pi3pi_err = 0.5*(up_bdt_3pi3pi - down_bdt_3pi3pi)
        eps_bdt_3pi3pi =  ufloat(N_bdt_3pi3pi/N_all_mass_vetoes_3pi3pi, eps_bdt_3pi3pi_err) 

        up_bdt_3pi3pipi0 = ROOT.TEfficiency.Wilson(N_all_mass_vetoes_3pi3pipi0, N_bdt_3pi3pipi0, 0.68, True)
        down_bdt_3pi3pipi0 = ROOT.TEfficiency.Wilson(N_all_mass_vetoes_3pi3pipi0, N_bdt_3pi3pipi0, 0.68, False)
        eps_bdt_3pi3pipi0_err = 0.5*(up_bdt_3pi3pipi0 - down_bdt_3pi3pipi0)
        eps_bdt_3pi3pipi0 =  ufloat(N_bdt_3pi3pipi0/N_all_mass_vetoes_3pi3pipi0, eps_bdt_3pi3pipi0_err) 

        up_bdt_3pi3pi2pi0 = ROOT.TEfficiency.Wilson(N_all_mass_vetoes_3pi3pi2pi0, N_bdt_3pi3pi2pi0, 0.68, True)
        down_bdt_3pi3pi2pi0 = ROOT.TEfficiency.Wilson(N_all_mass_vetoes_3pi3pi2pi0, N_bdt_3pi3pi2pi0, 0.68, False)
        eps_bdt_3pi3pi2pi0_err = 0.5*(up_bdt_3pi3pi2pi0 - down_bdt_3pi3pi2pi0)
        eps_bdt_3pi3pi2pi0 =  ufloat(N_bdt_3pi3pi2pi0/N_all_mass_vetoes_3pi3pi2pi0, eps_bdt_3pi3pi2pi0_err) 

        up_bdt_MC = ROOT.TEfficiency.Wilson(N_all_mass_vetoes_MC, N_bdt_MC, 0.68, True)
        down_bdt_MC = ROOT.TEfficiency.Wilson(N_all_mass_vetoes_MC, N_bdt_MC, 0.68, False)
        eps_bdt_MC_err = 0.5*(up_bdt_MC - down_bdt_MC)
        eps_bdt_MC =  ufloat(N_bdt_MC/N_all_mass_vetoes_MC, eps_bdt_MC_err) 
    
        up_bdt_RS = ROOT.TEfficiency.Wilson(N_all_mass_vetoes_RS, N_bdt_RS, 0.68, True)
        down_bdt_RS = ROOT.TEfficiency.Wilson(N_all_mass_vetoes_RS, N_bdt_RS, 0.68, False)
        eps_bdt_RS_err = 0.5*(up_bdt_RS - down_bdt_RS)
        eps_bdt_RS =  ufloat(N_bdt_RS/N_all_mass_vetoes_RS, eps_bdt_RS_err) 

        up_bdt_WS = ROOT.TEfficiency.Wilson(N_all_mass_vetoes_WS, N_bdt_WS, 0.68, True)
        down_bdt_WS = ROOT.TEfficiency.Wilson(N_all_mass_vetoes_WS, N_bdt_WS, 0.68, False)
        eps_bdt_WS_err = 0.5*(up_bdt_WS - down_bdt_WS)
        eps_bdt_WS =  ufloat(N_bdt_WS/N_all_mass_vetoes_WS, eps_bdt_WS_err) 

        bdt_columns = ["BDT efficiency"]
        bdt_rows = ["$3\\pi 3\\pi$", "$3\\pi 3\\pi \\pi^0$", "$3\\pi 3\\pi 2\\pi^0$", "All MC", "RS data", "WS data"]

        bdt_table = pd.DataFrame(index=bdt_rows, columns=bdt_columns)
        bdt_table.loc["$3\\pi 3\\pi$", "BDT efficiency"] = f"${eps_bdt_3pi3pi.nominal_value*100:.2f} \\pm {eps_bdt_3pi3pi.std_dev*100:.2f}$"
        bdt_table.loc["$3\\pi 3\\pi \\pi^0$", "BDT efficiency"] = f"${eps_bdt_3pi3pipi0.nominal_value*100:.2f} \\pm {eps_bdt_3pi3pipi0.std_dev*100:.2f}$"
        bdt_table.loc["$3\\pi 3\\pi 2\\pi^0$", "BDT efficiency"] = f"${eps_bdt_3pi3pi2pi0.nominal_value*100:.2f} \\pm {eps_bdt_3pi3pi2pi0.std_dev*100:.2f}$"
        bdt_table.loc["All MC", "BDT efficiency"] = f"${eps_bdt_MC.nominal_value*100:.2f} \\pm {eps_bdt_MC.std_dev*100:.2f}$"
        # bdt_table.loc["RS data", "BDT efficiency"] = f"${eps_bdt_RS.nominal_value*100:.2f} \\pm {eps_bdt_RS.std_dev*100:.2f}$"
        bdt_table.loc["RS data", "BDT efficiency"] = f"${-1} \\pm {-1}$"
        bdt_table.loc["WS data", "BDT efficiency"] = f"${eps_bdt_WS.nominal_value*100:.2f} \\pm {eps_bdt_WS.std_dev*100:.2f}$"
    
    elif(species == 100):
        f_BuDDKp_post_sel = ROOT.TFile("/panfs/felician/B2Ktautau/workflow/create_post_selection_tree/Species_100/post_sel_tree_bdt_0.0.root")
        t_BuDDKp_post_sel = f_BuDDKp_post_sel.Get("DecayTree")
        f_BuDDKp_best_cand = ROOT.TFile("/panfs/felician/B2Ktautau/workflow/candidate_selection/Species_100/candidate_bdt_0.999.root")
        t_BuDDKp_best_cand = f_BuDDKp_best_cand.Get("DecayTree")

        f_BdDDKp_post_sel = ROOT.TFile("/panfs/felician/B2Ktautau/workflow/create_post_selection_tree/Species_110/post_sel_tree_bdt_0.0.root")
        t_BdDDKp_post_sel = f_BdDDKp_post_sel.Get("DecayTree")
        f_BdDDKp_best_cand = ROOT.TFile("/panfs/felician/B2Ktautau/workflow/candidate_selection/Species_110/candidate_bdt_0.999.root")
        t_BdDDKp_best_cand = f_BdDDKp_best_cand.Get("DecayTree")

        f_BsDDKp_post_sel = ROOT.TFile("/panfs/felician/B2Ktautau/workflow/create_post_selection_tree/Species_120/post_sel_tree_bdt_0.0.root")
        t_BsDDKp_post_sel = f_BsDDKp_post_sel.Get("DecayTree")
        f_BsDDKp_best_cand = ROOT.TFile("/panfs/felician/B2Ktautau/workflow/candidate_selection/Species_120/candidate_bdt_0.999.root")
        t_BsDDKp_best_cand = f_BsDDKp_best_cand.Get("DecayTree")

        f_BuDDK0_post_sel = ROOT.TFile("/panfs/felician/B2Ktautau/workflow/create_post_selection_tree/Species_130/post_sel_tree_bdt_0.0.root")
        t_BuDDK0_post_sel = f_BuDDK0_post_sel.Get("DecayTree")
        f_BuDDK0_best_cand = ROOT.TFile("/panfs/felician/B2Ktautau/workflow/candidate_selection/Species_130/candidate_bdt_0.999.root")
        t_BuDDK0_best_cand = f_BuDDK0_best_cand.Get("DecayTree")

        f_BuDD_post_sel = ROOT.TFile("/panfs/felician/B2Ktautau/workflow/create_post_selection_tree/Species_150/post_sel_tree_bdt_0.0.root")
        t_BuDD_post_sel = f_BuDD_post_sel.Get("DecayTree")
        f_BuDD_best_cand = ROOT.TFile("/panfs/felician/B2Ktautau/workflow/candidate_selection/Species_150/candidate_bdt_0.999.root")
        t_BuDD_best_cand = f_BuDD_best_cand.Get("DecayTree")

        N_bdt_BuDDKp = t_BuDDKp_post_sel.GetEntries(trigger+" && "+all_rectangular_cuts+" && "+pass_mass_fit+" && "+fit_range+" && "+mass_vetoes+" && (BDT > 0.999)") 
        N_bdt_BdDDKp = t_BdDDKp_post_sel.GetEntries(trigger+" && "+all_rectangular_cuts+" && "+pass_mass_fit+" && "+fit_range+" && "+mass_vetoes+" && (BDT > 0.999)") 
        N_bdt_BsDDKp = t_BsDDKp_post_sel.GetEntries(trigger+" && "+all_rectangular_cuts+" && "+pass_mass_fit+" && "+fit_range+" && "+mass_vetoes+" && (BDT > 0.999)") 
        N_bdt_BuDDK0 = t_BuDDK0_post_sel.GetEntries(trigger+" && "+all_rectangular_cuts+" && "+pass_mass_fit+" && "+fit_range+" && "+mass_vetoes+" && (BDT > 0.999)") 
        N_bdt_BuDD = t_BuDD_post_sel.GetEntries(trigger+" && "+all_rectangular_cuts+" && "+pass_mass_fit+" && "+fit_range+" && "+mass_vetoes+" && (BDT > 0.999)") 

        up_bdt_BuDDKp = ROOT.TEfficiency.Wilson(N_mass_vetoes_BuDDKp, N_bdt_BuDDKp, 0.68, True)
        down_bdt_BuDDKp = ROOT.TEfficiency.Wilson(N_mass_vetoes_BuDDKp, N_bdt_BuDDKp, 0.68, False)
        eps_bdt_BuDDKp_err = 0.5*(up_bdt_BuDDKp - down_bdt_BuDDKp)
        eps_bdt_BuDDKp = ufloat(N_bdt_BuDDKp/N_mass_vetoes_BuDDKp, eps_bdt_BuDDKp_err)

        up_bdt_BdDDKp = ROOT.TEfficiency.Wilson(N_mass_vetoes_BdDDKp, N_bdt_BdDDKp, 0.68, True)
        down_bdt_BdDDKp = ROOT.TEfficiency.Wilson(N_mass_vetoes_BdDDKp, N_bdt_BdDDKp, 0.68, False)
        eps_bdt_BdDDKp_err = 0.5*(up_bdt_BdDDKp - down_bdt_BdDDKp)
        eps_bdt_BdDDKp = ufloat(N_bdt_BdDDKp/N_mass_vetoes_BdDDKp, eps_bdt_BdDDKp_err)

        up_bdt_BsDDKp = ROOT.TEfficiency.Wilson(N_mass_vetoes_BsDDKp, N_bdt_BsDDKp, 0.68, True)
        down_bdt_BsDDKp = ROOT.TEfficiency.Wilson(N_mass_vetoes_BsDDKp, N_bdt_BsDDKp, 0.68, False)
        eps_bdt_BsDDKp_err = 0.5*(up_bdt_BsDDKp - down_bdt_BsDDKp)
        eps_bdt_BsDDKp = ufloat(N_bdt_BsDDKp/N_mass_vetoes_BsDDKp, eps_bdt_BsDDKp_err)

        up_bdt_BuDDK0 = ROOT.TEfficiency.Wilson(N_mass_vetoes_BuDDK0, N_bdt_BuDDK0, 0.68, True)
        down_bdt_BuDDK0 = ROOT.TEfficiency.Wilson(N_mass_vetoes_BuDDK0, N_bdt_BuDDK0, 0.68, False)
        eps_bdt_BuDDK0_err = 0.5*(up_bdt_BuDDK0 - down_bdt_BuDDK0)
        eps_bdt_BuDDK0 = ufloat(N_bdt_BuDDK0/N_mass_vetoes_BuDDK0, eps_bdt_BuDDK0_err)

        up_bdt_BuDD = ROOT.TEfficiency.Wilson(N_mass_vetoes_BuDD, N_bdt_BuDD, 0.68, True)
        down_bdt_BuDD = ROOT.TEfficiency.Wilson(N_mass_vetoes_BuDD, N_bdt_BuDD, 0.68, False)
        eps_bdt_BuDD_err = 0.5*(up_bdt_BuDD - down_bdt_BuDD)
        eps_bdt_BuDD = ufloat(N_bdt_BuDD/N_mass_vetoes_BuDD, eps_bdt_BuDD_err)

        bdt_columns = ["BDT efficiency"]
        bdt_rows = ["$B^+ \\to D D K^+$", "$B^0 \\to D D K^+$", "$B^0_s \\to D D K^+$", "$B^+ \\to D D K^0$", "$B^+ \\to D D$"]

        bdt_table = pd.DataFrame(index=bdt_rows, columns=bdt_columns)
        bdt_table.loc["$B^+ \\to D D K^+$", "BDT efficiency"] = f"${eps_bdt_BuDDKp.nominal_value*100:.2f} \\pm {eps_bdt_BuDDKp.std_dev*100:.2f}$"
        bdt_table.loc["$B^0 \\to D D K^+$", "BDT efficiency"] = f"${eps_bdt_BdDDKp.nominal_value*100:.2f} \\pm {eps_bdt_BdDDKp.std_dev*100:.2f}$"
        bdt_table.loc["$B^0_s \\to D D K^+$", "BDT efficiency"] = f"${eps_bdt_BsDDKp.nominal_value*100:.2f} \\pm {eps_bdt_BsDDKp.std_dev*100:.2f}$"
        bdt_table.loc["$B^+ \\to D D K^0$", "BDT efficiency"] = f"${eps_bdt_BuDDK0.nominal_value*100:.2f} \\pm {eps_bdt_BuDDK0.std_dev*100:.2f}$"
        bdt_table.loc["$B^+ \\to D D$", "BDT efficiency"] = f"${eps_bdt_BuDD.nominal_value*100:.2f} \\pm {eps_bdt_BuDD.std_dev*100:.2f}$"

    if((species == 100) or (species == 1)):
        with open(f'/panfs/felician/B2Ktautau/workflow/selections_efficiency_tables/Species_{species}/bdt_efficiency_table.tex', 'w') as fout_bdt:
            fout_bdt.write(bdt_table.to_latex())


    # Best candidate
    print("Best candidate efficiency")
    if(species == 7):
        N_best_cand_MC = t_reco_mc_2016.GetEntries(truthMatch+" && "+trigger+" && "+all_rectangular_cuts_MC+" && "+pass_mass_fit+" && "+fit_range+" && "+best_cand) + t_reco_mc_2017.GetEntries(trigger+" && "+all_rectangular_cuts_MC+" && "+pass_mass_fit+" && "+fit_range+" && "+best_cand) + t_reco_mc_2018.GetEntries(trigger+" && "+all_rectangular_cuts_MC+" && "+pass_mass_fit+" && "+fit_range+" && "+best_cand)

        up_best_cand_MC = ROOT.TEfficiency.Wilson(N_fit_region_MC, N_best_cand_MC, 0.68, True)
        down_best_cand_MC = ROOT.TEfficiency.Wilson(N_fit_region_MC, N_best_cand_MC, 0.68, False)
        eps_best_cand_MC_err = 0.5*(up_best_cand_MC - down_best_cand_MC)
        eps_best_cand_MC = ufloat(N_best_cand_MC/N_fit_region_MC, eps_best_cand_MC_err)

        best_cand_columns = ["Best candidate efficiency (\\%)"]
        best_cand_rows  = ["MC"]

        best_cand_table = pd.DataFrame(index=best_cand_rows, columns=best_cand_columns)
        best_cand_table.loc["MC", "Best candidate efficiency (\\%)"] = f"${eps_best_cand_MC.nominal_value*100:.2f} \\pm {eps_best_cand_MC.std_dev*100:.2f}$"

        with open(f'/panfs/felician/B2Ktautau/workflow/selections_efficiency_tables/Species_{species}/best_candidate_efficiency_table.tex', 'w') as fout_best_cand:
            fout_best_cand.write(best_cand_table.to_latex())

    elif(species == 1):
        N_best_cand_MC = t_mc_best_cand.GetEntries("(is_best_cand == 1)")

        up_best_cand_MC = ROOT.TEfficiency.Wilson(N_bdt_MC, N_best_cand_MC, 0.68, True)
        down_best_cand_MC = ROOT.TEfficiency.Wilson(N_bdt_MC, N_best_cand_MC, 0.68, False)
        eps_best_cand_MC_err = 0.5*(up_best_cand_MC - down_best_cand_MC)
        eps_best_cand_MC =  ufloat(N_best_cand_MC/N_bdt_MC, eps_best_cand_MC_err) 

        best_cand_columns = ["Candidate selection efficiency"]
        best_cand_rows = ["All MC"]

        best_cand_table = pd.DataFrame(index=best_cand_rows, columns=best_cand_columns)
        best_cand_table.loc["All MC", "Candidate selection efficiency"] = f"${eps_best_cand_MC.nominal_value*100:.2f} \\pm {eps_best_cand_MC.std_dev*100:.2f}$"

    elif(species == 100):
        N_best_cand_BuDDKp = t_BuDDKp_best_cand.GetEntries("(is_best_cand == 1)") 
        N_best_cand_BdDDKp = t_BdDDKp_best_cand.GetEntries("(is_best_cand == 1)") 
        N_best_cand_BsDDKp = t_BsDDKp_best_cand.GetEntries("(is_best_cand == 1)") 
        N_best_cand_BuDDK0 = t_BuDDK0_best_cand.GetEntries("(is_best_cand == 1)") 
        N_best_cand_BuDD = t_BuDD_best_cand.GetEntries("(is_best_cand == 1)") 

        up_best_cand_BuDDKp = ROOT.TEfficiency.Wilson(N_bdt_BuDDKp, N_best_cand_BuDDKp, 0.68, True)
        down_best_cand_BuDDKp = ROOT.TEfficiency.Wilson(N_bdt_BuDDKp, N_best_cand_BuDDKp, 0.68, False)
        eps_best_cand_BuDDKp_err = 0.5*(up_best_cand_BuDDKp - down_best_cand_BuDDKp)
        eps_best_cand_BuDDKp = ufloat(N_best_cand_BuDDKp/N_bdt_BuDDKp, eps_best_cand_BuDDKp_err)

        up_best_cand_BdDDKp = ROOT.TEfficiency.Wilson(N_bdt_BdDDKp, N_best_cand_BdDDKp, 0.68, True)
        down_best_cand_BdDDKp = ROOT.TEfficiency.Wilson(N_bdt_BdDDKp, N_best_cand_BdDDKp, 0.68, False)
        eps_best_cand_BdDDKp_err = 0.5*(up_best_cand_BdDDKp - down_best_cand_BdDDKp)
        eps_best_cand_BdDDKp = ufloat(N_best_cand_BdDDKp/N_bdt_BdDDKp, eps_best_cand_BdDDKp_err)

        up_best_cand_BsDDKp = ROOT.TEfficiency.Wilson(N_bdt_BsDDKp, N_best_cand_BsDDKp, 0.68, True)
        down_best_cand_BsDDKp = ROOT.TEfficiency.Wilson(N_bdt_BsDDKp, N_best_cand_BsDDKp, 0.68, False)
        eps_best_cand_BsDDKp_err = 0.5*(up_best_cand_BsDDKp - down_best_cand_BsDDKp)
        eps_best_cand_BsDDKp = ufloat(N_best_cand_BsDDKp/N_bdt_BsDDKp, eps_best_cand_BsDDKp_err)

        up_best_cand_BuDDK0 = ROOT.TEfficiency.Wilson(N_bdt_BuDDK0, N_best_cand_BuDDK0, 0.68, True)
        down_best_cand_BuDDK0 = ROOT.TEfficiency.Wilson(N_bdt_BuDDK0, N_best_cand_BuDDK0, 0.68, False)
        eps_best_cand_BuDDK0_err = 0.5*(up_best_cand_BuDDK0 - down_best_cand_BuDDK0)
        if(N_bdt_BuDDK0 == 0):
            eps_best_cand_BuDDK0 = ufloat(0, eps_best_cand_BuDDK0_err)
        else:
            eps_best_cand_BuDDK0 = ufloat(N_best_cand_BuDDK0/N_bdt_BuDDK0, eps_best_cand_BuDDK0_err)

        up_best_cand_BuDD = ROOT.TEfficiency.Wilson(N_bdt_BuDD, N_best_cand_BuDD, 0.68, True)
        down_best_cand_BuDD = ROOT.TEfficiency.Wilson(N_bdt_BuDD, N_best_cand_BuDD, 0.68, False)
        eps_best_cand_BuDD_err = 0.5*(up_best_cand_BuDD - down_best_cand_BuDD)
        eps_best_cand_BuDD = ufloat(N_best_cand_BuDD/N_bdt_BuDD, eps_best_cand_BuDD_err)

        best_cand_columns = ["Candidate efficiency"]
        best_cand_rows = ["$B^+ \\to D D K^+$", "$B^0 \\to D D K^+$", "$B^0_s \\to D D K^+$", "$B^+ \\to D D K^0$", "$B^+ \\to D D$"]

        best_cand_table = pd.DataFrame(index=best_cand_rows, columns=best_cand_columns)
        best_cand_table.loc["$B^+ \\to D D K^+$", "Candidate efficiency"] = f"${eps_best_cand_BuDDKp.nominal_value*100:.2f} \\pm {eps_best_cand_BuDDKp.std_dev*100:.2f}$"
        best_cand_table.loc["$B^0 \\to D D K^+$", "Candidate efficiency"] = f"${eps_best_cand_BdDDKp.nominal_value*100:.2f} \\pm {eps_best_cand_BdDDKp.std_dev*100:.2f}$"
        best_cand_table.loc["$B^0_s \\to D D K^+$", "Candidate efficiency"] = f"${eps_best_cand_BsDDKp.nominal_value*100:.2f} \\pm {eps_best_cand_BsDDKp.std_dev*100:.2f}$"
        best_cand_table.loc["$B^+ \\to D D K^0$", "Candidate efficiency"] = f"${eps_best_cand_BuDDK0.nominal_value*100:.2f} \\pm {eps_best_cand_BuDDK0.std_dev*100:.2f}$"
        best_cand_table.loc["$B^+ \\to D D$", "Candidate efficiency"] = f"${eps_best_cand_BuDD.nominal_value*100:.2f} \\pm {eps_best_cand_BuDD.std_dev*100:.2f}$"

    if((species == 100) or (species == 1)):
        with open(f'/panfs/felician/B2Ktautau/workflow/selections_efficiency_tables/Species_{species}/best_cand_efficiency_table.tex', 'w') as fout_best_cand:
            fout_best_cand.write(best_cand_table.to_latex())        


    # Overall efficiency
    print("Overall efficiency")

    if(species == 7):
        up_overall_MC = ROOT.TEfficiency.Wilson(N_gen, N_best_cand_MC, 0.68, True)
        down_overall_MC = ROOT.TEfficiency.Wilson(N_gen, N_best_cand_MC, 0.68, False)
        eps_overall_MC_err = 0.5*(up_overall_MC - down_overall_MC)
        eps_overall_MC = ufloat(N_best_cand_MC/N_gen, eps_overall_MC_err)

        eps_overall_MC *= eps_acc

        overall_columns = ["Efficiency (\\%)"]
        overall_rows  = ["Acceptance", "Reconstruction", "Trigger", "Rectangular cuts", "Pass mass fit", "Fit range", "Overall"]

        overall_table = pd.DataFrame(index=overall_rows, columns=overall_columns)
        overall_table.loc["Acceptance", "Efficiency (\\%)"] = f"${eps_acc.nominal_value*100:.2f} \\pm {eps_acc.std_dev*100:.2f}$"
        overall_table.loc["Reconstruction", "Efficiency (\\%)"] = f"${eps_reco.nominal_value*100:.3f} \\pm {eps_reco.std_dev*100:.3f}$"
        overall_table.loc["Trigger", "Efficiency (\\%)"] = f"${eps_trigger.nominal_value*100:.2f} \\pm {eps_trigger.std_dev*100:.2f}$"
        overall_table.loc["Rectangular cuts", "Efficiency (\\%)"] = f"${eps_rect_cuts_MC_all.nominal_value*100:.1f} \\pm {eps_rect_cuts_MC_all.std_dev*100:.1f}$"
        overall_table.loc["Pass mass fit", "Efficiency (\\%)"] = f"${eps_pass_fit_MC.nominal_value*100:.2f} \\pm {eps_pass_fit_MC.std_dev*100:.2f}$"
        overall_table.loc["Fit range", "Efficiency (\\%)"] = f"${eps_fit_region_MC.nominal_value*100:.2f} \\pm {eps_fit_region_MC.std_dev*100:.2f}$"
        overall_table.loc["Overall", "Efficiency (\\%)"] = f"${eps_overall_MC.nominal_value*100:.4f} \\pm {eps_overall_MC.std_dev*100:.4f}$"

    elif(species == 1):
        up_overall_MC = ROOT.TEfficiency.Wilson(N_gen, N_best_cand_MC, 0.68, True)
        down_overall_MC = ROOT.TEfficiency.Wilson(N_gen, N_best_cand_MC, 0.68, False)
        eps_overall_MC_err = 0.5*(up_overall_MC - down_overall_MC)
        eps_overall_MC = ufloat(N_best_cand_MC/N_gen, eps_overall_MC_err)

        eps_overall_MC *= eps_acc*eps_strip

        overall_columns = ["Efficiency (\\%)"]
        overall_rows  = ["Acceptance", "Stripping", "Reconstruction", "Trigger", "Rectangular cuts", "Pass mass fit", "Fit range", "Mass vetoes", "BDT", "Overall"]

        overall_table = pd.DataFrame(index=overall_rows, columns=overall_columns)
        overall_table.loc["Acceptance", "Efficiency (\\%)"] = f"${eps_acc.nominal_value*100:.5f} \\pm {eps_acc.std_dev*100:.5f}$"
        overall_table.loc["Stripping", "Efficiency (\\%)"] = f"${eps_strip.nominal_value*100:.5f} \\pm {eps_strip.std_dev*100:.5f}$"
        overall_table.loc["Reconstruction", "Efficiency (\\%)"] = f"${eps_reco.nominal_value*100:.3f} \\pm {eps_reco.std_dev*100:.3f}$"
        overall_table.loc["Trigger", "Efficiency (\\%)"] = f"${eps_trigger.nominal_value*100:.2f} \\pm {eps_trigger.std_dev*100:.2f}$"
        overall_table.loc["Rectangular cuts", "Efficiency (\\%)"] = f"${eps_rect_cuts_MC_all.nominal_value*100:.1f} \\pm {eps_rect_cuts_MC_all.std_dev*100:.1f}$"
        overall_table.loc["Pass mass fit", "Efficiency (\\%)"] = f"${eps_pass_fit_MC.nominal_value*100:.2f} \\pm {eps_pass_fit_MC.std_dev*100:.2f}$"
        overall_table.loc["Fit range", "Efficiency (\\%)"] = f"${eps_fit_region_MC.nominal_value*100:.2f} \\pm {eps_fit_region_MC.std_dev*100:.2f}$"
        overall_table.loc["Mass vetoes", "Efficiency (\\%)"] = f"${eps_all_mass_vetoes_MC.nominal_value*100:.2f} \\pm {eps_all_mass_vetoes_MC.std_dev*100:.2f}$"
        overall_table.loc["BDT", "Efficiency (\\%)"] = f"${eps_bdt_MC.nominal_value*100:.2f} \\pm {eps_bdt_MC.std_dev*100:.2f}$"
        overall_table.loc["Overall", "Efficiency (\\%)"] = f"${eps_overall_MC.nominal_value*100:.6f} \\pm {eps_overall_MC.std_dev*100:.6f}$"
   
    elif(species == 100):
        up_overall_BuDDKp = ROOT.TEfficiency.Wilson(N_gen[0], N_best_cand_BuDDKp, 0.68, True)
        down_overall_BuDDKp = ROOT.TEfficiency.Wilson(N_gen[0], N_best_cand_BuDDKp, 0.68, False)
        eps_overall_BuDDKp_err = 0.5*(up_overall_BuDDKp - down_overall_BuDDKp)
        eps_overall_BuDDKp = ufloat(N_best_cand_BuDDKp/N_gen[0], eps_overall_BuDDKp_err)
        eps_overall_BuDDKp *= eps_acc[0]*eps_strip[0]*eps_reco1[0]

        up_overall_BdDDKp = ROOT.TEfficiency.Wilson(N_gen[1], N_best_cand_BdDDKp, 0.68, True)
        down_overall_BdDDKp = ROOT.TEfficiency.Wilson(N_gen[1], N_best_cand_BdDDKp, 0.68, False)
        eps_overall_BdDDKp_err = 0.5*(up_overall_BdDDKp - down_overall_BdDDKp)
        eps_overall_BdDDKp = ufloat(N_best_cand_BdDDKp/N_gen[1], eps_overall_BdDDKp_err)
        eps_overall_BdDDKp *= eps_acc[1]*eps_strip[1]*eps_reco1[1]

        up_overall_BsDDKp = ROOT.TEfficiency.Wilson(N_gen[2], N_best_cand_BsDDKp, 0.68, True)
        down_overall_BsDDKp = ROOT.TEfficiency.Wilson(N_gen[2], N_best_cand_BsDDKp, 0.68, False)
        eps_overall_BsDDKp_err = 0.5*(up_overall_BsDDKp - down_overall_BsDDKp)
        eps_overall_BsDDKp = ufloat(N_best_cand_BsDDKp/N_gen[2], eps_overall_BsDDKp_err)
        eps_overall_BsDDKp *= eps_acc[2]*eps_strip[2]*eps_reco1[2]

        up_overall_BuDDK0 = ROOT.TEfficiency.Wilson(N_gen[3], N_best_cand_BuDDK0, 0.68, True)
        down_overall_BuDDK0 = ROOT.TEfficiency.Wilson(N_gen[3], N_best_cand_BuDDK0, 0.68, False)
        eps_overall_BuDDK0_err = 0.5*(up_overall_BuDDK0 - down_overall_BuDDK0)
        eps_overall_BuDDK0 = ufloat(N_best_cand_BuDDK0/N_gen[3], eps_overall_BuDDK0_err)
        eps_overall_BuDDK0 *= eps_acc[3]*eps_strip[3]*eps_reco1[3]

        up_overall_BuDD = ROOT.TEfficiency.Wilson(N_gen[4], N_best_cand_BuDD, 0.68, True)
        down_overall_BuDD = ROOT.TEfficiency.Wilson(N_gen[4], N_best_cand_BuDD, 0.68, False)
        eps_overall_BuDD_err = 0.5*(up_overall_BuDD - down_overall_BuDD)
        eps_overall_BuDD = ufloat(N_best_cand_BuDD/N_gen[4], eps_overall_BuDD_err)
        eps_overall_BuDD *= eps_acc[4]*eps_strip[4]*eps_reco1[4]

        overall_columns = ["Overall efficiency (\\%)"]
        overall_rows = ["$B^+ \\to D D K^+$", "$B^0 \\to D D K^+$", "$B^0_s \\to D D K^+$", "$B^+ \\to D D K^0$", "$B^+ \\to D D$"]

        overall_table = pd.DataFrame(index=overall_rows, columns=overall_columns)
        overall_table.loc["$B^+ \\to D D K^+$", "Overall efficiency (\\%)"] = f"${eps_overall_BuDDKp.nominal_value*100:.8f} \\pm {eps_overall_BuDDKp.std_dev*100:.8f}$"
        overall_table.loc["$B^0 \\to D D K^+$", "Overall efficiency (\\%)"] = f"${eps_overall_BdDDKp.nominal_value*100:.8f} \\pm {eps_overall_BdDDKp.std_dev*100:.8f}$"
        overall_table.loc["$B^0_s \\to D D K^+$", "Overall efficiency (\\%)"] = f"${eps_overall_BsDDKp.nominal_value*100:.8f} \\pm {eps_overall_BsDDKp.std_dev*100:.8f}$"
        overall_table.loc["$B^+ \\to D D K^0$", "Overall efficiency (\\%)"] = f"${eps_overall_BuDDK0.nominal_value*100:.8f} \\pm {eps_overall_BuDDK0.std_dev*100:.8f}$"
        overall_table.loc["$B^+ \\to D D$", "Overall efficiency (\\%)"] = f"${eps_overall_BuDD.nominal_value*100:.8f} \\pm {eps_overall_BuDD.std_dev*100:.8f}$"

    with open(f'/panfs/felician/B2Ktautau/workflow/selections_efficiency_tables/Species_{species}/overall_efficiency_table.tex', 'w') as fout_overall:
        fout_overall.write(overall_table.to_latex())



def main(argv):
    species = argv[1]
    create_table = argv[2]

    species = int(species)
    create_table = create_table == "True"

    is_KtautauMC = False
    if((species == 1) or (species == 10)):
        is_KtautauMC = True
    
    is_KtautauData = False
    if((species == 2) or (species == 3)):
        is_KtautauData = True
    
    is_cocktailMC = False
    if((species == 100) or (species == 110) or (species == 120) or (species == 130) or (species == 140) or (species == 150) or (species == 160) or (species == 170) or (species == 1000) or (species == 1100) or (species == 1200) or (species == 1300) or (species == 1500)):
        is_cocktailMC = True

    ########################################### Definition of analysys selections ##################################################
    # Truth-match cuts
    truthMatch = ""
    if(species == 1): # B+ -> K+ tau+ tau- TM cuts
        truthMatch = "(abs(Kp_TRUEID) == 321) && (abs(taup_pi1_TRUEID) == 211) && (abs(taup_pi2_TRUEID) == 211) && (abs(taup_pi3_TRUEID) == 211) && (abs(taum_pi1_TRUEID) == 211) && (abs(taum_pi2_TRUEID) == 211) && (abs(taum_pi3_TRUEID) == 211) && (abs(taup_TRUEID) == 15) && (abs(taum_TRUEID) == 15) && (abs(Bp_TRUEID) == 521)"
    if(species == 4): # B+ -> D+ D- K+ TM cuts
        truthMatch = "(abs(Kp_TRUEID) == 321) && (abs(Dp_K_TRUEID) == 321) && (abs(Dp_pi1_TRUEID) == 211) && (abs(Dp_pi2_TRUEID) == 211) && (abs(Dm_K_TRUEID) == 321) && (abs(Dm_pi1_TRUEID) == 211) && (abs(Dm_pi2_TRUEID) == 211) && (abs(Dp_TRUEID) == 411) && (abs(Dm_TRUEID) == 411) && (abs(Bp_TRUEID) == 521)"
    if((species == 7) or (species == 71)): # B+ -> D0bar Ds+ TM cuts
        truthMatch = "(abs(D0bar_K_TRUEID) == 321) && (abs(Dsp_K1_TRUEID) == 321) && (abs(Dsp_K2_TRUEID) == 321) && (abs(D0bar_pi_TRUEID) == 211) && (abs(Dsp_pi_TRUEID) == 211) && (abs(D0bar_TRUEID) == 421) && (abs(Dsp_TRUEID) == 431) && (abs(Bp_TRUEID) == 521)"
    if(species == 9): # B+ -> D0 D0 K+ TM cuts
        truthMatch = "(abs(D0bar_K_TRUEID) == 321) && (abs(D0bar_pi1_TRUEID) == 211) && (abs(D0bar_pi2_TRUEID) == 211) && (abs(D0bar_pi3_TRUEID) == 211) && (abs(D0_K_TRUEID) == 321) && (abs(D0_pi1_TRUEID) == 211) && (abs(D0_pi2_TRUEID) == 211) && (abs(D0_pi3_TRUEID) == 211) && (abs(D0bar_TRUEID) == 421) && (abs(D0_TRUEID) == 421) && (abs(Kp_TRUEID) == 321) && (abs(Bp_TRUEID) == 521)"
    if(species == 1000):
        truthMatch = truthMatch_cocktailMC("BuD0D0Kp")+" || "+truthMatch_cocktailMC("BuD0starD0Kp")+" || "+truthMatch_cocktailMC("BuD0D0starKp")+" || "+truthMatch_cocktailMC("BuD0starD0starKp")+" || "+truthMatch_cocktailMC("BuDpDmKp")+" || "+truthMatch_cocktailMC("BuDpstarDmKp")+" || "+truthMatch_cocktailMC("BuDpDmstarKp")+" || "+truthMatch_cocktailMC("BuDpstarDmstarKp")+" || "+truthMatch_cocktailMC("BuDsDsKp")+" || "+truthMatch_cocktailMC("BuDsstarDsKp")+" || "+truthMatch_cocktailMC("BuDsDsstarKp")+" || "+truthMatch_cocktailMC("BuDsstarDsstarKp")
    if(species == 1100):
        truthMatch = truthMatch_cocktailMC("BdDmD0Kp")+" || "+truthMatch_cocktailMC("BdDmstarD0Kp")+" || "+truthMatch_cocktailMC("BdDmD0starKp")+" || "+truthMatch_cocktailMC("BdDmstarD0starKp")
    if(species == 1200):
        truthMatch = truthMatch_cocktailMC("BsDsD0Kp")+" || "+truthMatch_cocktailMC("BsDsstarD0Kp")+" || "+truthMatch_cocktailMC("BsDsD0starKp")+" || "+truthMatch_cocktailMC("BsDsstarD0starKp")
    if(species == 1300):
        truthMatch = truthMatch_cocktailMC("BuD0DpK0")+" || "+truthMatch_cocktailMC("BuD0starDpK0")+" || "+truthMatch_cocktailMC("BuD0DpstarK0")+" || "+truthMatch_cocktailMC("BuD0starDpstarK0")
    if(species == 1500):
        truthMatch = truthMatch_cocktailMC("BuD0Ds")+" || "+truthMatch_cocktailMC("BuD0starDs")+" || "+truthMatch_cocktailMC("BuD0Dsstar")+" || "+truthMatch_cocktailMC("BuD0starDsstar")+" || "+truthMatch_cocktailMC("BuD0Dp")+" || "+truthMatch_cocktailMC("BuD0starDp")+" || "+truthMatch_cocktailMC("BuD0Dpstar")+" || "+truthMatch_cocktailMC("BuD0starDpstar")

    # Trigger
    L0_trigger = "( (Bp_L0HadronDecision_TOS==1) || ((Bp_L0HadronDecision_TIS==1) || (Bp_L0MuonDecision_TIS==1) || (Bp_L0ElectronDecision_TIS==1) || (Bp_L0PhotonDecision_TIS==1)) )"
    HLT1_trigger = "( (Bp_Hlt1TrackMVADecision_TOS==1) || (Bp_Hlt1TwoTrackMVADecision_TOS==1) )"
    HLT2_trigger = "( (Bp_Hlt2Topo2BodyDecision_TOS==1) || (Bp_Hlt2Topo3BodyDecision_TOS==1) || (Bp_Hlt2Topo4BodyDecision_TOS==1) )"
    trigger = L0_trigger+" && "+HLT1_trigger+" && "+HLT2_trigger

    # Rectangular cuts
    rectangular_cuts = []
    rectangular_cuts_names = []

    if(is_KtautauMC or is_KtautauData or is_cocktailMC):
        rectangular_cuts.append("(taup_M > 750) && (taup_M < 1650)")
        rectangular_cuts.append("(taum_M > 750) && (taum_M < 1650)")
        rectangular_cuts.append("(Bp_VTXISODCHI2MASSONETRACK_B > 3600)")
        rectangular_cuts.append("(Bp_VTXISOBDTHARDFIRSTVALUE_B < 0)")
        rectangular_cuts.append("(Bp_BPVVD > 4)")
        rectangular_cuts.append("(TMath::Max(Bp_B2Ksttautau_ISOBDTTHIRDVALUE_taup,Bp_B2Ksttautau_ISOBDTTHIRDVALUE_taum) > -0.1)")
        rectangular_cuts.append("(TMath::Min(TMath::Log10(1 - TMath::Abs(taup_DIRA_ORIVX)) * TMath::Sign(1.0, taup_DIRA_ORIVX), TMath::Log10(1 - TMath::Abs(taum_DIRA_ORIVX)) * TMath::Sign(1.0, taum_DIRA_ORIVX)) < -1)")

        if(is_KtautauMC or is_cocktailMC): # stripping PID cuts (re-applied)
            rectangular_cuts.append("(Kp_ProbNNk_pidgen_default > 0.2)")
            rectangular_cuts.append("(taup_pi1_ProbNNpi_pidgen_default > 0.55) && (taup_pi2_ProbNNpi_pidgen_default > 0.55) && (taup_pi3_ProbNNpi_pidgen_default > 0.55) && (taum_pi1_ProbNNpi_pidgen_default > 0.55) && (taum_pi2_ProbNNpi_pidgen_default > 0.55) && (taum_pi3_ProbNNpi_pidgen_default > 0.55)")

        rectangular_cuts_names = ["$m_{\\tau^+} \\in ]750,1650[$\\, MeV", "$m_{\\tau^-} \\in ]750,1650[$\\, MeV", "$m_{B+1,track}^{vtx. iso} > 3600$\\, MeV", "$BDT_{1st}^{vtx.iso} < 0$", "$B^+$ FD from PV $> 4$\\,mm", "max($BDT_{3rd}^{tau. iso}$)$_{\\tau^+,\\tau^-}$ $> -0.1$", "min($\log10(1-|DIRA\\_BV|)\\times \\text{sign}(DIRA\\_BV)$)$_{\\tau^+,\\tau^-} < -1$", "Kp\\_ProbNNk $> 0.2$", "Tau\\_pions\\_ProbNNpi $> 0.55$"]

    if((species == 7) or (species == 71) or (species == 72) or  (species == 8)):
        if(species != 71):
            rectangular_cuts.append("(Bp_BPVVD > 4) && (Bp_BPVVD < 80)")
            rectangular_cuts.append("(D0bar_M > 1820) && (D0bar_M < 1910)")
            rectangular_cuts.append("(Dsp_M > 1930) && (Dsp_M < 2000)")

        if((species == 7) or (species == 72) or (species == 71)):
            if(species != 71):
                # offline cuts
                rectangular_cuts.append("(D0bar_K_ProbNNk_pidgen_default > 0.55)")
                rectangular_cuts.append("(Dsp_K1_ProbNNk_pidgen_default > 0.55)")
                rectangular_cuts.append("(Dsp_K2_ProbNNk_pidgen_default > 0.55)")

            # stripping PID cuts (re-applied)
            rectangular_cuts.append("(D0bar_pi_ProbNNpi_pidgen_default > 0.55)")
            rectangular_cuts.append("(Dsp_pi_ProbNNpi_pidgen_default > 0.55)")
        else:
            # offline cuts
            rectangular_cuts.append("(D0bar_K_ProbNNk > 0.55)")
            rectangular_cuts.append("(Dsp_K1_ProbNNk > 0.55)")
            rectangular_cuts.append("(Dsp_K2_ProbNNk > 0.55)")
        
        rectangular_cuts_names = ["$B^+$ FD from PV $\\in ]4,80[$\\,mm", "$m_{\\bar{D}^0} \\in ]1820, 1910[$\\, MeV", "$m_{D_s^+} \\in ]1930, 2000[$\\, MeV", "D0bar\\_K\\_ProbNNk $> 0.55$", "Dsp\\_K1\\_ProbNNk $> 0.55$", "Dsp\\_K2\\_ProbNNk $>0.55$", "D0bar\\_pi\\_ProbNNpi $> 0.55$", "Dsp\\_pi\\_ProbNNpi $> 0.55$"]

    if((species == 4) or (species == 5) or (species == 6)): 
        rectangular_cuts.append("(Dp_K_ProbNNk > 0.1) && (Dm_K_ProbNNk > 0.1)")
        rectangular_cuts.append("(TMath::Min( TMath::Log10(1-TMath::Abs(Dp_DIRA_ORIVX))*TMath::Sign(1,Dp_DIRA_ORIVX), TMath::Log10(1-TMath::Abs(Dm_DIRA_ORIVX))*TMath::Sign(1,Dm_DIRA_ORIVX) ) < -2.)")
        rectangular_cuts.append("(TMath::Log10(TMath::Abs(1-Bp_DIRA_OWNPV))*TMath::Sign(1,Bp_DIRA_OWNPV) < -4)")
        rectangular_cuts.append("(Bp_ENDVERTEX_CHI2 < 40)")
        rectangular_cuts.append("(abs(Dp_M-1869.66) < 50) && (abs(Dm_M < 1869.66) < 50)")
        rectangular_cuts.append("(Bp_FD_OWNPV > 2)")

    if((species == 9) or (species == 0) or (species == -1)):
        rectangular_cuts.append("(Kp_PT > 1000)")
        rectangular_cuts.append("(Kp_ProbNNk > 0.55)")
        rectangular_cuts.append("(D0bar_K_ProbNNk > 0.55)")
        rectangular_cuts.append("(D0_K_ProbNNk > 0.55)")
        rectangular_cuts.append("(D0bar_M > 1840) && (D0bar_M < 1890)")
        rectangular_cuts.append("(D0_M > 1840) && (D0_M < 1890)")

    all_rectangular_cuts = ""
    N = len(rectangular_cuts)
    for i in range(N):
        if i == N-1:
            all_rectangular_cuts += rectangular_cuts[i]
        else:
            all_rectangular_cuts += rectangular_cuts[i]+" && "

    # Pass mass reconstruction
    pass_mass_fit = ""
    if(is_KtautauMC or is_KtautauData or is_cocktailMC):
        pass_mass_fit = "(df_status == 0)"
    else:
        pass_mass_fit = "(Bp_dtf_status[0]==0)"
    
    # Fit range
    fit_range = ""
    if(is_KtautauMC or is_KtautauData or is_cocktailMC):
        fit_range = "(df_Bp_M > 4000) && (df_Bp_M < 9000)"
    if((species == 7) or (species == 71) or (species == 72) or  (species == 8) or (species == 81)):
        fit_range = "(Bp_dtf_M[0] > 5235) && (Bp_dtf_M[0] < 5355)"
    if((species == 4) or (species == 5) or (species == 6)):
        fit_range = "(Bp_dtf_M[0] > 5180) && (Bp_dtf_M[0] < 5600)"

    # Mass vetoes
    mass_vetoes_cuts = []
    mass_vetoes = ""
    if(is_KtautauMC or is_KtautauData or is_cocktailMC):
        mass_vetoes_cuts.append("(TMath::Abs(Bp_M02-1864.84) > 30) && (TMath::Abs(Bp_M04-1864.84) > 30) && (TMath::Abs(Bp_M06-1864.84) > 30)")
        mass_vetoes_cuts.append("(TMath::Abs(Bp_M046-1869.66) > 20)")
        mass_vetoes_cuts.append("(TMath::Abs(Bp_M0456-1864.84) > 30)")
        mass_vetoes_cuts.append("(TMath::Abs(Bp_M01246-2010.26) > 20) && (TMath::Abs(Bp_M02346-2010.26) > 20) && (TMath::Abs(Bp_M02456-2010.26) > 20)")

        Cx_taup = "(df_DV1y-df_BVy)*df_Kp_PZ - (df_DV1z-df_BVz)*df_Kp_PY"
        Cy_taup = "(df_DV1z-df_BVz)*df_Kp_PX - (df_DV1x-df_BVx)*df_Kp_PZ"
        Cz_taup = "(df_DV1x-df_BVx)*df_Kp_PY - (df_DV1y-df_BVy)*df_Kp_PX"
        C_taup = "TMath::Sqrt( pow("+Cx_taup+",2) + pow("+Cy_taup+",2) + pow("+Cz_taup+",2) )"
        IP_taup = "2*"+C_taup+"/(TMath::Sqrt( pow(df_Kp_PX,2) + pow(df_Kp_PY,2) + pow(df_Kp_PZ,2) ))"
        
        Cx_taum = "(df_DV2y-df_BVy)*df_Kp_PZ - (df_DV2z-df_BVz)*df_Kp_PY"
        Cy_taum = "(df_DV2z-df_BVz)*df_Kp_PX - (df_DV2x-df_BVx)*df_Kp_PZ"
        Cz_taum = "(df_DV2x-df_BVx)*df_Kp_PY - (df_DV2y-df_BVy)*df_Kp_PX"
        C_taum = "TMath::Sqrt( pow("+Cx_taum+",2) + pow("+Cy_taum+",2) + pow("+Cz_taum+",2) )"
        IP_taum = "2*"+C_taum+"/(TMath::Sqrt( pow(df_Kp_PX,2) + pow(df_Kp_PY,2) + pow(df_Kp_PZ,2) ))"

        IP_taup_cut = "("+IP_taup+" > 0.25) || (TMath::Abs(Bp_M02-892) > 50)"
        IP_taum_cut = "("+IP_taum+" > 0.25) || ((TMath::Abs(Bp_M04-892) > 50) && (TMath::Abs(Bp_M06-892) > 50))"
        
        IP_2_particles_cut = "("+IP_taup_cut+") &&  ("+IP_taum_cut+")"
        IP_4_particles_cut = "( (("+IP_taup+" > 0.15) && ("+IP_taum+" > 0.15)) || ((TMath::Abs(Bp_M0124-1864.84) > 30) && (TMath::Abs(Bp_M0126-1864.84) > 30) && (TMath::Abs(Bp_M0146-1864.84) > 30) && (TMath::Abs(Bp_M0234-1864.84) > 30) && (TMath::Abs(Bp_M0236-1864.84) > 30) && (TMath::Abs(Bp_M0245-1864.84) > 30) && (TMath::Abs(Bp_M0256-1864.84) > 30)) )"

        mass_vetoes_cuts.append(IP_2_particles_cut)
        mass_vetoes_cuts.append(IP_4_particles_cut)

        N1 = len(mass_vetoes_cuts)
        for i in range(N1):
            if i == N1-1:
                mass_vetoes += mass_vetoes_cuts[i]
            else:
                mass_vetoes += mass_vetoes_cuts[i]+" && "
        
    # Best candidate selection
    best_cand = ""
    if((species != 1) and (species != 7)):
        best_cand = "is_best_cand == 1"

    if(not create_table):
        # save selections in root file
        f = ROOT.TFile(f"/panfs/felician/B2Ktautau/workflow/selections_efficiency_tables/Species_{species}/selections.root", "RECREATE")
        f.cd()
        truthMatch_expression = ROOT.TObjString(truthMatch)
        trigger_expression = ROOT.TObjString(trigger)
        rectangular_cuts_expression = ROOT.TObjString(all_rectangular_cuts)
        pass_mass_fit_expression = ROOT.TObjString(pass_mass_fit)
        fit_range_expression = ROOT.TObjString(fit_range)
        mass_vetoes_expression = ROOT.TObjString(mass_vetoes)
        best_cand_expression = ROOT.TObjString(best_cand)

        truthMatch_expression.Write("truthMatch")
        trigger_expression.Write("trigger")
        rectangular_cuts_expression.Write("rectangular_cuts")
        pass_mass_fit_expression.Write("pass_mass_fit")
        fit_range_expression.Write("fit_range")
        mass_vetoes_expression.Write("mass_vetoes")
        best_cand_expression.Write("best_candidate")
        f.Close()
    else:
        if((species != 1) and (species != 7) and (species != 100)):
            print("Wrong species for creating tables. Use 1 for Ktautau, 7 for DDs or 100 for cocktail MCs.")
            quit()
        else:
            create_tables(species, truthMatch, L0_trigger, HLT1_trigger, HLT2_trigger, trigger, rectangular_cuts_names, rectangular_cuts, pass_mass_fit, fit_range, mass_vetoes, mass_vetoes_cuts, best_cand)

if __name__ == "__main__":
    main(sys.argv)