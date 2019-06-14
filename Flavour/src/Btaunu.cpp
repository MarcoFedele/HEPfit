/*
 * Copyright (C) 2015 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "Btaunu.h"
#include "StandardModel.h"
#include "std_make_vector.h"

Btaunu::Btaunu(const StandardModel& SM_i, QCD::meson meson_i)
: ThObservable(SM_i)
{
    meson = meson_i;
    SM.initializeMeson(meson);
    
    ChiralBasisflag = SM.getFlavour().getFlagChiralBasis();
    
    setParametersForObservable(make_vector<std::string>()
            << "regV" << "regA" << "regP" << "imgV" << "imgA" << "imgP"
            << "regVL" << "regSL" << "regSR" << "imgVL" << "imgSL" << "imgSR");
};

double Btaunu::computeThValue()
{
    //gslpp::vector<gslpp::complex> ** allcoeff = SM.getFlavour().ComputeCoeffbtaunu(meson);
    double mtau = SM.getLeptons(StandardModel::TAU).getMass();
    double mB = SM.getMesons(meson).getMass();
    double mb = SM.getQuarks(QCD::BOTTOM).getMass();
    double fact = 1.; /*factor introduced to scale the decay constant from that of the neutral B to the charged B.*/
    //double fact = 0.989;
    
    double Vckm2;
    
    if (meson == QCD::B_P) {
        Vckm2 = SM.getCKM().getV_ub().abs2();
        gV = 0.;
        gA = 0.;
        gP = 0.;
    }
    else if (meson == QCD::B_C) {
        Vckm2 = SM.getCKM().getV_cb().abs2();
        if (ChiralBasisflag) {
            gV = SM.getOptionalParameter("regVL") + gslpp::complex::i()*SM.getOptionalParameter("imgVL");
            gA = - SM.getOptionalParameter("regVL") - gslpp::complex::i()*SM.getOptionalParameter("imgVL");
            gP = SM.getOptionalParameter("regSR") + gslpp::complex::i()*SM.getOptionalParameter("imgSR")
                - SM.getOptionalParameter("regSL") - gslpp::complex::i()*SM.getOptionalParameter("imgSL");
        }
        else {
            gV = SM.getOptionalParameter("regV") + gslpp::complex::i()*SM.getOptionalParameter("imgV");
            gA = SM.getOptionalParameter("regA") + gslpp::complex::i()*SM.getOptionalParameter("imgA");
            gP = SM.getOptionalParameter("regP") + gslpp::complex::i()*SM.getOptionalParameter("imgP");
        }
    }
    
    return pow(fact * SM.getMesons(meson).getDecayconst() * SM.getGF(), 2.) * Vckm2 * mB / (8. * M_PI) * 
            mtau * mtau * pow(1. - mtau * mtau / mB / mB, 2.) / SM.getMesons(meson).computeWidth() * 
            (1. + (gV - gA)/2. + mB * mB / mb / mtau * gP).abs2(); // PLEASE NOTE THE DECAY CONST
    
    //return 1. / (64. * M_PI) * mtau * mtau * pow(fact * SM.getMesons(meson).getDecayconst(), 2.) * mB * pow(1. - mtau * mtau / mB / mB, 2.) / SM.getMesons(meson).computeWidth() * ((*(allcoeff[LO]))(0) 
    //           + mB * mB / mb / mtau * ((*(allcoeff[LO]))(1) + (*(allcoeff[LO]))(2))).abs2(); // PLEASE NOTE THE DECAY CONST
}