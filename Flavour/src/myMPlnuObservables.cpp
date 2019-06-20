/*
 * Copyright (C) 2018 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "myMPlnuObservables.h"
#include "myMPlnu.h"
#include "StandardModel.h"


/**
* Gamma
**/
Gamma_myMPlnu::Gamma_myMPlnu(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson pseudoscalar_i, QCD::lepton lep_i)
: ThObservable(SM_i)
{
    lep = lep_i;
    meson = meson_i;
    pseudoscalar = pseudoscalar_i;

    setParametersForObservable(SM.getFlavour().getmyMPlnu(meson, pseudoscalar, lep).initializemyMPlnuParameters());
}

double Gamma_myMPlnu::computeGamma(QCD::lepton lep)
{
    QCD::lepton lep_i = lep;

    return 2. * ( SM.getFlavour().getmyMPlnu(meson, pseudoscalar, lep_i).integrateI(0)
                + 1./3. * SM.getFlavour().getmyMPlnu(meson, pseudoscalar, lep_i).integrateI(2) );
}

double Gamma_myMPlnu::computeThValue()
{
    return computeGamma(lep);
}


/**
* A_FB
**/
A_FB_myMPlnu::A_FB_myMPlnu(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson pseudoscalar_i, QCD::lepton lep_i)
: Gamma_myMPlnu(SM_i, meson_i, pseudoscalar_i, lep_i)
{
    lep = lep_i;
    meson = meson_i;
    pseudoscalar = pseudoscalar_i;

    setParametersForObservable(SM.getFlavour().getmyMPlnu(meson, pseudoscalar, lep).initializemyMPlnuParameters());
}

double A_FB_myMPlnu::computeThValue()
{
    return 1. / computeGamma(lep) * SM.getFlavour().getmyMPlnu(meson, pseudoscalar, lep).integrateI(1);
}


/**
* A_lam
**/
A_lam_myMPlnu::A_lam_myMPlnu(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson pseudoscalar_i, QCD::lepton lep_i)
: Gamma_myMPlnu(SM_i, meson_i, pseudoscalar_i, lep_i)
{
    lep = lep_i;
    meson = meson_i;
    pseudoscalar = pseudoscalar_i;

    setParametersForObservable(SM.getFlavour().getmyMPlnu(meson, pseudoscalar, lep).initializemyMPlnuParameters());
}

double A_lam_myMPlnu::computeThValue()
{
    return 1. - 2./computeGamma(lep) * SM.getFlavour().getmyMPlnu(meson, pseudoscalar, lep).integrateI(3);
}



/*********************
 *                   *
 *    LFUV RATIOS    *
 *                   *
 *********************/



/**
* R_D
**/
R_D::R_D(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson pseudoscalar_i,
            QCD::lepton lep_1, QCD::lepton lep_2, QCD::lepton lep_3)
: Gamma_myMPlnu(SM_i, meson_i, pseudoscalar_i, lep_1)
{
    lep1 = lep_1;
    lep2 = lep_2;
    lep3 = lep_3;
    meson = meson_i;
    pseudoscalar = pseudoscalar_i;

    setParametersForObservable(SM.getFlavour().getmyMPlnu(meson, pseudoscalar, QCD::lepton::TAU).initializemyMPlnuParameters());
}

double R_D::computeThValue()
{
    return computeGamma(lep1) / (1./2. * (computeGamma(lep2) + computeGamma(lep3)) );
}


/**
* R_A_lam
**/
R_A_lam_myMPlnu::R_A_lam_myMPlnu(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson pseudoscalar_i,
            QCD::lepton lep_1, QCD::lepton lep_2, QCD::lepton lep_3)
: Gamma_myMPlnu(SM_i, meson_i, pseudoscalar_i, lep_1)
{
    lep1 = lep_1;
    lep2 = lep_2;
    lep3 = lep_3;
    meson = meson_i;
    pseudoscalar = pseudoscalar_i;

    setParametersForObservable(SM.getFlavour().getmyMPlnu(meson, pseudoscalar, QCD::lepton::TAU).initializemyMPlnuParameters());
}

double R_A_lam_myMPlnu::computeThValue()
{
    double Alamtau = 1. - 2./computeGamma(lep1) * SM.getFlavour().getmyMPlnu(meson, pseudoscalar, lep1).integrateI(3)
    double Alammu = 1. - 2./computeGamma(lep2) * SM.getFlavour().getmyMPlnu(meson, pseudoscalar, lep2).integrateI(3)
    double Alame = 1. - 2./computeGamma(lep3) * SM.getFlavour().getmyMPlnu(meson, pseudoscalar, lep3).integrateI(3)

    double res = 2. * Alamtau / (Alammu + Alame);

    return res;
}


/**
* R_A_FB
**/
R_A_FB_myMPlnu::R_A_FB_myMPlnu(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson pseudoscalar_i,
            QCD::lepton lep_1, QCD::lepton lep_2, QCD::lepton lep_3)
: Gamma_myMPlnu(SM_i, meson_i, pseudoscalar_i, lep_1)
{
    lep1 = lep_1;
    lep2 = lep_2;
    lep3 = lep_3;
    meson = meson_i;
    pseudoscalar = pseudoscalar_i;

    setParametersForObservable(SM.getFlavour().getmyMPlnu(meson, pseudoscalar, QCD::lepton::TAU).initializemyMPlnuParameters());
}

double R_A_FB_myMPlnu::computeThValue()
{
    double Mlep12 = SM.getLeptons(lep1).getMass() * SM.getLeptons(lep1).getMass();
    double Mlep22 = SM.getLeptons(lep2).getMass() * SM.getLeptons(lep2).getMass();
    double Mlep32 = SM.getLeptons(lep3).getMass() * SM.getLeptons(lep3).getMass();

    double AFBtau = SM.getFlavour().getmyMPlnu(meson, pseudoscalar, lep1).integrateI(1) / computeGamma(lep1);
    double AFBmu = SM.getFlavour().getmyMPlnu(meson, pseudoscalar, lep2).integrateI(1) / computeGamma(lep2);
    double AFBe = SM.getFlavour().getmyMPlnu(meson, pseudoscalar, lep3).integrateI(1) / computeGamma(lep3);

    double res = 2. * AFBtau/Mlep12 / (AFBmu/Mlep22 + AFBe/Mlep32);

    return res;
}
