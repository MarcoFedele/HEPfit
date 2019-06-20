/*
 * Copyright (C) 2018 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "myMVlnuObservables.h"
#include "myMVlnu.h"
#include "StandardModel.h"


/**
* Gamma
**/
Gamma_myMVlnu::Gamma_myMVlnu(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i, QCD::lepton lep_i)
: ThObservable(SM_i)
{
    lep = lep_i;
    meson = meson_i;
    vectorM = vector_i;

    setParametersForObservable(SM.getFlavour().getmyMVlnu(meson, vectorM, lep).initializemyMVlnuParameters());
}

double Gamma_myMVlnu::computeGamma(QCD::lepton lep)
{
    QCD::lepton lep_i = lep;

    return 3./4. * SM.getFlavour().getmyMVlnu(meson, vectorM, lep_i).integrateI(0)
            + 3./2. * SM.getFlavour().getmyMVlnu(meson, vectorM, lep_i).integrateI(1)
            - 1./4. * SM.getFlavour().getmyMVlnu(meson, vectorM, lep_i).integrateI(2)
            - 1./2. * SM.getFlavour().getmyMVlnu(meson, vectorM, lep_i).integrateI(3);
}

double Gamma_myMVlnu::computeThValue()
{
    return computeGamma(lep);
}


/**
* A_lam
**/
A_lam_myMVlnu::A_lam_myMVlnu(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i, QCD::lepton lep_i)
: Gamma_myMVlnu(SM_i, meson_i, vector_i, lep_i)
{
    lep = lep_i;
    meson = meson_i;
    vectorM = vector_i;

    setParametersForObservable(SM.getFlavour().getmyMVlnu(meson, vectorM, lep).initializemyMVlnuParameters());
}

double A_lam_myMVlnu::computeThValue()
{
    return 1. - 2. / computeGamma(lep) *
            ( 1./4. * SM.getFlavour().getmyMVlnu(meson, vectorM, lep).integrateI(0)
            + 1./2. * SM.getFlavour().getmyMVlnu(meson, vectorM, lep).integrateI(1)
            + 1./4. * SM.getFlavour().getmyMVlnu(meson, vectorM, lep).integrateI(2)
            - 3./2. * SM.getFlavour().getmyMVlnu(meson, vectorM, lep).integrateI(3)
            + SM.getFlavour().getmyMVlnu(meson, vectorM, lep).integrateI(12) );
}


/**
* R_LT
**/
R_LT_myMVlnu::R_LT_myMVlnu(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i, QCD::lepton lep_i)
: ThObservable(SM_i)
{
    lep = lep_i;
    meson = meson_i;
    vectorM = vector_i;

    setParametersForObservable(SM.getFlavour().getmyMVlnu(meson, vectorM, lep).initializemyMVlnuParameters());
}

double R_LT_myMVlnu::computeThValue()
{
    return ( 3./4. * SM.getFlavour().getmyMVlnu(meson, vectorM, lep).integrateI(0)
            - 1./4. * SM.getFlavour().getmyMVlnu(meson, vectorM, lep).integrateI(2) ) /
           ( 3./2. * SM.getFlavour().getmyMVlnu(meson, vectorM, lep).integrateI(1)
            - 1./2. * SM.getFlavour().getmyMVlnu(meson, vectorM, lep).integrateI(3) );
}


/**
* R_AB
**/
R_AB_myMVlnu::R_AB_myMVlnu(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i, QCD::lepton lep_i)
: ThObservable(SM_i)
{
    lep = lep_i;
    meson = meson_i;
    vectorM = vector_i;

    setParametersForObservable(SM.getFlavour().getmyMVlnu(meson, vectorM, lep).initializemyMVlnuParameters());
}

double R_AB_myMVlnu::computeThValue()
{
    return ( 1./4. * SM.getFlavour().getmyMVlnu(meson, vectorM, lep).integrateI(0)
            + 1./2. * SM.getFlavour().getmyMVlnu(meson, vectorM, lep).integrateI(1)
            - 3./4. * SM.getFlavour().getmyMVlnu(meson, vectorM, lep).integrateI(2)
            - 3./2. * SM.getFlavour().getmyMVlnu(meson, vectorM, lep).integrateI(3) ) /
           ( 1./2. * SM.getFlavour().getmyMVlnu(meson, vectorM, lep).integrateI(0)
            + SM.getFlavour().getmyMVlnu(meson, vectorM, lep).integrateI(1)
            + 1./2. * SM.getFlavour().getmyMVlnu(meson, vectorM, lep).integrateI(2)
            + SM.getFlavour().getmyMVlnu(meson, vectorM, lep).integrateI(3) );
}


/**
* A_FB
**/
A_FB_myMVlnu::A_FB_myMVlnu(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i, QCD::lepton lep_i)
: Gamma_myMVlnu(SM_i, meson_i, vector_i, lep_i)
{
    lep = lep_i;
    meson = meson_i;
    vectorM = vector_i;

    setParametersForObservable(SM.getFlavour().getmyMVlnu(meson, vectorM, lep).initializemyMVlnuParameters());
}

double A_FB_myMVlnu::computeThValue()
{
    return 1. / computeGamma(lep) *
            ( 3./8. * SM.getFlavour().getmyMVlnu(meson, vectorM, lep).integrateI(7)
            + 3./4. * SM.getFlavour().getmyMVlnu(meson, vectorM, lep).integrateI(8) );
}


/**
* A_3
**/
A_3_myMVlnu::A_3_myMVlnu(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i, QCD::lepton lep_i)
: Gamma_myMVlnu(SM_i, meson_i, vector_i, lep_i)
{
    lep = lep_i;
    meson = meson_i;
    vectorM = vector_i;

    setParametersForObservable(SM.getFlavour().getmyMVlnu(meson, vectorM, lep).initializemyMVlnuParameters());
}

double A_3_myMVlnu::computeA3(QCD::lepton lep)
{
    return 1. / computeGamma(lep) * SM.getFlavour().getmyMVlnu(meson, vectorM, lep).integrateI(4);
}

double A_3_myMVlnu::computeThValue()
{
    return computeA3(lep);
}


/**
* A_4
**/
A_4_myMVlnu::A_4_myMVlnu(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i, QCD::lepton lep_i)
: Gamma_myMVlnu(SM_i, meson_i, vector_i, lep_i)
{
    lep = lep_i;
    meson = meson_i;
    vectorM = vector_i;

    setParametersForObservable(SM.getFlavour().getmyMVlnu(meson, vectorM, lep).initializemyMVlnuParameters());
}

double A_4_myMVlnu::computeThValue()
{
    return 1. / computeGamma(lep) * (- 2./M_PI) * SM.getFlavour().getmyMVlnu(meson, vectorM, lep).integrateI(5);
}


/**
* A_5
**/
A_5_myMVlnu::A_5_myMVlnu(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i, QCD::lepton lep_i)
: Gamma_myMVlnu(SM_i, meson_i, vector_i, lep_i)
{
    lep = lep_i;
    meson = meson_i;
    vectorM = vector_i;

    setParametersForObservable(SM.getFlavour().getmyMVlnu(meson, vectorM, lep).initializemyMVlnuParameters());
}

double A_5_myMVlnu::computeThValue()
{
    return 1. / computeGamma(lep) * (- 3./4.) * SM.getFlavour().getmyMVlnu(meson, vectorM, lep).integrateI(6);
}


/**
* A_6
**/
A_6_myMVlnu::A_6_myMVlnu(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i, QCD::lepton lep_i)
: Gamma_myMVlnu(SM_i, meson_i, vector_i, lep_i)
{
    lep = lep_i;
    meson = meson_i;
    vectorM = vector_i;

    setParametersForObservable(SM.getFlavour().getmyMVlnu(meson, vectorM, lep).initializemyMVlnuParameters());
}

double A_6_myMVlnu::computeThValue()
{
    return 1. / computeGamma(lep) * (- 27./8.) * SM.getFlavour().getmyMVlnu(meson, vectorM, lep).integrateI(8);
}


/**
* A_7
**/
A_7_myMVlnu::A_7_myMVlnu(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i, QCD::lepton lep_i)
: Gamma_myMVlnu(SM_i, meson_i, vector_i, lep_i)
{
    lep = lep_i;
    meson = meson_i;
    vectorM = vector_i;

    setParametersForObservable(SM.getFlavour().getmyMVlnu(meson, vectorM, lep).initializemyMVlnuParameters());
}

double A_7_myMVlnu::computeThValue()
{
    return 1. / computeGamma(lep) * (- 3./4.) * SM.getFlavour().getmyMVlnu(meson, vectorM, lep).integrateI(9);
}


/**
* A_8
**/
A_8_myMVlnu::A_8_myMVlnu(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i, QCD::lepton lep_i)
: Gamma_myMVlnu(SM_i, meson_i, vector_i, lep_i)
{
    lep = lep_i;
    meson = meson_i;
    vectorM = vector_i;

    setParametersForObservable(SM.getFlavour().getmyMVlnu(meson, vectorM, lep).initializemyMVlnuParameters());
}

double A_8_myMVlnu::computeThValue()
{
    return 1. / computeGamma(lep) * 2./M_PI * SM.getFlavour().getmyMVlnu(meson, vectorM, lep).integrateI(10);
}


/**
* A_9
**/
A_9_myMVlnu::A_9_myMVlnu(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i, QCD::lepton lep_i)
: Gamma_myMVlnu(SM_i, meson_i, vector_i, lep_i)
{
    lep = lep_i;
    meson = meson_i;
    vectorM = vector_i;

    setParametersForObservable(SM.getFlavour().getmyMVlnu(meson, vectorM, lep).initializemyMVlnuParameters());
}

double A_9_myMVlnu::computeThValue()
{
    return 1. / computeGamma(lep) * SM.getFlavour().getmyMVlnu(meson, vectorM, lep).integrateI(11);
}



/*********************
 *                   *
 *    LFUV RATIOS    *
 *                   *
 *********************/



/**
* R_Dst
**/
R_Dst::R_Dst(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i,
            QCD::lepton lep_1, QCD::lepton lep_2, QCD::lepton lep_3)
: Gamma_myMVlnu(SM_i, meson_i, vector_i, lep_1)
{
    lep1 = lep_1;
    lep2 = lep_2;
    lep3 = lep_3;
    meson = meson_i;
    vectorM = vector_i;

    setParametersForObservable(SM.getFlavour().getmyMVlnu(meson, vectorM, QCD::lepton::TAU).initializemyMVlnuParameters());
}

double R_Dst::computeThValue()
{
    return computeGamma(lep1) / (1./2. * (computeGamma(lep2) + computeGamma(lep3)) );
}


/**
* R_A_lam
**/
R_A_lam_myMVlnu::R_A_lam_myMVlnu(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i,
            QCD::lepton lep_1, QCD::lepton lep_2, QCD::lepton lep_3)
: ThObservable(SM_i)
{
    lep1 = lep_1;
    lep2 = lep_2;
    lep3 = lep_3;
    meson = meson_i;
    vectorM = vector_i;

    setParametersForObservable(SM.getFlavour().getmyMVlnu(meson, vectorM, QCD::lepton::TAU).initializemyMVlnuParameters());
}

double R_A_lam_myMVlnu::computeThValue()
{
    double Alamtau = 1. - 2. / computeGamma(lep1) *
            ( 1./4. * SM.getFlavour().getmyMVlnu(meson, vectorM, lep1).integrateI(0)
            + 1./2. * SM.getFlavour().getmyMVlnu(meson, vectorM, lep1).integrateI(1)
            + 1./4. * SM.getFlavour().getmyMVlnu(meson, vectorM, lep1).integrateI(2)
            - 3./2. * SM.getFlavour().getmyMVlnu(meson, vectorM, lep1).integrateI(3)
            + SM.getFlavour().getmyMVlnu(meson, vectorM, lep1).integrateI(12) );
    double Alammu = 1. - 2. / computeGamma(lep2) *
            ( 1./4. * SM.getFlavour().getmyMVlnu(meson, vectorM, lep2).integrateI(0)
            + 1./2. * SM.getFlavour().getmyMVlnu(meson, vectorM, lep2).integrateI(1)
            + 1./4. * SM.getFlavour().getmyMVlnu(meson, vectorM, lep2).integrateI(2)
            - 3./2. * SM.getFlavour().getmyMVlnu(meson, vectorM, lep2).integrateI(3)
            + SM.getFlavour().getmyMVlnu(meson, vectorM, lep2).integrateI(12) );
    double Alame = 1. - 2. / computeGamma(lep3) *
            ( 1./4. * SM.getFlavour().getmyMVlnu(meson, vectorM, lep3).integrateI(0)
            + 1./2. * SM.getFlavour().getmyMVlnu(meson, vectorM, lep3).integrateI(1)
            + 1./4. * SM.getFlavour().getmyMVlnu(meson, vectorM, lep3).integrateI(2)
            - 3./2. * SM.getFlavour().getmyMVlnu(meson, vectorM, lep3).integrateI(3)
            + SM.getFlavour().getmyMVlnu(meson, vectorM, lep3).integrateI(12) );

    double res = 2. * Alamtau / (Alammu + Alame);

    return res;
}


/**
* R_R_LT
**/
R_R_LT_myMVlnu::R_R_LT_myMVlnu(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i,
            QCD::lepton lep_1, QCD::lepton lep_2, QCD::lepton lep_3)
: R_LT_myMVlnu(SM_i, meson_i, vector_i, lep_1)
{
    lep1 = lep_1;
    lep2 = lep_2;
    lep3 = lep_3;
    meson = meson_i;
    vectorM = vector_i;

    setParametersForObservable(SM.getFlavour().getmyMVlnu(meson, vectorM, QCD::lepton::TAU).initializemyMVlnuParameters());
}

double R_R_LT_myMVlnu::computeThValue()
{
    return R_LT_myMVlnu(SM,meson,vectorM,lep1).computeThValue() /
            (1./2. * (R_LT_myMVlnu(SM,meson,vectorM,lep2).computeThValue()
            + R_LT_myMVlnu(SM,meson,vectorM,lep3).computeThValue()) );
}


/**
* R_R_AB
**/
R_R_AB_myMVlnu::R_R_AB_myMVlnu(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i,
            QCD::lepton lep_1, QCD::lepton lep_2, QCD::lepton lep_3)
: R_AB_myMVlnu(SM_i, meson_i, vector_i, lep_1)
{
    lep1 = lep_1;
    lep2 = lep_2;
    lep3 = lep_3;
    meson = meson_i;
    vectorM = vector_i;

    setParametersForObservable(SM.getFlavour().getmyMVlnu(meson, vectorM, QCD::lepton::TAU).initializemyMVlnuParameters());
}

double R_R_AB_myMVlnu::computeThValue()
{
    return R_AB_myMVlnu(SM,meson,vectorM,lep1).computeThValue() /
            (1./2. * (R_AB_myMVlnu(SM,meson,vectorM,lep2).computeThValue()
            + R_AB_myMVlnu(SM,meson,vectorM,lep3).computeThValue()) );
}


/**
* R_A_FB
**/
R_A_FB_myMVlnu::R_A_FB_myMVlnu(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i,
            QCD::lepton lep_1, QCD::lepton lep_2, QCD::lepton lep_3)
: ThObservable(SM_i)
{
    lep1 = lep_1;
    lep2 = lep_2;
    lep3 = lep_3;
    meson = meson_i;
    vectorM = vector_i;

    setParametersForObservable(SM.getFlavour().getmyMVlnu(meson, vectorM, QCD::lepton::TAU).initializemyMVlnuParameters());
}

double R_A_FB_myMVlnu::computeThValue()
{
    double AFBtau = ( SM.getFlavour().getmyMVlnu(meson, vectorM, lep1).integrateI(7)
            + 2. * SM.getFlavour().getmyMVlnu(meson, vectorM, lep1).integrateI(8) ) / computeGamma(lep1);
    double AFBmu = ( SM.getFlavour().getmyMVlnu(meson, vectorM, lep2).integrateI(7)
            + 2. * SM.getFlavour().getmyMVlnu(meson, vectorM, lep2).integrateI(8) ) / computeGamma(lep2);
    double AFBe = ( SM.getFlavour().getmyMVlnu(meson, vectorM, lep3).integrateI(7)
            + 2. * SM.getFlavour().getmyMVlnu(meson, vectorM, lep3).integrateI(8) ) / computeGamma(lep3);

    double res = 2. * AFBtau / (AFBmu + AFBe);

    return res;
}


/**
* R_A_3
**/
R_A_3_myMVlnu::R_A_3_myMVlnu(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i,
            QCD::lepton lep_1, QCD::lepton lep_2, QCD::lepton lep_3)
: Gamma_myMVlnu(SM_i, meson_i, vector_i, lep_1)
{
    lep1 = lep_1;
    lep2 = lep_2;
    lep3 = lep_3;
    meson = meson_i;
    vectorM = vector_i;

    setParametersForObservable(SM.getFlavour().getmyMVlnu(meson, vectorM, QCD::lepton::TAU).initializemyMVlnuParameters());
}

double R_A_3_myMVlnu::computeThValue()
{
    double A3tau = SM.getFlavour().getmyMVlnu(meson, vectorM, lep1).integrateI(4) / computeGamma(lep1);
    double A3mu = SM.getFlavour().getmyMVlnu(meson, vectorM, lep2).integrateI(4) / computeGamma(lep2);
    double A3e = SM.getFlavour().getmyMVlnu(meson, vectorM, lep3).integrateI(4) / computeGamma(lep3);

    double res = 2. * A3tau / (A3mu + A3e);

    return res;
}


/**
* R_A_4
**/
R_A_4_myMVlnu::R_A_4_myMVlnu(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i,
            QCD::lepton lep_1, QCD::lepton lep_2, QCD::lepton lep_3)
: Gamma_myMVlnu(SM_i, meson_i, vector_i, lep_1)
{
    lep1 = lep_1;
    lep2 = lep_2;
    lep3 = lep_3;
    meson = meson_i;
    vectorM = vector_i;

    setParametersForObservable(SM.getFlavour().getmyMVlnu(meson, vectorM, QCD::lepton::TAU).initializemyMVlnuParameters());
}

double R_A_4_myMVlnu::computeThValue()
{
    double A4tau = SM.getFlavour().getmyMVlnu(meson, vectorM, lep1).integrateI(5) / computeGamma(lep1);
    double A4mu = SM.getFlavour().getmyMVlnu(meson, vectorM, lep2).integrateI(5) / computeGamma(lep2);
    double A4e = SM.getFlavour().getmyMVlnu(meson, vectorM, lep3).integrateI(5) / computeGamma(lep3);

    double res = 2. * A4tau / (A4mu + A4e);

    return res;
}


/**
* R_A_5
**/
R_A_5_myMVlnu::R_A_5_myMVlnu(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i,
            QCD::lepton lep_1, QCD::lepton lep_2, QCD::lepton lep_3)
: Gamma_myMVlnu(SM_i, meson_i, vector_i, lep_1)
{
    lep1 = lep_1;
    lep2 = lep_2;
    lep3 = lep_3;
    meson = meson_i;
    vectorM = vector_i;

    setParametersForObservable(SM.getFlavour().getmyMVlnu(meson, vectorM, QCD::lepton::TAU).initializemyMVlnuParameters());
}

double R_A_5_myMVlnu::computeThValue()
{
    double A5tau = SM.getFlavour().getmyMVlnu(meson, vectorM, lep1).integrateI(6) / computeGamma(lep1);
    double A5mu = SM.getFlavour().getmyMVlnu(meson, vectorM, lep2).integrateI(6) / computeGamma(lep2);
    double A5e = SM.getFlavour().getmyMVlnu(meson, vectorM, lep3).integrateI(6) / computeGamma(lep3);

    double res = 2. * A5tau / (A5mu + A5e);

    return res;
}


/**
* R_A_6
**/
R_A_6_myMVlnu::R_A_6_myMVlnu(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i,
            QCD::lepton lep_1, QCD::lepton lep_2, QCD::lepton lep_3)
: Gamma_myMVlnu(SM_i, meson_i, vector_i, lep_1)
{
    lep1 = lep_1;
    lep2 = lep_2;
    lep3 = lep_3;
    meson = meson_i;
    vectorM = vector_i;

    setParametersForObservable(SM.getFlavour().getmyMVlnu(meson, vectorM, QCD::lepton::TAU).initializemyMVlnuParameters());
}

double R_A_6_myMVlnu::computeThValue()
{
    double A6tau = SM.getFlavour().getmyMVlnu(meson, vectorM, lep1).integrateI(8) / computeGamma(lep1);
    double A6mu = SM.getFlavour().getmyMVlnu(meson, vectorM, lep2).integrateI(8) / computeGamma(lep2);
    double A6e = SM.getFlavour().getmyMVlnu(meson, vectorM, lep3).integrateI(8) / computeGamma(lep3);

    double res = 2. * A6tau / (A6mu + A6e);

    return res;
}


/**
* D_A_7
**/
D_A_7_myMVlnu::D_A_7_myMVlnu(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i,
            QCD::lepton lep_1, QCD::lepton lep_2, QCD::lepton lep_3)
: A_7_myMVlnu(SM_i, meson_i, vector_i, lep_1)
{
    lep1 = lep_1;
    lep2 = lep_2;
    lep3 = lep_3;
    meson = meson_i;
    vectorM = vector_i;

    setParametersForObservable(SM.getFlavour().getmyMVlnu(meson, vectorM, QCD::lepton::TAU).initializemyMVlnuParameters());
}

double D_A_7_myMVlnu::computeThValue()
{
    return A_7_myMVlnu(SM,meson,vectorM,lep1).computeThValue() -
            (1./2. * (A_7_myMVlnu(SM,meson,vectorM,lep2).computeThValue()
            + A_7_myMVlnu(SM,meson,vectorM,lep3).computeThValue()) );
}


/**
* D_A_8
**/
D_A_8_myMVlnu::D_A_8_myMVlnu(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i,
            QCD::lepton lep_1, QCD::lepton lep_2, QCD::lepton lep_3)
: A_8_myMVlnu(SM_i, meson_i, vector_i, lep_1)
{
    lep1 = lep_1;
    lep2 = lep_2;
    lep3 = lep_3;
    meson = meson_i;
    vectorM = vector_i;

    setParametersForObservable(SM.getFlavour().getmyMVlnu(meson, vectorM, QCD::lepton::TAU).initializemyMVlnuParameters());
}

double D_A_8_myMVlnu::computeThValue()
{
    return A_8_myMVlnu(SM,meson,vectorM,lep1).computeThValue() -
            (1./2. * (A_8_myMVlnu(SM,meson,vectorM,lep2).computeThValue()
            + A_8_myMVlnu(SM,meson,vectorM,lep3).computeThValue()) );
}


/**
* D_A_9
**/
D_A_9_myMVlnu::D_A_9_myMVlnu(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i,
            QCD::lepton lep_1, QCD::lepton lep_2, QCD::lepton lep_3)
: A_9_myMVlnu(SM_i, meson_i, vector_i, lep_1)
{
    lep1 = lep_1;
    lep2 = lep_2;
    lep3 = lep_3;
    meson = meson_i;
    vectorM = vector_i;

    setParametersForObservable(SM.getFlavour().getmyMVlnu(meson, vectorM, QCD::lepton::TAU).initializemyMVlnuParameters());
}

double D_A_9_myMVlnu::computeThValue()
{
    return A_9_myMVlnu(SM,meson,vectorM,lep1).computeThValue() -
            (1./2. * (A_9_myMVlnu(SM,meson,vectorM,lep2).computeThValue()
            + A_9_myMVlnu(SM,meson,vectorM,lep3).computeThValue()) );
}
