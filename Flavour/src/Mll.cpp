/*
 * Copyright (C) 2018 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "Mll.h"
#include "StandardModel.h"
#include "HeffDB1.h"

Mll::Mll(const StandardModel& SM_i, int obsFlag, QCD::meson meson_i, QCD::lepton lep_i)
: ThObservable(SM_i)
{
    meson = meson_i;
    lep = lep_i;
    if (obsFlag > 0 and obsFlag < 5) obs = obsFlag;
    else throw std::runtime_error("obsFlag in Bsmumu(myFlavour, obsFlag) called from ThFactory::ThFactory() can only be 1 (BR) or 2 (BRbar) or 3 (Amumu) or 4 (Smumu)");
    SM.initializeMeson(meson);
    FixedWCbtos = SM.getFlavour().getFlagFixedWCbtos();
    LoopModelDM = SM.getFlavour().getFlagLoopModelDM();
    if (FixedWCbtos) setParametersForObservable({ "C10_SM" });
    if (FixedWCbtos && LoopModelDM) setParametersForObservable({ "C10_SM", "QB", "ysybgD", "gD", "gmuV_NP", "rAV_NP", "mB_NP", "mchi_NP", "mV_NP"});
};

double Mll::computeThValue()
{
    computeObs(FULLNLO, FULLNLO_QED);
    double FBs = SM.getMesons(meson).getDecayconst();
    
//    double coupling = SM.getGF() * SM.getGF() * SM.Mw() * SM.Mw() /M_PI /M_PI ; /* Double GF for including EW corrections*/
    double coupling = SM.getGF() * SM.getAle() / 4. / M_PI; /* Single GF for excluding EW corrections*/
 
//    double PRF = pow(coupling, 2.) / M_PI /8. / SM.getMesons(meson).computeWidth() * pow(FBs, 2.) * pow(mlep, 2.) * mBs * beta; /* Double GF for including EW corrections*/
    double PRF = pow(coupling, 2.) / M_PI / SM.getMesons(meson).computeWidth() * pow(FBs, 2.) * pow(mlep, 2.) * mBs * beta; /* Single GF for excluding EW corrections*/
    timeInt = (1. + Amumu * ys) / (1. - ys * ys); // Note modification in form due to algorithm
     
    if (obs == 1) return( PRF * ampSq);
    if (obs == 2) return( PRF * ampSq * timeInt);
    if (obs == 3) return( Amumu );
    if (obs == 4) return( Smumu );
    
    throw std::runtime_error("Bsmumu::computeThValue(): Observable type not defined. Can be only any of (1,2,3,4)");
    return (EXIT_FAILURE);
}

void Mll::computeObs(orders order, orders_qed order_qed)
{   
    double mu = SM.getMub();  
    
    mlep = SM.getLeptons(lep).getMass();
    mBs = SM.getMesons(meson).getMass();
    mW = SM.Mw();
    mb = SM.getQuarks(QCD::BOTTOM).getMass();
    if (meson == QCD::B_S) {
        ms = SM.getQuarks(QCD::STRANGE).getMass();
        CKM_factor = SM.getCKM().computelamt_s();
        ys = SM.getMesons(meson).getDgamma_gamma()/2.; // For now. To be explicitly calculated.
    } else if (meson == QCD::B_D) {
        ms = SM.getQuarks(QCD::DOWN).getMass();
        CKM_factor = SM.getCKM().computelamt_d();
        ys = 0.;
    }
    chiral = pow(mBs, 2.) / 2. / mlep * mb / (mb + ms);
    beta = sqrt(1. - pow(2. * mlep / mBs, 2.));
    
    if (LoopModelDM) {
        ysybgD_logscale = SM.getFlavour().getFlagysybgD_logscale();
        gmu_logscale = SM.getFlavour().getFlaggmu_logscale();
        rAV_parametric = SM.getFlavour().getFlagrAV_parametric();

        gD = SM.getOptionalParameter("gD");
        QB = SM.getOptionalParameter("QB");
        mV_NP = SM.getOptionalParameter("mV_NP");
        if (ysybgD_logscale){
            ysybgD = pow(10.,SM.getOptionalParameter("ysybgD"));
        } else {
            ysybgD = SM.getOptionalParameter("ysybgD");
        }
        if (gmu_logscale){
            gmuV_NP = pow(10.,SM.getOptionalParameter("gmuV_NP"));
        } else {
            gmuV_NP = SM.getOptionalParameter("gmuV_NP");
        }
        if (!rAV_parametric) {
            rAV = SM.getOptionalParameter("rAV_NP");
        } else {
            rAV = -0.455 + 0.0244/mV_NP - 0.000652/mV_NP/mV_NP ;
        }
        gmuA_NP = rAV * gmuV_NP;
        mB_NP = SM.getOptionalParameter("mB_NP");
        mchi_NP = SM.getOptionalParameter("mchi_NP");

        mB2_NP = mB_NP * mB_NP;
        mchi2_NP = mchi_NP * mchi_NP;
        mV2_NP = mV_NP * mV_NP;
        y_NP = mchi2_NP / mB2_NP;
        Norm_NP = - sqrt(2.) / (4. * SM.getGF() * 0.0411494) / (4. * M_PI * SM.getAle()); // N.B. took the abs of CKM, hence changed overall sign
    }
    
    computeAmpSq(order, order_qed, mu);
    Amumu = (absP * absP * cos(2. * argP - phiNP) -  absS * absS * cos(2. * argS - phiNP)) / (absP * absP + absS * absS);
    Smumu = (absP * absP * sin(2. * argP - phiNP) -  absS * absS * sin(2. * argS - phiNP)) / (absP * absP + absS * absS);
}

double Mll::computeAmumu(orders order)
{
    computeObs(FULLNLO, FULLNLO_QED);
    return(Amumu);
}

double Mll::computeSmumu(orders order)
{
    computeObs(FULLNLO, FULLNLO_QED);
    return(Smumu);
}

double Mll::F9(double x)
{
    double xm14 = (x - 1.)*(x - 1.)*(x - 1.)*(x - 1.);

    if (x == 1.)
        return -1./24.;
    else
        return (- 2. + 9.*x - 18.*x*x + 11.*x*x*x - 6.*x*x*x*log(x)) / 36. / xm14;
}

double Mll::G9(double x)
{
    double xm14 = (x - 1.)*(x - 1.)*(x - 1.)*(x - 1.);

    if (x == 1.)
        return 1./8.;
    else
        return (- 16. + 45.*x - 36.*x*x + 7.*x*x*x + 6.*(3.*x - 2.)*log(x)) / 36. / xm14;
}

gslpp::complex Mll::C10_NP(double q2, double gmu_V, double gmu_A)
{
    double mchi2omV2 = mchi_NP*mchi_NP/mV_NP/mV_NP;
    double mmu2omV2 = mlep*mlep/mV_NP/mV_NP;

    double gDterm, gmuterm;

    if (4.*mchi2omV2 > 1.)
        gDterm = 0;
    else
        gDterm = sqrt(1.-4.*mchi2omV2) * (2.*mchi2omV2 + 1.);

    if (4.*mmu2omV2 > 1.)
        gmuterm = 0;
    else
        gmuterm = sqrt(1.-4.*mmu2omV2);

    double GammaV = (gD*gD * gDterm
                    + gmu_V*gmu_V * gmuterm * (2.*mmu2omV2 + 1.)
                    + gmu_A*gmu_A * gmuterm * (1.-4.*mmu2omV2))/12./M_PI;

    return Norm_NP * QB * ysybgD * gmu_A / mB2_NP * (F9(y_NP) + G9(y_NP)) * q2 /
            ( q2 - mV2_NP + gslpp::complex::i()*mV2_NP*GammaV );
}

void Mll::computeAmpSq(orders order, orders_qed order_qed, double mu)
{  
    if (SM.getFlavour().getHDB1().getCoeffsmumu().getOrder() < order % 3){
        std::stringstream out;
        out << order;
        throw std::runtime_error("Bsmumu::computeAmpSq(): required cofficient of "
                                 "order " + out.str() + " not computed");
    }
    /* Temporary usage of MVll class here below  */
//    gslpp::vector<gslpp::complex> ** allcoeffmumu; /* Double GF for including EW corrections*/
//    if (meson == QCD::B_S) allcoeffmumu = SM.getFlavour().ComputeCoeffsmumu(mu, NDR); /* Double GF for including EW corrections*/
//    if (meson == QCD::B_D) allcoeffmumu = SM.getFlavour().ComputeCoeffdmumu(mu, NDR); /* Double GF for including EW corrections*/
    
    allcoeff = SM.getFlavour().ComputeCoeffBMll(mu, lep); /* Single GF for excluding EW corrections*/
    allcoeffprime = SM.getFlavour().ComputeCoeffprimeBMll(mu, lep); /* Single GF for excluding EW corrections*/

//    double alsmu = SM.Als(mu, FULLNNLO, true)/4./M_PI; /* tilde */ /* Double GF for including EW corrections*/
//    double alemu = SM.Ale(mu, FULLNNLO)/4./M_PI; /* tilde */ /* Double GF for including EW corrections*/

//    double sw = sqrt( (M_PI * SM.getAle() ) / ( sqrt(2.) * SM.getGF() * SM.Mw() * SM.Mw()) ); /* Spurious sw */

    C_10p = (*(allcoeffprime[LO]))(9) + (*(allcoeffprime[NLO]))(9);
    C_S = (*(allcoeff[LO]))(10) + (*(allcoeff[NLO]))(10);
    C_Sp = (*(allcoeffprime[LO]))(10) + (*(allcoeffprime[NLO]))(10);
    C_P = (*(allcoeff[LO]))(11) + (*(allcoeff[NLO]))(11);
    C_Pp = (*(allcoeffprime[LO]))(11) + (*(allcoeffprime[NLO]))(11);
    
    if (FixedWCbtos) {
        allcoeff_noSM = SM.getFlavour().ComputeCoeffBMll(mu, lep, true); /* Single GF for excluding EW corrections*/
        C_10 = SM.getOptionalParameter("C10_SM") + ((*(allcoeff_noSM[LO]))(9) + (*(allcoeff_noSM[NLO]))(9));
    }
    else C_10 = ((*(allcoeff[LO]))(9) + (*(allcoeff[NLO]))(9));
    
    if (LoopModelDM && lep == QCD::MU) {
        C_10 += C10_NP(mBs*mBs, gmuV_NP, gmuA_NP);
        C_P += - 2. * mlep * (mb + ms) / (mV_NP * mV_NP) * C10_NP(mBs*mBs, gmuV_NP, gmuA_NP);
    }
    
    if ((order == FULLNLO) && (order_qed == FULLNLO_QED)) {

        switch (order_qed) {
            case FULLNLO_QED:
            {
                /* Implementation to be corrected and updated with new EVO: At present better to use MVll!*/ /* Double GF for including EW corrections*/
                //            gslpp::complex C10_SM = (*(allcoeffmumu[LO]))(7) /alemu  + (*(allcoeffmumu[NLO]))(7) * alsmu/alemu 
                //                    + (*(allcoeffmumu[NNLO]))(7) * alsmu * alsmu/alemu + (*(allcoeffmumu[LO_QED ]))(7) /alsmu
                //                    + (*(allcoeffmumu[NLO_QED11]))(7) + (*(allcoeffmumu[NLO_QED02]))(7) * alemu /alsmu /alsmu 
                //                    + (*(allcoeffmumu[NLO_QED21]))(7) * alsmu 
                //                    + (*(allcoeffmumu[NLO_QED12]))(7) * alemu /alsmu+ (*(allcoeffmumu[NLO_QED22]))(7) * alemu;

                /* Temporary usage of MVll result */
                //            std::cout << " C10_SM " << C10_SM << std::endl;
                //            gslpp::complex C10_SM_plus_NP = CKM_factor * sw * sw * ((*(allcoeff[LO]))(9) + (*(allcoeff[NLO]))(9));
                //            std::cout << " C10_SM_plus_NP " << C10_SM_plus_NP / sw / sw / CKM_factor << std::endl;

                //            gslpp::complex CC_P = C10_SM; /* Double GF for including EW corrections*/
                //            gslpp::complex CC_P = CKM_factor * sw * sw * ( C_10 - C_10p + mBs*mBs / ( 2.*mlep*(mb+ms) ) * (C_P - C_Pp) ); /* Spurious sw */
                gslpp::complex CC_P = CKM_factor * (C_10 - C_10p + mBs * mBs / (2. * mlep * (mb + ms)) * (C_P - C_Pp)); /* Single GF for excluding EW corrections*/

                absP = CC_P.abs(); //contains only SM contributions (P, P', S, S' not added)
                argP = CC_P.arg();

                //            gslpp::complex CC_S = CKM_factor * sw * sw * ( beta * mBs*mBs / ( 2.*mlep*(mb+ms) ) * (C_S - C_Sp) ); /* Spurious sw */
                gslpp::complex CC_S = CKM_factor * (beta * mBs * mBs / (2. * mlep * (mb + ms)) * (C_S - C_Sp)); /* Single GF for excluding EW corrections*/

                absS = CC_S.abs();
                argS = CC_S.arg();

                phiNP = 0.;

                ampSq = absP * absP + absS * absS;

            }
                break;
            default:
                std::stringstream out;
                out << order;
                throw std::runtime_error("Bsmumu::computeAmpSq(): order " + out.str() + " not implemented");
        }
    }
        
}


