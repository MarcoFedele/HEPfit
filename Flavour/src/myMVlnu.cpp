/*
 * Copyright (C) 2018 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "StandardModel.h"
#include "myMVlnu.h"
#include "std_make_vector.h"
#include "gslpp_function_adapter.h"
#include <boost/bind.hpp>

myMVlnu::myMVlnu(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i, QCD::lepton lep_i)
: mySM(SM_i),
meson_cache(2, 0.),
h_A1_cache(2, 0.),
h_A1_BGL_cache(3, 0.),
R_1_BGL_cache(3, 0.),
R_2_BGL_cache(2, 0.),
quark_cache(2, 0.),
N_cache(3, 0.)
{
    lep = lep_i;
    meson = meson_i;
    vectorM = vector_i;
    CLNflag = false;
    BGLflag = false;
    MSflag = false;
    ChiralBasisflag = false;

    I1c_updated = 0;
    I1s_updated = 0;
    I2c_updated = 0;
    I2s_updated = 0;
    I2n_updated = 0;
    I3_updated = 0;
    I4_updated = 0;
    I5_updated = 0;
    I6c_updated = 0;
    I6s_updated = 0;
    I7_updated = 0;
    I8_updated = 0;
    I9_updated = 0;

    w_I = gsl_integration_cquad_workspace_alloc (100);
}

myMVlnu::~myMVlnu()
{}


std::vector<std::string> myMVlnu::initializemyMVlnuParameters()
{
    CLNflag = mySM.getFlavour().getFlagCLN();
    BGLflag = mySM.getFlavour().getFlagBGL();
    MSflag = mySM.getFlavour().getFlagMS();
    ChiralBasisflag = mySM.getFlavour().getFlagChiralBasis();

    if (MSflag + CLNflag + BGLflag != true) throw std::runtime_error("myMVlnu: only one between MSflag, CLNflag and BGLflag can be true");

    if (MSflag) {
        if (vectorM == StandardModel::D_star_P && lep == StandardModel::TAU) mvlnuParameters = make_vector<std::string>()
            << "V0" << "A00" << "A10" << "A20" << "T10" << "T30"
            << "regV" << "regA" << "regP" << "regT" << "imgV" << "imgA" << "imgP" << "imgT"
            << "regVL" << "regVR" << "regSL" << "regSR" << "imgVL" << "imgVR" << "imgSL" << "imgSR" ;
        else if (vectorM == StandardModel::D_star_P && (lep == StandardModel::MU || lep == StandardModel::ELECTRON))
            mvlnuParameters = make_vector<std::string>() << "V0" << "A00" << "A10" << "A20" << "T10" << "T30" ;
        else {
            std::stringstream out;
            out << vectorM;
            throw std::runtime_error("myMVlnu: vector " + out.str() + " not implemented");
        }
    } else if (CLNflag) {
        if (vectorM == StandardModel::D_star_P && lep == StandardModel::TAU) mvlnuParameters = make_vector<std::string>()
            << "h_A1_1" << "rho2" << "R_1_1" << "R_2_1" << "A_0_err" << "T_1_err"  << "T_2_err"  << "T_3_err"
            << "regV" << "regA" << "regP" << "regT" << "imgV" << "imgA" << "imgP" << "imgT"
            << "regVL" << "regVR" << "regSL" << "regSR" << "imgVL" << "imgVR" << "imgSL" << "imgSR" ;
        else if (vectorM == StandardModel::D_star_P && (lep == StandardModel::MU || lep == StandardModel::ELECTRON))
            mvlnuParameters = make_vector<std::string>() << "h_A1_1" << "rho2" << "R_1_1" << "R_2_1"
                                                         << "A_0_err" << "T_1_err"  << "T_2_err"  << "T_3_err" ;
        else {
            std::stringstream out;
            out << vectorM;
            throw std::runtime_error("myMVlnu: vector " + out.str() + " not implemented");
        }
    } else if (BGLflag) {
        if (vectorM == StandardModel::D_star_P && lep == StandardModel::TAU) mvlnuParameters = make_vector<std::string>()
            << "af0" << "af1" << "af2" << "ag0" << "ag1" << "ag2" << "aF11" << "aF12" << "AbsVcb"
            << "A_0_err" << "T_1_err"  << "T_2_err"  << "T_3_err"
            << "regV" << "regA" << "regP" << "regT" << "imgV" << "imgA" << "imgP" << "imgT";
        else if (vectorM == StandardModel::D_star_P && (lep == StandardModel::MU || lep == StandardModel::ELECTRON))
            mvlnuParameters = make_vector<std::string>() << "af0" << "af1" << "af2" << "ag0" << "ag1" << "ag2" << "aF11" << "aF12" << "AbsVcb"
                                                         << "A_0_err" << "T_1_err"  << "T_2_err"  << "T_3_err" ;
        else {
            std::stringstream out;
            out << vectorM;
            throw std::runtime_error("myMVlnu: vector " + out.str() + " not implemented");
        }
    }

    mySM.initializeMeson(meson);
    mySM.initializeMeson(vectorM);
    return mvlnuParameters;
}


void myMVlnu::checkCache()
{

    if (gV == gV_cache) {
        gV_updated = 1;
    } else {
        gV_updated = 0;
        gV_cache = gV;
    }

    if (gA == gA_cache) {
        gA_updated = 1;
    } else {
        gA_updated = 0;
        gA_cache = gA;
    }

    if (gP == gP_cache) {
        gP_updated = 1;
    } else {
        gP_updated = 0;
        gP_cache = gP;
    }

    if (gT == gT_cache) {
        gT_updated = 1;
    } else {
        gT_updated = 0;
        gT_cache = gT;
    }

    if (MM == meson_cache(0) && MV == meson_cache(1)) {
        meson_updated = 1;
    } else {
        meson_updated = 0;
        meson_cache(0) = MM;
        meson_cache(1) = MV;
    }

    if (CLNflag) {
        if (h_A1_1 == h_A1_cache(0) && rho2 == h_A1_cache(1) && Vcb == h_A1c_cache) {
            h_A1_updated = meson_updated;
        } else {
            h_A1_updated = 0;
            h_A1_cache(0) = h_A1_1;
            h_A1_cache(1) = rho2;
            h_A1c_cache = Vcb;
        }

        if (R_1_1 == R_1_1_cache) {
            R_1_1_updated = 1;
        } else {
            R_1_1_updated = 0;
            R_1_1_cache = R_1_1;
        }

        if (R_2_1 == R_2_1_cache) {
            R_2_1_updated = 1;
        } else {
            R_2_1_updated = 0;
            R_2_1_cache = R_2_1;
        }

        if (A_0_err == A_0_err_cache) {
            A_0_err_updated = 1;
        } else {
            A_0_err_updated = 0;
            A_0_err_cache = A_0_err;
        }

        if (T_1_err == T_1_err_cache) {
            T_1_err_updated = 1;
        } else {
            T_1_err_updated = 0;
            T_1_err_cache = T_1_err;
        }

        if (T_2_err == T_2_err_cache) {
            T_2_err_updated = 1;
        } else {
            T_2_err_updated = 0;
            T_2_err_cache = T_2_err;
        }

        if (T_3_err == T_3_err_cache) {
            T_3_err_updated = 1;
        } else {
            T_3_err_updated = 0;
            T_3_err_cache = T_3_err;
        }

        V_updated = h_A1_updated*R_1_1_updated*meson_updated;
        A_0_updated = h_A1_updated*A_0_err_updated;
        A_1_updated = h_A1_updated;
        A_2_updated = h_A1_updated*R_2_1_updated*meson_updated;
        T_1_updated = h_A1_updated*T_1_err_updated;
        T_2_updated = h_A1_updated*T_2_err_updated;
        T_3_updated = h_A1_updated*T_3_err_updated;

    } else if (BGLflag) {
        if (af0 == h_A1_BGL_cache(0) && af1 == h_A1_BGL_cache(1) && af2 == h_A1_BGL_cache(2)) {
            h_A1_BGL_updated = meson_updated;
        } else {
            h_A1_BGL_updated = 0;
            h_A1_BGL_cache(0) = af0;
            h_A1_BGL_cache(1) = af1;
            h_A1_BGL_cache(2) = af2;
        }

        if (ag0 == R_1_BGL_cache(0) && ag1 == R_1_BGL_cache(1) && ag2 == R_1_BGL_cache(2)) {
            R_1_BGL_updated = meson_updated;
        } else {
            R_1_BGL_updated = 0;
            R_1_BGL_cache(0) = ag0;
            R_1_BGL_cache(1) = ag1;
            R_1_BGL_cache(2) = ag2;
        }

        if (aF11 == R_2_BGL_cache(0) && aF12 == R_2_BGL_cache(1)) {
            R_2_BGL_updated = meson_updated;
        } else {
            R_2_BGL_updated = 0;
            R_2_BGL_cache(0) = aF11;
            R_2_BGL_cache(1) = aF12;
        }

        if (A_0_err == A_0_err_cache) {
            A_0_err_updated = 1;
        } else {
            A_0_err_updated = 0;
            A_0_err_cache = A_0_err;
        }

        if (T_1_err == T_1_err_cache) {
            T_1_err_updated = 1;
        } else {
            T_1_err_updated = 0;
            T_1_err_cache = T_1_err;
        }

        if (T_2_err == T_2_err_cache) {
            T_2_err_updated = 1;
        } else {
            T_2_err_updated = 0;
            T_2_err_cache = T_2_err;
        }

        if (T_3_err == T_3_err_cache) {
            T_3_err_updated = 1;
        } else {
            T_3_err_updated = 0;
            T_3_err_cache = T_3_err;
        }

        V_updated = h_A1_BGL_updated*R_1_BGL_updated;
        A_0_updated = h_A1_BGL_updated*A_0_err_updated;
        A_1_updated = h_A1_BGL_updated;
        A_2_updated = h_A1_BGL_updated*R_2_BGL_updated;
        T_1_updated = h_A1_BGL_updated*T_1_err_updated;
        T_2_updated = h_A1_BGL_updated*T_2_err_updated;
        T_3_updated = h_A1_BGL_updated*T_3_err_updated;

    } else if (MSflag) {
        if (V0 == V_cache) {
            V_updated = 1;
        } else {
            V_updated = 0;
            V_cache = V0;
        }

        if (A00 == A_0_cache) {
            A_0_updated = 1;
        } else {
            A_0_updated = 0;
            A_0_cache = A00;
        }

        if (A10 == A_1_cache) {
            A_1_updated = 1;
        } else {
            A_1_updated = 0;
            A_1_cache = A10;
        }

        if (A20 == A_2_cache) {
            A_2_updated = 1;
        } else {
            A_2_updated = 0;
            A_2_cache = A20;
        }

        if (T10 == T_1_cache) {
            T_1_updated = 1;
        } else {
            T_1_updated = 0;
            T_1_cache = T10;
        }

        if (T20 == T_2_cache) {
            T_2_updated = 1;
        } else {
            T_2_updated = 0;
            T_2_cache = T20;
        }

        if (T30 == T_3_cache) {
            T_3_updated = 1;
        } else {
            T_3_updated = 0;
            T_3_cache = T30;
        }
    }

    if (Mb == quark_cache(0) && Mc == quark_cache(1)) {
        quark_updated = 1;
    } else {
        quark_updated = 0;
        quark_cache(0) = Mb;
        quark_cache(1) = Mc;
    }

    if (GF == N_cache(0) && Mlep == N_cache(1) && MM == N_cache(2) && Vcb == Nc_cache) {
        N_updated = 1;
    } else {
        N_updated = 0;
        N_cache(0) = GF;
        N_cache(1) = Mlep;
        N_cache(2) = MM;
        Nc_cache = Vcb;
    }

    H_V_p_updated = gV_updated*meson_updated*V_updated;
    H_A_t_updated = gA_updated*meson_updated*A_0_updated;
    H_A_0_updated = gA_updated*meson_updated*A_1_updated*A_2_updated;
    H_A_p_updated = gA_updated*meson_updated*A_1_updated;
    H_P_updated = gP_updated*meson_updated*quark_updated*A_0_updated;
    H_T_pt_updated = gT_updated*meson_updated*T_1_updated;
    H_T_p0_updated = gT_updated*meson_updated*T_2_updated;
    H_T_0_updated = gT_updated*meson_updated*T_2_updated*T_3_updated;

    H_pm_updated = H_V_p_updated*H_A_p_updated;
    H_0_updated = H_A_0_updated;
    H_t_updated = H_A_t_updated*H_P_updated;
    H_T_pm_updated = H_T_pt_updated*H_T_p0_updated;

    H_sum_pm_updated = H_pm_updated*H_T_pm_updated; /* Mlep is checked by N_updated */
    H_sum_0_updated = H_0_updated*H_T_0_updated; /* Mlep is checked by N_updated */

    I1c_updated = N_updated*H_sum_pm_updated;
    I1s_updated = N_updated*H_sum_pm_updated;
    I2c_updated = N_updated*H_sum_0_updated;
    I2s_updated = N_updated*H_sum_0_updated;
    I2n_updated = N_updated*H_t_updated;
    I3_updated = N_updated*H_pm_updated*H_T_pm_updated;
    I4_updated = N_updated*H_pm_updated*H_0_updated*H_T_pm_updated*H_T_0_updated;
    I5_updated = N_updated*H_sum_pm_updated*H_sum_0_updated*H_t_updated;
    I6c_updated = N_updated*H_sum_0_updated*H_t_updated;
    I6s_updated = N_updated*H_sum_pm_updated;
    I7_updated = N_updated*H_sum_pm_updated*H_sum_0_updated*H_t_updated;
    I8_updated = N_updated*H_pm_updated*H_0_updated*H_T_pm_updated*H_T_0_updated;
    I9_updated = N_updated*H_pm_updated*H_T_pm_updated;

    return;
}



void myMVlnu::updateParameters()
{
    if (!mySM.getFlavour().getUpdateFlag(meson, vectorM, lep)) return;

    GF = mySM.getGF();
    Mlep = mySM.getLeptons(lep).getMass();
    MM = mySM.getMesons(meson).getMass();
    MV = mySM.getMesons(vectorM).getMass();
    Mb = mySM.getQuarks(QCD::BOTTOM).getMass();
    Mc = mySM.getQuarks(QCD::CHARM).getMass();

    if (lep == StandardModel::TAU) {
        if (ChiralBasisflag) {
            gV = mySM.getOptionalParameter("regVR") + gslpp::complex::i()*mySM.getOptionalParameter("imgVR")
                + mySM.getOptionalParameter("regVL") + gslpp::complex::i()*mySM.getOptionalParameter("imgVL");
            gA = mySM.getOptionalParameter("regVR") + gslpp::complex::i()*mySM.getOptionalParameter("imgVR")
                - mySM.getOptionalParameter("regVL") - gslpp::complex::i()*mySM.getOptionalParameter("imgVL");
            gP = mySM.getOptionalParameter("regSR") + gslpp::complex::i()*mySM.getOptionalParameter("imgSR")
                - mySM.getOptionalParameter("regSL") - gslpp::complex::i()*mySM.getOptionalParameter("imgSL");
            gT = 2.*(mySM.getOptionalParameter("regT") + gslpp::complex::i()*mySM.getOptionalParameter("imgT"));
            //gT = 2.*(mySM.getOptionalParameter("regSL") + gslpp::complex::i()*mySM.getOptionalParameter("imgSL"))/8.14;
        }
        else {
            gV = mySM.getOptionalParameter("regV") + gslpp::complex::i()*mySM.getOptionalParameter("imgV");
            gA = mySM.getOptionalParameter("regA") + gslpp::complex::i()*mySM.getOptionalParameter("imgA");
            gP = mySM.getOptionalParameter("regP") + gslpp::complex::i()*mySM.getOptionalParameter("imgP");
            gT = mySM.getOptionalParameter("regT") + gslpp::complex::i()*mySM.getOptionalParameter("imgT");
        }
    } else {
        gV = 0;
        gA = 0;
        gP = 0;
        gT = 0;
    }


    if (MSflag) {
        switch (vectorM) {
        case StandardModel::D_star_P:
            V0 = mySM.getOptionalParameter("V0");
            A00 = mySM.getOptionalParameter("A00");
            A10 = mySM.getOptionalParameter("A10");
            A20 = mySM.getOptionalParameter("A20");
            T10 = mySM.getOptionalParameter("T10");
            T20 = T10;
            T30 = mySM.getOptionalParameter("T30");
            Vcb = mySM.getCKM().getV_cb();
            //Vcb = 0.0411;
            break;
        default:
            std::stringstream out;
            out << vectorM;
            throw std::runtime_error("myMVlnu: vector " + out.str() + " not implemented");
        }
    } else if (CLNflag) {
        switch (vectorM) {
        case StandardModel::D_star_P:
            h_A1_1 = mySM.getOptionalParameter("h_A1_1");
            rho2 = mySM.getOptionalParameter("rho2");
            R_1_1 = mySM.getOptionalParameter("R_1_1");
            R_2_1 = mySM.getOptionalParameter("R_2_1");
            A_0_err = mySM.getOptionalParameter("A_0_err");
            T_1_err = mySM.getOptionalParameter("T_1_err");
            T_2_err = mySM.getOptionalParameter("T_2_err");
            T_3_err = mySM.getOptionalParameter("T_3_err");
            Vcb = mySM.getCKM().getV_cb();
            //Vcb = 0.0411;
            break;
        default:
            std::stringstream out;
            out << vectorM;
            throw std::runtime_error("myMVlnu: vector " + out.str() + " not implemented");
        }
    } else if (BGLflag) {
        switch (vectorM) {
        case StandardModel::D_star_P:
            af0 = mySM.getOptionalParameter("af0");
            af1 = mySM.getOptionalParameter("af1");
            af2 = mySM.getOptionalParameter("af2");
            ag0 = mySM.getOptionalParameter("ag0");
            ag1 = mySM.getOptionalParameter("ag1");
            ag2 = mySM.getOptionalParameter("ag2");
            aF11 = mySM.getOptionalParameter("aF11");
            aF12 = mySM.getOptionalParameter("aF12");
            A_0_err = mySM.getOptionalParameter("A_0_err");
            T_1_err = mySM.getOptionalParameter("T_1_err");
            T_2_err = mySM.getOptionalParameter("T_2_err");
            T_3_err = mySM.getOptionalParameter("T_3_err");
            Vcb = mySM.getOptionalParameter("AbsVcb");
            break;
        default:
            std::stringstream out;
            out << vectorM;
            throw std::runtime_error("myMVlnu: vector " + out.str() + " not implemented");
        }
    }


    checkCache();

    MM2 = MM*MM;
    MV2 = MV*MV;
    Mlep2 = Mlep*Mlep;
    NN = (GF * GF * Vcb.abs2()) / (384. * M_PI*M_PI*M_PI * MM2*MM) ;

    mBcstV1 = 6.329;
    mBcstV2 = 6.920;
    mBcstV3 = 7.020;
    mBcstV4 = 7.280;
    mBcstA1 = 6.739;
    mBcstA2 = 6.750;
    mBcstA3 = 7.145;
    mBcstA4 = 7.150;

    zV1 = zP(mBcstV1*mBcstV1);
    zV2 = zP(mBcstV2*mBcstV2);
    zV3 = zP(mBcstV3*mBcstV3);
    zV4 = zP(mBcstV4*mBcstV4);

    zA1 = zP(mBcstA1*mBcstA1);
    zA2 = zP(mBcstA2*mBcstA2);
    zA3 = zP(mBcstA3*mBcstA3);
    zA4 = zP(mBcstA4*mBcstA4);

    etaI = 2.6;
    chiTV = 5.131e-4;
    chiTA = 3.894e-4;

    /*std::cout << "I1c(4.) = " << I_1c(4.) << std::endl;
    std::cout << "I1s(4.) = " << I_1s(4.) << std::endl;
    std::cout << "I2c(4.) = " << I_2c(4.) << std::endl;
    std::cout << "I2s(4.) = " << I_2s(4.) << std::endl;
    std::cout << "I3(4.) = " << I_3(4.) << std::endl;
    std::cout << "I4(4.) = " << I_4(4.) << std::endl;
    std::cout << "I5(4.) = " << I_5(4.) << std::endl;
    std::cout << "I6c(4.) = " << I_6c(4.) << std::endl;
    std::cout << "I6s(4.) = " << I_6s(4.) << std::endl;
    std::cout << "I7(4.) = " << I_7(4.) << std::endl;
    std::cout << "I8(4.) = " << I_8(4.) << std::endl;
    std::cout << "I9(4.) = " << I_9(4.) << std::endl << std::endl;


    std::cout << "R1 = " << R_1(4.) << std::endl;
    std::cout << "R2 = " << R_2(4.) << std::endl;
    std::cout << "V = " << V(4.) << std::endl;
    std::cout << "A0 = " << A_0(4.) << std::endl;
    std::cout << "A1 = " << A_1(4.) << std::endl;
    std::cout << "A2 = " << A_2(4.) << std::endl;
    std::cout << "T1 = " << T_1(4.) << std::endl;
    std::cout << "T2 = " << T_2(4.) << std::endl;
    std::cout << "T3 = " << T_3(4.) << std::endl << std::endl;*/
    /*std::cout << "V(1) = " << V(1.) << std::endl;
    std::cout << "V(2) = " << V(2.) << std::endl;
    std::cout << "V(3) = " << V(3.) << std::endl;
    std::cout << "V(4) = " << V(4.) << std::endl;
    std::cout << "V(5) = " << V(5.) << std::endl;
    std::cout << "V(6) = " << V(6.) << std::endl;
    std::cout << "V(7) = " << V(7.) << std::endl;
    std::cout << "V(8) = " << V(8.) << std::endl;
    std::cout << "V(9) = " << V(9.) << std::endl;
    std::cout << "V(10) = " << V(10.) << std::endl;*/
    mySM.getFlavour().setUpdateFlag(meson, vectorM, lep, false);

    return;
}

double myMVlnu::FF_fit1(double q2, double f_0, double s_1)
{
    double MR_2 = 6.4*6.4;
    return f_0 / (1. - q2 / MR_2) / (1. - s_1 * q2 / MR_2);
}

double myMVlnu::FF_fit2(double q2, double f_0, double s_1, double s_2)
{
    double MR_2 = 6.4*6.4;
    return f_0 / (1. - s_1 * q2 / MR_2 + s_2 * q2 * q2 / MR_2 / MR_2);
}

double myMVlnu::w(double q2)
{
    return (MM2 + MV2 - q2)/(2. * MM * MV);
}

double myMVlnu::z(double q2)
{
    return (sqrt(w(q2) + 1.) - sqrt(2.))/(sqrt(w(q2) + 1.) + sqrt(2.));
}

double myMVlnu::zP(double M2)
{
    double num = sqrt( (MM+MV)*(MM+MV) - M2 ) - sqrt( (MM+MV)*(MM+MV) - (MM-MV)*(MM-MV) );
    double den = sqrt( (MM+MV)*(MM+MV) - M2 ) + sqrt( (MM+MV)*(MM+MV) - (MM-MV)*(MM-MV) );

    return num/den;
}

double myMVlnu::phi_f(double z){
    double prefac = 4.*(MV/MM)/MM/MM*sqrt(etaI/(3.*M_PI*chiTA));
    double num = (1.+z)*sqrt((1.-z)*(1.-z)*(1.-z));
    double den = (1.+MV/MM)*(1.-z)+2.*sqrt(MV/MM)*(1.+z);
    double den4 = den*den*den*den;
    return prefac*num/den4;
}

double myMVlnu::f_BGL(double q2)
{
    double zz = z(q2);
    double Pfacf = (zz-zA1)/(1.-zz*zA1) * (zz-zA2)/(1.-zz*zA2) * (zz-zA3)/(1.-zz*zA3) * (zz-zA4)/(1.-zz*zA4);
    double phif = phi_f(zz);

    return (af0 + af1*zz + af2*zz*zz)/phif/Pfacf;
}

double myMVlnu::phi_g(double z){
    double prefac = sqrt(etaI/(3.*M_PI*chiTV));
    double num = 16.*(MV/MM)*(MV/MM)*(1.+z)*(1.+z)/sqrt(1.-z);
    double den = (1.+MV/MM)*(1.-z)+2.*sqrt(MV/MM)*(1.+z);
    double den4 = den*den*den*den;

    return prefac*num/den4;
}

double myMVlnu::g_BGL(double q2)
{
    double zz = z(q2);
    double Pfacg = (zz-zV1)/(1.-zz*zV1) * (zz-zV2)/(1.-zz*zV2) * (zz-zV3)/(1.-zz*zV3) * (zz-zV4)/(1.-zz*zV4);
    double phig = phi_g(zz);

    return (ag0 + ag1*zz + ag2*zz*zz)/phig/Pfacg;
}

double myMVlnu::phi_F1(double z){
    double prefac = 4.*(MV/MM)/MM/MM/MM*sqrt(etaI/(6.*M_PI*chiTA));
    double num = (1.+z)*sqrt((1.-z)*(1.-z)*(1.-z)*(1.-z)*(1.-z));
    double den = (1.+MV/MM)*(1.-z)+2.*sqrt(MV/MM)*(1.+z);
    double den5 = den*den*den*den*den;

    return prefac*num/den5;
}

double myMVlnu::F1_BGL(double q2)
{
    double zz = z(q2);
    double PfacF1 = (zz-zA1)/(1.-zz*zA1) * (zz-zA2)/(1.-zz*zA2) * (zz-zA3)/(1.-zz*zA3) * (zz-zA4)/(1.-zz*zA4);
    double phiF1 = phi_F1(zz);
    double aF10 = (MM-MV)*(phi_F1(0.)/phi_f(0.))*af0; // F1(z=0) = (MM-MV)*f(z=0)

    return (aF10 + aF11*zz + aF12*zz*zz)/phiF1/PfacF1;
}

double myMVlnu::h_A1(double q2)
{
    if (CLNflag) return h_A1_1/Vcb.abs() * (1. - 8.*rho2*z(q2) + (53.*rho2 - 15.)*z(q2)*z(q2) - (231.*rho2 - 91.)*z(q2)*z(q2)*z(q2));
    else return f_BGL(q2)/sqrt(MM*MV)/(1.+w(q2));
}

double myMVlnu::R_1(double q2)
{
    if (CLNflag) {
        double wm1 = w(q2) - 1.;
        return R_1_1 - 0.12*wm1 + 0.05*wm1*wm1;
    } else return (w(q2) + 1.) * MM*MV * g_BGL(q2)/f_BGL(q2);
}

double myMVlnu::R_2(double q2)
{
    double wm1 = w(q2) - 1.;

    if (CLNflag) return R_2_1 + 0.11*wm1 - 0.06*wm1*wm1;
    else return (w(q2) - (MV/MM))/wm1 - F1_BGL(q2)/f_BGL(q2)/MM/wm1;
}

double myMVlnu::R_0(double q2)
{
    return A_0_err * ( 1.1564489360511934
            + 0.03013646526546414*q2
            + 7.042838690876247e-4*q2*q2 );
}

double myMVlnu::R_T1(double q2)
{
    return T_1_err * ( 1.1263096014995009
            + 0.02635893266068262*q2
            + 6.100057859189564e-4*q2*q2
            + 1.4110100370762471e-5*q2*q2*q2 );
}

double myMVlnu::R_T2(double q2)
{
    return T_2_err * ( 1.126309601499501
            - 0.004002948915710529*q2
            - 9.702562452714152e-5*q2*q2
            - 2.3095776821592515e-6*q2*q2*q2 );
}

double myMVlnu::R_T3(double q2)
{
    return T_3_err * ( 0.543560340394165
            + 0.008838194731515387*q2
            + 1.3273876193195384e-4*q2*q2
            + 1.698295069136254e-6*q2*q2*q2 );
}

double myMVlnu::V(double q2)
{
    if (MSflag) return FF_fit1(q2, V0, 0.57);
    else return h_A1(q2) * R_1(q2) * (MM + MV)/( 2. * sqrt(MM*MV) );
}

double myMVlnu::A_0(double q2)
{
    if (MSflag) return FF_fit1(q2, A00, 0.58);
    else return A_1(q2) * R_0(q2);
}

double myMVlnu::A_1(double q2)
{
    if (MSflag) return FF_fit2(q2, A10, 0.78, 0);
    else return h_A1(q2) * ( w(q2) + 1. )/2. * ( 2. * sqrt(MM*MV) )/(MM + MV);
}

double myMVlnu::A_2(double q2)
{
    if (MSflag) return FF_fit2(q2, A20, 1.40, 0.41);
    else return h_A1(q2) * R_2(q2) * (MM + MV)/( 2. * sqrt(MM*MV) );
}

double myMVlnu::T_1(double q2)
{
    if (MSflag) return FF_fit1(q2, T10, 0.57);
    else return A_1(q2) * R_T1(q2);
}

double myMVlnu::T_2(double q2)
{
    if (MSflag) return FF_fit2(q2, T20, 0.64, 0);
    else return A_1(q2) * R_T2(q2);
}

double myMVlnu::T_3(double q2)
{
    if (MSflag) return FF_fit2(q2, T30, 1.46, 0.50);
    else return A_1(q2) * R_T3(q2);
}

double myMVlnu::lambda(double q2)
{
    return (MM2*MM2 + q2 * q2 + MV2*MV2 - 2. * (MV2 * q2 + MM2 * (q2 + MV2)));
}

gslpp::complex myMVlnu::H_V_p(double q2)
{
    return - (1. + gV) * sqrt(lambda(q2))/(MM + MV) * V(q2);
}

gslpp::complex myMVlnu::H_A_t(double q2)
{
    return (-1. + gA) * sqrt(lambda(q2)/q2) * A_0(q2);
}

gslpp::complex myMVlnu::H_A_0(double q2)
{
    return (-1. + gA)/(2. * MV * sqrt(q2)) * ( (MM + MV)*(MM2 - MV2 -q2)*A_1(q2) - lambda(q2)/(MM + MV)*A_2(q2) ) ;
}

gslpp::complex myMVlnu::H_A_p(double q2)
{
    return -(-1. + gA) * (MM + MV) * A_1(q2);
}

gslpp::complex myMVlnu::H_P(double q2)
{
    return -gP * sqrt(lambda(q2))/(Mb + Mc) * A_0(q2);
}

gslpp::complex myMVlnu::H_T_pt(double q2)
{
    return gT * sqrt(lambda(q2)/q2) * T_1(q2);
}

gslpp::complex myMVlnu::H_T_p0(double q2)
{
    return gT * (MM2 - MV2)/sqrt(q2) * T_2(q2);
}

gslpp::complex myMVlnu::H_T_0(double q2)
{
    return -gT/(2.*MV) * ( (MM2 + 3.*MV2 -q2)*T_2(q2) - lambda(q2)/(MM2 - MV2)*T_3(q2) ) ;
}

gslpp::complex myMVlnu::H_p(double q2)
{
    return H_V_p(q2) + H_A_p(q2);
}

gslpp::complex myMVlnu::H_m(double q2)
{
    return - H_V_p(q2) + H_A_p(q2);
}

gslpp::complex myMVlnu::H_0(double q2)
{
    return H_A_0(q2);
}

gslpp::complex myMVlnu::H_t(double q2)
{
    return H_A_t(q2) + sqrt(q2)/Mlep * H_P(q2);
}

gslpp::complex myMVlnu::H_T_p(double q2)
{
    return H_T_pt(q2) + H_T_p0(q2);
}

gslpp::complex myMVlnu::H_T_m(double q2)
{
    return - H_T_pt(q2) + H_T_p0(q2);
}

gslpp::complex myMVlnu::H_sum_p(double q2)
{
    return H_p(q2) - 2. * Mlep/sqrt(q2) * H_T_p(q2);
}

gslpp::complex myMVlnu::H_sum_m(double q2)
{
    return H_m(q2) - 2. * Mlep/sqrt(q2) * H_T_m(q2);
}

gslpp::complex myMVlnu::H_sum_0(double q2)
{
    return H_0(q2) - 2. * Mlep/sqrt(q2) * H_T_0(q2);
}

gslpp::complex myMVlnu::Ht_sum_p(double q2)
{
    return H_p(q2) - 2. * sqrt(q2)/Mlep * H_T_p(q2);
}

gslpp::complex myMVlnu::Ht_sum_m(double q2)
{
    return H_m(q2) - 2. * sqrt(q2)/Mlep * H_T_m(q2);
}

gslpp::complex myMVlnu::Ht_sum_0(double q2)
{
    return H_0(q2) - 2. * sqrt(q2)/Mlep * H_T_0(q2);
}

double myMVlnu::beta2(double q2)
{
    return 1. - Mlep2/q2;
}

double myMVlnu::Gamma_0(double q2)
{
    return NN * q2 * sqrt(lambda(q2)) * beta2(q2)*beta2(q2);
}

double myMVlnu::I_2n(double q2)
{
    return 2. * Gamma_0(q2) * Mlep2/q2 * H_t(q2).abs2();
}

double myMVlnu::I_1c(double q2)
{
    return 2. * Gamma_0(q2) * ( H_sum_0(q2).abs2() + Mlep2/q2 * Ht_sum_0(q2).abs2() + 2. * Mlep2/q2 * H_t(q2).abs2() );
}

double myMVlnu::I_1s(double q2)
{
    return Gamma_0(q2) / 2. *( 3. * ( H_sum_p(q2).abs2() + H_sum_m(q2).abs2() ) + Mlep2/q2 * ( Ht_sum_p(q2).abs2() + Ht_sum_m(q2).abs2() ) );
}

double myMVlnu::I_2c(double q2)
{
    return 2. * Gamma_0(q2) * ( - H_sum_0(q2).abs2() + Mlep2/q2 * Ht_sum_0(q2).abs2());
}

double myMVlnu::I_2s(double q2)
{
    return Gamma_0(q2) / 2. *( ( H_sum_p(q2).abs2() + H_sum_m(q2).abs2() ) - Mlep2/q2 * ( Ht_sum_p(q2).abs2() + Ht_sum_m(q2).abs2() ) );
}

double myMVlnu::I_3(double q2)
{
    return -2. * Gamma_0(q2) * beta2(q2) *
            ( H_p(q2)*H_m(q2).conjugate() - 4.*H_T_p(q2)*H_T_m(q2).conjugate() ).real();
}

double myMVlnu::I_4(double q2)
{
    return Gamma_0(q2) * beta2(q2) *
            ( (H_p(q2) + H_m(q2))*H_0(q2).conjugate() -
            4.*(H_T_p(q2) + H_T_m(q2))*H_T_0(q2).conjugate() ).real();
}

double myMVlnu::I_5(double q2)
{
    return 2. * Gamma_0(q2) *
            ( (H_sum_p(q2) - H_sum_m(q2))*H_sum_0(q2).conjugate() -
            Mlep2/q2 * (Ht_sum_p(q2) + Ht_sum_m(q2))*H_t(q2).conjugate() ).real();
}

double myMVlnu::I_6c(double q2)
{
    return 8. * Gamma_0(q2) * Mlep2/q2 * ( Ht_sum_0(q2)*H_t(q2).conjugate() ).real();
}

double myMVlnu::I_6s(double q2)
{
    return 2. * Gamma_0(q2) * ( H_sum_p(q2).abs2() - H_sum_m(q2).abs2() );
}

double myMVlnu::I_7(double q2)
{
    return 2. * Gamma_0(q2) *
            ( (H_sum_p(q2) + H_sum_m(q2))*H_sum_0(q2).conjugate() -
            Mlep2/q2 * (Ht_sum_p(q2) - Ht_sum_m(q2))*H_t(q2).conjugate() ).imag();
}

double myMVlnu::I_8(double q2)
{
    return Gamma_0(q2) * beta2(q2) *
            ( (H_p(q2) - H_m(q2))*H_0(q2).conjugate() -
            4.*(H_T_p(q2) - H_T_m(q2))*H_T_0(q2).conjugate() ).imag();
}

double myMVlnu::I_9(double q2)
{
    return -2. * Gamma_0(q2) * beta2(q2) *
            ( H_p(q2)*H_m(q2).conjugate() - 4.*H_T_p(q2)*H_T_m(q2).conjugate() ).imag();
}

double myMVlnu::integrateI(int i)
{
    updateParameters();

    double q_min = Mlep2;
    double q_max = (MM-MV)*(MM-MV);

    old_handler = gsl_set_error_handler_off();

    switch (i) {
        case 0:
            if (I1c_updated == 0) {
                FI = convertToGslFunction(boost::bind(&myMVlnu::I_1c, &(*this), _1));
                if (gsl_integration_cquad(&FI, q_min, q_max, 1.e-2, 1.e-1, w_I, &avaI, &errI, NULL) != 0) return std::numeric_limits<double>::quiet_NaN();
                cacheI1c = avaI;
                I1c_updated = 1;
            }
            return cacheI1c;
            break;
        case 1:
            if (I1s_updated == 0) {
                FI = convertToGslFunction(boost::bind(&myMVlnu::I_1s, &(*this), _1));
                if (gsl_integration_cquad(&FI, q_min, q_max, 1.e-2, 1.e-1, w_I, &avaI, &errI, NULL) != 0) return std::numeric_limits<double>::quiet_NaN();
                cacheI1s = avaI;
                I1s_updated = 1;
            }
            return cacheI1s;
            break;
        case 2:
            if (I2c_updated == 0) {
                FI = convertToGslFunction(boost::bind(&myMVlnu::I_2c, &(*this), _1));
                if (gsl_integration_cquad(&FI, q_min, q_max, 1.e-2, 1.e-1, w_I, &avaI, &errI, NULL) != 0) return std::numeric_limits<double>::quiet_NaN();
                cacheI2c = avaI;
                I2c_updated = 1;
            }
            return cacheI2c;
            break;
        case 3:
            if (I2s_updated == 0) {
                FI = convertToGslFunction(boost::bind(&myMVlnu::I_2s, &(*this), _1));
                if (gsl_integration_cquad(&FI, q_min, q_max, 1.e-2, 1.e-1, w_I, &avaI, &errI, NULL) != 0) return std::numeric_limits<double>::quiet_NaN();
                cacheI2s = avaI;
                I2s_updated = 1;
            }
            return cacheI2s;
            break;
        case 4:
            if (I3_updated == 0) {
                FI = convertToGslFunction(boost::bind(&myMVlnu::I_3, &(*this), _1));
                if (gsl_integration_cquad(&FI, q_min, q_max, 1.e-2, 1.e-1, w_I, &avaI, &errI, NULL) != 0) return std::numeric_limits<double>::quiet_NaN();
                cacheI3 = avaI;
                I3_updated = 1;
            }
            return cacheI3;
            break;
        case 5:
            if (I4_updated == 0) {
                FI = convertToGslFunction(boost::bind(&myMVlnu::I_4, &(*this), _1));
                if (gsl_integration_cquad(&FI, q_min, q_max, 1.e-2, 1.e-1, w_I, &avaI, &errI, NULL) != 0) return std::numeric_limits<double>::quiet_NaN();
                cacheI4 = avaI;
                I4_updated = 1;
            }
            return cacheI4;
            break;
        case 6:
            if (I5_updated == 0) {
                FI = convertToGslFunction(boost::bind(&myMVlnu::I_5, &(*this), _1));
                if (gsl_integration_cquad(&FI, q_min, q_max, 1.e-2, 1.e-1, w_I, &avaI, &errI, NULL) != 0) return std::numeric_limits<double>::quiet_NaN();
                cacheI5 = avaI;
                I5_updated = 1;
            }
            return cacheI5;
            break;
        case 7:
            if (I6c_updated == 0) {
                FI = convertToGslFunction(boost::bind(&myMVlnu::I_6c, &(*this), _1));
                if (gsl_integration_cquad(&FI, q_min, q_max, 1.e-2, 1.e-1, w_I, &avaI, &errI, NULL) != 0) return std::numeric_limits<double>::quiet_NaN();
                cacheI6c = avaI;
                I6c_updated = 1;
            }
            return cacheI6c;
            break;
        case 8:
            if (I6s_updated == 0) {
                FI = convertToGslFunction(boost::bind(&myMVlnu::I_6s, &(*this), _1));
                if (gsl_integration_cquad(&FI, q_min, q_max, 1.e-2, 1.e-1, w_I, &avaI, &errI, NULL) != 0) return std::numeric_limits<double>::quiet_NaN();
                cacheI6s = avaI;
                I6s_updated = 1;
            }
            return cacheI6s;
            break;
        case 9:
            if (I7_updated == 0) {
                FI = convertToGslFunction(boost::bind(&myMVlnu::I_7, &(*this), _1));
                if (gsl_integration_cquad(&FI, q_min, q_max, 1.e-2, 1.e-1, w_I, &avaI, &errI, NULL) != 0) return std::numeric_limits<double>::quiet_NaN();
                cacheI7 = avaI;
                I7_updated = 1;
            }
            return cacheI7;
            break;
        case 10:
            if (I8_updated == 0) {
                FI = convertToGslFunction(boost::bind(&myMVlnu::I_8, &(*this), _1));
                if (gsl_integration_cquad(&FI, q_min, q_max, 1.e-2, 1.e-1, w_I, &avaI, &errI, NULL) != 0) return std::numeric_limits<double>::quiet_NaN();
                cacheI8 = avaI;
                I8_updated = 1;
            }
            return cacheI8;
            break;
        case 11:
            if (I9_updated == 0) {
                FI = convertToGslFunction(boost::bind(&myMVlnu::I_9, &(*this), _1));
                if (gsl_integration_cquad(&FI, q_min, q_max, 1.e-2, 1.e-1, w_I, &avaI, &errI, NULL) != 0) return std::numeric_limits<double>::quiet_NaN();
                cacheI9 = avaI;
                I9_updated = 1;
            }
            return cacheI9;
            break;
        case 12:
            if (I2n_updated == 0) {
                FI = convertToGslFunction(boost::bind(&myMVlnu::I_2n, &(*this), _1));
                if (gsl_integration_cquad(&FI, q_min, q_max, 1.e-2, 1.e-1, w_I, &avaI, &errI, NULL) != 0) return std::numeric_limits<double>::quiet_NaN();
                cacheI2n = avaI;
                I2n_updated = 1;
            }
            return cacheI2n;
            break;
        default:
            std::stringstream out;
            out << i;
            throw std::runtime_error("myMVlnu::integrateI: index " + out.str() + " not implemented");
    }

    gsl_set_error_handler(old_handler);

}
