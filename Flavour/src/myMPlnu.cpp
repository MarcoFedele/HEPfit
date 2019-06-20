/*
 * Copyright (C) 2018 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "StandardModel.h"
#include "myMPlnu.h"
#include "std_make_vector.h"
#include "gslpp_function_adapter.h"
#include <boost/bind.hpp>

myMPlnu::myMPlnu(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson pseudoscalar_i, QCD::lepton lep_i)
: mySM(SM_i),
meson_cache(2, 0.),
f_p_Vcache(4, 0.),
f_0_Vcache(3, 0.),
quark_cache(2, 0.),
N_cache(3, 0.)
{
    lep = lep_i;
    meson = meson_i;
    pseudoscalar = pseudoscalar_i;
    MSflag = false;
    ChiralBasisflag = false;
    gSL814gTflag = false;
    gSLm814gTflag = false;

    a_updated = 0;
    b_updated = 0;
    c_updated = 0;
    Gamma_p_updated = 0;
    Gamma_m_updated = 0;

    w_I = gsl_integration_cquad_workspace_alloc (100);
}

myMPlnu::~myMPlnu()
{}


std::vector<std::string> myMPlnu::initializemyMPlnuParameters()
{
    MSflag = mySM.getFlavour().getFlagMS();
    ChiralBasisflag = mySM.getFlavour().getFlagChiralBasis();

    if (MSflag) {
        if (pseudoscalar == StandardModel::D_P && lep == StandardModel::TAU) mplnuParameters = make_vector<std::string>()
            << "fp0" << "fT0"
            << "regV" << "regS" << "regT" << "imgV" << "imgS" << "imgT"
            << "regVL" << "regVR" << "regSL" << "regSR" << "imgVL" << "imgVR" << "imgSL" << "imgSR" ;
        else if (pseudoscalar == StandardModel::D_P && (lep == StandardModel::MU || lep == StandardModel::ELECTRON))
            mplnuParameters = make_vector<std::string>() << "fp0" << "fT0" ;
        else {
            std::stringstream out;
            out << pseudoscalar;
            throw std::runtime_error("myMPlnu: pseudoscalar " + out.str() + " not implemented");
        }
    } else {
        if (pseudoscalar == StandardModel::D_P && lep == StandardModel::TAU) mplnuParameters = make_vector<std::string>()
            << "ap0" << "ap1" << "ap2" << "ap3" << "a01" << "a02" << "a03" << "RT_err"
            << "regV" << "regS" << "regT" << "imgV" << "imgS" << "imgT"
            << "regVL" << "regVR" << "regSL" << "regSR" << "imgVL" << "imgVR" << "imgSL" << "imgSR" ;
        else if (pseudoscalar == StandardModel::D_P && (lep == StandardModel::MU || lep == StandardModel::ELECTRON))
            mplnuParameters = make_vector<std::string>() << "ap0" << "ap1" << "ap2" << "ap3" << "a01" << "a02" << "a03" << "RT_err" ;
        else {
            std::stringstream out;
            out << pseudoscalar;
            throw std::runtime_error("myMPlnu: pseudoscalar " + out.str() + " not implemented");
        }
    }

    mySM.initializeMeson(meson);
    mySM.initializeMeson(pseudoscalar);
    return mplnuParameters;
}



void myMPlnu::checkCache()
{

    if (gV == gV_cache) {
        gV_updated = 1;
    } else {
        gV_updated = 0;
        gV_cache = gV;
    }

    if (gS == gS_cache) {
        gS_updated = 1;
    } else {
        gS_updated = 0;
        gS_cache = gS;
    }

    if (gT == gT_cache) {
        gT_updated = 1;
    } else {
        gT_updated = 0;
        gT_cache = gT;
    }

    if (MM == meson_cache(0) && MP == meson_cache(1)) {
        meson_cache = 1;
    } else {
        meson_updated = 0;
        meson_cache(0) = MM;
        meson_cache(1) = MP;
    }

    if (!MSflag) {
        if (ap0 == f_p_Vcache(0) && ap1 == f_p_Vcache(1) && ap2 == f_p_Vcache(2) && ap3 == f_p_Vcache(3)) {
            f_p_updated = meson_updated;
        } else {
            f_p_updated = 0;
            f_p_Vcache(0) = ap0;
            f_p_Vcache(1) = ap1;
            f_p_Vcache(2) = ap2;
            f_p_Vcache(3) = ap3;
        }

        if (a01 == f_0_Vcache(0) && a02 == f_0_Vcache(1) && a03 == f_0_Vcache(2)) {
            f_0_updated = f_p_updated;
        } else {
            f_0_updated = 0;
            f_0_Vcache(0) = a01;
            f_0_Vcache(1) = a02;
            f_0_Vcache(2) = a03;
        }

        if (RT_err == R_T_cache) {
            R_T_updated = 1;
        } else {
            R_T_updated = 0;
            R_T_cache = RT_err;
        }
        f_T_updated = f_p_updated*R_T_updated;

    } else {
        if (f00 == f_0_cache) {
            f_0_updated = 1;
        } else {
            f_0_updated = 0;
            f_0_cache = f00;
        }

        if (fp0 == f_p_cache) {
            f_p_updated = 1;
        } else {
            f_p_updated = 0;
            f_p_cache = fp0;
        }

        if (fT0 == f_T_cache) {
            f_T_updated = 1;
        } else {
            f_T_updated = 0;
            f_T_cache = fT0;
        }
    }

    if (Mb == quark_cache(0) && Mc == quark_cache(1)) {
        quark_cache = 1;
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

    H_V_t_updated = gV_updated*meson_updated*f_0_updated;
    H_V_updated = gV_updated*meson_updated*f_p_updated;
    H_S_updated = gS_updated*meson_updated*quark_updated*f_0_updated;
    H_T_updated = gT_updated*meson_updated*f_T_updated;

    H_t_updated = H_V_t_updated*H_S_updated;

    a_updated = N_updated*H_V_updated*H_T_updated*H_t_updated;
    b_updated = N_updated*H_V_updated*H_T_updated*H_t_updated;
    c_updated = N_updated*H_V_updated*H_T_updated;
    Gamma_p_updated = N_updated*H_V_updated*H_T_updated*H_t_updated;
    Gamma_m_updated = N_updated*H_V_updated*H_T_updated;

    return;
}



void myMPlnu::updateParameters()
{
    if (!mySM.getFlavour().getUpdateFlag(meson, pseudoscalar, lep)) return;

    GF = mySM.getGF();
    Mlep = mySM.getLeptons(lep).getMass();
    MM = mySM.getMesons(meson).getMass();
    MP = mySM.getMesons(pseudoscalar).getMass();
    Mb = mySM.getQuarks(QCD::BOTTOM).getMass();
    Mc = mySM.getQuarks(QCD::CHARM).getMass();
    Vcb = mySM.getCKM().getV_cb();
    //Vcb = 0.0411;

    if (lep == StandardModel::TAU) {
        if (ChiralBasisflag) {
            gV = mySM.getOptionalParameter("regVR") + gslpp::complex::i()*mySM.getOptionalParameter("imgVR")
                + mySM.getOptionalParameter("regVL") + gslpp::complex::i()*mySM.getOptionalParameter("imgVL");
            gS = mySM.getOptionalParameter("regSR") + gslpp::complex::i()*mySM.getOptionalParameter("imgSR")
                + mySM.getOptionalParameter("regSL") + gslpp::complex::i()*mySM.getOptionalParameter("imgSL");
            gT = 2.*(mySM.getOptionalParameter("regT") + gslpp::complex::i()*mySM.getOptionalParameter("imgT"));
            if (gSL814gTflag)
              gT = 2.*(mySM.getOptionalParameter("regSL") + gslpp::complex::i()*mySM.getOptionalParameter("imgSL"))/8.14;
            if (gSLm814gTflag)
              gT = -2.*(mySM.getOptionalParameter("regSL") + gslpp::complex::i()*mySM.getOptionalParameter("imgSL"))/8.14;
        }
        else {
            gV = mySM.getOptionalParameter("regV") + gslpp::complex::i()*mySM.getOptionalParameter("imgV");
            gS = mySM.getOptionalParameter("regS") + gslpp::complex::i()*mySM.getOptionalParameter("imgS");
            gT = mySM.getOptionalParameter("regT") + gslpp::complex::i()*mySM.getOptionalParameter("imgT");
        }
    } else {
        gV = 0;
        gS = 0;
        gT = 0;
    }

    if (MSflag) {
        switch (pseudoscalar) {
        case StandardModel::D_P:
            fp0 = mySM.getOptionalParameter("fp0");
            f00 = fp0;
            fT0 = mySM.getOptionalParameter("fT0");
            break;
        default:
            std::stringstream out;
            out << pseudoscalar;
            throw std::runtime_error("myMPlnu: pseudoscalar " + out.str() + " not implemented");
        }
    } else {
        switch (pseudoscalar) {
        case StandardModel::D_P:
            ap0 = mySM.getOptionalParameter("ap0");
            ap1 = mySM.getOptionalParameter("ap1");
            ap2 = mySM.getOptionalParameter("ap2");
            ap3 = mySM.getOptionalParameter("ap3");
            a01 = mySM.getOptionalParameter("a01");
            a02 = mySM.getOptionalParameter("a02");
            a03 = mySM.getOptionalParameter("a03");
            RT_err = mySM.getOptionalParameter("RT_err");
            break;
        default:
            std::stringstream out;
            out << pseudoscalar;
            throw std::runtime_error("myMPlnu: pseudoscalar " + out.str() + " not implemented");
        }
    }

    checkCache();

    MM2 = MM*MM;
    MP2 = MP*MP;
    Mlep2 = Mlep*Mlep;
    NN = (GF * GF * Vcb.abs2()) / (256. * M_PI*M_PI*M_PI * MM2*MM) ;

    /*
    std::cout << "f_p(4.) = " << f_p(4.) << std::endl;
    std::cout << "f_0(4.) = " << f_0(4.) << std::endl;

    std::cout << "lep = " << lep << std::endl;
    std::cout << "a_th(4.) = " << a_theta(4.) << std::endl;
    std::cout << "b_th(4.) = " << b_theta(4.) << std::endl;
    std::cout << "c_th(4.) = " << c_theta(4.) << std::endl;
    std::cout << "G_p(4.) = " << Gamma_p(4.) << std::endl;*/

    mySM.getFlavour().setUpdateFlag(meson, pseudoscalar, lep, false);

    return;
}

double myMPlnu::FF_fit1(double q2, double f_0, double s_1)
{
    double MR_2 = 6.4*6.4;
    return f_0 / (1. - q2 / MR_2) / (1. - s_1 * q2 / MR_2);
}

double myMPlnu::FF_fit2(double q2, double f_0, double s_1)
{
    double MR_2 = 6.4*6.4;
    return f_0 / (1. - s_1 * q2 / MR_2);
}

double myMPlnu::w(double q2)
{
    return (MM2 + MP2 - q2)/(2. * MM * MP);
}

double myMPlnu::z(double q2)
{
    return (sqrt(w(q2) + 1.) - sqrt(2.))/(sqrt(w(q2) + 1.) + sqrt(2.));
}

double myMPlnu::phip(double q2)
{
    double r = MP/MM;
    double Phip = 1.1213;

    return Phip * pow( 1. + z(q2), 2.) * pow( 1. - z(q2), 1./2.) *
            pow( (1. + r)*(1. - z(q2)) + 2.*sqrt(r)*(1. + z(q2)) , -5.);
}

double myMPlnu::phi0(double q2)
{
    double r = MP/MM;
    double Phi0 = 0.5299;

    return Phi0 * (1. + z(q2)) * pow( 1. - z(q2), 3./2.) *
            pow( (1. + r)*(1. - z(q2)) + 2.*sqrt(r)*(1. + z(q2)) , -4.);
}

double myMPlnu::R_T(double q2)
{
    return RT_err * ( 1.0645 - 2.6694e-5*q2 - 3.8882e-6*q2*q2 -6.4292e-8*q2*q2*q2 );
}

double myMPlnu::f_p(double q2)
{
    if (MSflag) return FF_fit1(q2, fp0, 0.57);
    else {
        double zz = z(q2);
        return 1./phip(q2) * ( ap0 + ap1*zz + ap2*zz*zz + ap3*zz*zz*zz );
    }
}

double myMPlnu::f_0(double q2)
{
    if (MSflag) return FF_fit2(q2, f00, 0.78);
    else {
        double zz = z(q2);
        double z0 = z(0.);

        double a00 = phi0(0.)*f_p(0.) - ( a01*z0 + a02*z0*z0 + a03*z0*z0*z0 );

        return 1./phi0(q2) * ( a00 + a01*zz + a02*zz*zz + a03*zz*zz*zz );
    }
}

double myMPlnu::f_T(double q2)
{
    if (MSflag) return FF_fit1(q2, fT0, 0.56);
    else return f_p(q2) * R_T(q2);
}

double myMPlnu::lambda(double q2)
{
    return (MM2*MM2 + q2 * q2 + MP2*MP2 - 2. * (MP2 * q2 + MM2 * (q2 + MP2)));
}

gslpp::complex myMPlnu::H_V(double q2)
{
    return (1. + gV) * sqrt(lambda(q2)/q2) * f_p(q2);
}

gslpp::complex myMPlnu::H_V_t(double q2)
{
    return (1. + gV) * (MM2 - MP2)/sqrt(q2) * f_0(q2);
}

gslpp::complex myMPlnu::H_S(double q2)
{
    return gS/(Mb - Mc) * (MM2 - MP2) * f_0(q2);
}

gslpp::complex myMPlnu::H_t(double q2)
{
    return H_V_t(q2) + sqrt(q2)/Mlep * H_S(q2);
}

gslpp::complex myMPlnu::H_T(double q2)
{
    return -gT * sqrt(lambda(q2))/(MM + MP) * f_T(q2);
}

double myMPlnu::beta2(double q2)
{
    return 1. - Mlep2/q2;
}

double myMPlnu::Gamma_0(double q2)
{
    return NN * q2 * sqrt(lambda(q2)) * beta2(q2)*beta2(q2);
}

double myMPlnu::a_theta(double q2)
{
    return Gamma_0(q2) * ( H_V(q2).abs2() + 4.*Mlep2/q2*H_T(q2).abs2()
            - 4.*Mlep/sqrt(q2)*(H_V(q2)*H_T(q2).conjugate()).real() + Mlep2/q2*H_t(q2).abs2() );
}

double myMPlnu::b_theta(double q2)
{
    return Gamma_0(q2) * ( 2.*Mlep2/q2*(H_V(q2)*H_t(q2).conjugate()).real()
            - 4.*Mlep/sqrt(q2)*(H_T(q2)*H_t(q2).conjugate()).real());
}

double myMPlnu::c_theta(double q2)
{
    return -Gamma_0(q2) * beta2(q2) * ( H_V(q2).abs2() - 4.*H_T(q2).abs2() );
}

double myMPlnu::Gamma_p(double q2)
{
    return Gamma_0(q2) * 4./6. * Mlep2/q2 * ( H_V(q2).abs2() + 4.*q2/Mlep2*H_T(q2).abs2()
            - 4.*sqrt(q2)/Mlep*(H_V(q2)*H_T(q2).conjugate()).real() + 3.*H_t(q2).abs2());
}

double myMPlnu::Gamma_m(double q2)
{
    return Gamma_0(q2) * 4./6. * ( H_V(q2).abs2() + 4.*Mlep2/q2*H_T(q2).abs2()
            - 4.*Mlep/sqrt(q2)*(H_V(q2)*H_T(q2).conjugate()).real() );
}

double myMPlnu::integrateI(int i)
{
    updateParameters();

    double q_min = Mlep2;
    double q_max = (MM-MP)*(MM-MP);

    old_handler = gsl_set_error_handler_off();

    switch (i) {
        case 0:
            if (a_updated == 0) {
                FI = convertToGslFunction(boost::bind(&myMPlnu::a_theta, &(*this), _1));
                if (gsl_integration_cquad(&FI, q_min, q_max, 1.e-2, 1.e-1, w_I, &avaI, &errI, NULL) != 0) return std::numeric_limits<double>::quiet_NaN();
                cachea = avaI;
                a_updated = 1;
            }
            return cachea;
            break;
        case 1:
            if (b_updated == 0) {
                FI = convertToGslFunction(boost::bind(&myMPlnu::b_theta, &(*this), _1));
                if (gsl_integration_cquad(&FI, q_min, q_max, 1.e-2, 1.e-1, w_I, &avaI, &errI, NULL) != 0) return std::numeric_limits<double>::quiet_NaN();
                cacheb = avaI;
                b_updated = 1;
            }
            return cacheb;
            break;
        case 2:
            if (c_updated == 0) {
                FI = convertToGslFunction(boost::bind(&myMPlnu::c_theta, &(*this), _1));
                if (gsl_integration_cquad(&FI, q_min, q_max, 1.e-2, 1.e-1, w_I, &avaI, &errI, NULL) != 0) return std::numeric_limits<double>::quiet_NaN();
                cachec = avaI;
                c_updated = 1;
            }
            return cachec;
            break;
        case 3:
            if (Gamma_p_updated == 0) {
                FI = convertToGslFunction(boost::bind(&myMPlnu::Gamma_p, &(*this), _1));
                if (gsl_integration_cquad(&FI, q_min, q_max, 1.e-2, 1.e-1, w_I, &avaI, &errI, NULL) != 0) return std::numeric_limits<double>::quiet_NaN();
                cachegp = avaI;
                Gamma_p_updated = 1;
            }
            return cachegp;
        case 4:
            if (Gamma_m_updated == 0) {
                FI = convertToGslFunction(boost::bind(&myMPlnu::Gamma_m, &(*this), _1));
                if (gsl_integration_cquad(&FI, q_min, q_max, 1.e-2, 1.e-1, w_I, &avaI, &errI, NULL) != 0) return std::numeric_limits<double>::quiet_NaN();
                cachegm = avaI;
                Gamma_m_updated = 1;
            }
            return cachegm;
            break;
        default:
            std::stringstream out;
            out << i;
            throw std::runtime_error("myMPlnu::integrateI: index " + out.str() + " not implemented");
    }

    gsl_set_error_handler(old_handler);

}
