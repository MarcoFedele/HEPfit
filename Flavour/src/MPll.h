/* 
 * Copyright (C) 2014 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef MPLL_H
#define	MPLL_H

#include <math.h>
#include "Flavour.h"
#include "MVll.h"
#include <StandardModel.h>
#include <ThObservable.h>
#include <gsl/gsl_integration.h>
#include <assert.h>


#define CUTOFF 10    //cutoff between LCSR and lattice values for Form Factors, in GeV^2

/*******************************************************************************
 * GSL Function Conversion BEGIN                                                  *
 * ****************************************************************************/

// Option 1. To be called with:
//gsl_function_pp Fp( boost::bind(&Class::member_function, &class, _1) );
//gsl_function *F = static_cast<gsl_function*>(&Fp);

//class gsl_function_pp : public gsl_function
//{
//public:
//    gsl_function_pp(std::function<double(double)> const& func) : _func(func){
//        function=&gsl_function_pp::invoke;
//        params=this;
//    }
//private:
//    std::function<double(double)> _func;
//    static double invoke(double x, void *params) {
//        return static_cast<gsl_function_pp*>(params)->_func(x);
//    }
//};


//Option 2. To be used with:
//gslFunction gslF = convertToGslFunction( boost::bind( &Class::member_function, &class, _1 ) );

/*template<class F>
static double gslFunctionAdapter( double x, void* p)
{
    // Here I do recover the "right" pointer, safer to use static_cast
    // than reinterpret_cast.
    F* function = static_cast<F*>( p );
    return (*function)( x );
}

template<class F>
gsl_function convertToGslFunction( const F& f )
{
    gsl_function gslFunction;
    
    const void* p = &f;
    assert (p != 0);
    
    gslFunction.function = &gslFunctionAdapter<F>;
    // Just to eliminate the const.
    gslFunction.params = const_cast<void*>( p );
    
    return gslFunction;
}*/

//Option 3. To be used with:
//Class* ptr2 = class;
//auto ptr = [=](double x)->double{return ptr2->member_function(x);};
//gsl_function_p<decltype(ptr)> Fp(ptr);
//gsl_function *F = static_cast<gsl_function*>(&Fp);

//template< typename F >
//class gsl_function_p : public gsl_function {
//public:
//    gsl_function_p(const F& func) : _func(func) {
//        function = &gsl_function_p::invoke;
//        params=this;
//    }
//private:
//    const F& _func;
//    static double invoke(double x, void *params) {
//        return static_cast<gsl_function_p*>(params)->_func(x);
//    }
//};

/*******************************************************************************
 * GSL Function conversion END                                                     *
 * ****************************************************************************/

/**
 * @class MPll
 * @ingroup flavour
 * @brief A class for the decay B -> K^*ll. 
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class MPll : public ThObservable {
public:
    MPll(const StandardModel& SM_i, StandardModel::lepton lep_i = StandardModel::MU);
    virtual ~MPll();
    void updateParameters();
    //void checkCache( double qmin, double qmax);
    virtual double computeThValue()=0;
    
    double GF;            //Fermi constant
    double ale;           //alpha electromagnetic
    double Mlep;          //muon mass
    double Me;            //electron mass
    double Mmu;           //muon mass
    double MB;            //B meson mass
    double MK;            //K star meson mass
    double Mb;            //b quark mass
    double mu_b;          //b mass scale
    double width_Bd;      //B meosn width
    double Ms;            //s quark mass
    double MW;            //W boson mass
    complex lambda_t;     //Vckm factor
    double b;             //BF of the decay K^* -> K pi
    complex h[3];         //parameter that contains the contribution from the hadronic hamiltonian  
    double q2;            //q^2 of the decay
    
    /*lattice fit parameters*/
    double a_0V, a_1V, dmV;
    double a_0A0, a_1A0, dmA0;
    double a_0A1, a_1A1, dmA1;
    double a_0A12, a_1A12, dmA12;
    double a_0T1, a_1T1, dmT1;
    double a_0T2, a_1T2, dmT2;
    double a_0T23, a_1T23, dmT23;
    
    /*LCSR fit parameters*/
    double r_1V, r_2V, m_RV, m_fit2V;
    double r_1A0, r_2A0, m_RA0, m_fit2A0;
    double r_2A1, m_fit2A1;
    double r_1A2, r_2A2, m_fit2A2;
    double r_1T1, r_2T1, m_RT1, m_fit2T1;
    double r_2T2, m_fit2T2;
    double r_1T3t, r_2T3t, m_fit2T3t;

    vector<complex> ** allcoeff;
    vector<complex> ** allcoeffprime;
    
    
    
    /**
    * @brief \f$ LCSR_fit1 \f$
    * @param[in] q2 q^2 of the decay
    * @param[in] r_1 fit parameter
    * @param[in] r_2 fit parameter
    * @param[in] m_R2 fit parameter
    * @param[in] m_fit2 fit parameter
    * @return return the first fit function from arXiv:hep-ph/0412079v1
    */
    double LCSR_fit1(double q2, double r_1, double r_2, double m_R2, double m_fit2);
    
    
    /**
    * @brief \f$ LCSR_fit2 \f$
    * @param[in] q2 q^2 of the decay
    * @param[in] r_1 fit parameter
    * @param[in] r_2 fit parameter
    * @param[in] m_fit2 fit parameter
    * @return return the second fit function from arXiv:hep-ph/0412079v1
    */
    double LCSR_fit2(double q2, double r_1, double r_2, double m_fit2);
    
    
    /**
    * @brief \f$ LCSR_fit3 \f$
    * @param[in] q2 q^2 of the decay
    * @param[in] r_2 fit parameter
    * @param[in] m_fit2 fit parameter
    * @return return the third fit function from arXiv:hep-ph/0412079v1
    */
    double LCSR_fit3(double q2, double r_2, double m_fit2);
    
    
    /**
    * @brief \f$ z \f$
    * @param[in] q2 q^2 of the decay
    * @return return the lattice parameter z from arXiv:1310.3722v3
    */
    double z(double q2);
    
    
    /**
    * @brief \f$ lat_fit(q^2) \f$
    * @param[in] q2 q^2 of the decay
    * @param[in] a_0 fit parameter
    * @param[in] a_1 fit parameter
    * @param[in] c_01 fit parameter
    * @param[in] c_01s fit parameter
    * @param[in] q2 q^2 of the decay
    * @return return the lattice pole factor P(q^2,dm) from arXiv:1310.3722v3
    */
    double lat_fit(double q2, double a_0, double a_1, double dm);
    
    
    /**
    * @brief \f$ V \f$
    * @param[in] q2 q^2 of the decay
    * @return return the transverse form factor V(q^2)
    */
    double V(double q2);

    
    /**
    * @brief \f$ A_0 \f$
    * @param[in] q2 q^2 of the decay
    * @return return the transverse form factor A_0(q^2)
    */
    double A_0(double q2);

    
    /**
    * @brief \f$ A_1 \f$
    * @param[in] q2 q^2 of the decay
    * @return return the transverse form factor A_1(q^2)
    */
    double A_1(double q2);

    
    /**
    * @brief \f$ A_2 \f$
    * @param[in] q2 q^2 of the decay
    * @return return the transverse form factor A_2(q^2)
    */
    double A_2(double q2);

    
    /**
    * @brief \f$ T_1 \f$
    * @param[in] q2 q^2 of the decay
    * @return return the transverse form factor T_1(q^2)
    */
    double T_1(double q2);

    
    /**
    * @brief \f$ V \f$
    * @param[in] q2 q^2 of the decay
    * @return return the transverse form factor V(q^2)
    */
    double T_2(double q2);

    
    /**
    * @brief \f$ T_3tilde \f$
    * @param[in] q2 q^2 of the decay
    * @return return the transverse form factor T_3tilde(q^2)
    */
    double T_3tilde(double q2);

    
    /**
    * @brief \f$ T_3 \f$
    * @param[in] q2 q^2 of the decay
    * @return return the transverse form factor T_3(q^2)
    */
    double T_3(double q2);
    
    
    /**
    * polarization   i
    *      0         0
    *      +         1
    *      -         2
    */

    
    /**
    * @brief \f$ V_L \f$
    * @param[in] i polarization
    * @param[in] q2 q^2 of the decay
    * @return return the helicity form factor V_L(lambda)
    */
    double V_L(int i, double q2);

    
    /**
    * @brief \f$ V_R \f$
    * @param[in] i polarization
    * @param[in] q2 q^2 of the decay
    * @return return the helicity form factor V_R(lambda)
    */
    double V_R(int i, double q2);


    /**
    * @brief \f$ T_L \f$
    * @param[in] i polarization
    * @param[in] q2 q^2 of the decay
    * @return return the helicity form factor T_L(lambda)
    */
    double T_L(int i, double q2);


    /**
    * @brief \f$ T_R \f$ 
    * @param[in] i polarization
    * @param[in] q2 q^2 of the decay
    * @return return the helicity form factor T_R(lambda)
    */
    double T_R(int i, double q2);


    /**
    * 
    * @brief \f$ S_L \f$
    * @param[in] q2 q^2 of the decay
    * @return return the helicity form factor S_L
    */
    double S_L(double q2);


    /**
    * 
    * @brief \f$ S_R \f$
    * @param[in] q2 q^2 of the decay
    * @return return the helicity form factor S_R
    */
    double S_R(double q2);
    
    
    /**
    * @brief \f$ N \f$ 
    * @return return the helicity amplitude normalization factor N
    */
    
    complex N();
    
    
    /**
    * @brief \f$ H_V(\lambda) \f$ 
    * @param[in] i polarization lambda
    * @param[in] q2 q^2 of the decay
    * @param[in] bar index to choose betwen regular coefficient (bar=0) and conjugated coefficient (bar=1)
    * @return return the helicity amplitude H_V(lambda)
    */
    gslpp::complex H_V(int i, double q2, int bar);


    /**
    * @brief \f$ H_A(\lambda) \f$ 
    * @param[in] i polarization lambda
    * @param[in] q2 q^2 of the decay
    * @param[in] bar index to choose betwen regular coefficient (bar=0) and conjugated coefficient (bar=1)
    * @return return the helicity amplitude H_A(lambda)
    */
    gslpp::complex H_A(int i, double q2, int bar);


    /**
    * @brief \f$ H_S \f$ 
    * @param[in] q2 q^2 of the decay
    * @param[in] bar index to choose betwen regular coefficient (bar=0) and conjugated coefficient (bar=1)
    * @return return the helicity amplitude H_S
    */
    gslpp::complex H_S(double q2, int bar);


    /**
    * @brief \f$ H_P \f$ 
    * @param[in] q2 q^2 of the decay
    * @param[in] bar index to choose betwen regular coefficient (bar=0) and conjugated coefficient (bar=1)
    * @return return the helicity amplitude H_P
    */
    gslpp::complex H_P(double q2, int bar);
    
    
    /**
    * @brief \f$ k^2 \f$ 
    * @param[in] q2 q^2 of the decay
    * @return return the square of the 3-momentum of the recoiling meson in the B rest frame
    */
    double k2 (double q2);
    
    
    /**
    * @brief \f$ beta \f$ 
    * @param[in] q2 q^2 of the decay
    * @return return the factor beta used in the angular coefficients I_i
    */
    double beta (double q2);
    
    
    /**
    * @brief \f$ lambda \f$ 
    * @param[in] q2 q^2 of the decay
    * @return return the factor lambda used in the angular coefficients I_i
    */
    double lambda(double q2);

    
    /**
    * @brief \f$ F \f$ 
    * @param[in] q2 q^2 of the decay
    * @param[in] b_i BF of the decay K* -> K pi
    * @return return the factor F used in the angular coefficients I_i
    */
    double F(double q2, double b_i);
    
    
    /**
    * i values:
    * 0 = 1c
    * 1 = 1s
    * 2 = 2c
    * 3 = 2s
    * 4 = 3
    * 5 = 4
    * 6 = 5
    * 7 = 6s
    * 8 = 6c
    * 9 = 7
    * 10 = 8
    * 11 = 9
    */
    
    
    /**
    * @brief \f$ I_{i} \f$ 
    * @param[in] i index of the angular coefficient
    * @param[in] q2 q^2 of the decay
    * @param[in] bar index to choose betwen regular coefficient (bar=0) and conjugated coefficient (bar=1)
    * @return return the angular coefficient I_i
    */
    double  I(int i, double q2, int bar);
    
    
    /**
    * @brief \f$ Sigma_{i} \f$ 
    * @param[in] i index of the angular coefficient I_i
    * @param[in] q2 q^2 of the decay
    * @return return the CP average Sigma_i
    */
    double Sigma(int i, double q2);
    
    
    /**
    * @brief \f$ Delta_{i} \f$ 
    * @param[in] i index of the angular coefficient I_i
    * @param[in] q2 q^2 of the decay
    * @return return the CP asymmetry Delta_i
    */
    double Delta(int i, double q2);
    
    /**
    * @brief \f$ Sigma_{1s} \f$ 
    * @param[in] q2 q^2 of the decay
    * @return return the CP average Sigma_1c
    */
    double getSigma0(double q2){
        return Sigma(0, q2);
    };
    
    /**
    * @brief \f$ Sigma_{1c} \f$ 
    * @param[in] q2 q^2 of the decay
    * @return return the CP average Sigma_1s
    */
    double getSigma1(double q2){
        return Sigma(1, q2);
    };
    
    /**
    * @brief \f$ Sigma_{2s} \f$ 
    * @param[in] q2 q^2 of the decay
    * @return return the CP average Sigma_2c
    */
    double getSigma2(double q2){
        return Sigma(2, q2);
    };
    
    /**
    * @brief \f$ Sigma_{2c} \f$ 
    * @param[in] q2 q^2 of the decay
    * @return return the CP average Sigma_2s
    */
    double getSigma3(double q2){
        return Sigma(3, q2);
    };
    
    /**
    * @brief \f$ Sigma_{3} \f$ 
    * @param[in] q2 q^2 of the decay
    * @return return the CP average Sigma_3
    */    
    double getSigma4(double q2){
        return Sigma(4, q2);
    };
    
    /**
    * @brief \f$ Sigma_{4} \f$ 
    * @param[in] q2 q^2 of the decay
    * @return return the CP average Sigma_4
    */
    double getSigma5(double q2){
        return Sigma(5, q2);
    };
    
    /**
    * @brief \f$ Sigma_{5} \f$ 
    * @param[in] q2 q^2 of the decay
    * @return return the CP average Sigma_5
    */
    double getSigma6(double q2){
        return Sigma(6, q2);
    };
    
    /**
    * @brief \f$ Sigma_{6s} \f$ 
    * @param[in] q2 q^2 of the decay
    * @return return the CP average Sigma_6s
    */
    double getSigma7(double q2){
        return Sigma(7, q2);
    };
    
    /**
    * @brief \f$ Sigma_{7} \f$ 
    * @param[in] q2 q^2 of the decay
    * @return return the CP average Sigma_7
    */
    double getSigma9(double q2){
        return Sigma(9, q2);
    };
    
    /**
    * @brief \f$ Sigma_{9} \f$ 
    * @param[in] q2 q^2 of the decay
    * @return return the CP average Sigma_9
    */
    double getSigma11(double q2){
        return Sigma(11, q2);
    };
    
    /**
    * @brief \f$ Delta_{1s} \f$ 
    * @param[in] q2 q^2 of the decay
    * @return return the CP asymmetry Delta_1s
    */
    double getDelta0(double q2){
        return Delta(0, q2);
    };
    
    /**
    * @brief \f$ Delta_{1c} \f$ 
    * @param[in] q2 q^2 of the decay
    * @return return the CP asymmetry Delta_1c
    */
    double getDelta1(double q2){
        return Delta(1, q2);
    };
    
    /**
    * @brief \f$ Delta_{2s} \f$ 
    * @param[in] q2 q^2 of the decay
    * @return return the CP asymmetry Delta_2s
    */    
    double getDelta2(double q2){
        return Delta(2, q2);
    };
    
    /**
    * @brief \f$ Delta_{2c} \f$ 
    * @param[in] q2 q^2 of the decay
    * @return return the CP asymmetry Delta_2c
    */
    double getDelta3(double q2){
        return Delta(3, q2);
    };
    
    /**
    * @brief \f$ Delta_{9} \f$ 
    * @param[in] q2 q^2 of the decay
    * @return return the CP asymmetry Delta_9
    */
    double getDelta11(double q2){
        return Delta(11, q2);
    };
    
    /**
    * @brief \f$ |H_V(0)|^2 \f$ 
    * @param[in] q2 q^2 of the decay
    * @return return the square of the absolute value of the helicity amplitude H_V(0)
    */
    double getHV0_abs2(double q2){
        return H_V(0, q2, 0).abs2();
    };
    
    /**
    * @brief \f$ |H_V(1)|^2 \f$ 
    * @param[in] q2 q^2 of the decay
    * @return return the square of the absolute value of the helicity amplitude H_V(1)
    */
    double getHV1_abs2(double q2){
        return H_V(1, q2, 0).abs2();
    };
    
    /**
    * @brief \f$ |H_V(2)|^2 \f$ 
    * @param[in] q2 q^2 of the decay
    * @return return the square of the absolute value of the helicity amplitude H_V(2)
    */
    double getHV2_abs2(double q2){
        return H_V(2, q2, 0).abs2();
    };
    
    /**
    * @brief \f$ |H_A(0)|^2 \f$ 
    * @param[in] q2 q^2 of the decay
    * @return return the square of the absolute value of the helicity amplitude H_A(0)
    */
    double getHA0_abs2(double q2){
        return H_A(0, q2, 0).abs2();
    };
    
    /**
    * @brief \f$ |H_A(1)|^2 \f$ 
    * @param[in] q2 q^2 of the decay
    * @return return the square of the absolute value of the helicity amplitude H_A(1)
    */
    double getHA1_abs2(double q2){
        return H_A(1, q2, 0).abs2();
    };
    
    /**
    * @brief \f$ |H_A(2)|^2 \f$ 
    * @param[in] q2 q^2 of the decay
    * @return return the square of the absolute value of the helicity amplitude H_A(2)
    */
    double getHA2_abs2(double q2){
        return H_A(2, q2, 0).abs2();
    };
    
    /**
    * @brief \f$ |H_S|^2 \f$ 
    * @param[in] q2 q^2 of the decay
    * @return return the square of the absolute value of the helicity amplitude H_S
    */
    double getHS_abs2(double q2){
        return H_S(q2, 0).abs2();
    }
    
    /**
    * @brief \f$ |H_P|^2 \f$ 
    * @param[in] q2 q^2 of the decay
    * @return return the square of the absolute value of the helicity amplitude H_P
    */
    double getHP_abs2(double q2){
        return H_P(q2, 0).abs2();
    }
    
    
    double getFactor(double q2){
        return q2/(2*Mlep*Mlep);
    }

private:
    const StandardModel& mySM;
    StandardModel::lepton lep;
/*    int it;
    
protected:
    
    unsigned int N_updated;
    vector<double> N_cache;
    complex Nc_cache;
    
    unsigned int V_updated;
    vector<double> V_cache;
    
    unsigned int A0_updated;
    vector<double> A0_cache;
    
    unsigned int A1_updated;
    vector<double> A1_cache;
    
    unsigned int A2_updated;
    vector<double> A2_cache;
    
    unsigned int T1_updated;
    vector<double> T1_cache;
    
    unsigned int T2_updated;
    vector<double> T2_cache;
    
    unsigned int T3t_updated;
    vector<double> T3t_cache;
    
    unsigned int T3_updated;
    
    unsigned int k2_updated;
    vector<double> k2_cache;
    
    unsigned int z_updated;
    
    unsigned int lambda_updated;
    double lambda_cache;
    
    unsigned int beta_updated;
    double beta_cache;
    
    unsigned int F_updated;
    
    unsigned int VL1_updated;
    unsigned int VL2_updated;
    
    unsigned int TL1_updated;
    unsigned int TL2_updated;

    unsigned int VR1_updated;
    unsigned int VR2_updated;
    
    unsigned int TR1_updated;
    unsigned int TR2_updated;
    
    unsigned int VL0_updated;
    vector<double> VL0_cache;
    
    unsigned int TL0_updated;
    vector<double> TL0_cache;
    
    unsigned int VR0_updated;
    
    unsigned int TR0_updated;
    
    unsigned int SL_updated;
    vector<double> SL_cache;
    
    unsigned int SR_updated;
    
    unsigned int C_7_updated;
    complex C_7_cache;

    unsigned int C_9_updated;
    complex C_9_cache;
    
    unsigned int C_10_updated;
    complex C_10_cache;
    
    unsigned int C_7p_updated;
    complex C_7p_cache;
    
    unsigned int C_9p_updated;
    complex C_9p_cache;
    
    unsigned int C_10p_updated;
    complex C_10p_cache;
    
    unsigned int C_S_updated;
    complex C_S_cache;
    
    unsigned int C_P_updated;
    complex C_P_cache;
    
    unsigned int C_Sp_updated;
    complex C_Sp_cache;
    
    unsigned int C_Pp_updated;
    complex C_Pp_cache;
    
    unsigned int H_V0updated;
    vector<double> H_V0cache;
    complex H_V0Ccache;
    
    unsigned int H_V1updated;
    vector<double> H_V1cache;
    complex H_V1Ccache;
    
    unsigned int H_V2updated;
    vector<double> H_V2cache;
    complex H_V2Ccache;
    
    unsigned int H_A0updated;
    unsigned int H_A1updated;
    unsigned int H_A2updated;
    
    unsigned int H_Supdated;
    vector<double> H_Scache;
    
    unsigned int H_Pupdated;
    vector<double> H_Pcache;
    
    unsigned int I0_updated;
    unsigned int I1_updated;
    unsigned int I2_updated;
    unsigned int I3_updated;
    unsigned int I4_updated;
    unsigned int I5_updated;
    unsigned int I6_updated;
    unsigned int I7_updated;
    unsigned int I8_updated;
    unsigned int I9_updated;
    unsigned int I10_updated;
    unsigned int I11_updated;*/
};



/**
 * @class GammaPrime
 * @ingroup flavour
 * @brief A class for the clean observable Gamma'. 
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class GammaPrime : public MPll{
public:
    
    /**
    * @brief \f$ Gamma' \f$ 
    */
    GammaPrime(const StandardModel& SM_i, StandardModel::lepton lep_i = StandardModel::MU);
    
    /**
    * @return return the clean observable Gamma'
    */
    double computeGammaPrime(double qmin, double qmax);
    double computeThValue ();
    
protected:
    
    
private:
    gsl_function F1, F2, F3, F4;
    double avaSigma0, errSigma0, avaSigma1, errSigma1, avaSigma2, errSigma2, avaSigma3, errSigma3;
};


/**
 * @class A_FB
 * @ingroup flavour
 * @brief A class for the clean observable A_{FB}. 
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class A_FB : public GammaPrime{
public:
    
    /**
    * @brief \f$ Gamma' \f$ 
    */
    A_FB(const StandardModel& SM_i, StandardModel::lepton lep_i = StandardModel::MU);
    
    /**
    * @return return the clean observable Gamma'
    */
    double computeThValue ();
protected:
    
    
private:
    gsl_function F5;
    double avaSigma7, errSigma7;
};


/**
 * @class BF
 * @ingroup flavour
 * @brief A class for the Branching Fraction. 
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class BR_MPll : public GammaPrime{
public:
    
    /**
    * @brief \f$B\to K^* l^+l^-\f$ 
    */
    BR_MPll(const StandardModel& SM_i, StandardModel::lepton lep_i = StandardModel::MU);
    
    /**
    * @return the branching fraction of \f$B\to K^* l^+l^-\f$
    */
    double computeThValue ();
};


/**
 * @class ACP
 * @ingroup flavour
 * @brief A class for the clean observable Gamma'. 
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class ACP : public GammaPrime{
public:
    
    /**
    * @brief \f$ A_{CP} \f$ 
    */
    ACP(const StandardModel& SM_i, StandardModel::lepton lep_i = StandardModel::MU);
    
    /**
    * @return return the clean observable ACP
    */
    double computeThValue ();
    
protected:
    
    
private:
    gsl_function F1, F2, F3, F4;
    double avaDelta0, errDelta0, avaDelta1, errDelta1, avaDelta2, errDelta2, avaDelta3, errDelta3;
};


/**
 * @class P3CP
 * @ingroup flavour
 * @brief A class for the clean observable Gamma'. 
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class P3CP : public MPll{
public:
    
    /**
    * @brief \f$ P_3^{CP} \f$ 
    */
    P3CP(const StandardModel& SM_i, StandardModel::lepton lep_i = StandardModel::MU);
    
    /**
    * @return return the clean observable P3CP
    */
    double computeThValue ();
    
protected:
    
    
private:
    gsl_function F1, F2;
    double avaDelta11, errDelta11, avaSigma3, errSigma3;
};


/**
 * @class F_L
 * @ingroup flavour
 * @brief A class for the clean observable F_L. 
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class F_L : public MPll{
public:
    
    /**
    * @brief \f$ F_L \f$ 
    */
    F_L(const StandardModel& SM_i, StandardModel::lepton lep_i = StandardModel::MU);

    
    /**
    * @return return the clean observable F_L
    */
    double computeThValue ();
    
protected:
    
    
private:
    gsl_function F1, F2, F3, F4;
    double avaSigma0, errSigma0, avaSigma1, errSigma1, avaSigma2, errSigma2, avaSigma3, errSigma3;
};



#endif	/* MPLL_H */

    