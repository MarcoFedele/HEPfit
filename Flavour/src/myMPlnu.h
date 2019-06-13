/*
 * Copyright (C) 2018 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef MYMPLNU_H
#define MYMPLNU_H

class StandardModel;
#include <gsl/gsl_integration.h>

class myMPlnu {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] meson_i initial meson of the decay
     * @param[in] pseudoscalar_i final vector meson of the decay
     * @param[in] lep_i final charged lepton of the decay
     */
    myMPlnu(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson pseudoscalar_i, QCD::lepton lep_i);

    /**
     * @brief Destructor.
     */
    virtual ~myMPlnu();

    /**
    * @brief A method for initializing the parameters necessary for MVlnu.
    * @return the vector of MVlnu specific parameters
    */
    std::vector<std::string> initializemyMPlnuParameters();

    /**
    * @brief The integral of \f$ I_{i} \f$ from \f$m_lep^2 \f$ to \f$(m_M - M_V)^2\f$
    * @param[in] i index of the angular coefficient \f$ I_{i} \f$
    * @return \f$ <\I_{i}> \f$
    */
    double integrateI(int i);

private:
    const StandardModel& mySM;/**< Model type */
    QCD::lepton lep;/**< Final charged lepton type */
    QCD::meson meson;/**< Initial meson type */
    QCD::meson pseudoscalar;/**< Final vector meson type */
    std::vector<std::string> mplnuParameters;/**< The string of mandatory myMPlnu parameters */
    bool MSflag;
    bool ChiralBasisflag;

    double GF;            /**<Fermi constant */
    double Mlep;          /**<Charged lepton mass */
    double MM;            /**<Initial meson mass */
    double MP;            /**<Final pseudoscalar meson mass */
    double Mb;            /**<b quark mass */
    double Mc;            /**<c quark mass */
    gslpp::complex Vcb;     /**<Vckm factor */

    double MM2;           /**< Cache variable */
    double MP2;           /**< Cache variable */
    double Mlep2;           /**< Cache variable */
    double NN;            /**< coupling including the CKM element */

    double fp0;
    double f00;
    double fT0;

    double ap0;
    double ap1;
    double ap2;
    double ap3;
    double a01;
    double a02;
    double a03;
    double RT_err;

    gslpp::complex gV;  /**<NP vector WC */
    gslpp::complex gS;  /**<NP scalar WC */
    gslpp::complex gT;  /**<NP tensor WC */

    double avaI;/**< GSL integral variable */
    double errI;/**< GSL integral variable */
    gsl_function FI;/**< GSL integral variable */
    gsl_integration_cquad_workspace * w_I;/**< GSL integral variable */
    gsl_error_handler_t * old_handler; /**< GSL error handler store */

    unsigned int gV_updated;/**< Cache variable */
    unsigned int gS_updated;/**< Cache variable */
    unsigned int gT_updated;/**< Cache variable */
    unsigned int meson_updated;/**< Cache variable */
    unsigned int R_T_updated;/**< Cache variable */
    unsigned int f_0_updated;/**< Cache variable */
    unsigned int f_p_updated;/**< Cache variable */
    unsigned int f_T_updated;/**< Cache variable */
    unsigned int quark_updated;/**< Cache variable */
    unsigned int N_updated;/**< Cache variable */

    gslpp::complex gV_cache;/**< Cache variable */
    gslpp::complex gS_cache;/**< Cache variable */
    gslpp::complex gT_cache;/**< Cache variable */
    gslpp::vector<double> meson_cache;/**< Cache variable */
    gslpp::vector<double> f_p_Vcache;/**< Cache variable */
    gslpp::vector<double> f_0_Vcache;/**< Cache variable */
    double R_T_cache;/**< Cache variable */
    double f_p_cache;/**< Cache variable */
    double f_0_cache;/**< Cache variable */
    double f_T_cache;/**< Cache variable */
    gslpp::vector<double> quark_cache;/**< Cache variable */
    gslpp::vector<double> N_cache;/**< Cache variable */
    gslpp::complex Nc_cache;/**< Cache variable */

    unsigned int H_V_t_updated;/**< Cache variable */
    unsigned int H_V_updated;/**< Cache variable */
    unsigned int H_S_updated;/**< Cache variable */
    unsigned int H_T_updated;/**< Cache variable */
    unsigned int H_t_updated;/**< Cache variable */

    unsigned int a_updated;/**< Cache variable */
    unsigned int b_updated;/**< Cache variable */
    unsigned int c_updated;/**< Cache variable */
    unsigned int Gamma_p_updated;/**< Cache variable */
    unsigned int Gamma_m_updated;/**< Cache variable */

    double cachea;/**< Cache variable */
    double cacheb;/**< Cache variable */
    double cachec;/**< Cache variable */
    double cachegp;/**< Cache variable */
    double cachegm;/**< Cache variable */


    /**
     * @brief The caching method for MVlnu.
     */
    void checkCache();


    /**
     * @brief The update parameter method for MVlnu.
     */
    void updateParameters();

    /**
    * @brief The fit function (9) from hep-ph/0001113, \f$ FF^{\rm fit} \f$.
    * @param[in] q2 \f$q^2\f$ of the decay
    * @param[in] f_0 fit parameter
    * @param[in] s_1 fit parameter
    * @return \f$ FF^{\rm fit} \f$
    */
    double FF_fit1(double q2, double f_0, double s_1);

    /**
    * @brief The fit function (10) from hep-ph/0001113, \f$ FF^{\rm fit} \f$.
    * @param[in] q2 \f$q^2\f$ of the decay
    * @param[in] f_0 fit parameter
    * @param[in] s_1 fit parameter
    * @return \f$ FF^{\rm fit} \f$
    */
    double FF_fit2(double q2, double f_0, double s_1);

    /**
    * @brief The \f$ w \f$ function used in the form factor definition from arXiv:1503.07237.
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ w \f$
    */
    double w(double q2);

    /**
    * @brief The \f$ z \f$ function used in the form factor definition from arXiv:1503.07237.
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ z \f$
    */
    double z(double q2);

    /**
    * @brief The \f$ \phi_+ \f$ function used in the form factor definition from arXiv:1503.07237.
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ \phi_+ \f$
    */
    double phip(double q2);

    /**
    * @brief The \f$ \phi_0 \f$ function used in the form factor definition from arXiv:1503.07237.
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ \phi_0 \f$
    */
    double phi0(double q2);

    /**
    * @brief The \f$ R_{T} \f$ function used in the form factor definition from arXiv:1703.05330.
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ R_{T} \f$
    */
    double R_T(double q2);

    /**
    * @brief The transverse form factor \f$ f_+ \f$.
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ f_+ \f$
    */
    double f_p(double q2);

    /**
    * @brief The transverse form factor \f$ f_0 \f$.
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ f_0 \f$
    */
    double f_0(double q2);

    /**
    * @brief The transverse form factor \f$ f_T \f$.
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ f_T \f$
    */
    double f_T(double q2);

    /**
    * @brief The factor \f$ \lambda \f$ used in the angular coefficients \f$I_i\f$.
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ \lambda \f$
    */
    double lambda(double q2);

    /**
    * @brief The helicity amplitude \f$ H_V \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ H_V \f$
    */
    gslpp::complex H_V(double q2);

    /**
    * @brief The helicity amplitude \f$ H_V^t \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ H_V^t \f$
    */
    gslpp::complex H_V_t(double q2);

    /**
    * @brief The helicity amplitude \f$ H_S \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ H_S \f$
    */
    gslpp::complex H_S(double q2);

    /**
    * @brief The helicity amplitude \f$ H_t \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ H_t \f$
    */
    gslpp::complex H_t(double q2);

    /**
    * @brief The helicity amplitude \f$ H_T \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ H_T \f$
    */
    gslpp::complex H_T(double q2);

    /**
    * @brief The factor \f$ \beta^2 \f$ used in the angular coefficients \f$I_i\f$.
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ \beta^2 \f$
    */
    double beta2 (double q2);

    /**
    * @brief The normalization of the angular coefficients \f$ \Gamma_0 \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ \Gamma_0 \f$
    */
    double Gamma_0(double q2);

    /**
    * @brief The angular coefficient \f$ a_{\theta} \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ a_{\theta} \f$
    */
    double  a_theta(double q2);

    /**
    * @brief The angular coefficient \f$ b_{\theta} \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ b_{\theta} \f$
    */
    double  b_theta(double q2);

    /**
    * @brief The angular coefficient \f$ c_{\theta} \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ c_{\theta} \f$
    */
    double  c_theta(double q2);

    /**
    * @brief The angular coefficient \f$ \Gamma_+ \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \Gamma_+ \f$
    */
    double  Gamma_p(double q2);

    /**
    * @brief The angular coefficient \f$ \Gamma_- \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \Gamma_- \f$
    */
    double  Gamma_m(double q2);

};

#endif /* MYMPLNU_H */
