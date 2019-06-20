/*
 * Copyright (C) 2018 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef MYMVLNU_H
#define MYMVLNU_H

class StandardModel;
#include <gsl/gsl_integration.h>

class myMVlnu {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] meson_i initial meson of the decay
     * @param[in] vector_i final vector meson of the decay
     * @param[in] lep_i final charged lepton of the decay
     */
    myMVlnu(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i, QCD::lepton lep_i);

    /**
     * @brief Destructor.
     */
    virtual ~myMVlnu();

    /**
    * @brief A method for initializing the parameters necessary for myMVlnu.
    * @return the vector of myMVlnu specific parameters
    */
    std::vector<std::string> initializemyMVlnuParameters();

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
    QCD::meson vectorM;/**< Final vector meson type */
    std::vector<std::string> mvlnuParameters;/**< The string of mandatory myMVlnu parameters */
    bool CLNflag;
    bool BGLflag;
    bool MSflag;
    bool ChiralBasisflag;
    bool gSL814gTflag;
    bool gSLm814gTflag;

    double GF;            /**<Fermi constant */
    double Mlep;          /**<Charged lepton mass */
    double MM;            /**<Initial meson mass */
    double MV;            /**<Final vector meson mass */
    double Mb;            /**<b quark mass */
    double Mc;            /**<c quark mass */
    gslpp::complex Vcb;     /**<Vckm factor */

    double MM2;           /**< Cache variable */
    double MV2;           /**< Cache variable */
    double Mlep2;           /**< Cache variable */
    double NN;            /**< coupling including the CKM element */

    double mBcstV1;
    double mBcstV2;
    double mBcstV3;
    double mBcstV4;
    double mBcstA1;
    double mBcstA2;
    double mBcstA3;
    double mBcstA4;

    double zV1;
    double zV2;
    double zV3;
    double zV4;
    double zA1;
    double zA2;
    double zA3;
    double zA4;

    double etaI;
    double chiTV;
    double chiTA;

    double af0;
    double af1;
    double af2;
    double ag0;
    double ag1;
    double ag2;
    double aF11;
    double aF12;

    double V0;
    double A00;
    double A10;
    double A20;
    double T10;
    double T20;
    double T30;

    double h_A1_1;
    double rho2;
    double R_1_1;
    double R_2_1;
    double A_0_err;
    double T_1_err;
    double T_2_err;
    double T_3_err;

    gslpp::complex gV;  /**< NP vector WC */
    gslpp::complex gA;  /**< NP axial WC */
    gslpp::complex gP;  /**< NP pseudoscalar WC */
    gslpp::complex gT;  /**< NP tensor WC */

    double avaI;/**< GSL integral variable */
    double errI;/**< GSL integral variable */
    gsl_function FI;/**< GSL integral variable */
    gsl_integration_cquad_workspace * w_I;/**< GSL integral variable */
    gsl_error_handler_t * old_handler; /**< GSL error handler store */

    unsigned int gV_updated;/**< Cache variable */
    unsigned int gA_updated;/**< Cache variable */
    unsigned int gP_updated;/**< Cache variable */
    unsigned int gT_updated;/**< Cache variable */
    unsigned int meson_updated;/**< Cache variable */
    unsigned int h_A1_updated;/**< Cache variable */
    unsigned int R_1_1_updated;/**< Cache variable */
    unsigned int R_2_1_updated;/**< Cache variable */
    unsigned int h_A1_BGL_updated;/**< Cache variable */
    unsigned int R_1_BGL_updated;/**< Cache variable */
    unsigned int R_2_BGL_updated;/**< Cache variable */
    unsigned int A_0_err_updated;/**< Cache variable */
    unsigned int T_1_err_updated;/**< Cache variable */
    unsigned int T_2_err_updated;/**< Cache variable */
    unsigned int T_3_err_updated;/**< Cache variable */
    unsigned int V_updated;/**< Cache variable */
    unsigned int A_0_updated;/**< Cache variable */
    unsigned int A_1_updated;/**< Cache variable */
    unsigned int A_2_updated;/**< Cache variable */
    unsigned int T_1_updated;/**< Cache variable */
    unsigned int T_2_updated;/**< Cache variable */
    unsigned int T_3_updated;/**< Cache variable */
    unsigned int quark_updated;/**< Cache variable */
    unsigned int N_updated;/**< Cache variable */

    gslpp::complex gV_cache;/**< Cache variable */
    gslpp::complex gA_cache;/**< Cache variable */
    gslpp::complex gP_cache;/**< Cache variable */
    gslpp::complex gT_cache;/**< Cache variable */
    gslpp::vector<double> meson_cache;/**< Cache variable */
    gslpp::vector<double> h_A1_cache;/**< Cache variable */
    gslpp::complex h_A1c_cache;/**< Cache variable */
    double R_1_1_cache;/**< Cache variable */
    double R_2_1_cache;/**< Cache variable */
    gslpp::vector<double> h_A1_BGL_cache;/**< Cache variable */
    gslpp::vector<double> R_1_BGL_cache;/**< Cache variable */
    gslpp::vector<double> R_2_BGL_cache;/**< Cache variable */
    double A_0_err_cache;/**< Cache variable */
    double T_1_err_cache;/**< Cache variable */
    double T_2_err_cache;/**< Cache variable */
    double T_3_err_cache;/**< Cache variable */
    double V_cache;/**< Cache variable */
    double A_0_cache;/**< Cache variable */
    double A_1_cache;/**< Cache variable */
    double A_2_cache;/**< Cache variable */
    double T_1_cache;/**< Cache variable */
    double T_2_cache;/**< Cache variable */
    double T_3_cache;/**< Cache variable */
    gslpp::vector<double> quark_cache;/**< Cache variable */
    gslpp::vector<double> N_cache;/**< Cache variable */
    gslpp::complex Nc_cache;/**< Cache variable */

    unsigned int H_V_p_updated;/**< Cache variable */
    unsigned int H_A_t_updated;/**< Cache variable */
    unsigned int H_A_0_updated;/**< Cache variable */
    unsigned int H_A_p_updated;/**< Cache variable */
    unsigned int H_P_updated;/**< Cache variable */
    unsigned int H_T_pt_updated;/**< Cache variable */
    unsigned int H_T_p0_updated;/**< Cache variable */
    unsigned int H_T_0_updated;/**< Cache variable */

    unsigned int H_pm_updated;/**< Cache variable */
    unsigned int H_0_updated;/**< Cache variable */
    unsigned int H_t_updated;/**< Cache variable */
    unsigned int H_T_pm_updated;/**< Cache variable */

    unsigned int H_sum_pm_updated;/**< Cache variable */
    unsigned int H_sum_0_updated;/**< Cache variable */

    unsigned int I1c_updated;/**< Cache variable */
    unsigned int I1s_updated;/**< Cache variable */
    unsigned int I2c_updated;/**< Cache variable */
    unsigned int I2s_updated;/**< Cache variable */
    unsigned int I2n_updated;/**< Cache variable */
    unsigned int I3_updated;/**< Cache variable */
    unsigned int I4_updated;/**< Cache variable */
    unsigned int I5_updated;/**< Cache variable */
    unsigned int I6c_updated;/**< Cache variable */
    unsigned int I6s_updated;/**< Cache variable */
    unsigned int I7_updated;/**< Cache variable */
    unsigned int I8_updated;/**< Cache variable */
    unsigned int I9_updated;/**< Cache variable */

    double cacheI1c;/**< Cache variable */
    double cacheI1s;/**< Cache variable */
    double cacheI2c;/**< Cache variable */
    double cacheI2s;/**< Cache variable */
    double cacheI2n;/**< Cache variable */
    double cacheI3;/**< Cache variable */
    double cacheI4;/**< Cache variable */
    double cacheI5;/**< Cache variable */
    double cacheI6c;/**< Cache variable */
    double cacheI6s;/**< Cache variable */
    double cacheI7;/**< Cache variable */
    double cacheI8;/**< Cache variable */
    double cacheI9;/**< Cache variable */


    /**
     * @brief The caching method for myMVlnu.
     */
    void checkCache();


    /**
     * @brief The update parameter method for myMVlnu.
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
    * @param[in] s_2 fit parameter
    * @return \f$ FF^{\rm fit} \f$
    */
    double FF_fit2(double q2, double f_0, double s_1, double s_2);

    /**
    * @brief The \f$ w \f$ function used in the form factor definition from arXiv:hep-ph/9712417.
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ w \f$
    */
    double w(double q2);

    /**
    * @brief The \f$ z \f$ function used in the form factor definition from arXiv:hep-ph/9712417.
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ z \f$
    */
    double z(double q2);

    /**
    * @brief The \f$ z \f$ function used in the form factor definition from arXiv:hep-ph/9712417.
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ z \f$
    */
    double zP(double M2);

    /**
    * @brief The \f$ z \f$ function used in the form factor definition from arXiv:hep-ph/9712417.
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ z \f$
    */
    double phi_f(double q2);

    /**
    * @brief The \f$ z \f$ function used in the form factor definition from arXiv:hep-ph/9712417.
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ z \f$
    */
    double f_BGL(double q2);

    /**
    * @brief The \f$ z \f$ function used in the form factor definition from arXiv:hep-ph/9712417.
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ z \f$
    */
    double phi_g(double q2);

    /**
    * @brief The \f$ z \f$ function used in the form factor definition from arXiv:hep-ph/9712417.
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ z \f$
    */
    double g_BGL(double q2);

    /**
    * @brief The \f$ z \f$ function used in the form factor definition from arXiv:hep-ph/9712417.
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ z \f$
    */
    double phi_F1(double q2);

    /**
    * @brief The \f$ z \f$ function used in the form factor definition from arXiv:hep-ph/9712417.
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ z \f$
    */
    double F1_BGL(double q2);

    /**
    * @brief The \f$ h_{A_1} \f$ function used in the form factor definition from arXiv:hep-ph/9712417.
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ h_{A_1} \f$
    */
    double h_A1(double q2);

    /**
    * @brief The \f$ R_1 \f$ function used in the form factor definition from arXiv:hep-ph/9712417.
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ R_1 \f$
    */
    double R_1(double q2);

    /**
    * @brief The \f$ R_2 \f$ function used in the form factor definition from arXiv:hep-ph/9712417.
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ R_2 \f$
    */
    double R_2(double q2);

    /**
    * @brief The \f$ R_0 \f$ function used in the form factor definition from arXiv:1703.05330.
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ R_2 \f$
    */
    double R_0(double q2);

    /**
    * @brief The \f$ R_{T_1} \f$ function used in the form factor definition from arXiv:1703.05330.
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ R_{T_1} \f$
    */
    double R_T1(double q2);

    /**
    * @brief The \f$ R_{T_2} \f$ function used in the form factor definition from arXiv:1703.05330.
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ R_{T_2} \f$
    */
    double R_T2(double q2);

    /**
    * @brief The \f$ R_{T_3} \f$ function used in the form factor definition from arXiv:1703.05330.
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ R_{T_3} \f$
    */
    double R_T3(double q2);

    /**
    * @brief The transverse form factor \f$ V \f$.
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ V \f$
    */
    double V(double q2);

    /**
    * @brief The transverse form factor \f$ A_0 \f$.
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ A_0 \f$
    */
    double A_0(double q2);

    /**
    * @brief The transverse form factor \f$ A_1 \f$.
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ A_1 \f$
    */
    double A_1(double q2);

    /**
    * @brief The transverse form factor \f$ A_2 \f$.
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ A_2 \f$
    */
    double A_2(double q2);

    /**
    * @brief The transverse form factor \f$ T_1 \f$.
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ T_1 \f$
    */
    double T_1(double q2);

    /**
    * @brief The transverse form factor \f$ T_2 \f$.
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ T_2 \f$
    */
    double T_2(double q2);

    /**
    * @brief The transverse form factor \f$ T_3 \f$.
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ T_3 \f$
    */
    double T_3(double q2);

    /**
    * @brief The factor \f$ \lambda \f$ used in the angular coefficients \f$I_i\f$.
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ \lambda \f$
    */
    double lambda(double q2);

    /**
    * @brief The helicity amplitude \f$ H_V^+ \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ H_V^+ \f$
    */
    gslpp::complex H_V_p(double q2);

    /**
    * @brief The helicity amplitude \f$ H_A^t \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ H_A^t \f$
    */
    gslpp::complex H_A_t(double q2);

    /**
    * @brief The helicity amplitude \f$ H_A^0 \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ H_A^0 \f$
    */
    gslpp::complex H_A_0(double q2);

    /**
    * @brief The helicity amplitude \f$ H_A^+ \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ H_A^+ \f$
    */
    gslpp::complex H_A_p(double q2);

    /**
    * @brief The helicity amplitude \f$ H_P \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ H_P \f$
    */
    gslpp::complex H_P(double q2);

    /**
    * @brief The helicity amplitude \f$ H_T^{+t} \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ H_T^{+t} \f$
    */
    gslpp::complex H_T_pt(double q2);

    /**
    * @brief The helicity amplitude \f$ H_T^{+0} \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ H_T^{+0} \f$
    */
    gslpp::complex H_T_p0(double q2);

    /**
    * @brief The helicity amplitude \f$ H_T^{+-} \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ H_T^{+-} \f$
    */
    gslpp::complex H_T_0(double q2);

    /**
    * @brief The helicity amplitude \f$ H_+ \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ H_+ \f$
    */
    gslpp::complex H_p(double q2);

    /**
    * @brief The helicity amplitude \f$ H_- \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ H_- \f$
    */
    gslpp::complex H_m(double q2);

    /**
    * @brief The helicity amplitude \f$ H_0 \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ H_0 \f$
    */
    gslpp::complex H_0(double q2);

    /**
    * @brief The helicity amplitude \f$ H_t \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ H_t \f$
    */
    gslpp::complex H_t(double q2);

    /**
    * @brief The helicity amplitude \f$ H_T^+ \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ H_T^+ \f$
    */
    gslpp::complex H_T_p(double q2);

    /**
    * @brief The helicity amplitude \f$ H_T^- \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ H_T^- \f$
    */
    gslpp::complex H_T_m(double q2);

    /**
    * @brief The helicity amplitude \f$ \mathcal{H}_+ \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ \mathcal{H}_+ \f$
    */
    gslpp::complex H_sum_p(double q2);

    /**
    * @brief The helicity amplitude \f$ \mathcal{H}_- \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ \mathcal{H}_- \f$
    */
    gslpp::complex H_sum_m(double q2);

    /**
    * @brief The helicity amplitude \f$ \mathcal{H}_0 \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ \mathcal{H}_0 \f$
    */
    gslpp::complex H_sum_0(double q2);

    /**
    * @brief The helicity amplitude \f$ \tilde{\mathcal{H}}_+ \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ \tilde{\mathcal{H}}_+ \f$
    */
    gslpp::complex Ht_sum_p(double q2);

    /**
    * @brief The helicity amplitude \f$ \tilde{\mathcal{H}}_- \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ \tilde{\mathcal{H}}_- \f$
    */
    gslpp::complex Ht_sum_m(double q2);

    /**
    * @brief The helicity amplitude \f$ \tilde{\mathcal{H}}_0 \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ \tilde{\mathcal{H}}_0 \f$
    */
    gslpp::complex Ht_sum_0(double q2);

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
    * @brief The angular coefficient \f$ I_{2n} \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ I_{2n} \f$
    */
    double  I_2n(double q2);

    /**
    * @brief The angular coefficient \f$ I_{1c} \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ I_{1c} \f$
    */
    double  I_1c(double q2);

    /**
    * @brief The angular coefficient \f$ I_{1s} \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ I_{1s} \f$
    */
    double  I_1s(double q2);

    /**
    * @brief The angular coefficient \f$ I_{1c} \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ I_{1c} \f$
    */
    double  I_2c(double q2);

    /**
    * @brief The angular coefficient \f$ I_{2s} \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ I_{2s} \f$
    */
    double  I_2s(double q2);

    /**
    * @brief The angular coefficient \f$ I_3 \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ I_3 \f$
    */
    double  I_3(double q2);

    /**
    * @brief The angular coefficient \f$ I_4 \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ I_4 \f$
    */
    double  I_4(double q2);

    /**
    * @brief The angular coefficient \f$ I_5 \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ I_5 \f$
    */
    double  I_5(double q2);

    /**
    * @brief The angular coefficient \f$ I_{6c} \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ I_{6c} \f$
    */
    double  I_6c(double q2);

    /**
    * @brief The angular coefficient \f$ I_{6s} \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ I_{6s} \f$
    */
    double  I_6s(double q2);

    /**
    * @brief The angular coefficient \f$ I_7 \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ I_7 \f$
    */
    double  I_7(double q2);

    /**
    * @brief The angular coefficient \f$ I_8 \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ I_8 \f$
    */
    double  I_8(double q2);

    /**
    * @brief The angular coefficient \f$ I_9 \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ I_9 \f$
    */
    double  I_9(double q2);

};

#endif /* MYMVLNU_H */
