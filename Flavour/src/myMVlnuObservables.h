/* 
 * Copyright (C) 2018 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef MYMVLNUOBSERVABLES_H
#define	MYMVLNUOBSERVABLES_H

#include "QCD.h"
#include "ThObservable.h"

/**
 * @class Gamma_myMVlnu
 * @ingroup Flavour
 * @brief A class for the binned observable @f$<\Gamma>@f$ in @f$M \to V l^-\nu@f$. 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the integrated observable @f$<\Gamma>@f$ in 
 * @f$M \to V l^-\nu@f$ in terms of the integrated helicity coefficients 
 * @f$<\I_i>@f$, computed in the myMVlnu class:
 * @f[
 * <\Gamma>=\,.
 * @f]
 */
class Gamma_myMVlnu : public ThObservable{
public:
    
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] meson_i initial meson of the decay
     * @param[in] vector_i final vector meson of the decay
     * @param[in] lep_i final charged lepton of the decay
     */
    Gamma_myMVlnu(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i, QCD::lepton lep_i);
    
    /**
     * @brief A method to compute the binned observable @f$<\Gamma>@f$ in @f$M \to V l^-\nu@f$.
     * @param[in] lep final charged lepton of the decay
     * @return @f$<\Gamma'>_{[qmin,qmax]}@f$
     */
    double computeGamma(QCD::lepton lep);
    
    /**
    * @brief The binned observable @f$<\Gamma>@f$ in @f$M \to V l^-\nu@f$.
    * @return @f$<\Gamma'>@f$
    */
    double computeThValue ();
   
private:
    QCD::lepton lep; /**< Final leptons type. */
    QCD::meson meson; /**< Initial meson type. */
    QCD::meson vectorM; /**< Final vector meson type. */

};


/**
 * @class A_lam_myMVlnu
 * @ingroup Flavour
 * @brief A class for the binned observable @f$<A_{\lambda}>@f$ in @f$M \to V l^-\nu@f$. 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the integrated observable @f$<A_{\lambda}>@f$ in 
 * @f$M \to V l^-\nu@f$ in terms of the integrated helicity coefficients 
 * @f$<\I_i>@f$, computed in the myMVlnu class:
 * @f[
 * <A_{\lambda}>=\,.
 * @f]
 */

class A_lam_myMVlnu : public Gamma_myMVlnu{
public:
    
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] meson_i initial meson of the decay
     * @param[in] vector_i final vector meson of the decay
     * @param[in] lep_i final leptons of the decay
     */
    A_lam_myMVlnu(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i, QCD::lepton lep_i);
    
    /**
    * @brief The binned observable @f$<A_{FB}>@f$ in @f$M \to V l^-\nu@f$.
    * @return @f$<A_{FB}>@f$
    */
    double computeThValue ();
   
private:
    QCD::lepton lep; /**< Final leptons type. */
    QCD::meson meson; /**< Initial meson type. */
    QCD::meson vectorM; /**< Final vector meson type. */

};


/**
 * @class 
 * @ingroup Flavour
 * @brief A class for the binned observable @f$< >@f$ in @f$M \to V l^-\nu@f$. 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the integrated observable @f$< >@f$ in 
 * @f$M \to V l^-\nu@f$ in terms of the integrated helicity coefficients 
 * @f$<\I_i>@f$, computed in the myMVlnu class:
 * @f[
 * < >=\,.
 * @f]
 */

class R_LT_myMVlnu : public ThObservable{
public:
    
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] meson_i initial meson of the decay
     * @param[in] vector_i final vector meson of the decay
     * @param[in] lep_i final leptons of the decay
     */
    R_LT_myMVlnu(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i, QCD::lepton lep_i);
    
    /**
    * @brief The binned observable @f$<A_{FB}>@f$ in @f$M \to V l^-\nu@f$.
    * @return @f$<A_{FB}>@f$
    */
    double computeThValue ();
   
private:
    QCD::lepton lep; /**< Final leptons type. */
    QCD::meson meson; /**< Initial meson type. */
    QCD::meson vectorM; /**< Final vector meson type. */

};


/**
 * @class 
 * @ingroup Flavour
 * @brief A class for the binned observable @f$< >@f$ in @f$M \to V l^-\nu@f$. 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the integrated observable @f$< >@f$ in 
 * @f$M \to V l^-\nu@f$ in terms of the integrated helicity coefficients 
 * @f$<\I_i>@f$, computed in the myMVlnu class:
 * @f[
 * < >=\,.
 * @f]
 */

class R_AB_myMVlnu : public ThObservable{
public:
    
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] meson_i initial meson of the decay
     * @param[in] vector_i final vector meson of the decay
     * @param[in] lep_i final leptons of the decay
     */
    R_AB_myMVlnu(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i, QCD::lepton lep_i);
    
    /**
    * @brief The binned observable @f$<A_{FB}>@f$ in @f$M \to V l^-\nu@f$.
    * @return @f$<A_{FB}>@f$
    */
    double computeThValue ();
   
private:
    QCD::lepton lep; /**< Final leptons type. */
    QCD::meson meson; /**< Initial meson type. */
    QCD::meson vectorM; /**< Final vector meson type. */

};


/**
 * @class A_FB_myMVlnu
 * @ingroup Flavour
 * @brief A class for the binned observable @f$<A_{FB}>@f$ in @f$M \to V l^-\nu@f$. 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the integrated observable @f$<A_{FB}>@f$ in 
 * @f$M \to V l^-\nu@f$ in terms of the integrated helicity coefficients 
 * @f$<\I_i>@f$, computed in the myMVlnu class:
 * @f[
 * <A_{FB}>=\,.
 * @f]
 */

class A_FB_myMVlnu : public Gamma_myMVlnu{
public:
    
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] meson_i initial meson of the decay
     * @param[in] vector_i final vector meson of the decay
     * @param[in] lep_i final leptons of the decay
     */
    A_FB_myMVlnu(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i, QCD::lepton lep_i);
    
    /**
    * @brief The binned observable @f$<A_{FB}>@f$ in @f$M \to V l^-\nu@f$.
    * @return @f$<A_{FB}>@f$
    */
    double computeThValue ();
   
private:
    QCD::lepton lep; /**< Final leptons type. */
    QCD::meson meson; /**< Initial meson type. */
    QCD::meson vectorM; /**< Final vector meson type. */

};


/**
 * @class A_3_myMVlnu
 * @ingroup Flavour
 * @brief A class for the binned observable @f$<A_3>@f$ in @f$M \to V l^-\nu@f$. 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the integrated observable @f$<A_3>@f$ in 
 * @f$M \to V l^-\nu@f$ in terms of the integrated helicity coefficients 
 * @f$<\I_i>@f$, computed in the myMVlnu class:
 * @f[
 * <A_3>=\,.
 * @f]
 */

class A_3_myMVlnu : public Gamma_myMVlnu{
public:
    
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] meson_i initial meson of the decay
     * @param[in] vector_i final vector meson of the decay
     * @param[in] lep_i final leptons of the decay
     */
    A_3_myMVlnu(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i, QCD::lepton lep_i);
    
    double computeA3(QCD::lepton lep);
    
    /**
    * @brief The binned observable @f$<A_{FB}>@f$ in @f$M \to V l^-\nu@f$.
    * @return @f$<A_{FB}>@f$
    */
    double computeThValue ();
   
private:
    QCD::lepton lep; /**< Final leptons type. */
    QCD::meson meson; /**< Initial meson type. */
    QCD::meson vectorM; /**< Final vector meson type. */

};


/**
 * @class A_4_myMVlnu
 * @ingroup Flavour
 * @brief A class for the binned observable @f$<A_4>@f$ in @f$M \to V l^-\nu@f$. 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the integrated observable @f$<A_4>@f$ in 
 * @f$M \to V l^-\nu@f$ in terms of the integrated helicity coefficients 
 * @f$<\I_i>@f$, computed in the myMVlnu class:
 * @f[
 * <A_4>=\,.
 * @f]
 */

class A_4_myMVlnu : public Gamma_myMVlnu{
public:
    
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] meson_i initial meson of the decay
     * @param[in] vector_i final vector meson of the decay
     * @param[in] lep_i final leptons of the decay
     */
    A_4_myMVlnu(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i, QCD::lepton lep_i);
    
    /**
    * @brief The binned observable @f$<A_{FB}>@f$ in @f$M \to V l^-\nu@f$.
    * @return @f$<A_{FB}>@f$
    */
    double computeThValue ();
   
private:
    QCD::lepton lep; /**< Final leptons type. */
    QCD::meson meson; /**< Initial meson type. */
    QCD::meson vectorM; /**< Final vector meson type. */

};


/**
 * @class A_5_myMVlnu
 * @ingroup Flavour
 * @brief A class for the binned observable @f$<A_5>@f$ in @f$M \to V l^-\nu@f$. 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the integrated observable @f$<A_5>@f$ in 
 * @f$M \to V l^-\nu@f$ in terms of the integrated helicity coefficients 
 * @f$<\I_i>@f$, computed in the myMVlnu class:
 * @f[
 * <A_5>=\,.
 * @f]
 */

class A_5_myMVlnu : public Gamma_myMVlnu{
public:
    
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] meson_i initial meson of the decay
     * @param[in] vector_i final vector meson of the decay
     * @param[in] lep_i final leptons of the decay
     */
    A_5_myMVlnu(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i, QCD::lepton lep_i);
    
    /**
    * @brief The binned observable @f$<A_{FB}>@f$ in @f$M \to V l^-\nu@f$.
    * @return @f$<A_{FB}>@f$
    */
    double computeThValue ();
   
private:
    QCD::lepton lep; /**< Final leptons type. */
    QCD::meson meson; /**< Initial meson type. */
    QCD::meson vectorM; /**< Final vector meson type. */

};


/**
 * @class A_6_myMVlnu
 * @ingroup Flavour
 * @brief A class for the binned observable @f$<A_6>@f$ in @f$M \to V l^-\nu@f$. 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the integrated observable @f$<A_6>@f$ in 
 * @f$M \to V l^-\nu@f$ in terms of the integrated helicity coefficients 
 * @f$<\I_i>@f$, computed in the myMVlnu class:
 * @f[
 * <A_6>=\,.
 * @f]
 */

class A_6_myMVlnu : public Gamma_myMVlnu{
public:
    
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] meson_i initial meson of the decay
     * @param[in] vector_i final vector meson of the decay
     * @param[in] lep_i final leptons of the decay
     */
    A_6_myMVlnu(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i, QCD::lepton lep_i);
    
    /**
    * @brief The binned observable @f$<A_{FB}>@f$ in @f$M \to V l^-\nu@f$.
    * @return @f$<A_{FB}>@f$
    */
    double computeThValue ();
   
private:
    QCD::lepton lep; /**< Final leptons type. */
    QCD::meson meson; /**< Initial meson type. */
    QCD::meson vectorM; /**< Final vector meson type. */

};


/**
 * @class A_7_myMVlnu
 * @ingroup Flavour
 * @brief A class for the binned observable @f$<A_7>@f$ in @f$M \to V l^-\nu@f$. 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the integrated observable @f$<A_7>@f$ in 
 * @f$M \to V l^-\nu@f$ in terms of the integrated helicity coefficients 
 * @f$<\I_i>@f$, computed in the myMVlnu class:
 * @f[
 * <A_7>=\,.
 * @f]
 */

class A_7_myMVlnu : public Gamma_myMVlnu{
public:
    
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] meson_i initial meson of the decay
     * @param[in] vector_i final vector meson of the decay
     * @param[in] lep_i final leptons of the decay
     */
    A_7_myMVlnu(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i, QCD::lepton lep_i);
    
    /**
    * @brief The binned observable @f$<A_{FB}>@f$ in @f$M \to V l^-\nu@f$.
    * @return @f$<A_{FB}>@f$
    */
    double computeThValue ();
   
private:
    QCD::lepton lep; /**< Final leptons type. */
    QCD::meson meson; /**< Initial meson type. */
    QCD::meson vectorM; /**< Final vector meson type. */

};


/**
 * @class A_8_myMVlnu
 * @ingroup Flavour
 * @brief A class for the binned observable @f$<A_8>@f$ in @f$M \to V l^-\nu@f$. 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the integrated observable @f$<A_8>@f$ in 
 * @f$M \to V l^-\nu@f$ in terms of the integrated helicity coefficients 
 * @f$<\I_i>@f$, computed in the myMVlnu class:
 * @f[
 * <A_8>=\,.
 * @f]
 */

class A_8_myMVlnu : public Gamma_myMVlnu{
public:
    
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] meson_i initial meson of the decay
     * @param[in] vector_i final vector meson of the decay
     * @param[in] lep_i final leptons of the decay
     */
    A_8_myMVlnu(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i, QCD::lepton lep_i);
    
    /**
    * @brief The binned observable @f$<A_{FB}>@f$ in @f$M \to V l^-\nu@f$.
    * @return @f$<A_{FB}>@f$
    */
    double computeThValue ();
   
private:
    QCD::lepton lep; /**< Final leptons type. */
    QCD::meson meson; /**< Initial meson type. */
    QCD::meson vectorM; /**< Final vector meson type. */

};


/**
 * @class A_9_myMVlnu
 * @ingroup Flavour
 * @brief A class for the binned observable @f$<A_9>@f$ in @f$M \to V l^-\nu@f$. 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the integrated observable @f$<A_9>@f$ in 
 * @f$M \to V l^-\nu@f$ in terms of the integrated helicity coefficients 
 * @f$<\I_i>@f$, computed in the myMVlnu class:
 * @f[
 * <A_9>=\,.
 * @f]
 */

class A_9_myMVlnu : public Gamma_myMVlnu{
public:
    
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] meson_i initial meson of the decay
     * @param[in] vector_i final vector meson of the decay
     * @param[in] lep_i final leptons of the decay
     */
    A_9_myMVlnu(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i, QCD::lepton lep_i);
    
    /**
    * @brief The binned observable @f$<A_{FB}>@f$ in @f$M \to V l^-\nu@f$.
    * @return @f$<A_{FB}>@f$
    */
    double computeThValue ();
   
private:
    QCD::lepton lep; /**< Final leptons type. */
    QCD::meson meson; /**< Initial meson type. */
    QCD::meson vectorM; /**< Final vector meson type. */

};



/*********************
 *                   *
 *    LFUV RATIOS    *
 *                   *
 *********************/



/**
 * @class R_D*
 * @ingroup Flavour
 * @brief A class for the observable @f$(D*)@f$ in @f$M \to V l^-\nu@f$. 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the observable @f$R(D*)@f$ in 
 * @f$M \to V l^-\nu@f$ in terms of the integrated helicity coefficients 
 * @f$<\I_i>@f$, computed in the myMVlnu class:
 * @f[
 * R(D*)=\,.
 * @f]
 */

class R_Dst : public Gamma_myMVlnu{
public:
    
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] meson_i initial meson of the decay
     * @param[in] vector_i final vector meson of the decay
     * @param[in] lep_1 final leptons of the decay
     * @param[in] lep_2 final leptons of the decay
     * @param[in] lep_3 final leptons of the decay
     */
    R_Dst(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i, 
            QCD::lepton lep_1, QCD::lepton lep_2, QCD::lepton lep_3);
    
    /**
    * @brief The binned observable @f$<A_{FB}>@f$ in @f$M \to V l^-\nu@f$.
    * @return @f$<A_{FB}>@f$
    */
    double computeThValue ();
   
private:
    QCD::lepton lep1; /**< Final leptons type. */
    QCD::lepton lep2; /**< Final leptons type. */
    QCD::lepton lep3; /**< Final leptons type. */
    QCD::meson meson; /**< Initial meson type. */
    QCD::meson vectorM; /**< Final vector meson type. */

};


/**
 * @class A_lam_myMVlnu
 * @ingroup Flavour
 * @brief A class for the binned observable @f$<A_{\lambda}>@f$ in @f$M \to V l^-\nu@f$. 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the integrated observable @f$<A_{\lambda}>@f$ in 
 * @f$M \to V l^-\nu@f$ in terms of the integrated helicity coefficients 
 * @f$<\I_i>@f$, computed in the myMVlnu class:
 * @f[
 * <A_{\lambda}>=\,.
 * @f]
 */

class R_A_lam_myMVlnu : public ThObservable{
public:
    
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] meson_i initial meson of the decay
     * @param[in] vector_i final vector meson of the decay
     * @param[in] lep_i final leptons of the decay
     */
    R_A_lam_myMVlnu(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i, 
            QCD::lepton lep_1, QCD::lepton lep_2, QCD::lepton lep_3);
    
    /**
    * @brief The binned observable @f$<A_{FB}>@f$ in @f$M \to V l^-\nu@f$.
    * @return @f$<A_{FB}>@f$
    */
    double computeThValue ();
   
private:
    QCD::lepton lep1; /**< Final leptons type. */
    QCD::lepton lep2; /**< Final leptons type. */
    QCD::lepton lep3; /**< Final leptons type. */
    QCD::meson meson; /**< Initial meson type. */
    QCD::meson vectorM; /**< Final vector meson type. */

};


/**
 * @class 
 * @ingroup Flavour
 * @brief A class for the binned observable @f$< >@f$ in @f$M \to V l^-\nu@f$. 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the integrated observable @f$< >@f$ in 
 * @f$M \to V l^-\nu@f$ in terms of the integrated helicity coefficients 
 * @f$<\I_i>@f$, computed in the myMVlnu class:
 * @f[
 * < >=\,.
 * @f]
 */

class R_R_LT_myMVlnu : public R_LT_myMVlnu{
public:
    
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] meson_i initial meson of the decay
     * @param[in] vector_i final vector meson of the decay
     * @param[in] lep_i final leptons of the decay
     */
    R_R_LT_myMVlnu(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i, 
            QCD::lepton lep_1, QCD::lepton lep_2, QCD::lepton lep_3);
    
    /**
    * @brief The binned observable @f$<A_{FB}>@f$ in @f$M \to V l^-\nu@f$.
    * @return @f$<A_{FB}>@f$
    */
    double computeThValue ();
   
private:
    QCD::lepton lep1; /**< Final leptons type. */
    QCD::lepton lep2; /**< Final leptons type. */
    QCD::lepton lep3; /**< Final leptons type. */
    QCD::meson meson; /**< Initial meson type. */
    QCD::meson vectorM; /**< Final vector meson type. */

};


/**
 * @class 
 * @ingroup Flavour
 * @brief A class for the binned observable @f$< >@f$ in @f$M \to V l^-\nu@f$. 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the integrated observable @f$< >@f$ in 
 * @f$M \to V l^-\nu@f$ in terms of the integrated helicity coefficients 
 * @f$<\I_i>@f$, computed in the myMVlnu class:
 * @f[
 * < >=\,.
 * @f]
 */

class R_R_AB_myMVlnu : public R_AB_myMVlnu{
public:
    
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] meson_i initial meson of the decay
     * @param[in] vector_i final vector meson of the decay
     * @param[in] lep_i final leptons of the decay
     */
    R_R_AB_myMVlnu(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i, 
            QCD::lepton lep_1, QCD::lepton lep_2, QCD::lepton lep_3);
    
    /**
    * @brief The binned observable @f$<A_{FB}>@f$ in @f$M \to V l^-\nu@f$.
    * @return @f$<A_{FB}>@f$
    */
    double computeThValue ();
   
private:
    QCD::lepton lep1; /**< Final leptons type. */
    QCD::lepton lep2; /**< Final leptons type. */
    QCD::lepton lep3; /**< Final leptons type. */
    QCD::meson meson; /**< Initial meson type. */
    QCD::meson vectorM; /**< Final vector meson type. */

};


/**
 * @class A_FB_myMVlnu
 * @ingroup Flavour
 * @brief A class for the binned observable @f$<A_{FB}>@f$ in @f$M \to V l^-\nu@f$. 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the integrated observable @f$<A_{FB}>@f$ in 
 * @f$M \to V l^-\nu@f$ in terms of the integrated helicity coefficients 
 * @f$<\I_i>@f$, computed in the myMVlnu class:
 * @f[
 * <A_{FB}>=\,.
 * @f]
 */

class R_A_FB_myMVlnu : public ThObservable{
public:
    
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] meson_i initial meson of the decay
     * @param[in] vector_i final vector meson of the decay
     * @param[in] lep_i final leptons of the decay
     */
    R_A_FB_myMVlnu(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i, 
            QCD::lepton lep_1, QCD::lepton lep_2, QCD::lepton lep_3);
    
    /**
    * @brief The binned observable @f$<A_{FB}>@f$ in @f$M \to V l^-\nu@f$.
    * @return @f$<A_{FB}>@f$
    */
    double computeThValue ();
   
private:
    QCD::lepton lep1; /**< Final leptons type. */
    QCD::lepton lep2; /**< Final leptons type. */
    QCD::lepton lep3; /**< Final leptons type. */
    QCD::meson meson; /**< Initial meson type. */
    QCD::meson vectorM; /**< Final vector meson type. */

};


/**
 * @class A_3_myMVlnu
 * @ingroup Flavour
 * @brief A class for the binned observable @f$<A_3>@f$ in @f$M \to V l^-\nu@f$. 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the integrated observable @f$<A_3>@f$ in 
 * @f$M \to V l^-\nu@f$ in terms of the integrated helicity coefficients 
 * @f$<\I_i>@f$, computed in the myMVlnu class:
 * @f[
 * <A_3>=\,.
 * @f]
 */

class R_A_3_myMVlnu : public ThObservable{
public:
    
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] meson_i initial meson of the decay
     * @param[in] vector_i final vector meson of the decay
     * @param[in] lep_i final leptons of the decay
     */
    R_A_3_myMVlnu(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i, 
            QCD::lepton lep_1, QCD::lepton lep_2, QCD::lepton lep_3);
    
    /**
    * @brief The binned observable @f$<A_{FB}>@f$ in @f$M \to V l^-\nu@f$.
    * @return @f$<A_{FB}>@f$
    */
    double computeThValue ();
   
private:
    QCD::lepton lep1; /**< Final leptons type. */
    QCD::lepton lep2; /**< Final leptons type. */
    QCD::lepton lep3; /**< Final leptons type. */
    QCD::meson meson; /**< Initial meson type. */
    QCD::meson vectorM; /**< Final vector meson type. */

};


/**
 * @class A_4_myMVlnu
 * @ingroup Flavour
 * @brief A class for the binned observable @f$<A_4>@f$ in @f$M \to V l^-\nu@f$. 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the integrated observable @f$<A_4>@f$ in 
 * @f$M \to V l^-\nu@f$ in terms of the integrated helicity coefficients 
 * @f$<\I_i>@f$, computed in the myMVlnu class:
 * @f[
 * <A_4>=\,.
 * @f]
 */

class R_A_4_myMVlnu : public ThObservable{
public:
    
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] meson_i initial meson of the decay
     * @param[in] vector_i final vector meson of the decay
     * @param[in] lep_i final leptons of the decay
     */
    R_A_4_myMVlnu(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i, 
            QCD::lepton lep_1, QCD::lepton lep_2, QCD::lepton lep_3);
    
    /**
    * @brief The binned observable @f$<A_{FB}>@f$ in @f$M \to V l^-\nu@f$.
    * @return @f$<A_{FB}>@f$
    */
    double computeThValue ();
   
private:
    QCD::lepton lep1; /**< Final leptons type. */
    QCD::lepton lep2; /**< Final leptons type. */
    QCD::lepton lep3; /**< Final leptons type. */
    QCD::meson meson; /**< Initial meson type. */
    QCD::meson vectorM; /**< Final vector meson type. */

};


/**
 * @class A_5_myMVlnu
 * @ingroup Flavour
 * @brief A class for the binned observable @f$<A_5>@f$ in @f$M \to V l^-\nu@f$. 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the integrated observable @f$<A_5>@f$ in 
 * @f$M \to V l^-\nu@f$ in terms of the integrated helicity coefficients 
 * @f$<\I_i>@f$, computed in the myMVlnu class:
 * @f[
 * <A_5>=\,.
 * @f]
 */

class R_A_5_myMVlnu : public ThObservable{
public:
    
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] meson_i initial meson of the decay
     * @param[in] vector_i final vector meson of the decay
     * @param[in] lep_i final leptons of the decay
     */
    R_A_5_myMVlnu(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i, 
            QCD::lepton lep_1, QCD::lepton lep_2, QCD::lepton lep_3);
    
    /**
    * @brief The binned observable @f$<A_{FB}>@f$ in @f$M \to V l^-\nu@f$.
    * @return @f$<A_{FB}>@f$
    */
    double computeThValue ();
   
private:
    QCD::lepton lep1; /**< Final leptons type. */
    QCD::lepton lep2; /**< Final leptons type. */
    QCD::lepton lep3; /**< Final leptons type. */
    QCD::meson meson; /**< Initial meson type. */
    QCD::meson vectorM; /**< Final vector meson type. */

};


/**
 * @class A_6_myMVlnu
 * @ingroup Flavour
 * @brief A class for the binned observable @f$<A_6>@f$ in @f$M \to V l^-\nu@f$. 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the integrated observable @f$<A_6>@f$ in 
 * @f$M \to V l^-\nu@f$ in terms of the integrated helicity coefficients 
 * @f$<\I_i>@f$, computed in the myMVlnu class:
 * @f[
 * <A_6>=\,.
 * @f]
 */

class R_A_6_myMVlnu : public ThObservable{
public:
    
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] meson_i initial meson of the decay
     * @param[in] vector_i final vector meson of the decay
     * @param[in] lep_i final leptons of the decay
     */
    R_A_6_myMVlnu(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i, 
            QCD::lepton lep_1, QCD::lepton lep_2, QCD::lepton lep_3);
    
    /**
    * @brief The binned observable @f$<A_{FB}>@f$ in @f$M \to V l^-\nu@f$.
    * @return @f$<A_{FB}>@f$
    */
    double computeThValue ();
   
private:
    QCD::lepton lep1; /**< Final leptons type. */
    QCD::lepton lep2; /**< Final leptons type. */
    QCD::lepton lep3; /**< Final leptons type. */
    QCD::meson meson; /**< Initial meson type. */
    QCD::meson vectorM; /**< Final vector meson type. */

};


/**
 * @class A_7_myMVlnu
 * @ingroup Flavour
 * @brief A class for the binned observable @f$<A_7>@f$ in @f$M \to V l^-\nu@f$. 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the integrated observable @f$<A_7>@f$ in 
 * @f$M \to V l^-\nu@f$ in terms of the integrated helicity coefficients 
 * @f$<\I_i>@f$, computed in the myMVlnu class:
 * @f[
 * <A_7>=\,.
 * @f]
 */

class D_A_7_myMVlnu : public A_7_myMVlnu{
public:
    
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] meson_i initial meson of the decay
     * @param[in] vector_i final vector meson of the decay
     * @param[in] lep_i final leptons of the decay
     */
    D_A_7_myMVlnu(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i, 
            QCD::lepton lep_1, QCD::lepton lep_2, QCD::lepton lep_3);
    
    /**
    * @brief The binned observable @f$<A_{FB}>@f$ in @f$M \to V l^-\nu@f$.
    * @return @f$<A_{FB}>@f$
    */
    double computeThValue ();
   
private:
    QCD::lepton lep1; /**< Final leptons type. */
    QCD::lepton lep2; /**< Final leptons type. */
    QCD::lepton lep3; /**< Final leptons type. */
    QCD::meson meson; /**< Initial meson type. */
    QCD::meson vectorM; /**< Final vector meson type. */

};


/**
 * @class A_8_myMVlnu
 * @ingroup Flavour
 * @brief A class for the binned observable @f$<A_8>@f$ in @f$M \to V l^-\nu@f$. 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the integrated observable @f$<A_8>@f$ in 
 * @f$M \to V l^-\nu@f$ in terms of the integrated helicity coefficients 
 * @f$<\I_i>@f$, computed in the myMVlnu class:
 * @f[
 * <A_8>=\,.
 * @f]
 */

class D_A_8_myMVlnu : public A_8_myMVlnu{
public:
    
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] meson_i initial meson of the decay
     * @param[in] vector_i final vector meson of the decay
     * @param[in] lep_i final leptons of the decay
     */
    D_A_8_myMVlnu(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i, 
            QCD::lepton lep_1, QCD::lepton lep_2, QCD::lepton lep_3);
    
    /**
    * @brief The binned ovector_i final vector meson of the decaybservable @f$<A_{FB}>@f$ in @f$M \to V l^-\nu@f$.
    * @return @f$<A_{FB}>@f$
    */
    double computeThValue ();
   
private:
    QCD::lepton lep1; /**< Final leptons type. */
    QCD::lepton lep2; /**< Final leptons type. */
    QCD::lepton lep3; /**< Final leptons type. */
    QCD::meson meson; /**< Initial meson type. */
    QCD::meson vectorM; /**< Final vector meson type. */

};


/**
 * @class A_9_myMVlnu
 * @ingroup Flavour
 * @brief A class for the binned observable @f$<A_9>@f$ in @f$M \to V l^-\nu@f$. 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the integrated observable @f$<A_9>@f$ in 
 * @f$M \to V l^-\nu@f$ in terms of the integrated helicity coefficients 
 * @f$<\I_i>@f$, computed in the myMVlnu class:
 * @f[
 * <A_9>=\,.
 * @f]
 */

class D_A_9_myMVlnu : public A_9_myMVlnu{
public:
    
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] meson_i initial meson of the decay
     * @param[in] vector_i final vector meson of the decay
     * @param[in] lep_i final leptons of the decay
     */
    D_A_9_myMVlnu(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i, 
            QCD::lepton lep_1, QCD::lepton lep_2, QCD::lepton lep_3);
    
    /**
    * @brief The binned observable @f$<A_{FB}>@f$ in @f$M \to V l^-\nu@f$.
    * @return @f$<A_{FB}>@f$
    */
    double computeThValue ();
   
private:
    QCD::lepton lep1; /**< Final leptons type. */
    QCD::lepton lep2; /**< Final leptons type. */
    QCD::lepton lep3; /**< Final leptons type. */
     QCD::meson meson; /**< Initial meson type. */
    QCD::meson vectorM; /**< Final vector meson type. */

};


#endif /* MYMVLNUOBSERVABLES_H */

