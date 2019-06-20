/*
 * Copyright (C) 2014 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef MYMPLNUOBSERVABLES_H
#define	MYMPLNUOBSERVABLES_H

#include "QCD.h"
#include "ThObservable.h"

/**
 * @class Gamma_myMPlnu
 * @ingroup Flavour
 * @brief A class for the binned observable @f$<\Gamma>@f$ in @f$M \to P l^-\nu@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the integrated observable @f$<\Gamma>@f$ in
 * @f$M \to P l^-\nu@f$ in terms of the integrated helicity coefficients
 * @f$<\I_i>@f$, computed in the MVlnu class:
 * @f[
 * <\Gamma>=\,.
 * @f]
 */
class Gamma_myMPlnu : public ThObservable{
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] meson_i initial meson of the decay
     * @param[in] pseudoscalar_i final vector meson of the decay
     * @param[in] lep_i final charged lepton of the decay
     */
    Gamma_myMPlnu(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson pseudoscalar_i, QCD::lepton lep_i);

    /**
     * @brief A method to compute the binned observable @f$<\Gamma>@f$ in @f$M \to P l^-\nu@f$.
     * @param[in] lep final charged lepton of the decay
     * @return @f$<\Gamma'>_{[qmin,qmax]}@f$
     */
    double computeGamma(QCD::lepton lep);

    /**
    * @brief The binned observable @f$<\Gamma>@f$ in @f$M \to P l^-\nu@f$.
    * @return @f$<\Gamma'>@f$
    */
    double computeThValue ();

private:
    QCD::lepton lep; /**< Final leptons type. */
    QCD::meson meson; /**< Initial meson type. */
    QCD::meson pseudoscalar; /**< Final vector meson type. */

};


/**
 * @class A_FB_myMPlnu
 * @ingroup Flavour
 * @brief A class for the binned observable @f$<A_{FB}>@f$ in @f$M \to P l^-\nu@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the integrated observable @f$<A_{FB}>@f$ in
 * @f$M \to P l^-\nu@f$ in terms of the integrated helicity coefficients
 * @f$<\I_i>@f$, computed in the MVlnu class:
 * @f[
 * <A_{FB}>=\,.
 * @f]
 */

class A_FB_myMPlnu : public Gamma_myMPlnu{
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] meson_i initial meson of the decay
     * @param[in] pseudoscalar_i final vector meson of the decay
     * @param[in] lep_i final leptons of the decay
     */
    A_FB_myMPlnu(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson pseudoscalar_i, QCD::lepton lep_i);

    /**
    * @brief The binned observable @f$<A_{FB}>@f$ in @f$M \to P l^-\nu@f$.
    * @return @f$<A_{FB}>@f$
    */
    double computeThValue ();

private:
    QCD::lepton lep; /**< Final leptons type. */
    QCD::meson meson; /**< Initial meson type. */
    QCD::meson pseudoscalar; /**< Final vector meson type. */

};


/**
 * @class A_lam_MVlnu
 * @ingroup Flavour
 * @brief A class for the binned observable @f$<A_{\lambda}>@f$ in @f$M \to P l^-\nu@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the integrated observable @f$<A_{\lambda}>@f$ in
 * @f$M \to P l^-\nu@f$ in terms of the integrated helicity coefficients
 * @f$<\I_i>@f$, computed in the MVlnu class:
 * @f[
 * <A_{\lambda}>=\,.
 * @f]
 */

class A_lam_myMPlnu : public Gamma_myMPlnu{
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] meson_i initial meson of the decay
     * @param[in] pseudoscalar_i final vector meson of the decay
     * @param[in] lep_i final leptons of the decay
     */
    A_lam_myMPlnu(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson pseudoscalar_i, QCD::lepton lep_i);

    /**
    * @brief The binned observable @f$<A_{FB}>@f$ in @f$M \to P l^-\nu@f$.
    * @return @f$<A_{FB}>@f$
    */
    double computeThValue ();

private:
    QCD::lepton lep; /**< Final leptons type. */
    QCD::meson meson; /**< Initial meson type. */
    QCD::meson pseudoscalar; /**< Final vector meson type. */

};



/*********************
 *                   *
 *    LFUV RATIOS    *
 *                   *
 *********************/



/**
 * @class R_D
 * @ingroup Flavour
 * @brief A class for the observable @f$(D*)@f$ in @f$M \to P l^+\nu@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the observable @f$R(D)@f$ in
 * @f$M \to P l^-\nu@f$ in terms of the integrated helicity coefficients
 * @f$<\I_i>@f$, computed in the MVlnu class:
 * @f[
 * R(D)=\,.
 * @f]
 */

class R_D : public Gamma_myMPlnu{
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] meson_i initial meson of the decay
     * @param[in] pseudoscalar_i final vector meson of the decay
     * @param[in] lep_1 final leptons of the decay
     * @param[in] lep_2 final leptons of the decay
     * @param[in] lep_3 final leptons of the decay
     */
    R_D(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson pseudoscalar_i,
            QCD::lepton lep_1, QCD::lepton lep_2, QCD::lepton lep_3);

    /**
    * @brief The binned observable @f$<A_{FB}>@f$ in @f$M \to P l^-\nu@f$.
    * @return @f$<A_{FB}>@f$
    */
    double computeThValue ();

private:
    QCD::lepton lep1; /**< Final leptons type. */
    QCD::lepton lep2; /**< Final leptons type. */
    QCD::lepton lep3; /**< Final leptons type. */
    QCD::meson meson; /**< Initial meson type. */
    QCD::meson pseudoscalar; /**< Final vector meson type. */

};


/**
 * @class A_lam_myMPlnu
 * @ingroup Flavour
 * @brief A class for the binned observable @f$<A_{\lambda}>@f$ in @f$M \to P l^-\nu@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the integrated observable @f$<A_{\lambda}>@f$ in
 * @f$M \to P l^-\nu@f$ in terms of the integrated helicity coefficients
 * @f$<\I_i>@f$, computed in the MVlnu class:
 * @f[
 * <A_{\lambda}>=\,.
 * @f]
 */

class R_A_lam_myMPlnu : public Gamma_myMPlnu{
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] meson_i initial meson of the decay
     * @param[in] pseudoscalar_i final vector meson of the decay
     * @param[in] lep_i final leptons of the decay
     */
    R_A_lam_myMPlnu(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson pseudoscalar_i,
            QCD::lepton lep_1, QCD::lepton lep_2, QCD::lepton lep_3);

    /**
    * @brief The binned observable @f$<A_{FB}>@f$ in @f$M \to P l^-\nu@f$.
    * @return @f$<A_{FB}>@f$
    */
    double computeThValue ();

private:
    QCD::lepton lep1; /**< Final leptons type. */
    QCD::lepton lep2; /**< Final leptons type. */
    QCD::lepton lep3; /**< Final leptons type. */
    QCD::meson meson; /**< Initial meson type. */
    QCD::meson pseudoscalar; /**< Final vector meson type. */

};


/**
 * @class A_FB_myMPlnu
 * @ingroup Flavour
 * @brief A class for the binned observable @f$<A_{FB}>@f$ in @f$M \to P l^+\nu@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the integrated observable @f$<A_{FB}>@f$ in
 * @f$M \to P l^-\nu@f$ in terms of the integrated helicity coefficients
 * @f$<\I_i>@f$, computed in the MVlnu class:
 * @f[
 * <A_{FB}>=\,.
 * @f]
 */

class R_A_FB_myMPlnu : public Gamma_myMPlnu{
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] meson_i initial meson of the decay
     * @param[in] pseudoscalar_i final vector meson of the decay
     * @param[in] lep_i final leptons of the decay
     */
    R_A_FB_myMPlnu(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson pseudoscalar_i,
            QCD::lepton lep_1, QCD::lepton lep_2, QCD::lepton lep_3);

    /**
    * @brief The binned observable @f$<A_{FB}>@f$ in @f$M \to P l^-\nu@f$.
    * @return @f$<A_{FB}>@f$
    */
    double computeThValue ();

private:
    QCD::lepton lep1; /**< Final leptons type. */
    QCD::lepton lep2; /**< Final leptons type. */
    QCD::lepton lep3; /**< Final leptons type. */
    QCD::meson meson; /**< Initial meson type. */
    QCD::meson pseudoscalar; /**< Final vector meson type. */

};


#endif /* MYMPLNUOBSERVABLES_H */
