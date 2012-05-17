/* 
 * File:   EWModel.h
 * Author: mishima
 */

#ifndef EWMODEL_H
#define	EWMODEL_H

#include <gslpp.h>
#include <StandardModel.h>
using namespace gslpp;


class EWModel {
public:

    // The virtual functions below have to be defined in child classes 
    // inherited from the current class for each (new-physics) model. 
    
    /**
     * @return the W boson mass
     */
    virtual double Mw() const {
        throw "EWModel::Mw() is undefined.";
    };

    /**
     * @return Mw^2/Mz^2
     */
    virtual double cW2() const {
        throw "EWModel::cW2() is undefined.";
    };
    
    /**
     * @return 1-Mw^2/Mz^2
     */
    virtual double sW2() const {
        throw "EWModel::sW2() is undefined.";
    };
    
    /**
     * @brief effective coupling rho_Z^l
     * @param[in] l name of a lepton 
     * @return rho_Z^l
     */
    virtual complex rhoZ_l(const StandardModel::lepton l) const {
        throw "EWModel::rhoZ_l() is undefined.";
    };

    /**
     * @brief 
     * @param[in] q name of a quark
     * @return rho_Z^q
     */
    virtual complex rhoZ_q(const StandardModel::quark q) const {
        throw "EWModel::rhoZ_q() is undefined.";
    };
    
    /**
     * @brief the ratio of the effective couplings for neutral-current interactions
     * @param[in] l name of a lepton 
     * @return g_V^l/g_A^l
     */
    virtual complex gZl_over_gAl(const StandardModel::lepton l) const {
        throw "EWModel::gZl_over_gAl() is undefined.";
    };
    
    /**
     * @brief the ratio of the effective couplings for neutral-current interactions
     * @param[in] q name of a quark
     * @return g_V^q/g_A^q
     */
    virtual complex gZq_over_gAq(const StandardModel::quark q) const {
        throw "EWModel::gZq_over_gAq() is undefined.";
    };

    /**
     * @return the total width of the W boson
     */
    virtual double GammaW() const {
        throw "EWModel::GammaW() is undefined.";
    };
    
    /**
     * @return NP contribution to oblique parameter S
     */
    virtual double obliqueS() const {
        throw "EWModel::obliqueS() is undefined.";
    };
        
    /**
     * @return NP contribution to oblique parameter T
     */
    virtual double obliqueT() const {
        throw "EWModel::obliqueT() is undefined.";
    };
    
    /**
     * @return NP contribution to oblique parameter U
     */
    virtual double obliqueU() const {
        throw "EWModel::obliqueU() is undefined.";
    };
    
    
};

#endif	/* EWMODEL_H */
