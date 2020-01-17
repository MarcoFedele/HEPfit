/* 
 * Copyright (C) 2018 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef MLL_H
#define MLL_H

class StandardModel;
#include "ThObservable.h"
#include "QCD.h"
#include "OrderScheme.h"

class Mll : public ThObservable {
public:   
    /**
     * constructor
     * @param Flavour
     */
    Mll(const StandardModel& SM_i, int obsFlag, QCD::meson meson_i, QCD::lepton lep_i);
    
    /*double F9(double x);
    double G9(double x);
    
    gslpp::complex C10_NP(double q2);*/
    
    /**
     * 
     * @brief hep-ph/9512380v2
     * @return theoretical value of |\f$ BR(B_s \rightarrow \mu \bar{\mu}) \f$|
     */
    double computeThValue();
    double computeAmumu(orders order);
    double computeSmumu(orders order);
    
    
protected:
    
    /**
     * 
     * @param order
     * @param order_qed
     * @return the short distance contribution to the 
     * |\f$ BR(B_s \rightarrow \mu \bar{\mu}) \f$|
     */
    void computeAmpSq(orders order, orders_qed order_qed, double mu);
    void computeObs(orders order, orders_qed order_qed);
    
private:
    
    QCD::lepton lep; /**< Final leptons type. */
    QCD::meson meson;
    double ys;
    gslpp::complex CKM_factor;
    double beta;
    double mBs;
    double mW;
    double mlep;
    double mb;
    double ms;
    double chiral;
    double absP;
    double argP;
    double absS;
    double argS;
    double ampSq;
    double Amumu;
    double Smumu;
    double phiNP;
    double timeInt;
    int obs;
    gslpp::complex C_10;
    gslpp::complex C_10p;
    gslpp::complex C_S;
    gslpp::complex C_Sp;
    gslpp::complex C_P;
    gslpp::complex C_Pp;
    bool FixedWCbtos;
    //bool LoopModelDM;
    
    gslpp::vector<gslpp::complex> ** allcoeff;
    gslpp::vector<gslpp::complex> ** allcoeffprime;
    gslpp::vector<gslpp::complex> ** allcoeff_noSM;
    
    /*double ysybgD;
    double rVA;
    double QB;
    double mB_NP;
    double mchi_NP;
    double mV_NP;
    double mB2_NP;
    double mchi2_NP;
    double mV2_NP;
    double y_NP;
    double gmuV_NP;
    double gmuA_NP;
    double Norm_NP;*/

};

#endif /* MLL_H */

