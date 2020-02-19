/*
 * Copyright (C) 2018 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef MYOBS_H
#define MYOBS_H

class StandardModel;
#include <gsl/gsl_integration.h>
#include "ThObservable.h"
#include "QCD.h"
#include "OrderScheme.h"

class gm2_Zprime : public ThObservable {
public:
    /**
     * constructor
     * @param Flavour
     */
    gm2_Zprime(const StandardModel& SM_i);

    /**
     *
     */
    double computeThValue();


protected:
    double IntFunct(double z);

private:
    bool LoopModelDM;
    bool gmu_logscale;
    bool rAV_parametric;

    gsl_function FInt;/**< GSL integral variable */
    double Int;/**< GSL integral variable */
    double errInt;/**< GSL integral variable */
    gsl_integration_cquad_workspace * w_Int;/**< GSL integral variable */

    double mV;
    double gmuV;
    double rAV;
    double gmuA;

    double mmu;
    double x;

};

#endif /* MYOBS_H */
