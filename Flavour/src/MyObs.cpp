/*
 * Copyright (C) 2018 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "MyObs.h"
#include "StandardModel.h"
#include "gslpp_function_adapter.h"
#include <boost/bind.hpp>

gm2_Zprime::gm2_Zprime(const StandardModel& SM_i)
: ThObservable(SM_i)
{
    LoopModelDM = false;
    
    LoopModelDM = SM.getFlavour().getFlagLoopModelDM();
    
    if (LoopModelDM) {
        setParametersForObservable({ "mV_NP", "gmuV_NP", "rAV_NP" });
        w_Int = gsl_integration_cquad_workspace_alloc(100);
    }
};


double gm2_Zprime::IntFunct(double z)
{
    double x2 = x*x;
    
    return (2.*gmuV*gmuV*(1.-z)*z*z + 
            gmuA*gmuA*(2.*(1.-z)*(z-4.)*z - 4.*x2*z*z*z)) /
            (x2*z + (1.-z)*(1-x2*z));
}


double gm2_Zprime::computeThValue()
{
    
    if (LoopModelDM) { 
        mV = SM.getOptionalParameter("mV_NP");
        gmuV = SM.getOptionalParameter("gmuV_NP");
        gmuA = SM.getOptionalParameter("rAV_NP") * gmuV;
        
        mmu = SM.getLeptons(QCD::MU).getMass();
        x = mmu/mV;
        
        FInt = convertToGslFunction(boost::bind(&gm2_Zprime::IntFunct, &(*this), _1));
        if (gsl_integration_cquad(&FInt, 0., 1., 1.e-2, 1.e-1, w_Int, &Int, &errInt, NULL) != 0) return std::numeric_limits<double>::quiet_NaN();
                
        
        return Int * x*x / 8. / M_PI;
        }
    else throw std::runtime_error("gm2_Zprime::computeThValue(): Observable only defined if LoopModelDM flag is set to true");
    return (EXIT_FAILURE);
}
