/* 
 * File:   LEP2sigmaMu.cpp
 * Author: mishima
 */

#include "LEP2sigmaMu.h"


double LEP2sigmaMu::getThValue() { 
    double s = sqrt_s*sqrt_s;
    double Mw = SM.Mw(); 
    double GammaZ = myEW.Gamma_Z();

    if (!SM.getEWSM()->checkForLEP2(SMparams_cache, bRCs_cache, s, Mw, GammaZ, bRCs))
        SMresult_cache = SM.sigma_l_LEP2(StandardModel::MU, s, Mw, GammaZ, bRCs);
    double sigma_mu = SMresult_cache;
    
    if ( myEW.checkModelForSTU() ) {
        sigma_mu += myLEP2oblique.sigma_l_LEP2_NP(StandardModel::MU, s);      
    }
    
    return ( sigma_mu*GeVminus2_to_nb*1000. );
}
        

