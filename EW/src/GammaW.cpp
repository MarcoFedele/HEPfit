/* 
 * Copyright (C) 2012-2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "GammaW.h"


double GammaW::getThValue() 
{  
    double Gamma_W;
    EW::EWTYPE myEWTYPE = myEW.getEWTYPE();

    if (myEWTYPE==EW::EWCHMN)  
        Gamma_W = myEW.getMyEW_CHMN().GammaW();
    else if (myEWTYPE==EW::EWABC || myEWTYPE==EW::EWABC2) 
        throw std::runtime_error("GammaW::getThValue() is not implemented for EW::EWABC");  
    else {
        Gamma_W = SM.GammaW();
        
        /* Oblique NP */
        if ( myEW.checkSTU() && !SM.IsFlagNotLinearizedNP() ) {
            double Wbar = 0.0;        
            if (SM.ModelName()=="NPSTUVWXY") {
                Wbar = (myEW.V() - myEW.W())/SM.alphaMz();
            }

            if(myEWTYPE==EW::EWBURGESS) {
                // TEST: the fit result by Gfitter in arXiv:1209.2716, 
                //       corresponding to MH=125.7 and Mt=173.52 
                //Gamma_W = 2.091; 
                
                Gamma_W *= 1.0 - 0.00723*myEW.S() + 0.0111*myEW.T() + 0.00849*myEW.U() + 0.00781*Wbar;
            } else {
                double alpha = myEW.alpha();  
                double c = sqrt(myEW.cW2_SM());
                double c2 = myEW.cW2_SM();
                double s2 = myEW.sW2_SM();
            
                // TEST!!
                //std::cout << 3.0*alpha*alpha*c*SM.getMz()/8.0/s2/(c2-s2) << " "
                //          << 3.0*alpha*alpha*c*SM.getMz()/8.0/s2/(c2-s2)*2.0*c2 << " "
                //          << 3.0*alpha*alpha*c*SM.getMz()/8.0/s2/2.0/s2 << std::endl;

                // TEST!!
                //Gamma_W *= 1.0 - 0.00723*myEW.S() + 0.0111*myEW.T() + 0.00849*myEW.U() + 0.00781*Wbar;
                
                Gamma_W -= 3.0*alpha*alpha*c*SM.getMz()/8.0/s2/(c2-s2)
                           *( myEW.S() - 2.0*c2*myEW.T() - (c2-s2)*myEW.U()/2.0/s2 
                              - 2.0*(c2 - s2)*Wbar );
            }
        }
    }
 
    return Gamma_W;
}
