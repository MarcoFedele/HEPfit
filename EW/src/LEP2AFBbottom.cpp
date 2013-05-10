/* 
 * Copyright (C) 2012-2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include <Math/Functor.h>
#include <Math/Integrator.h>
#include <Math/AllIntegrationTypes.h>
#include "LEP2AFBbottom.h"


double LEP2AFBbottom::getThValue()
{ 
    Mw = SM.Mw(); 
    GammaZ = myEW.Gamma_Z();

    if (!checkSMparams(s, Mw, GammaZ)) {
        mq_cache = m_q(SM.BOTTOM, sqrt_s);
        
        double AFB_noBox, sigma;
        if (!flag[ISR])
            AFB_noBox = AFB_NoISR_q();
        else {
            ROOT::Math::Functor1D wf(this, &LEP2AFBbottom::Integrand_AFBnumeratorWithISR_q);
            ROOT::Math::Integrator ig(wf, ROOT::Math::IntegrationOneDim::kADAPTIVESINGULAR);
            ig.SetAbsTolerance(1.E-14); // desired absolute error
            ig.SetRelTolerance(1.E-4); // desired relative error
            double numerator = ig.Integral(0.0, 1.0-0.85*0.85); // interval
            //std::cout << numerator << " +- " << ig.Error() << std::endl;
            // results: 4.0e-9 -- 1.8e-8
            
            // denominator
            myLEP2sigmaBottom.setFlags(flag);
            sigma = myLEP2sigmaBottom.getThValue()/GeVminus2_to_nb/1000.0;
            
            AFB_noBox = numerator/sigma;
        }
        SMresult_cache = AFB_noBox;
        
        if (flag[WeakBox]) {
            // numerator
            ROOT::Math::Functor1D wf_F(this, &LEP2AFBbottom::Integrand_dsigmaBox_q);
            ROOT::Math::Integrator ig_F(wf_F, ROOT::Math::IntegrationOneDim::kADAPTIVESINGULAR);
            ig_F.SetAbsTolerance(1.E-17); // desired absolute error
            ig_F.SetRelTolerance(1.E-4); // desired relative error
            double sigma_box_F = ig_F.Integral(0.0, 1.0); // interval
            //std::cout << sigma_box_F << " +- " << ig_F.Error()  << std::endl;
            // results: 3.7e-12 -- 1.5e-10
            //
            ROOT::Math::Functor1D wf_B(this, &LEP2AFBbottom::Integrand_dsigmaBox_q);
            ROOT::Math::Integrator ig_B(wf_B, ROOT::Math::IntegrationOneDim::kADAPTIVESINGULAR);
            ig_B.SetAbsTolerance(1.E-17); // desired absolute error
            ig_B.SetRelTolerance(1.E-4); // desired relative error
            double sigma_box_B = ig_B.Integral(-1.0, 0.0); // interval
            //std::cout << sigma_box_B << " +- " << ig_B.Error()  << std::endl;
            // results: 2.2e-12 -- 5.0e-11
                        
            // denominator
            if (!flag[ISR]) {
                myLEP2sigmaBottom.setFlags(flag);
                sigma = myLEP2sigmaBottom.getThValue()/GeVminus2_to_nb/1000.0;
            }
            
            SMresult_cache += (sigma_box_F - sigma_box_B)/sigma;
        }

        if ( myEW.checkSTUVWXY() && SM.IsFlagFixedAllSMparams() ) {
            double ObParam[7];
            for (int i=0; i<7; i++) {
                SetObParam((LEP2oblique::Oblique)i, ObParam);
                Coeff_cache[i] 
                    = myLEP2oblique.AFB_q_LEP2_NP(StandardModel::BOTTOM, s, mq_cache, ObParam);
            }
        }
    }
    double AFB_b = SMresult_cache;
    
    #ifdef LEP2TEST
    AFB_b = myTEST.AFBbottomTEST(sqrt_s);
    #endif
    
    if ( myEW.checkSTUVWXY() ) {
        if ( SM.IsFlagFixedAllSMparams() ) {
            AFB_b += Coeff_cache[myLEP2oblique.Shat]*myEW.Shat()
                   + Coeff_cache[myLEP2oblique.That]*myEW.That()
                   + Coeff_cache[myLEP2oblique.Uhat]*myEW.Uhat()
                   + Coeff_cache[myLEP2oblique.V]*myEW.V()
                   + Coeff_cache[myLEP2oblique.W]*myEW.W()
                   + Coeff_cache[myLEP2oblique.X]*myEW.X()
                   + Coeff_cache[myLEP2oblique.Y]*myEW.Y();
        } else {
            double ObParam[7] = {myEW.Shat(), myEW.That(), myEW.Uhat(),
                                 myEW.V(), myEW.W(), myEW.X(), myEW.Y()};     
            AFB_b += myLEP2oblique.AFB_q_LEP2_NP(StandardModel::BOTTOM, s, mq_cache, ObParam);
        }
    }
    
    return AFB_b;
}
        
