/* 
 * Copyright (C) 2012-2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "sigmaHadron.h"


double sigmaHadron::getThValue() 
{ 
    double sigma_had;
    EW::EWTYPE myEWTYPE = myEW.getEWTYPE();

    if (myEWTYPE==EW::EWCHMN)  
        sigma_had = myEW.getMyEW_CHMN().sigma0_had();
    else if (myEWTYPE==EW::EWABC) 
        sigma_had = myEW.getMyEW_ABC().sigma0_had(SM.epsilon1(),SM.epsilon3(),SM.epsilonb());
    else if (myEWTYPE==EW::EWABC2) {
        double delta_als = (SM.Als(SM.getMz(),FULLNNLO) - 0.119)/M_PI;
        double delta_alpha = (SM.alphaMz() - 1.0/128.90)/SM.getAle();
        double sigma_h0 = 41.420*(1.0 - 0.41*delta_als + 0.03*delta_alpha)/GeVminus2_to_nb;
        sigma_had = sigma_h0*(1.0 - 0.03*SM.epsilon1() + 0.04*SM.epsilon3() - 0.20*SM.epsilonb());
    } else {   
        sigma_had = myEW.sigma0_had();
        
        /* Oblique NP */
        if ( myEW.checkSTU() && !SM.IsFlagNotLinearizedNP() ) {
            if(myEWTYPE==EW::EWBURGESS) {
                // TEST: the fit result by Gfitter in arXiv:1209.2716, 
                //       corresponding to MH=125.7 and Mt=173.52 
                //sigma_had = 41.479/GeVminus2_to_nb;                
                
                double delta_l = - 0.000192*myEW.S() + 0.000790*myEW.T();
                double delta_had = - 0.00901*myEW.S() + 0.0200*myEW.T();
                double delta_Z = - 0.00961*myEW.S() + 0.0263*myEW.T();
                sigma_had *= 1.0 + delta_l/myEW.Gamma_l(SM.ELECTRON)
                             + delta_had/myEW.Gamma_had()
                             - 2.0*delta_Z/myEW.Gamma_Z();
            } else {
                double alpha = myEW.alpha();  
                double Mz = SM.getMz();
                double c2 = myEW.cW2_SM();
                double s2 = myEW.sW2_SM();
                double s4 = s2*s2;
                double s6 = s4*s2;        
                double s8 = s6*s2;
                
                sigma_had -= 72.0*M_PI*alpha
                             *(729.0-4788.0*s2+8352.0*s4-6176.0*s6+640.0*s8)
                             /Mz/Mz/pow(63.0-120.0*s2+160.0*s4, 3.0)/(c2-s2)
                             *( myEW.S() - 4.0*c2*s2*myEW.T() );         
            }
        }

        /* NP contribution to the Zff vertex */
        if ( !SM.IsFlagNotLinearizedNP() ) {
            bool nonZeroNP = false;

            double delGVl[6], delGAl[6], delGVq[6], delGAq[6];
            for (int p=0; p<6; ++p) {
                delGVl[p] = SM.deltaGVl((StandardModel::lepton)p);
                delGAl[p] = SM.deltaGAl((StandardModel::lepton)p);
                delGVq[p] = SM.deltaGVq((StandardModel::quark)p);
                delGAq[p] = SM.deltaGAq((StandardModel::quark)p);
                if (delGVl[p]!=0.0 || delGAl[p]!=0.0
                        || delGVq[p]!=0.0 || delGAq[p]!=0.0)
                    nonZeroNP = true;
            }

            if (nonZeroNP) {
                double gVf, gAf;
                double Gl[6], deltaGl[6], Gq[6], deltaGq[6];
                double Gq_sum = 0.0, delGq_sum = 0.0;
                double Gf_sum = 0.0, delGf_sum = 0.0;
                for (int p=0; p<6; ++p) {
                    gVf = SM.StandardModel::gVl((StandardModel::lepton)p).real();
                    gAf = SM.StandardModel::gAl((StandardModel::lepton)p).real();
                    Gl[p] = gVf*gVf + gAf*gAf;
                    deltaGl[p] = 2.0*(gVf*delGVl[p] + gAf*delGAl[p]);

                    gVf = SM.StandardModel::gVq((StandardModel::quark)p).real();
                    gAf = SM.StandardModel::gAq((StandardModel::quark)p).real();
                    Gq[p] = gVf*gVf + gAf*gAf;
                    deltaGq[p] = 2.0*(gVf*delGVq[p] + gAf*delGAq[p]);

                    Gq_sum += 3.0*Gq[p];
                    Gf_sum += Gl[p] + 3.0*Gq[p];
                    delGq_sum += 3.0*deltaGq[p];
                    delGf_sum += deltaGl[p] + 3.0*deltaGq[p];
                }

                sigma_had += 12.0*M_PI/SM.getMz()/SM.getMz()
                             *Gl[(int)SM.ELECTRON]*Gq_sum/Gf_sum/Gf_sum
                             *( deltaGl[(int)SM.ELECTRON]/Gl[(int)SM.ELECTRON]
                                + delGq_sum/Gq_sum - 2.0*delGf_sum/Gf_sum );
            }
        }
        
        /* TEST */
        //sigma_had -= myEW.sigma0_had();
    }
    
    return ( sigma_had*GeVminus2_to_nb );
}
        


