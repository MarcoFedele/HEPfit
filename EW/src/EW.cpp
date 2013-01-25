/* 
 * Copyright (C) 2012 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "EW.h"
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <stdexcept>


EW::EW(const StandardModel& SM_i) 
: ThObsType(SM_i), SM(SM_i), myEW_CHMN(SM_i), myEW_ABC(SM_i) 
{
    bDebug = SM_i.isBDebug();
}


////////////////////////////////////////////////////////////////////////

EW::EWTYPE EW::getEWTYPE() const 
{
    if ( SM.IsFlagEWABC() && SM.IsFlagEWABC2() ) 
        throw std::runtime_error("Flags EWABC and EWABC2 cannot be set to true simultaneously");
    if ( SM.IsFlagEWBURGESS() && SM.IsFlagEWCHMN() ) 
        throw std::runtime_error("Flags EWBURGESS and EWCHMN cannot be set to true simultaneously");
    if ( SM.IsFlagEWABC() && SM.IsFlagEWCHMN() ) 
        throw std::runtime_error("Flags EWABC and EWCHMN cannot be set to true simultaneously");
    if ( SM.IsFlagEWABC2() && SM.IsFlagEWCHMN() ) 
        throw std::runtime_error("Flags EWABC2 and EWCHMN cannot be set to true simultaneously");
    if ( SM.IsFlagEWABC() && SM.IsFlagEWBURGESS() ) 
        throw std::runtime_error("Flags EWABC and EWBURGESS cannot be set to true simultaneously");
    if ( SM.IsFlagEWABC2() && SM.IsFlagEWBURGESS() ) 
        throw std::runtime_error("Flags EWABC2 and EWBURGESS cannot be set to true simultaneously");
            
    if ( SM.IsFlagEWBURGESS() ) return EWBURGESS;
    else if ( SM.IsFlagEWCHMN() ) return EWCHMN;
    else if ( SM.IsFlagEWABC() ) return EWABC;
    else if ( SM.IsFlagEWABC2() ) return EWABC2;
    else return EWDEFAULT;
}


bool EW::checkModelForSTU() const 
{
    std::string Model = SM.ModelName();
    //std::cout << "Model in EW: " << Model << std::endl; // TEST

    if (Model=="StandardModel" || Model=="SUSY" || Model=="NPEpsilons" 
            || Model=="NPHiggs" || Model=="NPZbbbar" )
        return false;
    else if (Model=="NPSTU" || Model=="NPSTUVWXY" || Model=="THDM")
        return true;
    else 
        throw std::runtime_error("Error in EW::checkModelForSTU()");  
}


double EW::Ql(const StandardModel::lepton l) const 
{
    return ( SM.getLeptons(l).getCharge() );
}        
 

double EW::Qq(const StandardModel::quark q) const 
{
    return ( SM.getQuarks(q).getCharge() );
}


double EW::alpha0() const 
{
    return ( SM.alphaMz() );    

//    return SM.getAle(); // TEST!!
}


double EW::Mw0() const 
{
    return ( SM.Mw0() );

//    return SM.Mw_tree(); // TEST!!
//    return SM.Mw(); // TEST!!
}


double EW::s02() const 
{
    return ( SM.s02() );

//    return ( 1.0 - c02() ); // TEST!!
//    return SM.sW2(); // TEST?!!
}


double EW::c02() const 
{
    return ( SM.c02() );

//    return ( SM.Mw_tree()*SM.Mw_tree()/SM.getMz()/SM.getMz() ); // TEST!!
//    return SM.cW2(); // TEST?!!
}


////////////////////////////////////////////////////////////////////////

double EW::sin2thetaEff(const StandardModel::lepton l) const 
{
    complex gV_over_gA = SM.gVl(l)/SM.gAl(l);
    double absQf = fabs(Ql(l));
    return ( 1.0/4.0/absQf*(1.0 - gV_over_gA.real()) );
}


double EW::sin2thetaEff(const StandardModel::quark q) const 
{
    complex gV_over_gA = SM.gVq(q)/SM.gAq(q);
    double absQf = fabs(Qq(q));
    return ( 1.0/4.0/absQf*(1.0 - gV_over_gA.real()) );
}


double EW::Gamma_l(const StandardModel::lepton l) const 
{
    complex rhoZ_l = SM.rhoZ_l(l);
    complex gV_over_gA = SM.gVl(l)/SM.gAl(l);
    double alphaMz = SM.alphaMz();
    double Q = Ql(l);
    double xl = pow(SM.getLeptons(l).getMass()/SM.getMz(), 2.0);
    double G0 = SM.getGF()*pow(SM.getMz(),3.0)/24.0/sqrt(2.0)/M_PI;
    double Gamma = G0*rhoZ_l.abs()*sqrt(1.0 - 4.0*xl)
                   * ( (1.0 + 2.0*xl)*(gV_over_gA.abs2() + 1.0) - 6.0*xl )
                   * ( 1.0 + 3.0/4.0*alphaMz/M_PI*pow(Q,2.0) );
    return Gamma;
}


double EW::Gamma_q(const StandardModel::quark q) const 
{
    if (q==StandardModel::TOP) return 0.0;

    complex rhoZ_q = SM.rhoZ_q(q);
    complex gV_over_gA = SM.gVq(q)/SM.gAq(q);
    
    double G0 = SM.getGF()*pow(SM.getMz(),3.0)/24.0/sqrt(2.0)/M_PI;    
    double Gamma = 3.0*G0*rhoZ_q.abs()
                   * ( gV_over_gA.abs2()*SM.RVq(q) + SM.RAq(q) ) + SM.Delta_EWQCD(q);
    return Gamma;
}


double EW::Gamma_inv() const 
{
    return ( Gamma_l(SM.NEUTRINO_1) + Gamma_l(SM.NEUTRINO_2) 
             + Gamma_l(SM.NEUTRINO_3) );
}


double EW::Gamma_had() const 
{
    double Gamma_had_tmp = Gamma_q(SM.UP) + Gamma_q(SM.DOWN) + Gamma_q(SM.CHARM)
                           + Gamma_q(SM.STRANGE) + Gamma_q(SM.BOTTOM);

    /* We neglect the singlet contribution, since it makes the connection between 
       rho_Z^b and the two-loop approximate formula of R_b^0 complicated, and is 
       negligible compared to the experimental uncertainty. */
    //double G0 = SM.getGF()*pow(SM.getMz(),3.0)/24.0/sqrt(2.0)/M_PI; 
    //Gamma_had_tmp += 4.0*3.0*G0*SM.RVh();

    return Gamma_had_tmp;    
}


double EW::Gamma_Z() const 
{
    return ( Gamma_l(SM.ELECTRON) + Gamma_l(SM.MU) + Gamma_l(SM.TAU) 
             + Gamma_inv() + Gamma_had() );
}


double EW::sigma0_l(const StandardModel::lepton l) const 
{
    return ( 12.0*M_PI*Gamma_l(SM.ELECTRON)*Gamma_l(l)
             /SM.getMz()/SM.getMz()/Gamma_Z()/Gamma_Z() );
}


double EW::sigma0_had() const 
{
     return (12.0*M_PI*Gamma_l(SM.ELECTRON)*Gamma_had()
            /SM.getMz()/SM.getMz()/Gamma_Z()/Gamma_Z());
}


double EW::A_l(const StandardModel::lepton l) const 
{
    double Re_gV_over_gA = (SM.gVl(l)/SM.gAl(l)).real();
    return ( 2.0*Re_gV_over_gA/(1.0+pow(Re_gV_over_gA,2.0)) );
}


double EW::A_q(const StandardModel::quark q) const 
{
    double Re_gV_over_gA = (SM.gVq(q)/SM.gAq(q)).real();    
    return ( 2.0*Re_gV_over_gA/(1.0+pow(Re_gV_over_gA,2.0)) );
}


