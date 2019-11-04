/*
 * Copyright (C) 2018 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "LoopMediatorsA.h"

const std::string LoopMediatorsA::LoopMediatorsAvars[NLoopMediatorsAvars] = {
        "GammaQ", "Gammamu", "mpsi", "mphiQ", "mphiL", "charge", "chi", "chiBB",
        "chi_M", "chiBB_M", "eta", "etaBB", "eta_M", "etaBB_M", "WCscale"};

LoopMediatorsA::LoopMediatorsA() : StandardModel(), LoopMediatorsAM(*this)
{

    SMM.setObj((StandardModelMatching&) LoopMediatorsAM.getObj());
    ModelParamMap.insert(std::make_pair("GammaQ", std::cref(GammaQ)));
    ModelParamMap.insert(std::make_pair("Gammamu", std::cref(Gammamu)));
    ModelParamMap.insert(std::make_pair("mpsi", std::cref(mpsi)));
    ModelParamMap.insert(std::make_pair("mphiQ", std::cref(mphiQ)));
    ModelParamMap.insert(std::make_pair("mphiL", std::cref(mphiL)));
    ModelParamMap.insert(std::make_pair("charge", std::cref(charge)));

    ModelParamMap.insert(std::make_pair("chi", std::cref(chi)));
    ModelParamMap.insert(std::make_pair("chiBB", std::cref(chiBB)));
    ModelParamMap.insert(std::make_pair("chi_M", std::cref(chi)));
    ModelParamMap.insert(std::make_pair("chiBB_M", std::cref(chiBB)));

    ModelParamMap.insert(std::make_pair("eta", std::cref(eta)));
    ModelParamMap.insert(std::make_pair("etaBB", std::cref(etaBB)));
    ModelParamMap.insert(std::make_pair("eta_M", std::cref(eta)));
    ModelParamMap.insert(std::make_pair("etaBB_M", std::cref(etaBB)));

    ModelParamMap.insert(std::make_pair("WCscale", std::cref(WCscale)));
}

LoopMediatorsA::~LoopMediatorsA()
{
    if (IsModelInitialized()) {
        }
}

///////////////////////////////////////////////////////////////////////////
// Initialization

bool LoopMediatorsA::InitializeModel()
{
    setModelInitialized(StandardModel::InitializeModel());
    return(true);
}

bool LoopMediatorsA::Init(const std::map<std::string, double>& DPars)
{
    return(StandardModel::Init(DPars));
}

bool LoopMediatorsA::PreUpdate()
{
    if(!StandardModel::PreUpdate()) return (false);

    return (true);
}

bool LoopMediatorsA::Update(const std::map<std::string, double>& DPars)
{

    if(!PreUpdate()) return (false);

    UpdateError = false;

    for (std::map<std::string, double>::const_iterator it = DPars.begin(); it != DPars.end(); it++)
        setParameter(it->first, it->second);

    if (UpdateError) return (false);

    if(!PostUpdate()) return (false);

    return (true);
}

bool LoopMediatorsA::PostUpdate()
{
    if(!StandardModel::PostUpdate()) return (false);

    double GammaQ2 = GammaQ*GammaQ;
    double Gammamu2 = Gammamu * Gammamu;
    double mpsi2 = mpsi * mpsi;
    double mphiQ2 = mphiQ * mphiQ;
    double mphiL2 = mphiL * mphiL;
    double yQ = mphiQ2 / mpsi2;
    double yL = mphiL2 / mpsi2;
    double qphi = -charge;
    double qpsi = charge - 1.;

    double M_PI2 = M_PI*M_PI;
    double mmu = getLeptons(MU).getMass();
    double Norm = - sqrt(2.) / (4. * GF * 0.0411494) / (32. * M_PI * ale); // N.B. took the abs of CKM, hence changed overall sign

    C1 = 1. / (128. * M_PI2 * mpsi2) * GammaQ2 * (chiBB * etaBB * F9(yQ,yL) + 2. * chiBB_M * etaBB_M * G9(yQ,yL));
    C2 = 0.;
    C3 = 0.;
    C4 = 0.;
    C5 = 0.;

    C1p = 0.;
    C2p = 0.;
    C3p = 0.;

    C7 = 0.;
    C8 = 0.;
    C9 =    Norm * GammaQ * Gammamu2 / mpsi2 * (chi * eta * F9(yQ,yL) + 2. * chi_M * eta_M * G9(yQ,yL));
    C10 = - Norm * GammaQ * Gammamu2 / mpsi2 * (chi * eta * F9(yQ,yL) + 2. * chi_M * eta_M * G9(yQ,yL));
    CS = 0.;
    CP = 0.;

    C7p = 0.;
    C8p = 0.;
    C9p = 0.;
    C10p = 0.;
    CSp = 0.;
    CPp = 0.;

    /*std::cout << "C1 " << C1 << std::endl;
    std::cout << "C9 " << C9.abs() << std::endl;*/

    /*Deltaamu = mmu*mmu / ( 8. * M_PI2 * mphi2) * (
            (GammamuL2 + GammamuR2) * (qphi*F7t(yE) - qpsi*F7(yE))
            + 4./sqrt(2.) * v() * lambdaE / mmu * (2.*GammamuL*GammamuR) * (qphi*G7t(yE) - qpsi*G7(yE))
          );*/

    /* Necessary for updating StandardModel parameters in StandardModelMatching,
     * and LoopMediatorsA and LoopMediatorsA-derived parameters in LoopMediatorsAMatching */
    LoopMediatorsAM.getObj().updateLoopMediatorsAParameters();

    return (true);
}

void LoopMediatorsA::setParameter(const std::string name, const double& value)
{
    if(name.compare("GammaQ") == 0)
        GammaQ = value;
    else if(name.compare("Gammamu") == 0)
        Gammamu = value;
    else if(name.compare("mpsi") == 0)
        mpsi = value;
    else if(name.compare("mphiQ") == 0)
        mphiQ = value;
    else if(name.compare("mphiL") == 0)
        mphiL = value;
    else if(name.compare("charge") == 0)
        charge = value;
    else if(name.compare("chi") == 0)
        chi = value;
    else if(name.compare("chiBB") == 0)
        chiBB = value;
    else if(name.compare("chi_M") == 0)
        chi_M = value;
    else if(name.compare("chiBB_M") == 0)
        chiBB_M = value;
    else if(name.compare("eta") == 0)
        eta = value;
    else if(name.compare("etaBB") == 0)
        etaBB = value;
    else if(name.compare("eta_M") == 0)
        eta_M = value;
    else if(name.compare("etaBB_M") == 0)
        etaBB_M = value;
    else if(name.compare("WCscale") == 0)
        WCscale = value;
    else
        StandardModel::setParameter(name,value);
}

bool LoopMediatorsA::CheckParameters(const std::map<std::string, double>& DPars)
{
    for (int i = 0; i < NLoopMediatorsAvars; i++) {
        if (DPars.find(LoopMediatorsAvars[i]) == DPars.end()) {
            std::cout << "ERROR: missing mandatory LoopMediatorsA parameter " << LoopMediatorsAvars[i] << std::endl;
            //raiseMissingModelParameterCount();
            //addMissingModelParameter(LoopMediatorsAvars[i]);
        }
    }
    return(StandardModel::CheckParameters(DPars));
}

///////////////////////////////////////////////////////////////////////////
// Flags

bool LoopMediatorsA::setFlag(const std::string name, const bool value)
{
    bool res = false;

    res = StandardModel::setFlag(name,value);

    return(res);
}

///////////////////////////////////////////////////////////////////////////
// Loop Functions

double LoopMediatorsA::F9(double x, double y)
{
    double ym1 = y - 1.;
    double xm1 = x - 1.;

    if (x == 1. && y == 1.)
        return 1./3.;
    else if (x == 1.)
        return (-3.*y*y + 2.*y*y*log(y) + 4.*y - 1.)/2./ym1/ym1/ym1;
    else if (y == 1.)
        return (-3.*x*x + 2.*x*x*log(x) + 4.*x - 1.)/2./xm1/xm1/xm1;
    else if (x == y)
        return (y*y - 2.*y*log(y) - 1.)/ym1/ym1/ym1;
    else
        return 1./xm1/ym1 + x*x*log(x)/xm1/xm1/(x-y) + y*y*log(y)/ym1/ym1/(y-x);
}

double LoopMediatorsA::G9(double x, double y)
{
    double ym1 = y - 1.;
    double xm1 = x - 1.;

    if (x == 1. && y == 1.)
        return -1./6.;
    else if (x == 1.)
        return (-y*y + 2.*y*log(y) + 1.)/2./ym1/ym1/ym1;
    else if (y == 1.)
        return (-x*x + 2.*x*log(x) + 1.)/2./xm1/xm1/xm1;
    else if (x == y)
        return (2.*ym1 - (y + 1.)*log(y))/ym1/ym1/ym1;
    else
        return 1./xm1/ym1 + x*log(x)/xm1/xm1/(x-y) + y*log(y)/ym1/ym1/(y-x);
}

double LoopMediatorsA::F7(double x)
{
    double xm1 = x - 1.;

    if (x == 1.)
        return 1./24.;
    else
        return (x*x*x - 6.*x*x + 3.*x + 2. + 6.*x*log(x)) / (12.*xm1*xm1*xm1*xm1);
}

double LoopMediatorsA::F7t(double x)
{
    if (x == 1.)
        return 1./24.;
    else
        return F7(1./x) / x;
}

double LoopMediatorsA::G7(double x)
{
    double xm1 = x - 1.;

    if (x == 1.)
        return 1./12.;
    else
        return (x*x - 4.*x + 3. + 2.*log(x)) / (8.*xm1*xm1*xm1);
}

double LoopMediatorsA::G7t(double x)
{
    double xm1 = x - 1.;

    if (x == 1.)
        return 1./24.;
    else
        return (x*x - 1. - 2.*x*log(x)) / (8.*xm1*xm1*xm1);
}




/////////////////////////////////////////////////////////////////////////



Deltaamu_LMA::Deltaamu_LMA(const StandardModel& SM_i)
: ThObservable(SM_i), myLM(static_cast<const LoopMediatorsA*> (&SM_i))
{
};

Deltaamu_LMA::~Deltaamu_LMA()
{
};

double Deltaamu_LMA::computeThValue()
{
    return myLM->getDeltaamu();
}



C9mC10_LMA::C9mC10_LMA(const StandardModel& SM_i)
: ThObservable(SM_i), myLM(static_cast<const LoopMediatorsA*> (&SM_i))
{
};

C9mC10_LMA::~C9mC10_LMA()
{
};

double C9mC10_LMA::computeThValue()
{
    return myLM->getC9();
}



mpsi_mphiQ::mpsi_mphiQ(const StandardModel& SM_i)
: ThObservable(SM_i), myLM(static_cast<const LoopMediatorsA*> (&SM_i))
{
};

mpsi_mphiQ::~mpsi_mphiQ()
{
};

double mpsi_mphiQ::computeThValue()
{
    return myLM->getmpsi() - myLM->getmphiQ();
}



mpsi_mphiL::mpsi_mphiL(const StandardModel& SM_i)
: ThObservable(SM_i), myLM(static_cast<const LoopMediatorsA*> (&SM_i))
{
};

mpsi_mphiL::~mpsi_mphiL()
{
};

double mpsi_mphiL::computeThValue()
{
    return myLM->getmpsi() - myLM->getmphiL();
}



mphiQ_mphiL::mphiQ_mphiL(const StandardModel& SM_i)
: ThObservable(SM_i), myLM(static_cast<const LoopMediatorsA*> (&SM_i))
{
};

mphiQ_mphiL::~mphiQ_mphiL()
{
};

double mphiQ_mphiL::computeThValue()
{
    return myLM->getmphiQ() - myLM->getmphiL();
}
