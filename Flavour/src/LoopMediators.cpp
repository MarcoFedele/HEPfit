/*
 * Copyright (C) 2018 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "LoopMediators.h"

const std::string LoopMediators::LoopMediatorsvars[NLoopMediatorsvars] = {
        "GammaQ", "Gammamu", "mphi", "mpsiQ", "mpsiL", "charge", "chi", "chiBB",
        "chi_M", "chiBB_M", "eta", "etaBB", "eta_M", "etaBB_M", "WCscale"};

LoopMediators::LoopMediators() : StandardModel(), LoopMediatorsM(*this)
{

    SMM.setObj((StandardModelMatching&) LoopMediatorsM.getObj());
    ModelParamMap.insert(std::make_pair("GammaQ", std::cref(GammaQ)));
    ModelParamMap.insert(std::make_pair("Gammamu", std::cref(Gammamu)));
    ModelParamMap.insert(std::make_pair("mphi", std::cref(mphi)));
    ModelParamMap.insert(std::make_pair("mpsiQ", std::cref(mpsiQ)));
    ModelParamMap.insert(std::make_pair("mpsiL", std::cref(mpsiL)));
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

LoopMediators::~LoopMediators()
{
    if (IsModelInitialized()) {
        }
}

///////////////////////////////////////////////////////////////////////////
// Initialization

bool LoopMediators::InitializeModel()
{
    setModelInitialized(StandardModel::InitializeModel());
    return(true);
}

bool LoopMediators::Init(const std::map<std::string, double>& DPars)
{
    return(StandardModel::Init(DPars));
}

bool LoopMediators::PreUpdate()
{
    if(!StandardModel::PreUpdate()) return (false);

    return (true);
}

bool LoopMediators::Update(const std::map<std::string, double>& DPars)
{

    if(!PreUpdate()) return (false);

    UpdateError = false;

    for (std::map<std::string, double>::const_iterator it = DPars.begin(); it != DPars.end(); it++)
        setParameter(it->first, it->second);

    if (UpdateError) return (false);

    if(!PostUpdate()) return (false);

    return (true);
}

bool LoopMediators::PostUpdate()
{
    if(!StandardModel::PostUpdate()) return (false);

    double GammaQ2 = GammaQ*GammaQ;
    double Gammamu2 = Gammamu * Gammamu;
    double mphi2 = mphi * mphi;
    double mpsiQ2 = mpsiQ * mpsiQ;
    double mpsiL2 = mpsiL * mpsiL;
    double yQ = mpsiQ2 / mphi2;
    double yL = mpsiL2 / mphi2;
    double qphi = -charge;
    double qpsi = charge - 1.;

    double M_PI2 = M_PI*M_PI;
    double mmu = getLeptons(MU).getMass();
    double Norm = - sqrt(2.) / (4. * GF * 0.0411494) / (32. * M_PI * ale); // N.B. took the abs of CKM, hence changed overall sign

    C1 = 1. / (128. * M_PI2 * mphi2) * GammaQ2 * (chiBB * etaBB - chiBB_M * etaBB_M) * F9(yQ,yQ);
    C2 = 0.;
    C3 = 0.;
    C4 = 0.;
    C5 = 0.;

    C1p = 0.;
    C2p = 0.;
    C3p = 0.;

    C7 = 0.;
    C8 = 0.;
    C9 = - Norm * GammaQ * Gammamu2 / mphi2 * (chi * eta - chi_M * eta_M) * F9(yQ,yL);
    C10 =  Norm * GammaQ * Gammamu2 / mphi2 * (chi * eta - chi_M * eta_M) * F9(yQ,yL);
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
     * and LoopMediators and LoopMediators-derived parameters in LoopMediatorsMatching */
    LoopMediatorsM.getObj().updateLoopMediatorsParameters();

    return (true);
}

void LoopMediators::setParameter(const std::string name, const double& value)
{
    if(name.compare("GammaQ") == 0)
        GammaQ = value;
    else if(name.compare("Gammamu") == 0)
        Gammamu = value;
    else if(name.compare("mphi") == 0)
        mphi = value;
    else if(name.compare("mpsiQ") == 0)
        mpsiQ = value;
    else if(name.compare("mpsiL") == 0)
        mpsiL = value;
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

bool LoopMediators::CheckParameters(const std::map<std::string, double>& DPars)
{
    for (int i = 0; i < NLoopMediatorsvars; i++) {
        if (DPars.find(LoopMediatorsvars[i]) == DPars.end()) {
            std::cout << "ERROR: missing mandatory LoopMediators parameter " << LoopMediatorsvars[i] << std::endl;
            //raiseMissingModelParameterCount();
            //addMissingModelParameter(LoopMediatorsvars[i]);
        }
    }
    return(StandardModel::CheckParameters(DPars));
}

///////////////////////////////////////////////////////////////////////////
// Flags

bool LoopMediators::setFlag(const std::string name, const bool value)
{
    bool res = false;

    res = StandardModel::setFlag(name,value);

    return(res);
}

///////////////////////////////////////////////////////////////////////////
// Loop Functions

double LoopMediators::F9(double x, double y)
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

double LoopMediators::F7(double x)
{
    double xm1 = x - 1.;

    if (x == 1.)
        return 1./24.;
    else
        return (x*x*x - 6.*x*x + 3.*x + 2. + 6.*x*log(x)) / (12.*xm1*xm1*xm1*xm1);
}

double LoopMediators::F7t(double x)
{
    if (x == 1.)
        return 1./24.;
    else
        return F7(1./x) / x;
}

double LoopMediators::G7(double x)
{
    double xm1 = x - 1.;

    if (x == 1.)
        return 1./12.;
    else
        return (x*x - 4.*x + 3. + 2.*log(x)) / (8.*xm1*xm1*xm1);
}

double LoopMediators::G7t(double x)
{
    double xm1 = x - 1.;

    if (x == 1.)
        return 1./24.;
    else
        return (x*x - 1. - 2.*x*log(x)) / (8.*xm1*xm1*xm1);
}




/////////////////////////////////////////////////////////////////////////



Deltaamu::Deltaamu(const StandardModel& SM_i)
: ThObservable(SM_i), myLM(static_cast<const LoopMediators*> (&SM_i))
{
};

Deltaamu::~Deltaamu()
{
};

double Deltaamu::computeThValue()
{
    return myLM->getDeltaamu();
}



C9mC10::C9mC10(const StandardModel& SM_i)
: ThObservable(SM_i), myLM(static_cast<const LoopMediators*> (&SM_i))
{
};

C9mC10::~C9mC10()
{
};

double C9mC10::computeThValue()
{
    return myLM->getC9();
}
