/*
 * Copyright (C) 2018 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "LoopMediatorsDM.h"

const std::string LoopMediatorsDM::LoopMediatorsDMvars[NLoopMediatorsDMvars] = {
        "ysybgD", "gmuV", "gmuA", "QB", "mB_NP", "mchi_NP", "mV_NP"};

LoopMediatorsDM::LoopMediatorsDM() : StandardModel(), LoopMediatorsDMM(*this)
{

    SMM.setObj((StandardModelMatching&) LoopMediatorsDMM.getObj());
    ModelParamMap.insert(std::make_pair("ysybgD", std::cref(ysybgD)));
    ModelParamMap.insert(std::make_pair("gmuV", std::cref(gmuV)));
    ModelParamMap.insert(std::make_pair("gmuA", std::cref(gmuA)));
    ModelParamMap.insert(std::make_pair("QB", std::cref(QB)));
    ModelParamMap.insert(std::make_pair("mB_NP", std::cref(mB_NP)));
    ModelParamMap.insert(std::make_pair("mchi_NP", std::cref(mchi_NP)));
    ModelParamMap.insert(std::make_pair("mV_NP", std::cref(mV_NP)));
}

LoopMediatorsDM::~LoopMediatorsDM()
{
    if (IsModelInitialized()) {
        }
}

///////////////////////////////////////////////////////////////////////////
// Initialization

bool LoopMediatorsDM::InitializeModel()
{
    setModelInitialized(StandardModel::InitializeModel());
    return(true);
}

bool LoopMediatorsDM::Init(const std::map<std::string, double>& DPars)
{
    return(StandardModel::Init(DPars));
}

bool LoopMediatorsDM::PreUpdate()
{
    if(!StandardModel::PreUpdate()) return (false);

    return (true);
}

bool LoopMediatorsDM::Update(const std::map<std::string, double>& DPars)
{

    if(!PreUpdate()) return (false);

    UpdateError = false;

    for (std::map<std::string, double>::const_iterator it = DPars.begin(); it != DPars.end(); it++)
        setParameter(it->first, it->second);

    if (UpdateError) return (false);

    if(!PostUpdate()) return (false);

    return (true);
}

bool LoopMediatorsDM::PostUpdate()
{
    if(!StandardModel::PostUpdate()) return (false);

    C9 =  0.;
    C10 = 0.;

    /*std::cout << "C1 " << C1 << std::endl;
    std::cout << "C9 " << C9.abs() << std::endl;*/

    /*Deltaamu = mmu*mmu / ( 8. * M_PI2 * mphi2) * (
            (GammamuL2 + GammamuR2) * (qphi*F7t(yE) - qpsi*F7(yE))
            + 4./sqrt(2.) * v() * lambdaE / mmu * (2.*GammamuL*GammamuR) * (qphi*G7t(yE) - qpsi*G7(yE))
          );*/

    /* Necessary for updating StandardModel parameters in StandardModelMatching,
     * and LoopMediatorsDM and LoopMediatorsDM-derived parameters in LoopMediatorsDMMatching */
    LoopMediatorsDMM.getObj().updateLoopMediatorsDMParameters();

    return (true);
}

void LoopMediatorsDM::setParameter(const std::string name, const double& value) 
{
    if(name.compare("ysybgD") == 0)
        ysybgD = value;
    else if(name.compare("gmuV") == 0)
        gmuV = value;
    else if(name.compare("gmuA") == 0)
        gmuA = value;
    else if(name.compare("QB") == 0)
        QB = value;
    else if(name.compare("mB_NP") == 0)
        mB_NP = value;
    else if(name.compare("mchi_NP") == 0)
        mchi_NP = value;
    else if(name.compare("mV_NP") == 0)
        mV_NP = value;
    else
        StandardModel::setParameter(name,value);
}

bool LoopMediatorsDM::CheckParameters(const std::map<std::string, double>& DPars)
{
    for (int i = 0; i < NLoopMediatorsDMvars; i++) {
        if (DPars.find(LoopMediatorsDMvars[i]) == DPars.end()) {
            std::cout << "ERROR: missing mandatory LoopMediatorsDM parameter " << LoopMediatorsDMvars[i] << std::endl;
            //raiseMissingModelParameterCount();
            //addMissingModelParameter(LoopMediatorsDMvars[i]);
        }
    }
    return(StandardModel::CheckParameters(DPars));
}

///////////////////////////////////////////////////////////////////////////
// Flags

bool LoopMediatorsDM::setFlag(const std::string name, const bool value)
{
    bool res = false;

    res = StandardModel::setFlag(name,value);

    return(res);
}

///////////////////////////////////////////////////////////////////////////
// Loop Functions

double LoopMediatorsDM::F9(double x, double y)
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

double LoopMediatorsDM::G9(double x, double y)
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




/////////////////////////////////////////////////////////////////////////



C9_LMDM::C9_LMDM(const StandardModel& SM_i)
: ThObservable(SM_i), myLM(static_cast<const LoopMediatorsDM*> (&SM_i))
{
};

C9_LMDM::~C9_LMDM()
{
};

double C9_LMDM::computeThValue()
{
    return myLM->getC9();
}



C10_LMDM::C10_LMDM(const StandardModel& SM_i)
: ThObservable(SM_i), myLM(static_cast<const LoopMediatorsDM*> (&SM_i))
{
};

C10_LMDM::~C10_LMDM()
{
};

double C10_LMDM::computeThValue()
{
    return myLM->getC10();
}