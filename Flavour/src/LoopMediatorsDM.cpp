/*
 * Copyright (C) 2018 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "LoopMediatorsDM.h"

const std::string LoopMediatorsDM::LoopMediatorsDMvars[NLoopMediatorsDMvars] = {
        "ysybgD", "rVA", "QB", "mB_NP", "mchi_NP", "mV_NP"};

LoopMediatorsDM::LoopMediatorsDM() : StandardModel(), LoopMediatorsDMM(*this)
{

    SMM.setObj((StandardModelMatching&) LoopMediatorsDMM.getObj());
    ModelParamMap.insert(std::make_pair("ysybgD", std::cref(ysybgD)));
    ModelParamMap.insert(std::make_pair("rVA", std::cref(rVA)));
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
    
    double mB2 = mB_NP * mB_NP;
    double mchi2 = mchi_NP * mchi_NP;
    double y = mB2 / mchi2;
    double gmuV = 0.007 * mV_NP;
    double gmuA = rVA * gmuV;
    
    double Norm = - sqrt(2.) / (4. * GF * 0.0411494) / (8. * M_PI * ale); // N.B. took the abs of CKM, hence changed overall sign


    C9 =  Norm * ysybgD * gmuV / mB2 * QB * (F9(y) + G9(y)) ;
    C10 = Norm * ysybgD * gmuA / mB2 * QB * (F9(y) + G9(y)) ;

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
    else if(name.compare("rVA") == 0)
        rVA = value;
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

double LoopMediatorsDM::F9(double x)
{
    double xm14 = pow(x - 1.,4);

    if (x == 1.)
        return -1./24.;
    else
        return (- 2. + 9.*x - 18.*x*x + 11.*x*x*x - 6.*x*x*x*log(x)) / 36. / xm14;
}

double LoopMediatorsDM::G9(double x)
{
    double xm14 = pow(x - 1.,4);

    if (x == 1.)
        return 1./8.;
    else
        return (- 16. + 45.*x - 36.*x*x + 7.*x*x*x + 6.*(3.*x - 2.)*log(x)) / 36. / xm14;
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