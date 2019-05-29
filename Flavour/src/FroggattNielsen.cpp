/* 
 * Copyright (C) 2018 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "FroggattNielsen.h"

const std::string FroggattNielsen::FroggattNielsenvars[NFroggattNielsenvars] = {"logma"};

FroggattNielsen::FroggattNielsen() : StandardModel() {   
    ModelParamMap.insert(std::make_pair("logma", std::cref(ma)));
}

FroggattNielsen::~FroggattNielsen(){
    if (IsModelInitialized()) {
        }
}

///////////////////////////////////////////////////////////////////////////
// Initialization

bool FroggattNielsen::InitializeModel()
{
    setModelInitialized(StandardModel::InitializeModel());
    return(true);
}
    
bool FroggattNielsen::Init(const std::map<std::string, double>& DPars) {
    return(StandardModel::Init(DPars));
}

bool FroggattNielsen::PreUpdate()
{    
    if(!StandardModel::PreUpdate()) return (false);

    return (true);
}

bool FroggattNielsen::Update(const std::map<std::string, double>& DPars) {
    
    if(!PreUpdate()) return (false);
    
    UpdateError = false;

    for (std::map<std::string, double>::const_iterator it = DPars.begin(); it != DPars.end(); it++)
        setParameter(it->first, it->second);

    if (UpdateError) return (false);

    if(!PostUpdate()) return (false);

    return (true);
}

bool FroggattNielsen::PostUpdate()
{
    if(!StandardModel::PostUpdate()) return (false);

    /* Necessary for updating StandardModel parameters in StandardModelMatching,
     * and Axions and Axions-derived parameters in AxionsMatching */
    //LoopMediatorsM.getObj().updateLoopMediatorsParameters();

    return (true);
}

void FroggattNielsen::setParameter(const std::string name, const double& value){    
    if(name.compare("logma") == 0) {
        ma = pow(10., value);
    }
    else
        StandardModel::setParameter(name,value);
}

bool FroggattNielsen::CheckParameters(const std::map<std::string, double>& DPars) {
    for (int i = 0; i < NFroggattNielsenvars; i++) {
        if (DPars.find(FroggattNielsenvars[i]) == DPars.end()) {
            std::cout << "missing mandatory Axions parameter " << FroggattNielsenvars[i] << std::endl;
            return false;
        }
    }
    return(StandardModel::CheckParameters(DPars));
}

///////////////////////////////////////////////////////////////////////////
// Flags

bool FroggattNielsen::setFlag(const std::string name, const bool value)
{
    bool res = false;
    
    res = StandardModel::setFlag(name,value);

    return(res);
}


///////////////////////////////////////////////////////////////////////////
// Observables


test::test(const StandardModel& SM_i)
: ThObservable(SM_i), myFroggattNielsen(static_cast<const FroggattNielsen*> (&SM_i))
{
};

test::~test()
{
};

double test::computeThValue()
{
    return myFroggattNielsen->getma();

}