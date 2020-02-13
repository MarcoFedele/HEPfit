/* 
 * Copyright (C) 2018 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "Axions.h"

const std::string Axions::Axionsvars[NAxionsvars] = {"gae", "gag",
                                "a_G117B15A", "b_G117B15A", "c_G117B15A", "d_G117B15A", 
                                "a_R548", "b_R548", "c_R548", "d_R548", 
                                "a_PG1351489", "b_PG1351489", "c_PG1351489", "d_PG1351489", 
                                "a_L192", "b_L192", "c_L192", "d_L192", 
                                "a_TRGB", "b_TRGB",
                                "Y_HBR"};

Axions::Axions() : StandardModel() {   

    ModelParamMap.insert(std::make_pair("gae", std::cref(gae)));
    ModelParamMap.insert(std::make_pair("gag", std::cref(gag)));
    
    ModelParamMap.insert(std::make_pair("a_G117B15A", std::cref(a_G117B15A)));
    ModelParamMap.insert(std::make_pair("b_G117B15A", std::cref(b_G117B15A)));
    ModelParamMap.insert(std::make_pair("c_G117B15A", std::cref(c_G117B15A)));
    ModelParamMap.insert(std::make_pair("d_G117B15A", std::cref(d_G117B15A)));
    
    ModelParamMap.insert(std::make_pair("a_R548", std::cref(a_R548)));
    ModelParamMap.insert(std::make_pair("b_R548", std::cref(b_R548)));
    ModelParamMap.insert(std::make_pair("c_R548", std::cref(c_R548)));
    ModelParamMap.insert(std::make_pair("d_R548", std::cref(d_R548)));
    
    ModelParamMap.insert(std::make_pair("a_PG1351489", std::cref(a_PG1351489)));
    ModelParamMap.insert(std::make_pair("b_PG1351489", std::cref(b_PG1351489)));
    ModelParamMap.insert(std::make_pair("c_PG1351489", std::cref(c_PG1351489)));
    ModelParamMap.insert(std::make_pair("d_PG1351489", std::cref(d_PG1351489)));
    
    ModelParamMap.insert(std::make_pair("a_L192", std::cref(a_L192)));
    ModelParamMap.insert(std::make_pair("b_L192", std::cref(b_L192)));
    ModelParamMap.insert(std::make_pair("c_L192", std::cref(c_L192)));
    ModelParamMap.insert(std::make_pair("d_L192", std::cref(d_L192)));
    
    ModelParamMap.insert(std::make_pair("a_TRGB", std::cref(a_TRGB)));
    ModelParamMap.insert(std::make_pair("b_TRGB", std::cref(b_TRGB)));
    
    ModelParamMap.insert(std::make_pair("Y_HBR", std::cref(Y_HBR)));
}

Axions::~Axions(){
    if (IsModelInitialized()) {
        }
}

///////////////////////////////////////////////////////////////////////////
// Initialization

bool Axions::InitializeModel()
{
    setModelInitialized(StandardModel::InitializeModel());
    return(true);
}
    
bool Axions::Init(const std::map<std::string, double>& DPars) {
    return(StandardModel::Init(DPars));
}

bool Axions::PreUpdate()
{    
    if(!StandardModel::PreUpdate()) return (false);

    return (true);
}

bool Axions::Update(const std::map<std::string, double>& DPars) {
    
    if(!PreUpdate()) return (false);
    
    UpdateError = false;

    for (std::map<std::string, double>::const_iterator it = DPars.begin(); it != DPars.end(); it++)
        setParameter(it->first, it->second);

    if (UpdateError) return (false);

    if(!PostUpdate()) return (false);

    return (true);
}

bool Axions::PostUpdate()
{
    if(!StandardModel::PostUpdate()) return (false);

    /* Necessary for updating StandardModel parameters in StandardModelMatching,
     * and Axions and Axions-derived parameters in AxionsMatching */
    //LoopMediatorsM.getObj().updateLoopMediatorsParameters();

    return (true);
}

void Axions::setParameter(const std::string name, const double& value){    
    if(name.compare("gae") == 0) 
        gae = value;  
    else if(name.compare("gag") == 0) 
        gag = value;
    else if(name.compare("a_G117B15A") == 0)
        a_G117B15A = value;
    else if(name.compare("b_G117B15A") == 0)
        b_G117B15A = value;
    else if(name.compare("c_G117B15A") == 0)
        c_G117B15A = value;
    else if(name.compare("d_G117B15A") == 0)
        d_G117B15A = value;
    else if(name.compare("a_R548") == 0)
        a_R548 = value;
    else if(name.compare("b_R548") == 0)
        b_R548 = value;
    else if(name.compare("c_R548") == 0)
        c_R548 = value;
    else if(name.compare("d_R548") == 0)
        d_R548 = value;
    else if(name.compare("a_PG1351489") == 0)
        a_PG1351489 = value;
    else if(name.compare("b_PG1351489") == 0)
        b_PG1351489 = value;
    else if(name.compare("c_PG1351489") == 0)
        c_PG1351489 = value;
    else if(name.compare("d_PG1351489") == 0)
        d_PG1351489 = value;
    else if(name.compare("a_L192") == 0)
        a_L192 = value;
    else if(name.compare("b_L192") == 0)
        b_L192 = value;
    else if(name.compare("c_L192") == 0)
        c_L192 = value;
    else if(name.compare("d_L192") == 0)
        d_L192 = value;
    else if(name.compare("a_TRGB") == 0)
        a_TRGB = value;
    else if(name.compare("b_TRGB") == 0)
        b_TRGB = value;
    else if(name.compare("Y_HBR") == 0)
        Y_HBR = value;
    else
        StandardModel::setParameter(name,value);
}

bool Axions::CheckParameters(const std::map<std::string, double>& DPars) {
    for (int i = 0; i < NAxionsvars; i++) {
        if (DPars.find(Axionsvars[i]) == DPars.end()) {
            std::cout << "missing mandatory Axions parameter " << Axionsvars[i] << std::endl;
            return false;
        }
    }
    return(StandardModel::CheckParameters(DPars));
}

///////////////////////////////////////////////////////////////////////////
// Flags

bool Axions::setFlag(const std::string name, const bool value)
{
    bool res = false;
    
    res = StandardModel::setFlag(name,value);

    return(res);
}



///////////////////////////////////////////////////////////////////////////
// Observables


mac2::mac2(const StandardModel& SM_i)
: ThObservable(SM_i), myAxions(static_cast<const Axions*> (&SM_i))
{
};

mac2::~mac2()
{
};

double mac2::computeThValue()
{
    double gae=myAxions->getgae();
    
    return gae/0.28;

}


G117B15A::G117B15A(const StandardModel& SM_i)
: ThObservable(SM_i), myAxions(static_cast<const Axions*> (&SM_i))
{
};

G117B15A::~G117B15A()
{
};

double G117B15A::computeThValue()
{
    double a=myAxions->geta_G117B15A();
    double b=myAxions->getb_G117B15A();
    double c=myAxions->getc_G117B15A();
    double d=myAxions->getd_G117B15A();
    
    double mac2=myAxions->getgae()/0.28;
    
    return a + mac2*b + mac2*mac2*c + mac2*mac2*mac2*d;

}


R548::R548(const StandardModel& SM_i)
: ThObservable(SM_i), myAxions(static_cast<const Axions*> (&SM_i))
{
};

R548::~R548()
{
};

double R548::computeThValue()
{
    double a=myAxions->geta_R548();
    double b=myAxions->getb_R548();
    double c=myAxions->getc_R548();
    double d=myAxions->getd_R548();
    
    double mac2=myAxions->getgae()/0.28;
    
    return a + mac2*b + mac2*mac2*c + mac2*mac2*mac2*d;

}


PG1351489::PG1351489(const StandardModel& SM_i)
: ThObservable(SM_i), myAxions(static_cast<const Axions*> (&SM_i))
{
};

PG1351489::~PG1351489()
{
};

double PG1351489::computeThValue()
{
    double a=myAxions->geta_PG1351489();
    double b=myAxions->getb_PG1351489();
    double c=myAxions->getc_PG1351489();
    double d=myAxions->getd_PG1351489();
    
    double mac2=myAxions->getgae()/0.28;
    
    return a + mac2*b + mac2*mac2*c + mac2*mac2*mac2*d;

}


L192::L192(const StandardModel& SM_i)
: ThObservable(SM_i), myAxions(static_cast<const Axions*> (&SM_i))
{
};

L192::~L192()
{
};

double L192::computeThValue()
{
    double a=myAxions->geta_L192();
    double b=myAxions->getb_L192();
    double c=myAxions->getc_L192();
    double d=myAxions->getd_L192();
    
    double mac2=myAxions->getgae()/0.28;
    
    return a + mac2*b + mac2*mac2*c + mac2*mac2*mac2*d;

}


TRGB::TRGB(const StandardModel& SM_i)
: ThObservable(SM_i), myAxions(static_cast<const Axions*> (&SM_i))
{
};

TRGB::~TRGB()
{
};

double TRGB::computeThValue()
{
    double a=myAxions->geta_TRGB();
    double b=myAxions->getb_TRGB();
    
    double gae=myAxions->getgae();
    
    return - 4.03 - 0.25*(sqrt(gae*gae + 0.93) - 0.96 - 0.17*pow(gae,1.5)) 
            + 0.039 + a + gae*b;

}


HBR::HBR(const StandardModel& SM_i)
: ThObservable(SM_i), myAxions(static_cast<const Axions*> (&SM_i))
{
};

HBR::~HBR()
{
};

double HBR::computeThValue()
{
    double Y=myAxions->getY_HBR();
    
    gslpp::complex gag=myAxions->getgag();
    double gae=myAxions->getgae();
    
    double alpha = gae*gae/4./M_PI;
    double dMc = 0.024*(sqrt(gae*gae + 1.23*1.23) - 1.23 -0.921*pow(alpha,0.75));
    
    return 7.33*Y + 0.02 - 0.095*sqrt(21.86 + 21.08*gag.abs()) - 1.61*dMc - 0.067*alpha;  // 1512.08108, Eq. (7.6)
    
    return 6.26*Y - 0.12 - 0.14*gag.abs2() - 1.61*dMc - 0.067*alpha;  // 1512.08108, Eq. (7.5)
    
    return 6.26*Y - 0.41*gag.abs2() - 0.12;  // 1406.6053, Eq. (1)

}