/* 
 * Copyright (C) 2018 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "Axions.h"

const std::string Axions::Axionsvars[NAxionsvars] = {"gae", "gag",
                                "a_G117B15A", "b_G117B15A",
                                "a_R548", "b_R548",
                                "a_PG1351489", "b_PG1351489",
                                "a_L113", "b_L113", 
                                "a_L192", "b_L192",
                                "a_TRGB", "b_TRGB",
                                "Y_HBR"};

Axions::Axions() : StandardModel() {   

    ModelParamMap.insert(std::make_pair("gae", std::cref(gae)));
    ModelParamMap.insert(std::make_pair("gag", std::cref(gag)));
    
    ModelParamMap.insert(std::make_pair("a_G117B15A", std::cref(a_G117B15A)));
    ModelParamMap.insert(std::make_pair("b_G117B15A", std::cref(b_G117B15A)));

    ModelParamMap.insert(std::make_pair("a_R548", std::cref(a_R548)));
    ModelParamMap.insert(std::make_pair("b_R548", std::cref(b_R548)));

    ModelParamMap.insert(std::make_pair("a_PG1351489", std::cref(a_PG1351489)));
    ModelParamMap.insert(std::make_pair("b_PG1351489", std::cref(b_PG1351489)));

    ModelParamMap.insert(std::make_pair("a_L113", std::cref(a_L113)));
    ModelParamMap.insert(std::make_pair("b_L113", std::cref(b_L113)));

    ModelParamMap.insert(std::make_pair("a_L192", std::cref(a_L192)));
    ModelParamMap.insert(std::make_pair("b_L192", std::cref(b_L192)));
    
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
        gag = value; //pow(10., value);
    else if(name.compare("a_G117B15A") == 0)
        a_G117B15A = value;
    else if(name.compare("b_G117B15A") == 0)
        b_G117B15A = value;
    else if(name.compare("a_R548") == 0)
        a_R548 = value;
    else if(name.compare("b_R548") == 0)
        b_R548 = value;
    else if(name.compare("a_PG1351489") == 0)
        a_PG1351489 = value;
    else if(name.compare("b_PG1351489") == 0)
        b_PG1351489 = value;
    else if(name.compare("a_L113") == 0)
        a_L113 = value;
    else if(name.compare("b_L113") == 0)
        b_L113 = value;
    else if(name.compare("a_L192") == 0)
        a_L192 = value;
    else if(name.compare("b_L192") == 0)
        b_L192 = value;
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



loggag::loggag(const StandardModel& SM_i)
: ThObservable(SM_i), myAxions(static_cast<const Axions*> (&SM_i))
{
};

loggag::~loggag()
{
};

double loggag::computeThValue()
{
    return log10(myAxions->getgag());

}


loggae::loggae(const StandardModel& SM_i)
: ThObservable(SM_i), myAxions(static_cast<const Axions*> (&SM_i))
{
};

loggae::~loggae()
{
};

double loggae::computeThValue()
{
    return log10(myAxions->getgae());

}


mac2::mac2(const StandardModel& SM_i)
: ThObservable(SM_i), myAxions(static_cast<const Axions*> (&SM_i))
{
};

mac2::~mac2()
{
};

double mac2::computeThValue()
{
    return myAxions->getgae()/0.28e-13;

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

    double a26=myAxions->getgae()*myAxions->getgae()/(4.*M_PI)/1.e-26;

    return a + a26*b;

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

    double a26=myAxions->getgae()*myAxions->getgae()/(4.*M_PI)/1.e-26;

    return a + a26*b;

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

    double a26=myAxions->getgae()*myAxions->getgae()/(4.*M_PI)/1.e-26;

    return a + a26*b;

}


L113::L113(const StandardModel& SM_i)
: ThObservable(SM_i), myAxions(static_cast<const Axions*> (&SM_i))
{
};

L113::~L113()
{
};

double L113::computeThValue()
{
    double a=myAxions->geta_L113();
    double b=myAxions->getb_L113();

    double a26=myAxions->getgae()*myAxions->getgae()/(4.*M_PI)/1.e-26;

    return a + a26*b;

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

    double a26=myAxions->getgae()*myAxions->getgae()/(4.*M_PI)/1.e-26;

    return a + a26*b;

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

    double gae=myAxions->getgae()/1.e-13;

    return - 4.08 - 0.25*(sqrt(gae*gae + 0.93) - 0.96 - 0.17*pow(gae,1.5))
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
    double gae=myAxions->getgae()/1.e-13;

    double alpha = gae*gae/4./M_PI;
    double g10 = gag.abs()*1.e10;

    double a = 6.26*Y - 0.12;
    double b = 0.41;

    double dMc = 0.024*(sqrt(gae*gae + 1.23*1.23) - 1.23 - 0.921*pow(alpha,0.75));
    
    // dMc = 0.024*(sqrt(9.*9. + 1.23*1.23) - 1.23 - 0.921*pow(9.*9./4./M_PI,0.75));  // FISSATO GAE=9 QUI!!!!!
    
    // return a / (1 + g10*g10*b/a) ;  // Versione per gag grande
    
    return 7.33*Y + 0.02 - 0.095*sqrt(21.86 + 21.08*g10) - 1.61*dMc - 0.067*alpha;  // 1512.08108, Eq. (7.6)

    return 6.26*Y - 0.12 - 0.14*g10*g10 - 1.61*dMc - 0.067*alpha;  // 1512.08108, Eq. (7.5)

    return 6.26*Y - 0.41*g10*g10 - 0.12;  // 1406.6053, Eq. (1)

}


Xenon::Xenon(const StandardModel& SM_i)
: ThObservable(SM_i), myAxions(static_cast<const Axions*> (&SM_i))
{
};

Xenon::~Xenon()
{
};

double Xenon::computeThValue()
{
    double gg10=myAxions->getgag()/1.e-10;
    double ge12=myAxions->getgae()/1.e-12;
    
    return pow( ge12*ge12 * (ge12*ge12 + 2.*gg10*gg10) ,0.25);  // 2003.01100, Eq. (247)

}


CaeoCag::CaeoCag(const StandardModel& SM_i)
: ThObservable(SM_i), myAxions(static_cast<const Axions*> (&SM_i))
{
};

CaeoCag::~CaeoCag()
{
};

double CaeoCag::computeThValue()
{
    double gag=myAxions->getgag();
    double gae=myAxions->getgae();
    
    double res = 7.2973525698e-3/(2 * M_PI * 5.109989e-4) * ( gae / gag);
    
    if (res < 0. || res > 1.) return 0.;
    else return res;

}