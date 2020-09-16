/*
 * Copyright (C) 2018 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

/*
 * Models:
 * 
 * -1: SM
 * 0: KSVZ
 * 1: DFSZ1
 * 2: DFSZ2
 * 3: M1
 * 4: M2
 * 5: M3
 * 6: M4
 * 7: 3HDM
 * 8: T2d
 * 9: T2u
 * 11: DFSZ1 general
 * 12: DFSZ2 general
 *
 */

#include "AxionsTHDM.h"

const std::string AxionsTHDM::Axionsvars[NAxionsvars] = {"logtanb", "logma", "model", "eps_L", "EoN", "logcgamma", "Chi3",
                                "C_pn_err_0", "C_pn_err_1", "C_pn_err_2",
                                "a_G117B15A", "b_G117B15A",
                                "a_R548", "b_R548",
                                "a_PG1351489", "b_PG1351489",
                                "a_L113", "b_L113", 
                                "a_L192", "b_L192",
                                "a_TRGB", "b_TRGB",
                                "Y_HBR"};

AxionsTHDM::AxionsTHDM() : StandardModel() {
    ModelParamMap.insert(std::make_pair("logtanb", std::cref(tanb)));
    ModelParamMap.insert(std::make_pair("logma", std::cref(ma)));
    ModelParamMap.insert(std::make_pair("model", std::cref(model)));
    ModelParamMap.insert(std::make_pair("eps_L", std::cref(eps_L)));
    ModelParamMap.insert(std::make_pair("EoN", std::cref(EoN)));
    ModelParamMap.insert(std::make_pair("logcgamma", std::cref(cgamma)));
    ModelParamMap.insert(std::make_pair("Chi3", std::cref(Chi3)));
    ModelParamMap.insert(std::make_pair("C_pn_err_0", std::cref(C_pn_err_0)));
    ModelParamMap.insert(std::make_pair("C_pn_err_1", std::cref(C_pn_err_1)));
    ModelParamMap.insert(std::make_pair("C_pn_err_2", std::cref(C_pn_err_2)));


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

AxionsTHDM::~AxionsTHDM(){
    if (IsModelInitialized()) {
        }
}

///////////////////////////////////////////////////////////////////////////
// Initialization

bool AxionsTHDM::InitializeModel()
{
    setModelInitialized(StandardModel::InitializeModel());
    return(true);
}

bool AxionsTHDM::Init(const std::map<std::string, double>& DPars) {
    return(StandardModel::Init(DPars));
}

bool AxionsTHDM::PreUpdate()
{
    if(!StandardModel::PreUpdate()) return (false);

    return (true);
}

bool AxionsTHDM::Update(const std::map<std::string, double>& DPars) {

    if(!PreUpdate()) return (false);

    UpdateError = false;

    for (std::map<std::string, double>::const_iterator it = DPars.begin(); it != DPars.end(); it++)
        setParameter(it->first, it->second);

    if (UpdateError) return (false);

    if(!PostUpdate()) return (false);

    return (true);
}

bool AxionsTHDM::PostUpdate()
{
    if(!StandardModel::PostUpdate()) return (false);

    /* Necessary for updating StandardModel parameters in StandardModelMatching,
     * and Axions and Axions-derived parameters in AxionsMatching */
    //LoopMediatorsM.getObj().updateLoopMediatorsParameters();

    return (true);
}

void AxionsTHDM::setParameter(const std::string name, const double& value){
    if(name.compare("logtanb") == 0) {
        tanb = pow(10., value);
        sinb = tanb / sqrt(1. + tanb*tanb);
        cosb = 1. / sqrt(1. + tanb*tanb);
        }
    else if(name.compare("logma") == 0) {
        ma = pow(10., value);
    }
    else if(name.compare("model") == 0)
        model = value;
    else if(name.compare("eps_L") == 0)
        eps_L = value;
    else if(name.compare("EoN") == 0)
        EoN = value;
    else if(name.compare("logcgamma") == 0)
        cgamma = pow(10., value);
    else if(name.compare("Chi3") == 0)
        Chi3 = value;
    else if(name.compare("C_pn_err_0") == 0)
        C_pn_err_0 = value;
    else if(name.compare("C_pn_err_1") == 0)
        C_pn_err_1 = value;
    else if(name.compare("C_pn_err_2") == 0)
        C_pn_err_2 = value;
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

bool AxionsTHDM::CheckParameters(const std::map<std::string, double>& DPars) {
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

bool AxionsTHDM::setFlag(const std::string name, const bool value)
{
    bool res = false;

    res = StandardModel::setFlag(name,value);

    return(res);
}


double AxionsTHDM::gag() const
{
    double Cag;

    if (model == - 1.)
        return 0.;
    else if (model == 0.)
        Cag = std::abs(EoN - 1.92);
    else if (model == 1.)
        Cag = 8./3. - 1.92;
    else if (model == 2.)
        Cag = -2./3. + 1.92;
    else if (model == 3.)
        Cag = -2./3. + 1.92;
    else if (model == 4.)
        Cag = 8./3. - 1.92;
    else if (model == 5.)
        Cag = 4./3. + 1.92;
    else if (model == 6.)
        Cag = 14./3. - 1.92;
    else if (model == 7.)
        Cag = 8./3. - 1.92;
    else if (model == 8.)
        Cag = 20./3. - 1.92;
    else if (model == 9.)
        Cag = 10./3. + 1.92;
    else if ((model == 11.) || (model == 12.))
        Cag = cgamma;
    else
        throw std::runtime_error("error in AzionsTHDM::gag, model can only be an integer between -1. and 9. !");

    return getAle()/2./M_PI * ma/(5.7e9) * Cag;
}

double AxionsTHDM::Cae() const
{
    double Cae;

    if (model == - 1.)
        return 0.;
    else if (model == 0.)
    {
        double ale = 1./137.;
        double me = 5.109989e-4;
        Cae = 3.*ale*ale/4/M_PI/M_PI * std::abs( EoN * log(2.e11/me) - 1.92 * log(1./me) );
    }
    else if ((model == 1.) || (model == 11.))
        Cae = sinb*sinb/3.;
    else if ((model == 2.) || (model == 12.))
        Cae = -cosb*cosb/3.;
    else if (model == 3.)
        Cae = -sinb*sinb + eps_L;
    else if (model == 4.)
        Cae = -sinb*sinb + eps_L;
    else if (model == 5.)
        Cae = sinb*sinb + eps_L;
    else if (model == 6.)
        Cae = -sinb*sinb + eps_L;
    else if (model == 7.)
        Cae = Chi3/3.;
    else if (model == 8.)
        Cae = cosb*cosb;
    else if (model == 9.)
        Cae = -sinb*sinb;
    else
        throw std::runtime_error("error in AzionsTHDM::gae, model can only be an integer between -1. and 9. !");

    return Cae;
}

double AxionsTHDM::gae() const
{
    return 0.84e-13 * ma * Cae();
}

double AxionsTHDM::gap() const
{
    double Cad, Cas, Cab;
    double Cau, Cac, Cat;

    if (model == -1.) {
        return 0.;
    }
    else if (model == 0.) {
        Cad = 0.;
        Cas = 0.;
        Cab = 0.;
        Cau = 0.;
        Cac = 0.;
        Cat = 0.;
    }
    else if ((model == 1.) || (model == 11.)){
        Cad = sinb*sinb/3.;
        Cas = sinb*sinb/3.;
        Cab = sinb*sinb/3.;
        Cau = cosb*cosb/3.;
        Cac = cosb*cosb/3.;
        Cat = cosb*cosb/3.;
    }
    else if ((model == 2.) || (model == 12.)){
        Cad = sinb*sinb/3.;
        Cas = sinb*sinb/3.;
        Cab = sinb*sinb/3.;
        Cau = cosb*cosb/3.;
        Cac = cosb*cosb/3.;
        Cat = cosb*cosb/3.;
    }
    else if (model == 3.){
        Cad = cosb*cosb;
        Cas = cosb*cosb;
        Cab = -sinb*sinb;
        Cau = sinb*sinb;
        Cac = sinb*sinb;
        Cat = -cosb*cosb;
    }
    else if (model == 4.){
        Cad = cosb*cosb;
        Cas = -sinb*sinb;
        Cab = -sinb*sinb;
        Cau = sinb*sinb;
        Cac = sinb*sinb;
        Cat = sinb*sinb;
    }
    else if (model == 5.){
        Cad = sinb*sinb;
        Cas = sinb*sinb;
        Cab = sinb*sinb;
        Cau = cosb*cosb;
        Cac = -sinb*sinb;
        Cat = -sinb*sinb;
    }
    else if (model == 6.){
        Cad = cosb*cosb;
        Cas = -sinb*sinb;
        Cab = -sinb*sinb;
        Cau = sinb*sinb;
        Cac = sinb*sinb;
        Cat = sinb*sinb;
    }
    else if (model == 7.){
        Cad = 1./3. + Chi3/3.;
        Cas = 1./3. + Chi3/3.;
        Cab = -2./3. + Chi3/3.;
        Cau = 2./3. - Chi3/3.;
        Cac = 2./3. - Chi3/3.;
        Cat = -1./3. - Chi3/3.;
    }
    else if (model == 8.){
        Cad = cosb*cosb;
        Cas = -sinb*sinb;
        Cab = -sinb*sinb;
        Cau = 0.1 + sinb*sinb;
        Cac = -0.09 + sinb*sinb;
        Cat = -0.01 + sinb*sinb;
    }
    else if (model == 9.){
        Cad = 0.1 + cosb*cosb;
        Cas = 0.55 + cosb*cosb;
        Cab = -0.66 + cosb*cosb;
        Cau = sinb*sinb;
        Cac = -0.0036 - cosb*cosb;
        Cat = 0.0036 - cosb*cosb;
    }
    else
        throw std::runtime_error("error in AzionsTHDM::gap, model can only be an integer between -1. and 9. !");

    double Cap = (-0.47+C_pn_err_0)
            + (0.88+C_pn_err_1)*Cau + (-0.39+C_pn_err_2)*Cad - 0.038*Cas - 0.012*Cac - 0.009*Cab - 0.0035*Cat;

    //Cap = - 0.435*sinb*sinb - 0.182;

    return 1.56e-10 * ma * Cap;
}

double AxionsTHDM::gan() const
{
    double Cad, Cas, Cab;
    double Cau, Cac, Cat;

    if (model == -1.) {
        return 0.;
    }
    else if (model == 0.) {
        Cad = 0.;
        Cas = 0.;
        Cab = 0.;
        Cau = 0.;
        Cac = 0.;
        Cat = 0.;
    }
    else if ((model == 1.) || (model == 11.)){
        Cad = sinb*sinb/3.;
        Cas = sinb*sinb/3.;
        Cab = sinb*sinb/3.;
        Cau = cosb*cosb/3.;
        Cac = cosb*cosb/3.;
        Cat = cosb*cosb/3.;
    }
    else if ((model == 2.) || (model == 12.)){
        Cad = sinb*sinb/3.;
        Cas = sinb*sinb/3.;
        Cab = sinb*sinb/3.;
        Cau = cosb*cosb/3.;
        Cac = cosb*cosb/3.;
        Cat = cosb*cosb/3.;
    }
    else if (model == 3.){
        Cad = cosb*cosb;
        Cas = cosb*cosb;
        Cab = -sinb*sinb;
        Cau = sinb*sinb;
        Cac = sinb*sinb;
        Cat = -cosb*cosb;
    }
    else if (model == 4.){
        Cad = cosb*cosb;
        Cas = -sinb*sinb;
        Cab = -sinb*sinb;
        Cau = sinb*sinb;
        Cac = sinb*sinb;
        Cat = sinb*sinb;
    }
    else if (model == 5.){
        Cad = sinb*sinb;
        Cas = sinb*sinb;
        Cab = sinb*sinb;
        Cau = cosb*cosb;
        Cac = -sinb*sinb;
        Cat = -sinb*sinb;
    }
    else if (model == 6.){
        Cad = cosb*cosb;
        Cas = -sinb*sinb;
        Cab = -sinb*sinb;
        Cau = sinb*sinb;
        Cac = sinb*sinb;
        Cat = sinb*sinb;
    }
    else if (model == 7.){
        Cad = 1./3. + Chi3/3.;
        Cas = 1./3. + Chi3/3.;
        Cab = -2./3. + Chi3/3.;
        Cau = 2./3. - Chi3/3.;
        Cac = 2./3. - Chi3/3.;
        Cat = -1./3. - Chi3/3.;
    }
    else if (model == 8.){
        Cad = cosb*cosb;
        Cas = -sinb*sinb;
        Cab = -sinb*sinb;
        Cau = 0.1 + sinb*sinb;
        Cac = -0.09 + sinb*sinb;
        Cat = -0.01 + sinb*sinb;
    }
    else if (model == 9.){
        Cad = 0.1 + cosb*cosb;
        Cas = 0.55 + cosb*cosb;
        Cab = -0.66 + cosb*cosb;
        Cau = sinb*sinb;
        Cac = -0.0036 - cosb*cosb;
        Cat = 0.0036 - cosb*cosb;
    }
    else
        throw std::runtime_error("error in AzionsTHDM::gap, model can only be an integer between -1. and 9. !");

    double Can = (-0.02+C_pn_err_0)
            + (0.88+C_pn_err_1)*Cad + (-0.39+C_pn_err_2)*Cau - 0.038*Cas - 0.012*Cac - 0.009*Cab - 0.0035*Cat;

    //Can = 0.414*sinb*sinb - 0.160;

    return 1.57e-10 * ma * Can;
}


///////////////////////////////////////////////////////////////////////////
// Observables


gagTHDM::gagTHDM(const StandardModel& SM_i)
: ThObservable(SM_i), myAxions(static_cast<const AxionsTHDM*> (&SM_i))
{
};

gagTHDM::~gagTHDM()
{
};

double gagTHDM::computeThValue()
{
    return myAxions->gag();

}


gaeTHDM::gaeTHDM(const StandardModel& SM_i)
: ThObservable(SM_i), myAxions(static_cast<const AxionsTHDM*> (&SM_i))
{
};

gaeTHDM::~gaeTHDM()
{
};

double gaeTHDM::computeThValue()
{
    return myAxions->gae();

}


logtbTHDM::logtbTHDM(const StandardModel& SM_i)
: ThObservable(SM_i), myAxions(static_cast<const AxionsTHDM*> (&SM_i))
{
};

logtbTHDM::~logtbTHDM()
{
};

double logtbTHDM::computeThValue()
{
    if (myAxions->getmodel() == -1.) {
        return 0.;
    }
    else return log10(myAxions->gettanb());

}


logmaTHDM::logmaTHDM(const StandardModel& SM_i)
: ThObservable(SM_i), myAxions(static_cast<const AxionsTHDM*> (&SM_i))
{
};

logmaTHDM::~logmaTHDM()
{
};

double logmaTHDM::computeThValue()
{
    if (myAxions->getmodel() == -1.) {
        return 0.;
    }
    else return log10(myAxions->getma());

}


loggagTHDM::loggagTHDM(const StandardModel& SM_i)
: ThObservable(SM_i), myAxions(static_cast<const AxionsTHDM*> (&SM_i))
{
};

loggagTHDM::~loggagTHDM()
{
};

double loggagTHDM::computeThValue()
{
    if (myAxions->getmodel() == -1.) {
        return 0.;
    }
    else return log10(myAxions->gag());

}


loggaeTHDM::loggaeTHDM(const StandardModel& SM_i)
: ThObservable(SM_i), myAxions(static_cast<const AxionsTHDM*> (&SM_i))
{
};

loggaeTHDM::~loggaeTHDM()
{
};

double loggaeTHDM::computeThValue()
{
    if (myAxions->getmodel() == -1.) {
        return 0.;
    }
    else return log10(myAxions->gae());

}


mac2THDM::mac2THDM(const StandardModel& SM_i)
: ThObservable(SM_i), myAxions(static_cast<const AxionsTHDM*> (&SM_i))
{
};

mac2THDM::~mac2THDM()
{
};

double mac2THDM::computeThValue()
{
    return myAxions->gae()/0.28e-13;

}


G117B15ATHDM::G117B15ATHDM(const StandardModel& SM_i)
: ThObservable(SM_i), myAxions(static_cast<const AxionsTHDM*> (&SM_i))
{
};

G117B15ATHDM::~G117B15ATHDM()
{
};

double G117B15ATHDM::computeThValue()
{
    double a=myAxions->geta_G117B15A();
    double b=myAxions->getb_G117B15A();

    double a26=myAxions->gae()*myAxions->gae()/(4.*M_PI)/1.e-26;

    return a + a26*b;

}


R548THDM::R548THDM(const StandardModel& SM_i)
: ThObservable(SM_i), myAxions(static_cast<const AxionsTHDM*> (&SM_i))
{
};

R548THDM::~R548THDM()
{
};

double R548THDM::computeThValue()
{
    double a=myAxions->geta_R548();
    double b=myAxions->getb_R548();

    double a26=myAxions->gae()*myAxions->gae()/(4.*M_PI)/1.e-26;

    return a + a26*b;

}


PG1351489THDM::PG1351489THDM(const StandardModel& SM_i)
: ThObservable(SM_i), myAxions(static_cast<const AxionsTHDM*> (&SM_i))
{
};

PG1351489THDM::~PG1351489THDM()
{
};

double PG1351489THDM::computeThValue()
{
    double a=myAxions->geta_PG1351489();
    double b=myAxions->getb_PG1351489();

    double a26=myAxions->gae()*myAxions->gae()/(4.*M_PI)/1.e-26;

    return a + a26*b;

}


L113THDM::L113THDM(const StandardModel& SM_i)
: ThObservable(SM_i), myAxions(static_cast<const AxionsTHDM*> (&SM_i))
{
};

L113THDM::~L113THDM()
{
};

double L113THDM::computeThValue()
{
    double a=myAxions->geta_L113();
    double b=myAxions->getb_L113();

    double a26=myAxions->gae()*myAxions->gae()/(4.*M_PI)/1.e-26;

    return a + a26*b;

}


L192THDM::L192THDM(const StandardModel& SM_i)
: ThObservable(SM_i), myAxions(static_cast<const AxionsTHDM*> (&SM_i))
{
};

L192THDM::~L192THDM()
{
};

double L192THDM::computeThValue()
{
    double a=myAxions->geta_L192();
    double b=myAxions->getb_L192();

    double a26=myAxions->gae()*myAxions->gae()/(4.*M_PI)/1.e-26;

    return a + a26*b;

}


TRGBTHDM::TRGBTHDM(const StandardModel& SM_i)
: ThObservable(SM_i), myAxions(static_cast<const AxionsTHDM*> (&SM_i))
{
};

TRGBTHDM::~TRGBTHDM()
{
};

double TRGBTHDM::computeThValue()
{
    double a=myAxions->geta_TRGB();
    double b=myAxions->getb_TRGB();

    double gae=myAxions->gae()/1.e-13;

    return - 4.08 - 0.25*(sqrt(gae*gae + 0.93) - 0.96 - 0.17*pow(gae,1.5))
            + 0.039 + a + gae*b;

}


HBRTHDM::HBRTHDM(const StandardModel& SM_i)
: ThObservable(SM_i), myAxions(static_cast<const AxionsTHDM*> (&SM_i))
{
};

HBRTHDM::~HBRTHDM()
{
};

double HBRTHDM::computeThValue()
{
    double Y=myAxions->getY_HBR();

    gslpp::complex gag=myAxions->gag();
    double gae=myAxions->gae()/1.e-13;

    double alpha = gae*gae/4./M_PI;
    double g10 = gag.abs()*1.e10;

    double dMc = 0.024*(sqrt(gae*gae + 1.23*1.23) - 1.23 - 0.921*pow(alpha,0.75));

    return 7.33*Y + 0.02 - 0.095*sqrt(21.86 + 21.08*g10) - 1.61*dMc - 0.067*alpha;  // 1512.08108, Eq. (7.6)

    return 6.26*Y - 0.12 - 0.14*g10*g10 - 1.61*dMc - 0.067*alpha;  // 1512.08108, Eq. (7.5)

    return 6.26*Y - 0.41*g10*g10 - 0.12;  // 1406.6053, Eq. (1)

}


GaNTHDM::GaNTHDM(const StandardModel& SM_i)
: ThObservable(SM_i), myAxions(static_cast<const AxionsTHDM*> (&SM_i))
{
};

GaNTHDM::~GaNTHDM()
{
};

double GaNTHDM::computeThValue()
{
    double gan=myAxions->gan();
    double gap=myAxions->gap();

    return gan*gan + 0.61*gap*gap + 0.53*gan*gap;
    
    return gan*gan + 0.29*gap*gap + 0.27*gan*gap; // old, published bound

}


ganTHDM::ganTHDM(const StandardModel& SM_i)
: ThObservable(SM_i), myAxions(static_cast<const AxionsTHDM*> (&SM_i))
{
};

ganTHDM::~ganTHDM()
{
};

double ganTHDM::computeThValue()
{
    if (myAxions->getmodel() == -1.) {
        return 0.;
    }
    else return std::abs(myAxions->gan());

}


CaeTHDM::CaeTHDM(const StandardModel& SM_i)
: ThObservable(SM_i), myAxions(static_cast<const AxionsTHDM*> (&SM_i))
{
};

CaeTHDM::~CaeTHDM()
{
};

double CaeTHDM::computeThValue()
{
    if (myAxions->getmodel() == -1.) {
        return 0.;
    }
    else return std::abs(myAxions->Cae());

}


logcgTHDM::logcgTHDM(const StandardModel& SM_i)
: ThObservable(SM_i), myAxions(static_cast<const AxionsTHDM*> (&SM_i))
{
};

logcgTHDM::~logcgTHDM()
{
};

double logcgTHDM::computeThValue()
{
    if ((myAxions->getmodel() == 11.) || (myAxions->getmodel() == 12.)) {
        return log10(myAxions->getcgamma());
    }
    else return 0.;

}
