/*
 * Copyright (C) 2018 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "AxionsTHDM.h"

const std::string AxionsTHDM::Axionsvars[NAxionsvars] = {"logtanb", "logma", "model", "eps_L", "C_pn_err",
                                "a_G117B15A", "b_G117B15A", "c_G117B15A", "d_G117B15A",
                                "a_R548", "b_R548", "c_R548", "d_R548",
                                "a_PG1351489", "b_PG1351489", "c_PG1351489", "d_PG1351489",
                                "a_L192", "b_L192", "c_L192", "d_L192",
                                "a_TRGB", "b_TRGB",
                                "Y_HBR"};

AxionsTHDM::AxionsTHDM() : StandardModel() {
    ModelParamMap.insert(std::make_pair("logtanb", std::cref(tanb)));
    ModelParamMap.insert(std::make_pair("logma", std::cref(ma)));
    ModelParamMap.insert(std::make_pair("model", std::cref(model)));
    ModelParamMap.insert(std::make_pair("eps_L", std::cref(eps_L)));
    ModelParamMap.insert(std::make_pair("C_pn_err", std::cref(C_pn_err)));

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
    else if(name.compare("C_pn_err") == 0)
        C_pn_err = value;
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

    if (model == 0.)
        Cag = - 1.92;
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
    else
        throw std::runtime_error("error in AzionsTHDM::gag, model can only be an integer between 0. and 6. !");

    return getAle()/2./M_PI * ma/(5.7e9) * Cag;
}

double AxionsTHDM::gae() const
{
    double Cae;

    if (model == 0.)
        Cae = 0.;
    else if (model == 1.)
        Cae = sinb*sinb/3.;
    else if (model == 2.)
        Cae = cosb*cosb/3.;
    else if (model == 3.)
        Cae = cosb*cosb + eps_L;
    else if (model == 4.)
        Cae = sinb*sinb + eps_L;
    else if (model == 5.)
        Cae = cosb*cosb + eps_L;
    else if (model == 6.)
        Cae = sinb*sinb + eps_L;
    else
        throw std::runtime_error("error in AzionsTHDM::gae, model can only be an integer between 0. and 6. !");

    return 0.84e-13 * ma * Cae;
}

double AxionsTHDM::gap() const
{
    double Cad, Cas, Cab;
    double Cau, Cac, Cat;

    if (model == 0.) {
        Cad = 0.;
        Cas = 0.;
        Cab = 0.;
        Cau = 0.;
        Cac = 0.;
        Cat = 0.;
    }
    else if (model == 1.){
        Cad = sinb*sinb/3.;
        Cas = sinb*sinb/3.;
        Cab = sinb*sinb/3.;
        Cau = cosb*cosb/3.;
        Cac = cosb*cosb/3.;
        Cat = cosb*cosb/3.;
    }
    else if (model == 2.){
        Cad = sinb*sinb/3.;
        Cas = sinb*sinb/3.;
        Cab = sinb*sinb/3.;
        Cau = cosb*cosb/3.;
        Cac = cosb*cosb/3.;
        Cat = cosb*cosb/3.;
    }
    else if (model == 3.){
        Cad = sinb*sinb;
        Cas = sinb*sinb;
        Cab = -cosb*cosb;
        Cau = cosb*cosb;
        Cac = cosb*cosb;
        Cat = -sinb*sinb;
    }
    else if (model == 4.){
        Cad = sinb*sinb;
        Cas = 0.;
        Cab = 0.;
        Cau = cosb*cosb;
        Cac = 0.;
        Cat = 0.;
    }
    else if (model == 5.){
        Cad = sinb*sinb;
        Cas = 0.;
        Cab = 0.;
        Cau = cosb*cosb;
        Cac = 0.;
        Cat = 0.;
    }
    else if (model == 6.){
        Cad = sinb*sinb;
        Cas = 0.;
        Cab = 0.;
        Cau = cosb*cosb;
        Cac = 0.;
        Cat = 0.;
    }
    else
        throw std::runtime_error("error in AzionsTHDM::gap, model can only be an integer between 0. and 6. !");
    
    double Cap = (-0.47+C_pn_err)
            + (0.88+C_pn_err)*Cau + (-0.39+C_pn_err)*Cad - 0.038*Cas - 0.012*Cac - 0.009*Cab - 0.0035*Cat;

    //Cap = - 0.435*sinb*sinb - 0.182;

    return 1.56e-10 * ma * Cap;
}

double AxionsTHDM::gan() const
{
    double Cad, Cas, Cab;
    double Cau, Cac, Cat;

    if (model == 0.) {
        Cad = 0.;
        Cas = 0.;
        Cab = 0.;
        Cau = 0.;
        Cac = 0.;
        Cat = 0.;
    }
    else if (model == 1.){
        Cad = sinb*sinb/3.;
        Cas = sinb*sinb/3.;
        Cab = sinb*sinb/3.;
        Cau = cosb*cosb/3.;
        Cac = cosb*cosb/3.;
        Cat = cosb*cosb/3.;
    }
    else if (model == 2.){
        Cad = sinb*sinb/3.;
        Cas = sinb*sinb/3.;
        Cab = sinb*sinb/3.;
        Cau = cosb*cosb/3.;
        Cac = cosb*cosb/3.;
        Cat = cosb*cosb/3.;
    }
    else if (model == 3.){
        Cad = sinb*sinb;
        Cas = sinb*sinb;
        Cab = -cosb*cosb;
        Cau = cosb*cosb;
        Cac = cosb*cosb;
        Cat = -sinb*sinb;
    }
    else if (model == 4.){
        Cad = sinb*sinb;
        Cas = 0.;
        Cab = 0.;
        Cau = cosb*cosb;
        Cac = 0.;
        Cat = 0.;
    }
    else if (model == 5.){
        Cad = sinb*sinb;
        Cas = 0.;
        Cab = 0.;
        Cau = cosb*cosb;
        Cac = 0.;
        Cat = 0.;
    }
    else if (model == 6.){
        Cad = sinb*sinb;
        Cas = 0.;
        Cab = 0.;
        Cau = cosb*cosb;
        Cac = 0.;
        Cat = 0.;
    }
    else
        throw std::runtime_error("error in AzionsTHDM::gan, model can only be an integer between 0. and 6. !");

    double Can = (-0.02+C_pn_err)
            + (0.88+C_pn_err)*Cad + (-0.39+C_pn_err)*Cau - 0.038*Cas - 0.012*Cac - 0.009*Cab - 0.0035*Cat;

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
    return log10(myAxions->gettanb());

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
    return log10(myAxions->getma());

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
    return log10(myAxions->gag());

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
    return log10(myAxions->gae());

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
    double c=myAxions->getc_G117B15A();
    double d=myAxions->getd_G117B15A();

    double mac2=myAxions->gae()/0.28e-13;

    return a + mac2*b + mac2*mac2*c + mac2*mac2*mac2*d;

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
    double c=myAxions->getc_R548();
    double d=myAxions->getd_R548();

    double mac2=myAxions->gae()/0.28e-13;

    return a + mac2*b + mac2*mac2*c + mac2*mac2*mac2*d;

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
    double c=myAxions->getc_PG1351489();
    double d=myAxions->getd_PG1351489();

    double mac2=myAxions->gae()/0.28e-13;

    return a + mac2*b + mac2*mac2*c + mac2*mac2*mac2*d;

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
    double c=myAxions->getc_L192();
    double d=myAxions->getd_L192();

    double mac2=myAxions->gae()/0.28e-13;

    return a + mac2*b + mac2*mac2*c + mac2*mac2*mac2*d;

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

    double gae=myAxions->gae();

    return - 4.03 - 0.25*(sqrt(gae*gae + 0.93) - 0.96 - 0.17*pow(gae,1.5))
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
    double gae=myAxions->gae();

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

    return gap*gap + gan*gan;

}
