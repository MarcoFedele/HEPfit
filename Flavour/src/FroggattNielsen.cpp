/* 
 * Copyright (C) 2018 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "FroggattNielsen.h"

const std::string FroggattNielsen::FroggattNielsenvars[NFroggattNielsenvars] = {
    "eps", "mtopEW",
    "s12CKM", "s13CKM", "s23CKM", "deltaCKM",
    "QL1", "QL2", "QL3",
    "QRu1", "QRu2", "QRu3",
    "QRd1", "QRd2", "QRd3",
    "thL1", "thL2", "thL3", "phiL1", "phiL2", "phiL3", "phiL4", "phiL5",
    "thRu1", "thRu2", "thRu3", "phiRu1", "phiRu2", "phiRu3", "phiRu4", "phiRu5",
    "thRd1", "thRd2", "thRd3", "phiRd1", "phiRd2", "phiRd3", "phiRd4", "phiRd5"};

FroggattNielsen::FroggattNielsen() : StandardModel() {   
    ModelParamMap.insert(std::make_pair("eps", std::cref(eps)));
    ModelParamMap.insert(std::make_pair("mtopEW", std::cref(mtopEW)));
    ModelParamMap.insert(std::make_pair("s12CKM", std::cref(s12CKM)));
    ModelParamMap.insert(std::make_pair("s13CKM", std::cref(s13CKM)));
    ModelParamMap.insert(std::make_pair("s23CKM", std::cref(s23CKM)));
    ModelParamMap.insert(std::make_pair("deltaCKM", std::cref(deltaCKM)));
    ModelParamMap.insert(std::make_pair("QL1", std::cref(QL1)));
    ModelParamMap.insert(std::make_pair("QL2", std::cref(QL2)));
    ModelParamMap.insert(std::make_pair("QL3", std::cref(QL3)));
    ModelParamMap.insert(std::make_pair("QRu1", std::cref(QRu1)));
    ModelParamMap.insert(std::make_pair("QRu2", std::cref(QRu2)));
    ModelParamMap.insert(std::make_pair("QRu3", std::cref(QRu3)));
    ModelParamMap.insert(std::make_pair("QRd1", std::cref(QRd1)));
    ModelParamMap.insert(std::make_pair("QRd2", std::cref(QRd2)));
    ModelParamMap.insert(std::make_pair("QRd3", std::cref(QRd3)));
    ModelParamMap.insert(std::make_pair("thL1", std::cref(thL1)));
    ModelParamMap.insert(std::make_pair("thL2", std::cref(thL2)));
    ModelParamMap.insert(std::make_pair("thL3", std::cref(thL3)));
    ModelParamMap.insert(std::make_pair("phiL1", std::cref(phiL1)));
    ModelParamMap.insert(std::make_pair("phiL2", std::cref(phiL2)));
    ModelParamMap.insert(std::make_pair("phiL3", std::cref(phiL3)));
    ModelParamMap.insert(std::make_pair("phiL4", std::cref(phiL4)));
    ModelParamMap.insert(std::make_pair("phiL5", std::cref(phiL5)));
    ModelParamMap.insert(std::make_pair("thRu1", std::cref(thRu1)));
    ModelParamMap.insert(std::make_pair("thRu2", std::cref(thRu2)));
    ModelParamMap.insert(std::make_pair("thRu3", std::cref(thRu3)));
    ModelParamMap.insert(std::make_pair("phiRu1", std::cref(phiRu1)));
    ModelParamMap.insert(std::make_pair("phiRu2", std::cref(phiRu2)));
    ModelParamMap.insert(std::make_pair("phiRu3", std::cref(phiRu3)));
    ModelParamMap.insert(std::make_pair("phiRu4", std::cref(phiRu4)));
    ModelParamMap.insert(std::make_pair("phiRu5", std::cref(phiRu5)));
    ModelParamMap.insert(std::make_pair("thRd1", std::cref(thRd1)));
    ModelParamMap.insert(std::make_pair("thRd2", std::cref(thRd2)));
    ModelParamMap.insert(std::make_pair("thRd3", std::cref(thRd3)));
    ModelParamMap.insert(std::make_pair("phiRd1", std::cref(phiRd1)));
    ModelParamMap.insert(std::make_pair("phiRd2", std::cref(phiRd2)));
    ModelParamMap.insert(std::make_pair("phiRd3", std::cref(phiRd3)));
    ModelParamMap.insert(std::make_pair("phiRd4", std::cref(phiRd4)));
    ModelParamMap.insert(std::make_pair("phiRd5", std::cref(phiRd5)));

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

    return (true);
}

void FroggattNielsen::setParameter(const std::string name, const double& value){    
    if(name.compare("eps") == 0) {
        eps = value;
    }
    else if(name.compare("mtopEW") == 0) {
        mtopEW = value;
    }
    else if(name.compare("s12CKM") == 0) {
        s12CKM = value;
        c12CKM = sqrt(1. - s12CKM*s12CKM);
    }
    else if(name.compare("s13CKM") == 0) {
        s13CKM = value;
        c13CKM = sqrt(1. - s13CKM*s13CKM);
    }
    else if(name.compare("s23CKM") == 0) {
        s23CKM = value;
        c23CKM = sqrt(1. - s23CKM*s23CKM);
    }
    else if(name.compare("deltaCKM") == 0) {
        deltaCKM = value;
    }
    else if(name.compare("QL1") == 0) {
        QL1 = value;
    }
    else if(name.compare("QL2") == 0) {
        QL2 = value;
    }
    else if(name.compare("QL3") == 0) {
        QL3 = value;
    }
    else if(name.compare("QRu1") == 0) {
        QRu1 = value;
    }
    else if(name.compare("QRu2") == 0) {
        QRu2 = value;
    }
    else if(name.compare("QRu3") == 0) {
        QRu3 = value;
    }
    else if(name.compare("QRd1") == 0) {
        QRd1 = value;
    }
    else if(name.compare("QRd2") == 0) {
        QRd2 = value;
    }
    else if(name.compare("QRd3") == 0) {
        QRd3 = value;
    }
    else if(name.compare("thL1") == 0) {
        thL1 = value;
    }
    else if(name.compare("thL2") == 0) {
        thL2 = value;
    }
    else if(name.compare("thL3") == 0) {
        thL3 = value;
    }
    else if(name.compare("phiL1") == 0) {
        phiL1 = value;
    }
    else if(name.compare("phiL2") == 0) {
        phiL2 = value;
    }
    else if(name.compare("phiL3") == 0) {
        phiL3 = value;
    }
    else if(name.compare("phiL4") == 0) {
        phiL4 = value;
    }
    else if(name.compare("phiL5") == 0) {
        phiL5 = value;
    }
    else if(name.compare("thRu1") == 0) {
        thRu1 = value;
    }
    else if(name.compare("thRu2") == 0) {
        thRu2 = value;
    }
    else if(name.compare("thRu3") == 0) {
        thRu3 = value;
    }
    else if(name.compare("phiRu1") == 0) {
        phiRu1 = value;
    }
    else if(name.compare("phiRu2") == 0) {
        phiRu2 = value;
    }
    else if(name.compare("phiRu3") == 0) {
        phiRu3 = value;
    }
    else if(name.compare("phiRu4") == 0) {
        phiRu4 = value;
    }
    else if(name.compare("phiRu5") == 0) {
        phiRu5 = value;
    }
    else if(name.compare("thRd1") == 0) {
        thRd1 = value;
    }
    else if(name.compare("thRd2") == 0) {
        thRd2 = value;
    }
    else if(name.compare("thRd3") == 0) {
        thRd3 = value;
    }
    else if(name.compare("phiRd1") == 0) {
        phiRd1 = value;
    }
    else if(name.compare("phiRd2") == 0) {
        phiRd2 = value;
    }
    else if(name.compare("phiRd3") == 0) {
        phiRd3 = value;
    }
    else if(name.compare("phiRd4") == 0) {
        phiRd4 = value;
    }
    else if(name.compare("phiRd5") == 0) {
        phiRd5 = value;
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

gslpp::matrix<double> FroggattNielsen::absCKM() const
{   
    gslpp::matrix<double> V(3,3);
    
    V.assign(0, 0, c12CKM*c13CKM);
    V.assign(0, 1, s12CKM*c13CKM);
    V.assign(0, 2, gslpp::complex(s13CKM, -deltaCKM, true).abs());

    V.assign(1, 0, (-s12CKM * c23CKM - gslpp::complex(c12CKM * s23CKM*s13CKM, deltaCKM, true)).abs());
    V.assign(1, 1, (c12CKM * c23CKM - gslpp::complex(s12CKM * s23CKM*s13CKM, deltaCKM, true)).abs());
    V.assign(1, 2, s23CKM*c13CKM);

    V.assign(2, 0, (s12CKM * s23CKM - gslpp::complex(c12CKM * c23CKM*s13CKM, deltaCKM, true)).abs());
    V.assign(2, 1, (-c12CKM * s23CKM - gslpp::complex(s12CKM * c23CKM*s13CKM, deltaCKM, true)).abs());
    V.assign(2, 2, c23CKM*c13CKM);  
    
    return V;
}

gslpp::matrix<double> FroggattNielsen::lambda_u() const
{   
    gslpp::matrix<double> lambdau(3,3,0);
    
    lambdau.assign(0, 0, getQuarks(QCD::UP).getMass()*sqrt(2.)/v());

    lambdau.assign(1, 1, getQuarks(QCD::CHARM).getMass()*sqrt(2.)/v());

    lambdau.assign(2, 2, getmtopEW()*sqrt(2.)/v());
    
    return lambdau;
}

gslpp::matrix<double> FroggattNielsen::lambda_d() const
{   
    gslpp::matrix<double> lambdad(3,3,0);
    
    lambdad.assign(0, 0, getQuarks(QCD::DOWN).getMass()*sqrt(2.)/v());

    lambdad.assign(1, 1, getQuarks(QCD::STRANGE).getMass()*sqrt(2.)/v());

    lambdad.assign(2, 2, getQuarks(QCD::BOTTOM).getMass()*sqrt(2.)/v());
    
    lambdad = absCKM()*lambdad;
    
    return lambdad;
}

gslpp::matrix<double> FroggattNielsen::create_V(double th1, double th2, double th3, double phi1, double phi2, double phi3, double phi4, double phi5) const
{   
    double c1 = cos(th1);
    double s1 = sin(th1);
    double c2 = cos(th2);
    double s2 = sin(th2);
    double c3 = cos(th3);
    double s3 = sin(th3);
    
    gslpp::matrix<double> V(3,3,0);
    
    V.assign(0, 0, gslpp::complex(c1*c2, phi1, true).abs());
    V.assign(0, 1, gslpp::complex(s1, phi3, true).abs());
    V.assign(0, 2, gslpp::complex(c1*s2, phi4, true).abs());
    
    V.assign(1, 0, (gslpp::complex(s2*s3, -phi4-phi5, true)-gslpp::complex(c2*c3*s1, phi1+phi2-phi3, true)).abs());
    V.assign(1, 1, gslpp::complex(c1*c3, phi2, true).abs());
    V.assign(1, 2, (gslpp::complex(-c2*s3, -phi1-phi5, true)-gslpp::complex(c3*s1*s2, phi2-phi3+phi4, true)).abs());
    
    V.assign(2, 0, (-gslpp::complex(c3*s2, -phi2-phi4, true)-gslpp::complex(c2*s1*s3, phi1-phi3+phi5, true)).abs());
    V.assign(2, 1, gslpp::complex(c1*s3, phi5, true).abs());
    V.assign(2, 2, (gslpp::complex(c2*c3, -phi1-phi2, true)-gslpp::complex(s1*s2*s3, -phi3+phi4+phi5, true)).abs());
    
    return V;
}

gslpp::matrix<double> FroggattNielsen::y_u() const
{   
    gslpp::matrix<double> VL(3,3,0);
    gslpp::matrix<double> VRu(3,3,0);
    
    VL = create_V(thL1, thL2, thL3, phiL1, phiL2, phiL3, phiL4, phiL5);
    VRu = create_V(thRu1, thRu2, thRu3, phiRu1, phiRu2, phiRu3, phiRu4, phiRu5);
    
    gslpp::matrix<double> yu(3,3,0);
    
    yu = VL.transpose()*lambda_u()*VRu;
    
    return yu;
}

gslpp::matrix<double> FroggattNielsen::y_d() const
{   
    gslpp::matrix<double> VL(3,3,0);
    gslpp::matrix<double> VRd(3,3,0);
    
    VL = create_V(thL1, thL2, thL3, phiL1, phiL2, phiL3, phiL4, phiL5);
    VRd = create_V(thRd1, thRd2, thRd3, phiRd1, phiRd2, phiRd3, phiRd4, phiRd5);
    
    gslpp::matrix<double> yd(3,3,0);
    
    yd = VL.transpose()*lambda_d()*VRd;
    
    return yd;
}

gslpp::matrix<double> FroggattNielsen::F_u() const
{   
    double eps = geteps();
    double QL1 = getQL1();
    double QL2 = getQL2();
    double QL3 = getQL3();
    double QRu1 = getQRu1();
    double QRu2 = getQRu2();
    double QRu3 = getQRu3();
    
    gslpp::matrix<double> Fu(3,3,0);
    
    Fu.assign(0, 0, pow(eps, QL1 - QRu1));
    Fu.assign(0, 1, pow(eps, QL1 - QRu2));
    Fu.assign(0, 2, pow(eps, QL1 - QRu3));
    
    Fu.assign(1, 0, pow(eps, QL2 - QRu1));
    Fu.assign(1, 1, pow(eps, QL2 - QRu2));
    Fu.assign(1, 2, pow(eps, QL2 - QRu3));
    
    Fu.assign(2, 0, pow(eps, QL3 - QRu1));
    Fu.assign(2, 1, pow(eps, QL3 - QRu2));
    Fu.assign(2, 2, pow(eps, QL3 - QRu3));
    
    return Fu;
}

gslpp::matrix<double> FroggattNielsen::F_d() const
{   
    double eps = geteps();
    double QL1 = getQL1();
    double QL2 = getQL2();
    double QL3 = getQL3();
    double QRd1 = getQRd1();
    double QRd2 = getQRd2();
    double QRd3 = getQRd3();
    
    gslpp::matrix<double> Fd(3,3,0);
    
    Fd.assign(0, 0, pow(eps, QL1 - QRd1));
    Fd.assign(0, 1, pow(eps, QL1 - QRd2));
    Fd.assign(0, 2, pow(eps, QL1 - QRd3));
    
    Fd.assign(1, 0, pow(eps, QL2 - QRd1));
    Fd.assign(1, 1, pow(eps, QL2 - QRd2));
    Fd.assign(1, 2, pow(eps, QL2 - QRd3));
    
    Fd.assign(2, 0, pow(eps, QL3 - QRd1));
    Fd.assign(2, 1, pow(eps, QL3 - QRd2));
    Fd.assign(2, 2, pow(eps, QL3 - QRd3));
    
    return Fd;
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
   /* std::cout << "md : " << SM.getQuarks(QCD::DOWN).getMass() << std::endl;
    std::cout << "mu : " << SM.getQuarks(QCD::UP).getMass() << std::endl;
    std::cout << "ms : " << SM.getQuarks(QCD::STRANGE).getMass() << std::endl;
    std::cout << "mc : " << SM.getQuarks(QCD::CHARM).getMass() << std::endl;
    std::cout << "mb : " << SM.getQuarks(QCD::BOTTOM).getMass() << std::endl;
    std::cout << "mt : " << myFroggattNielsen->getmtopEW() << std::endl << std::endl;
    std::cout << "vev : " << SM.v()/sqrt(2.) << std::endl << std::endl;
    std::cout << "ckm : " << myFroggattNielsen->absCKM() << std::endl << std::endl;
    std::cout << "QL1 : " << myFroggattNielsen->getQL1() << std::endl;
    std::cout << "QL2 : " << myFroggattNielsen->getQL2() << std::endl;
    std::cout << "QL3 : " << myFroggattNielsen->getQL3() << std::endl;
    std::cout << "QRu1 : " << myFroggattNielsen->getQRu1() << std::endl;
    std::cout << "QRu2 : " << myFroggattNielsen->getQRu2() << std::endl;
    std::cout << "QRu3 : " << myFroggattNielsen->getQRu3() << std::endl;
    std::cout << "QRd1 : " << myFroggattNielsen->getQRd1() << std::endl;
    std::cout << "QRd2 : " << myFroggattNielsen->getQRd2() << std::endl;
    std::cout << "QRd3 : " << myFroggattNielsen->getQRd3() << std::endl << std::endl;*/
    
    std::cout << "lambda_u : " << myFroggattNielsen->lambda_u() << std::endl;
    std::cout << "y_u : " << myFroggattNielsen->y_u() << std::endl;
    std::cout << "F_u : " << myFroggattNielsen->F_u() << std::endl;
    std::cout << "lambda_d : " << myFroggattNielsen->lambda_d() << std::endl;
    std::cout << "y_d : " << myFroggattNielsen->y_d() << std::endl;
    std::cout << "F_d : " << myFroggattNielsen->F_d() << std::endl << std::endl;
    
    return myFroggattNielsen->geteps();

}