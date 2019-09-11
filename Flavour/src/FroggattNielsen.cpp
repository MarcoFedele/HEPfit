/* 
 * Copyright (C) 2018 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "FroggattNielsen.h"

const std::string FroggattNielsen::FroggattNielsenvars[NFroggattNielsenvars] = {
    "eps", "mupFN", "mcharmFN", "mtopFN", "mdownFN", "mstrangeFN", "mbottomFN",
    "s12CKM", "s13CKM", "s23CKM", "deltaCKM",
    "QL1", "QL2", "QL3",
    "QRu1", "QRu2", "QRu3",
    "QRd1", "QRd2", "QRd3",
    "thL1", "thL2", "thL3", "phiL1", "phiL2", "phiL3",
    "thRu1", "thRu2", "thRu3", "phiRu1", "phiRu2", "phiRu3",
    "thRd1", "thRd2", "thRd3", "phiRd1", "phiRd2", "phiRd3"};

FroggattNielsen::FroggattNielsen() : StandardModel() {   
    ModelParamMap.insert(std::make_pair("eps", std::cref(eps)));
    ModelParamMap.insert(std::make_pair("mupFN", std::cref(mup)));
    ModelParamMap.insert(std::make_pair("mcharmFN", std::cref(mcharm)));
    ModelParamMap.insert(std::make_pair("mtopFN", std::cref(mtop)));
    ModelParamMap.insert(std::make_pair("mdownFN", std::cref(mdown)));
    ModelParamMap.insert(std::make_pair("mstrangeFN", std::cref(mstrange)));
    ModelParamMap.insert(std::make_pair("mbottomFN", std::cref(mbottom)));
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
    ModelParamMap.insert(std::make_pair("thRu1", std::cref(thRu1)));
    ModelParamMap.insert(std::make_pair("thRu2", std::cref(thRu2)));
    ModelParamMap.insert(std::make_pair("thRu3", std::cref(thRu3)));
    ModelParamMap.insert(std::make_pair("phiRu1", std::cref(phiRu1)));
    ModelParamMap.insert(std::make_pair("phiRu2", std::cref(phiRu2)));
    ModelParamMap.insert(std::make_pair("phiRu3", std::cref(phiRu3)));
    ModelParamMap.insert(std::make_pair("thRd1", std::cref(thRd1)));
    ModelParamMap.insert(std::make_pair("thRd2", std::cref(thRd2)));
    ModelParamMap.insert(std::make_pair("thRd3", std::cref(thRd3)));
    ModelParamMap.insert(std::make_pair("phiRd1", std::cref(phiRd1)));
    ModelParamMap.insert(std::make_pair("phiRd2", std::cref(phiRd2)));
    ModelParamMap.insert(std::make_pair("phiRd3", std::cref(phiRd3)));

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
    
    JCKM = c12CKM*c13CKM*c13CKM*s12CKM*s13CKM*s23CKM*deltaCKM;

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
    else if(name.compare("mupFN") == 0) {
        mup = value;
    }
    else if(name.compare("mcharmFN") == 0) {
        mcharm = value;
    }
    else if(name.compare("mtopFN") == 0) {
        mtop = value;
    }
    else if(name.compare("mdownFN") == 0) {
        mdown = value;
    }
    else if(name.compare("mstrangeFN") == 0) {
        mstrange = value;
    }
    else if(name.compare("mbottomFN") == 0) {
        mbottom = value;
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

gslpp::matrix<gslpp::complex> FroggattNielsen::CKM() const
{   
    gslpp::matrix<gslpp::complex> V(3,3);
    
    V.assign(0, 0, c12CKM*c13CKM);
    V.assign(0, 1, s12CKM*c13CKM);
    V.assign(0, 2, gslpp::complex(s13CKM, -deltaCKM, true));

    V.assign(1, 0, (-s12CKM*c23CKM - gslpp::complex(c12CKM*s23CKM*s13CKM, deltaCKM, true)));
    V.assign(1, 1, (c12CKM*c23CKM - gslpp::complex(s12CKM*s23CKM*s13CKM, deltaCKM, true)));
    V.assign(1, 2, s23CKM*c13CKM);

    V.assign(2, 0, (s12CKM*s23CKM - gslpp::complex(c12CKM*c23CKM*s13CKM, deltaCKM, true)));
    V.assign(2, 1, (-c12CKM*s23CKM - gslpp::complex(s12CKM*c23CKM*s13CKM, deltaCKM, true)));
    V.assign(2, 2, c23CKM*c13CKM);  
    
    return V;
}

gslpp::matrix<gslpp::complex> FroggattNielsen::lambda_u() const
{   
    gslpp::matrix<gslpp::complex> lambdau(3,3,0);
    
    lambdau.assign(0, 0, getmup()*sqrt(2.)/v());

    lambdau.assign(1, 1, getmcharm()*sqrt(2.)/v());

    lambdau.assign(2, 2, getmtop()*sqrt(2.)/v());
    
    return lambdau;
}

gslpp::matrix<gslpp::complex> FroggattNielsen::lambda_d() const
{   
    gslpp::matrix<gslpp::complex> lambdad(3,3,0);
    
    lambdad.assign(0, 0, getmdown()*sqrt(2.)/v());

    lambdad.assign(1, 1, getmstrange()*sqrt(2.)/v());

    lambdad.assign(2, 2, getmbottom()*sqrt(2.)/v());
    
    lambdad = CKM()*lambdad;
    
    return lambdad;
}

gslpp::matrix<gslpp::complex> FroggattNielsen::create_V(double th1, double th2, double th3, double phi1, double phi2, double phi3) const
{   
    double c1 = cos(th1);
    double s1 = sin(th1);
    double c2 = cos(th2);
    double s2 = sin(th2);
    double c3 = cos(th3);
    double s3 = sin(th3);
    
    gslpp::matrix<gslpp::complex> V(3,3,0);
    
    V.assign(0, 0, c1*c2);
    V.assign(0, 1, c1*s2);
    V.assign(0, 2, gslpp::complex(s1, -phi1, true));
    
    V.assign(1, 0, (-gslpp::complex(c3*s2, phi2, true)-gslpp::complex(c2*s1*s3, phi1+phi2, true)));
    V.assign(1, 1, (gslpp::complex(c2*c3, phi2, true)-gslpp::complex(s1*s2*s3, phi1+phi2, true)));
    V.assign(1, 2, gslpp::complex(c1*s3, phi2, true).abs());
    
    V.assign(2, 0, (gslpp::complex(s2*s3, phi3, true)-gslpp::complex(c2*c3*s1, phi1+phi3, true)));
    V.assign(2, 1, (-gslpp::complex(c2*s3, phi3, true)-gslpp::complex(c3*s1*s2, phi1+phi3, true)));
    V.assign(2, 2, gslpp::complex(c1*c3, phi3, true).abs());
    
    return V;
}

gslpp::matrix<gslpp::complex> FroggattNielsen::y_u() const
{   
    gslpp::matrix<gslpp::complex> VL(3,3,0);
    gslpp::matrix<gslpp::complex> VRu(3,3,0);
    
    VL = create_V(thL1, thL2, thL3, phiL1, phiL2, phiL3);
    VRu = create_V(thRu1, thRu2, thRu3, phiRu1, phiRu2, phiRu3);
    
    gslpp::matrix<gslpp::complex> yu(3,3,0);
    
    yu = VL.hconjugate()*lambda_u()*VRu;
    
    return yu;
}

gslpp::matrix<gslpp::complex> FroggattNielsen::y_d() const
{   
    gslpp::matrix<gslpp::complex> VL(3,3,0);
    gslpp::matrix<gslpp::complex> VRd(3,3,0);
    
    VL = create_V(thL1, thL2, thL3, phiL1, phiL2, phiL3);
    VRd = create_V(thRd1, thRd2, thRd3, phiRd1, phiRd2, phiRd3);
    
    gslpp::matrix<gslpp::complex> yd(3,3,0);
    
    yd = VL.hconjugate()*lambda_d()*VRd;
    
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

double FroggattNielsen::delta() const
{
    gslpp::matrix<gslpp::complex> delta_yu(3,3,0), m_u(3,3,0), m_u_2(3,3,0), V_u(3,3,0);
    gslpp::matrix<gslpp::complex> delta_yd(3,3,0), m_d(3,3,0), m_d_2(3,3,0), V_d(3,3,0);
    gslpp::matrix<gslpp::complex> ckm(3,3,0), C(3,3,0);
    vector<double> m_u_2_diag(3,0), m_d_2_diag(3,0);
    double s_12, c_12, s_13, c_13, s_23, c_23, T, B, JCPV;
    gslpp::complex Cdet;
    double delta_mu, delta_mc, delta_mt, delta_md, delta_ms, delta_mb, delta_s12, delta_s13, delta_s23, delta_J, c_over_delta_c;
    double delta = 0.;
    double delta_c = 0.001;
    
    /* real up */
    for (int i=0; i<3; i++)
    {
        for (int j=0; j<3; j++)
        {
            c_over_delta_c = ( y_u()(i,j)/F_u()(i,j) ).real() / delta_c;

            delta_yu = y_u();
            delta_yu.assign(i,j, delta_yu(i,j) + delta_c * delta_yu(i,j).real() );
            m_u = delta_yu * v() / sqrt(2.);
            m_u_2 = m_u * m_u.hconjugate();
            m_u_2.eigensystem(V_u, m_u_2_diag);

            m_d = y_d() * v() / sqrt(2.);
            m_d_2 = m_d * m_d.hconjugate();
            m_d_2.eigensystem(V_d, m_d_2_diag);
            
            ckm = V_u.hconjugate()*V_d;
            
            s_13 = ckm(0,2).abs();
            c_13 = sqrt(1. - s_13*s_13);
            c_12 = ckm(0,0).abs() / c_13;
            s_12 = sqrt(1. - c_12*c_12);
            s_23 = ckm(1,2).abs() / c_13;
            c_23 = sqrt(1. - s_23*s_23);
            
            C = m_u_2*m_d_2 - m_d_2*m_u_2;
            T = (m_u_2_diag(2) - m_u_2_diag(1)) * (m_u_2_diag(2) - m_u_2_diag(0)) * (m_u_2_diag(1) - m_u_2_diag(0));
            B = (m_d_2_diag(2) - m_d_2_diag(1)) * (m_d_2_diag(2) - m_d_2_diag(0)) * (m_d_2_diag(1) - m_d_2_diag(0));
            Cdet = C(0,0)*(C(1,1)*C(2,2) - C(1,2)*C(2,1))
                    - C(0,1)*(C(1,0)*C(2,2) - C(1,2)*C(2,0))
                    + C(0,2)*(C(1,0)*C(2,1) - C(1,1)*C(2,0));
            JCPV = ( Cdet / (2. * T * B) ).imag();

            delta_mu = abs(c_over_delta_c * (sqrt(m_u_2_diag(0)) - mup)) / mup;
            delta_mc = abs(c_over_delta_c * (sqrt(m_u_2_diag(1)) - mcharm)) / mcharm;
            delta_mt = abs(c_over_delta_c * (sqrt(m_u_2_diag(2)) - mtop)) / mtop;
            delta_md = abs(c_over_delta_c * (sqrt(m_d_2_diag(0)) - mdown)) / mdown;
            delta_ms = abs(c_over_delta_c * (sqrt(m_d_2_diag(1)) - mstrange)) / mstrange;
            delta_mb = abs(c_over_delta_c * (sqrt(m_d_2_diag(2)) - mbottom)) / mbottom;
            delta_s12 = abs(c_over_delta_c * (s_12 - s12CKM)) / s12CKM;
            delta_s13 = abs(c_over_delta_c * (s_13 - s13CKM)) / s13CKM;
            delta_s23 = abs(c_over_delta_c * (s_23 - s23CKM)) / s23CKM;
            delta_J = abs(c_over_delta_c * (JCPV - JCKM)) / JCKM;
            
            if (delta_mu > delta)
                delta = delta_mu;
            if (delta_mc > delta)
                delta = delta_mc;
            if (delta_mt > delta)
                delta = delta_mt;
            if (delta_md > delta)
                delta = delta_md;
            if (delta_ms > delta)
                delta = delta_ms;
            if (delta_mb > delta)
                delta = delta_mb;
            if (delta_s12 > delta)
                delta = delta_s12;
            if (delta_s13 > delta)
                delta = delta_s13;
            if (delta_s23 > delta)
                delta = delta_s23;
            if (delta_J > delta)
                delta = delta_J;
            
            /*std::cout << "real u" << std::endl;
            std::cout << "c_over_delta_c = " << c_over_delta_c << std::endl;
            std::cout << "m_u_2_diag = " << m_u_2_diag << std::endl;
            std::cout << "m_d_2_diag = " << m_d_2_diag << std::endl;
            std::cout << "JCKM = " << JCKM << std::endl << std::endl;
            std::cout << "delta_mu = " << delta_mu << std::endl;
            std::cout << "delta_mc = " << delta_mc << std::endl;
            std::cout << "delta_mt = " << delta_mt << std::endl;
            std::cout << "delta_md = " << delta_md << std::endl;
            std::cout << "delta_ms = " << delta_ms << std::endl;
            std::cout << "delta_mb = " << delta_mb << std::endl;
            std::cout << "delta_s12 = " << delta_s12 << std::endl;
            std::cout << "delta_s13 = " << delta_s13 << std::endl;
            std::cout << "delta_s23 = " << delta_s23 << std::endl;
            std::cout << "delta_J = " << delta_J << std::endl;
            std::cout << "delta = " << delta << std::endl << std::endl;*/
        }
    }
    
    /* imag up */
    for (int i=0; i<3; i++)
    {
        for (int j=0; j<3; j++)
        {
            c_over_delta_c = ( y_u()(i,j)/F_u()(i,j) ).imag() / delta_c;

            delta_yu = y_u();
            delta_yu.assign(i,j, delta_yu(i,j) + delta_c * gslpp::complex::i() * delta_yu(i,j).imag() );
            m_u = delta_yu * v() / sqrt(2.);
            m_u_2 = m_u * m_u.hconjugate();
            m_u_2.eigensystem(V_u, m_u_2_diag);

            m_d = y_d() * v() / sqrt(2.);
            m_d_2 = m_d * m_d.hconjugate();
            m_d_2.eigensystem(V_d, m_d_2_diag);
            
            ckm = V_u.hconjugate()*V_d;
            
            s_13 = ckm(0,2).abs();
            c_13 = sqrt(1. - s_13*s_13);
            c_12 = ckm(0,0).abs() / c_13;
            s_12 = sqrt(1. - c_12*c_12);
            s_23 = ckm(1,2).abs() / c_13;
            c_23 = sqrt(1. - s_23*s_23);
            
            C = m_u_2*m_d_2 - m_d_2*m_u_2;
            T = (m_u_2_diag(2) - m_u_2_diag(1)) * (m_u_2_diag(2) - m_u_2_diag(0)) * (m_u_2_diag(1) - m_u_2_diag(0));
            B = (m_d_2_diag(2) - m_d_2_diag(1)) * (m_d_2_diag(2) - m_d_2_diag(0)) * (m_d_2_diag(1) - m_d_2_diag(0));
            Cdet = C(0,0)*(C(1,1)*C(2,2) - C(1,2)*C(2,1))
                    - C(0,1)*(C(1,0)*C(2,2) - C(1,2)*C(2,0))
                    + C(0,2)*(C(1,0)*C(2,1) - C(1,1)*C(2,0));
            JCPV = ( Cdet / (2. * T * B) ).imag();

            delta_mu = abs(c_over_delta_c * (sqrt(m_u_2_diag(0)) - mup)) / mup;
            delta_mc = abs(c_over_delta_c * (sqrt(m_u_2_diag(1)) - mcharm)) / mcharm;
            delta_mt = abs(c_over_delta_c * (sqrt(m_u_2_diag(2)) - mtop)) / mtop;
            delta_md = abs(c_over_delta_c * (sqrt(m_d_2_diag(0)) - mdown)) / mdown;
            delta_ms = abs(c_over_delta_c * (sqrt(m_d_2_diag(1)) - mstrange)) / mstrange;
            delta_mb = abs(c_over_delta_c * (sqrt(m_d_2_diag(2)) - mbottom)) / mbottom;
            delta_s12 = abs(c_over_delta_c * (s_12 - s12CKM)) / s12CKM;
            delta_s13 = abs(c_over_delta_c * (s_13 - s13CKM)) / s13CKM;
            delta_s23 = abs(c_over_delta_c * (s_23 - s23CKM)) / s23CKM;
            delta_J = abs(c_over_delta_c * (JCPV - JCKM)) / JCKM;
            
            if (delta_mu > delta)
                delta = delta_mu;
            if (delta_mc > delta)
                delta = delta_mc;
            if (delta_mt > delta)
                delta = delta_mt;
            if (delta_md > delta)
                delta = delta_md;
            if (delta_ms > delta)
                delta = delta_ms;
            if (delta_mb > delta)
                delta = delta_mb;
            if (delta_s12 > delta)
                delta = delta_s12;
            if (delta_s13 > delta)
                delta = delta_s13;
            if (delta_s23 > delta)
                delta = delta_s23;
            if (delta_J > delta)
                delta = delta_J;
            
            /*std::cout << "imag u" << std::endl;
            std::cout << "c_over_delta_c = " << c_over_delta_c << std::endl;
            std::cout << "m_u_2_diag = " << m_u_2_diag << std::endl;
            std::cout << "m_d_2_diag = " << m_d_2_diag << std::endl;
            std::cout << "JCKM = " << JCKM << std::endl << std::endl;
            std::cout << "delta_mu = " << delta_mu << std::endl;
            std::cout << "delta_mc = " << delta_mc << std::endl;
            std::cout << "delta_mt = " << delta_mt << std::endl;
            std::cout << "delta_md = " << delta_md << std::endl;
            std::cout << "delta_ms = " << delta_ms << std::endl;
            std::cout << "delta_mb = " << delta_mb << std::endl;
            std::cout << "delta_s12 = " << delta_s12 << std::endl;
            std::cout << "delta_s13 = " << delta_s13 << std::endl;
            std::cout << "delta_s23 = " << delta_s23 << std::endl;
            std::cout << "delta_J = " << delta_J << std::endl;
            std::cout << "delta = " << delta << std::endl << std::endl;*/
        }
    }
    
    /* real down */
    for (int i=0; i<3; i++)
    {
        for (int j=0; j<3; j++)
        {
            c_over_delta_c = ( y_d()(i,j)/F_d()(i,j) ).real() / delta_c;

            m_u = y_u() * v() / sqrt(2.);
            m_u_2 = m_u * m_u.hconjugate();
            m_u_2.eigensystem(V_u, m_u_2_diag);

            delta_yd = y_d();
            delta_yd.assign(i,j, delta_yd(i,j) + delta_c * delta_yd(i,j).real() );
            m_d = delta_yd * v() / sqrt(2.);
            m_d_2 = m_d * m_d.hconjugate();
            m_d_2.eigensystem(V_d, m_d_2_diag);
            
            ckm = V_u.hconjugate()*V_d;
            
            s_13 = ckm(0,2).abs();
            c_13 = sqrt(1. - s_13*s_13);
            c_12 = ckm(0,0).abs() / c_13;
            s_12 = sqrt(1. - c_12*c_12);
            s_23 = ckm(1,2).abs() / c_13;
            c_23 = sqrt(1. - s_23*s_23);
            
            C = m_u_2*m_d_2 - m_d_2*m_u_2;
            T = (m_u_2_diag(2) - m_u_2_diag(1)) * (m_u_2_diag(2) - m_u_2_diag(0)) * (m_u_2_diag(1) - m_u_2_diag(0));
            B = (m_d_2_diag(2) - m_d_2_diag(1)) * (m_d_2_diag(2) - m_d_2_diag(0)) * (m_d_2_diag(1) - m_d_2_diag(0));
            Cdet = C(0,0)*(C(1,1)*C(2,2) - C(1,2)*C(2,1))
                    - C(0,1)*(C(1,0)*C(2,2) - C(1,2)*C(2,0))
                    + C(0,2)*(C(1,0)*C(2,1) - C(1,1)*C(2,0));
            JCPV = ( Cdet / (2. * T * B) ).imag();

            delta_mu = abs(c_over_delta_c * (sqrt(m_u_2_diag(0)) - mup)) / mup;
            delta_mc = abs(c_over_delta_c * (sqrt(m_u_2_diag(1)) - mcharm)) / mcharm;
            delta_mt = abs(c_over_delta_c * (sqrt(m_u_2_diag(2)) - mtop)) / mtop;
            delta_md = abs(c_over_delta_c * (sqrt(m_d_2_diag(0)) - mdown)) / mdown;
            delta_ms = abs(c_over_delta_c * (sqrt(m_d_2_diag(1)) - mstrange)) / mstrange;
            delta_mb = abs(c_over_delta_c * (sqrt(m_d_2_diag(2)) - mbottom)) / mbottom;
            delta_s12 = abs(c_over_delta_c * (s_12 - s12CKM)) / s12CKM;
            delta_s13 = abs(c_over_delta_c * (s_13 - s13CKM)) / s13CKM;
            delta_s23 = abs(c_over_delta_c * (s_23 - s23CKM)) / s23CKM;
            delta_J = abs(c_over_delta_c * (JCPV - JCKM)) / JCKM;
            
            if (delta_mu > delta)
                delta = delta_mu;
            if (delta_mc > delta)
                delta = delta_mc;
            if (delta_mt > delta)
                delta = delta_mt;
            if (delta_md > delta)
                delta = delta_md;
            if (delta_ms > delta)
                delta = delta_ms;
            if (delta_mb > delta)
                delta = delta_mb;
            if (delta_s12 > delta)
                delta = delta_s12;
            if (delta_s13 > delta)
                delta = delta_s13;
            if (delta_s23 > delta)
                delta = delta_s23;
            if (delta_J > delta)
                delta = delta_J;
            
            /*std::cout << "real d" << std::endl;
            std::cout << "c_over_delta_c = " << c_over_delta_c << std::endl;
            std::cout << "m_u_2_diag = " << m_u_2_diag << std::endl;
            std::cout << "m_d_2_diag = " << m_d_2_diag << std::endl;
            std::cout << "JCKM = " << JCKM << std::endl << std::endl;
            std::cout << "delta_mu = " << delta_mu << std::endl;
            std::cout << "delta_mc = " << delta_mc << std::endl;
            std::cout << "delta_mt = " << delta_mt << std::endl;
            std::cout << "delta_md = " << delta_md << std::endl;
            std::cout << "delta_ms = " << delta_ms << std::endl;
            std::cout << "delta_mb = " << delta_mb << std::endl;
            std::cout << "delta_s12 = " << delta_s12 << std::endl;
            std::cout << "delta_s13 = " << delta_s13 << std::endl;
            std::cout << "delta_s23 = " << delta_s23 << std::endl;
            std::cout << "delta_J = " << delta_J << std::endl;
            std::cout << "delta = " << delta << std::endl << std::endl;*/
        }
    }
    
    /* imag down */
    for (int i=0; i<3; i++)
    {
        for (int j=0; j<3; j++)
        {
            c_over_delta_c = ( y_d()(i,j)/F_d()(i,j) ).imag() / delta_c;

            m_u = y_u() * v() / sqrt(2.);
            m_u_2 = m_u * m_u.hconjugate();
            m_u_2.eigensystem(V_u, m_u_2_diag);

            delta_yd = y_d();
            delta_yd.assign(i,j, delta_yd(i,j) + delta_c * gslpp::complex::i() * delta_yd(i,j).imag() );
            m_d = delta_yd * v() / sqrt(2.);
            m_d_2 = m_d * m_d.hconjugate();
            m_d_2.eigensystem(V_d, m_d_2_diag);
            
            ckm = V_u.hconjugate()*V_d;
            
            s_13 = ckm(0,2).abs();
            c_13 = sqrt(1. - s_13*s_13);
            c_12 = ckm(0,0).abs() / c_13;
            s_12 = sqrt(1. - c_12*c_12);
            s_23 = ckm(1,2).abs() / c_13;
            c_23 = sqrt(1. - s_23*s_23);
            
            C = m_u_2*m_d_2 - m_d_2*m_u_2;
            T = (m_u_2_diag(2) - m_u_2_diag(1)) * (m_u_2_diag(2) - m_u_2_diag(0)) * (m_u_2_diag(1) - m_u_2_diag(0));
            B = (m_d_2_diag(2) - m_d_2_diag(1)) * (m_d_2_diag(2) - m_d_2_diag(0)) * (m_d_2_diag(1) - m_d_2_diag(0));
            Cdet = C(0,0)*(C(1,1)*C(2,2) - C(1,2)*C(2,1))
                    - C(0,1)*(C(1,0)*C(2,2) - C(1,2)*C(2,0))
                    + C(0,2)*(C(1,0)*C(2,1) - C(1,1)*C(2,0));
            JCPV = ( Cdet / (2. * T * B) ).imag();

            delta_mu = abs(c_over_delta_c * (sqrt(m_u_2_diag(0)) - mup)) / mup;
            delta_mc = abs(c_over_delta_c * (sqrt(m_u_2_diag(1)) - mcharm)) / mcharm;
            delta_mt = abs(c_over_delta_c * (sqrt(m_u_2_diag(2)) - mtop)) / mtop;
            delta_md = abs(c_over_delta_c * (sqrt(m_d_2_diag(0)) - mdown)) / mdown;
            delta_ms = abs(c_over_delta_c * (sqrt(m_d_2_diag(1)) - mstrange)) / mstrange;
            delta_mb = abs(c_over_delta_c * (sqrt(m_d_2_diag(2)) - mbottom)) / mbottom;
            delta_s12 = abs(c_over_delta_c * (s_12 - s12CKM)) / s12CKM;
            delta_s13 = abs(c_over_delta_c * (s_13 - s13CKM)) / s13CKM;
            delta_s23 = abs(c_over_delta_c * (s_23 - s23CKM)) / s23CKM;
            delta_J = abs(c_over_delta_c * (JCPV - JCKM)) / JCKM;
            
            if (delta_mu > delta)
                delta = delta_mu;
            if (delta_mc > delta)
                delta = delta_mc;
            if (delta_mt > delta)
                delta = delta_mt;
            if (delta_md > delta)
                delta = delta_md;
            if (delta_ms > delta)
                delta = delta_ms;
            if (delta_mb > delta)
                delta = delta_mb;
            if (delta_s12 > delta)
                delta = delta_s12;
            if (delta_s13 > delta)
                delta = delta_s13;
            if (delta_s23 > delta)
                delta = delta_s23;
            if (delta_J > delta)
                delta = delta_J;
            
            /*std::cout << "imag d" << std::endl;
            std::cout << "c_over_delta_c = " << c_over_delta_c << std::endl;
            std::cout << "m_u_2_diag = " << m_u_2_diag << std::endl;
            std::cout << "m_d_2_diag = " << m_d_2_diag << std::endl;
            std::cout << "JCKM = " << JCKM << std::endl << std::endl;
            std::cout << "delta_mu = " << delta_mu << std::endl;
            std::cout << "delta_mc = " << delta_mc << std::endl;
            std::cout << "delta_mt = " << delta_mt << std::endl;
            std::cout << "delta_md = " << delta_md << std::endl;
            std::cout << "delta_ms = " << delta_ms << std::endl;
            std::cout << "delta_mb = " << delta_mb << std::endl;
            std::cout << "delta_s12 = " << delta_s12 << std::endl;
            std::cout << "delta_s13 = " << delta_s13 << std::endl;
            std::cout << "delta_s23 = " << delta_s23 << std::endl;
            std::cout << "delta_J = " << delta_J << std::endl;
            std::cout << "delta = " << delta << std::endl << std::endl;*/
        }
    }
    
    return delta;
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
    std::cout << "md : " << myFroggattNielsen->getmdown() << std::endl;
    std::cout << "mu : " << myFroggattNielsen->getmup() << std::endl;
    std::cout << "ms : " << myFroggattNielsen->getmstrange() << std::endl;
    std::cout << "mc : " << myFroggattNielsen->getmcharm() << std::endl;
    std::cout << "mb : " << myFroggattNielsen->getmbottom() << std::endl;
    std::cout << "mt : " << myFroggattNielsen->getmtop() << std::endl << std::endl;
    std::cout << "vev : " << SM.v()/sqrt(2.) << std::endl << std::endl;
    std::cout << "ckm : " << myFroggattNielsen->CKM() << std::endl << std::endl;
    std::cout << "QL1 : " << myFroggattNielsen->getQL1() << std::endl;
    std::cout << "QL2 : " << myFroggattNielsen->getQL2() << std::endl;
    std::cout << "QL3 : " << myFroggattNielsen->getQL3() << std::endl;
    std::cout << "QRu1 : " << myFroggattNielsen->getQRu1() << std::endl;
    std::cout << "QRu2 : " << myFroggattNielsen->getQRu2() << std::endl;
    std::cout << "QRu3 : " << myFroggattNielsen->getQRu3() << std::endl;
    std::cout << "QRd1 : " << myFroggattNielsen->getQRd1() << std::endl;
    std::cout << "QRd2 : " << myFroggattNielsen->getQRd2() << std::endl;
    std::cout << "QRd3 : " << myFroggattNielsen->getQRd3() << std::endl << std::endl;
    
    std::cout << "lambda_u : " << myFroggattNielsen->lambda_u() << std::endl;
    std::cout << "y_u : " << myFroggattNielsen->y_u() << std::endl;
    std::cout << "F_u : " << myFroggattNielsen->F_u() << std::endl;
    std::cout << "lambda_d : " << myFroggattNielsen->lambda_d() << std::endl;
    std::cout << "y_d : " << myFroggattNielsen->y_d() << std::endl;
    std::cout << "F_d : " << myFroggattNielsen->F_d() << std::endl << std::endl;/*
    
    std::cout << "delta : " << myFroggattNielsen->delta() << std::endl << std::endl;*/
    
    return myFroggattNielsen->geteps();

}

Ru11::Ru11(const StandardModel& SM_i)
: ThObservable(SM_i), myFroggattNielsen(static_cast<const FroggattNielsen*> (&SM_i))
{
};

Ru11::~Ru11()
{
};

double Ru11::computeThValue()
{
    gslpp::matrix<gslpp::complex> yu = myFroggattNielsen->y_u();
    gslpp::matrix<double> Fu = myFroggattNielsen->F_u();
    
    return yu(0,0).abs()/Fu(0,0);
}

Ru12::Ru12(const StandardModel& SM_i)
: ThObservable(SM_i), myFroggattNielsen(static_cast<const FroggattNielsen*> (&SM_i))
{
};

Ru12::~Ru12()
{
};

double Ru12::computeThValue()
{
    gslpp::matrix<gslpp::complex> yu = myFroggattNielsen->y_u();
    gslpp::matrix<double> Fu = myFroggattNielsen->F_u();
    
    return yu(0,1).abs()/Fu(0,1);
}

Ru13::Ru13(const StandardModel& SM_i)
: ThObservable(SM_i), myFroggattNielsen(static_cast<const FroggattNielsen*> (&SM_i))
{
};

Ru13::~Ru13()
{
};

double Ru13::computeThValue()
{
    gslpp::matrix<gslpp::complex> yu = myFroggattNielsen->y_u();
    gslpp::matrix<double> Fu = myFroggattNielsen->F_u();
    
    return yu(0,2).abs()/Fu(0,2);
}

Ru21::Ru21(const StandardModel& SM_i)
: ThObservable(SM_i), myFroggattNielsen(static_cast<const FroggattNielsen*> (&SM_i))
{
};

Ru21::~Ru21()
{
};

double Ru21::computeThValue()
{
    gslpp::matrix<gslpp::complex> yu = myFroggattNielsen->y_u();
    gslpp::matrix<double> Fu = myFroggattNielsen->F_u();
    
    return yu(1,0).abs()/Fu(1,0);
}

Ru22::Ru22(const StandardModel& SM_i)
: ThObservable(SM_i), myFroggattNielsen(static_cast<const FroggattNielsen*> (&SM_i))
{
};

Ru22::~Ru22()
{
};

double Ru22::computeThValue()
{
    gslpp::matrix<gslpp::complex> yu = myFroggattNielsen->y_u();
    gslpp::matrix<double> Fu = myFroggattNielsen->F_u();
    
    return yu(1,1).abs()/Fu(1,1);
}

Ru23::Ru23(const StandardModel& SM_i)
: ThObservable(SM_i), myFroggattNielsen(static_cast<const FroggattNielsen*> (&SM_i))
{
};

Ru23::~Ru23()
{
};

double Ru23::computeThValue()
{
    gslpp::matrix<gslpp::complex> yu = myFroggattNielsen->y_u();
    gslpp::matrix<double> Fu = myFroggattNielsen->F_u();
    
    return yu(1,2).abs()/Fu(1,2);
}

Ru31::Ru31(const StandardModel& SM_i)
: ThObservable(SM_i), myFroggattNielsen(static_cast<const FroggattNielsen*> (&SM_i))
{
};

Ru31::~Ru31()
{
};

double Ru31::computeThValue()
{
    gslpp::matrix<gslpp::complex> yu = myFroggattNielsen->y_u();
    gslpp::matrix<double> Fu = myFroggattNielsen->F_u();
    
    return yu(2,0).abs()/Fu(2,0);
}

Ru32::Ru32(const StandardModel& SM_i)
: ThObservable(SM_i), myFroggattNielsen(static_cast<const FroggattNielsen*> (&SM_i))
{
};

Ru32::~Ru32()
{
};

double Ru32::computeThValue()
{
    gslpp::matrix<gslpp::complex> yu = myFroggattNielsen->y_u();
    gslpp::matrix<double> Fu = myFroggattNielsen->F_u();
    
    return yu(2,1).abs()/Fu(2,1);
}

Ru33::Ru33(const StandardModel& SM_i)
: ThObservable(SM_i), myFroggattNielsen(static_cast<const FroggattNielsen*> (&SM_i))
{
};

Ru33::~Ru33()
{
};

double Ru33::computeThValue()
{
    gslpp::matrix<gslpp::complex> yu = myFroggattNielsen->y_u();
    gslpp::matrix<double> Fu = myFroggattNielsen->F_u();
    
    return yu(2,2).abs()/Fu(2,2);
}

Rd11::Rd11(const StandardModel& SM_i)
: ThObservable(SM_i), myFroggattNielsen(static_cast<const FroggattNielsen*> (&SM_i))
{
};

Rd11::~Rd11()
{
};

double Rd11::computeThValue()
{
    gslpp::matrix<gslpp::complex> yd = myFroggattNielsen->y_d();
    gslpp::matrix<double> Fd = myFroggattNielsen->F_d();
    
    return yd(0,0).abs()/Fd(0,0);
}

Rd12::Rd12(const StandardModel& SM_i)
: ThObservable(SM_i), myFroggattNielsen(static_cast<const FroggattNielsen*> (&SM_i))
{
};

Rd12::~Rd12()
{
};

double Rd12::computeThValue()
{
    gslpp::matrix<gslpp::complex> yd = myFroggattNielsen->y_d();
    gslpp::matrix<double> Fd = myFroggattNielsen->F_d();
    
    return yd(0,1).abs()/Fd(0,1);
}

Rd13::Rd13(const StandardModel& SM_i)
: ThObservable(SM_i), myFroggattNielsen(static_cast<const FroggattNielsen*> (&SM_i))
{
};

Rd13::~Rd13()
{
};

double Rd13::computeThValue()
{
    gslpp::matrix<gslpp::complex> yd = myFroggattNielsen->y_d();
    gslpp::matrix<double> Fd = myFroggattNielsen->F_d();
    
    return yd(0,2).abs()/Fd(0,2);
}

Rd21::Rd21(const StandardModel& SM_i)
: ThObservable(SM_i), myFroggattNielsen(static_cast<const FroggattNielsen*> (&SM_i))
{
};

Rd21::~Rd21()
{
};

double Rd21::computeThValue()
{
    gslpp::matrix<gslpp::complex> yd = myFroggattNielsen->y_d();
    gslpp::matrix<double> Fd = myFroggattNielsen->F_d();
    
    return yd(1,0).abs()/Fd(1,0);
}

Rd22::Rd22(const StandardModel& SM_i)
: ThObservable(SM_i), myFroggattNielsen(static_cast<const FroggattNielsen*> (&SM_i))
{
};

Rd22::~Rd22()
{
};

double Rd22::computeThValue()
{
    gslpp::matrix<gslpp::complex> yd = myFroggattNielsen->y_d();
    gslpp::matrix<double> Fd = myFroggattNielsen->F_d();
    
    return yd(1,1).abs()/Fd(1,1);
}

Rd23::Rd23(const StandardModel& SM_i)
: ThObservable(SM_i), myFroggattNielsen(static_cast<const FroggattNielsen*> (&SM_i))
{
};

Rd23::~Rd23()
{
};

double Rd23::computeThValue()
{
    gslpp::matrix<gslpp::complex> yd = myFroggattNielsen->y_d();
    gslpp::matrix<double> Fd = myFroggattNielsen->F_d();
    
    return yd(1,2).abs()/Fd(1,2);
}

Rd31::Rd31(const StandardModel& SM_i)
: ThObservable(SM_i), myFroggattNielsen(static_cast<const FroggattNielsen*> (&SM_i))
{
};

Rd31::~Rd31()
{
};

double Rd31::computeThValue()
{
    gslpp::matrix<gslpp::complex> yd = myFroggattNielsen->y_d();
    gslpp::matrix<double> Fd = myFroggattNielsen->F_d();
    
    return yd(2,0).abs()/Fd(2,0);
}

Rd32::Rd32(const StandardModel& SM_i)
: ThObservable(SM_i), myFroggattNielsen(static_cast<const FroggattNielsen*> (&SM_i))
{
};

Rd32::~Rd32()
{
};

double Rd32::computeThValue()
{
    gslpp::matrix<gslpp::complex> yd = myFroggattNielsen->y_d();
    gslpp::matrix<double> Fd = myFroggattNielsen->F_d();
    
    return yd(2,1).abs()/Fd(2,1);
}

Rd33::Rd33(const StandardModel& SM_i)
: ThObservable(SM_i), myFroggattNielsen(static_cast<const FroggattNielsen*> (&SM_i))
{
};

Rd33::~Rd33()
{
};

double Rd33::computeThValue()
{
    gslpp::matrix<gslpp::complex> yd = myFroggattNielsen->y_d();
    gslpp::matrix<double> Fd = myFroggattNielsen->F_d();
    
    return yd(2,2).abs()/Fd(2,2);
}

FN_delta::FN_delta(const StandardModel& SM_i)
: ThObservable(SM_i), myFroggattNielsen(static_cast<const FroggattNielsen*> (&SM_i))
{
};

FN_delta::~FN_delta()
{
};

double FN_delta::computeThValue()
{
    return myFroggattNielsen->delta();
}