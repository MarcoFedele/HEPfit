/* 
 * Copyright (C) 2018 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "LoopMediatorsAMatching.h"
#include "LoopMediatorsA.h"
#include <stdexcept>

LoopMediatorsAMatching::LoopMediatorsAMatching(const LoopMediatorsA & LoopMediatorsA_i) :

    StandardModelMatching(LoopMediatorsA_i),
    myLoopMediatorsA(LoopMediatorsA_i),
    mcbsg(8, NDR, NNLO),
    mcprimebsg(8, NDR, NNLO),
    mcBMll(13, NDR, NLO),
    mcprimeBMll(13, NDR, NLO),
    mcdbs2(5, NDR, NLO)
{}

void LoopMediatorsAMatching::updateLoopMediatorsAParameters()
{   
    
    C1NP = myLoopMediatorsA.getC1();
    C2NP = myLoopMediatorsA.getC2();
    C3NP = myLoopMediatorsA.getC3();
    C4NP = myLoopMediatorsA.getC4();
    C5NP = myLoopMediatorsA.getC5();
    
    C1pNP = myLoopMediatorsA.getC1p();
    C2pNP = myLoopMediatorsA.getC2p();
    C3pNP = myLoopMediatorsA.getC3p();
    
    C7NP = myLoopMediatorsA.getC7();
    C8NP = myLoopMediatorsA.getC8();
    C9NPmu = myLoopMediatorsA.getC9();
    C10NPmu = myLoopMediatorsA.getC10();
    CSNPmu = myLoopMediatorsA.getCS();
    CPNPmu = myLoopMediatorsA.getCP();
    
    C7pNP = myLoopMediatorsA.getC7p();
    C8pNP = myLoopMediatorsA.getC8p();
    C9pNPmu = myLoopMediatorsA.getC9p();
    C10pNPmu = myLoopMediatorsA.getC10p();
    CSpNPmu = myLoopMediatorsA.getCSp();
    CPpNPmu = myLoopMediatorsA.getCPp();
    
    WCscale = myLoopMediatorsA.getWCscale();

    StandardModelMatching::updateSMParameters();
}

LoopMediatorsAMatching::~LoopMediatorsAMatching()
{}


std::vector<WilsonCoefficient>& LoopMediatorsAMatching::CMbsg()
{
    vmcbsg.clear();
    for (std::vector<WilsonCoefficient>::iterator it = StandardModelMatching::CMbsg().begin(); it != StandardModelMatching::CMbsg().end(); it++ ) vmcbsg.push_back(*it);

    switch (mcbsg.getScheme()) {
        case NDR:

            break;
        default:
            std::stringstream out;
            out << mcbsg.getScheme();
            throw std::runtime_error("LoopMediatorsAMatching::CMbsg(): scheme " + out.str() + "not implemented"); 
    }

    mcbsg.setMu(WCscale);

    switch (mcbsg.getOrder()) {
        case NNLO:
            mcbsg.setCoeff(6, 0., NNLO);
        case NLO:
            mcbsg.setCoeff(6, 0., NLO);
        case LO:
            mcbsg.setCoeff(6, C7NP, LO);
            break;
        default:
            std::stringstream out;
            out << mcbsg.getOrder();
            throw std::runtime_error("LoopMediatorsAMatching::CMbsg(): order " + out.str() + "not implemented"); 
    }

    vmcbsg.push_back(mcbsg);
    return (vmcbsg);
}
 
std::vector<WilsonCoefficient>& LoopMediatorsAMatching::CMprimebsg()
{
    vmcprimebsg.clear();
    for (std::vector<WilsonCoefficient>::iterator it = StandardModelMatching::CMprimebsg().begin(); it != StandardModelMatching::CMprimebsg().end(); it++ ) vmcprimebsg.push_back(*it);

    switch (mcprimebsg.getScheme()) {
        case NDR:

            break;
        default:
            std::stringstream out;
            out << mcprimebsg.getScheme();
            throw std::runtime_error("LoopMediatorsAMatching::CMprimebsg(): scheme " + out.str() + "not implemented"); 
    }

    mcprimebsg.setMu(WCscale);

    switch (mcprimebsg.getOrder()) {
        case NNLO:
            mcprimebsg.setCoeff(6, 0., NNLO);
        case NLO:
            mcprimebsg.setCoeff(6, 0., NLO);
        case LO:
            mcprimebsg.setCoeff(6, C7pNP, LO);
            break;
        default:
            std::stringstream out;
            out << mcprimebsg.getOrder();
            throw std::runtime_error("LoopMediatorsAMatching::CMprimebsg(): order " + out.str() + "not implemented"); 
    }

    vmcprimebsg.push_back(mcprimebsg);
    return (vmcprimebsg);
}

std::vector<WilsonCoefficient>& LoopMediatorsAMatching::CMBMll(QCD::lepton lepton)
{
    vmcBMll.clear();
    for (std::vector<WilsonCoefficient>::iterator it = StandardModelMatching::CMBMll(lepton).begin(); it != StandardModelMatching::CMBMll(lepton).end(); it++ ) vmcBMll.push_back(*it);

    switch (mcBMll.getScheme()) {
        case NDR:

            break;
        default:
            std::stringstream out;
            out << mcBMll.getScheme();
            throw std::runtime_error("LoopMediatorsAMatching::CMBMll(): scheme " + out.str() + "not implemented"); 
    }

    mcBMll.setMu(WCscale);

    switch (mcBMll.getOrder()) {
        case NNLO:
        case NLO:
                mcBMll.setCoeff(6, 0., NLO);
                mcBMll.setCoeff(8, 0., NLO);
                mcBMll.setCoeff(9, 0., NLO);
                mcBMll.setCoeff(10, 0., NLO);
                mcBMll.setCoeff(11, 0., NLO);
        case LO:
            mcBMll.setCoeff(6, C7NP, LO);
            if(lepton == LoopMediatorsA::MU){
                mcBMll.setCoeff(8, C9NPmu, LO);
                mcBMll.setCoeff(9, C10NPmu, LO);
                mcBMll.setCoeff(10, CSNPmu, LO);
                mcBMll.setCoeff(11, CPNPmu, LO);               
            }
            break;
        default:
            std::stringstream out;
            out << mcBMll.getOrder();
            throw std::runtime_error("LoopMediatorsAMatching::CMBMll(): order " + out.str() + "not implemented"); 
    }

    vmcBMll.push_back(mcBMll);
    return (vmcBMll);
}

std::vector<WilsonCoefficient>& LoopMediatorsAMatching::CMprimeBMll(QCD::lepton lepton)
{
    vmcprimeBMll.clear();
    for (std::vector<WilsonCoefficient>::iterator it = StandardModelMatching::CMprimeBMll(lepton).begin(); it != StandardModelMatching::CMprimeBMll(lepton).end(); it++ ) vmcprimeBMll.push_back(*it);

    switch (mcprimeBMll.getScheme()) {
        case NDR:

            break;
        default:
            std::stringstream out;
            out << mcprimeBMll.getScheme();
            throw std::runtime_error("LoopMediatorsAMatching::CMprimeBMll(): scheme " + out.str() + "not implemented"); 
    }

    mcprimeBMll.setMu(WCscale);

    switch (mcprimeBMll.getOrder()) {
        case NNLO:
        case NLO:
                mcprimeBMll.setCoeff(6, 0., NLO);
                mcprimeBMll.setCoeff(8, 0., NLO);
                mcprimeBMll.setCoeff(9, 0., NLO);
                mcprimeBMll.setCoeff(10, 0., NLO);
                mcprimeBMll.setCoeff(11, 0., NLO);
        case LO:
            mcprimeBMll.setCoeff(6, C7pNP, LO);
            if(lepton == LoopMediatorsA::MU){
                mcprimeBMll.setCoeff(8, C9pNPmu, LO);
                mcprimeBMll.setCoeff(9, C10pNPmu, LO);
                mcprimeBMll.setCoeff(10, CSpNPmu, LO);
                mcprimeBMll.setCoeff(11, CPpNPmu, LO);              
            }
            break;
        default:
            std::stringstream out;
            out << mcprimeBMll.getOrder();
            throw std::runtime_error("LoopMediatorsAMatching::CMprimeBMll(): order " + out.str() + "not implemented"); 
    }

    vmcprimeBMll.push_back(mcprimeBMll);
    return (vmcprimeBMll);
}
 
std::vector<WilsonCoefficient>& LoopMediatorsAMatching::CMdbs2()
{
    vmcdbs2.clear();
    for (std::vector<WilsonCoefficient>::iterator it = StandardModelMatching::CMdbs2().begin(); it != StandardModelMatching::CMdbs2().end(); it++ ) vmcdbs2.push_back(*it);

    switch (mcdbs2.getScheme()) {
        case NDR:

            break;
        default:
            std::stringstream out;
            out << mcdbs2.getScheme();
            throw std::runtime_error("LoopMediatorsAMatching::CMdbs2(): scheme " + out.str() + "not implemented"); 
    }

    mcdbs2.setMu(WCscale);
    
    /*std::cout << "high scale :" << std::endl;
    std::cout << "C1 :" << C1NP + C1pNP << std::endl;
    std::cout << "C2 :" << C2NP + C2pNP << std::endl;
    std::cout << "C3 :" << C3NP + C3pNP << std::endl;
    std::cout << "C4 :" << C4NP << std::endl;
    std::cout << "C5 :" << C5NP << std::endl << std::endl;*/

    switch (mcdbs2.getOrder()) {
        case NNLO:
        case NLO:
            mcdbs2.setCoeff(0, 0., NLO);
            mcdbs2.setCoeff(1, 0., NLO);
            mcdbs2.setCoeff(2, 0., NLO);
            mcdbs2.setCoeff(3, 0., NLO);
            mcdbs2.setCoeff(4, 0., NLO);
        case LO:
            mcdbs2.setCoeff(0, C1NP + C1pNP, LO);
            mcdbs2.setCoeff(1, C2NP + C2pNP, LO);
            mcdbs2.setCoeff(2, C3NP + C3pNP, LO);
            mcdbs2.setCoeff(3, C4NP, LO);
            mcdbs2.setCoeff(4, C5NP, LO);
            break;
        default:
            std::stringstream out;
            out << mcdbs2.getOrder();
            throw std::runtime_error("LoopMediatorsAMatching::CMdbs2(): order " + out.str() + "not implemented"); 
    }

    vmcdbs2.push_back(mcdbs2);
    return (vmcdbs2);
}

