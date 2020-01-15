/* 
 * Copyright (C) 2018 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "LoopMediatorsDMMatching.h"
#include "LoopMediatorsDM.h"
#include <stdexcept>

LoopMediatorsDMMatching::LoopMediatorsDMMatching(const LoopMediatorsDM & LoopMediatorsDM_i) :

    StandardModelMatching(LoopMediatorsDM_i),
    myLoopMediatorsDM(LoopMediatorsDM_i),
    mcBMll(13, NDR, NLO)
{}

void LoopMediatorsDMMatching::updateLoopMediatorsDMParameters()
{   
    
    C9NPmu = myLoopMediatorsDM.getC9();
    C10NPmu = myLoopMediatorsDM.getC10();

    StandardModelMatching::updateSMParameters();
}

LoopMediatorsDMMatching::~LoopMediatorsDMMatching()
{}


std::vector<WilsonCoefficient>& LoopMediatorsDMMatching::CMBMll(QCD::lepton lepton)
{
    vmcBMll.clear();
    for (std::vector<WilsonCoefficient>::iterator it = StandardModelMatching::CMBMll(lepton).begin(); it != StandardModelMatching::CMBMll(lepton).end(); it++ ) vmcBMll.push_back(*it);

    switch (mcBMll.getScheme()) {
        case NDR:

            break;
        default:
            std::stringstream out;
            out << mcBMll.getScheme();
            throw std::runtime_error("LoopMediatorsDMMatching::CMBMll(): scheme " + out.str() + "not implemented"); 
    }

    mcBMll.setMu(1000.);

    switch (mcBMll.getOrder()) {
        case NNLO:
        case NLO:
                mcBMll.setCoeff(6, 0., NLO);
                mcBMll.setCoeff(8, 0., NLO);
                mcBMll.setCoeff(9, 0., NLO);
                mcBMll.setCoeff(10, 0., NLO);
                mcBMll.setCoeff(11, 0., NLO);
        case LO:
            mcBMll.setCoeff(6, 0., LO);
            if(lepton == LoopMediatorsDM::MU){
                mcBMll.setCoeff(8, C9NPmu, LO);
                mcBMll.setCoeff(9, C10NPmu, LO);
                mcBMll.setCoeff(10, 0., LO);
                mcBMll.setCoeff(11, 0., LO);               
            }
            break;
        default:
            std::stringstream out;
            out << mcBMll.getOrder();
            throw std::runtime_error("LoopMediatorsDMMatching::CMBMll(): order " + out.str() + "not implemented"); 
    }

    vmcBMll.push_back(mcBMll);
    return (vmcBMll);
}
