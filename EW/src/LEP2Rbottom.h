/* 
 * File:   LEP2Rbottom.h
 * Author: mishima
 */

#ifndef LEP2RBOTTOM_H
#define	LEP2RBOTTOM_H

#include "LEP2ThObservable.h"
#include "LEP2sigmaHadron.h"


class LEP2Rbottom : public LEP2ThObservable {
public:

    /**
     * @brief LEP2Rbottom constructor
     * @param[in] EW_i an object of EW class
     * @param[in] sqrt_s_i the CM energy of the e^+ e^- pair
     */
    LEP2Rbottom(const EW& EW_i, const double sqrt_s_i) : LEP2ThObservable(EW_i, sqrt_s_i),
            myLEP2sigmaHadron(EW_i, sqrt_s_i) {
        q_flavor = StandardModel::BOTTOM;
    }

    /**
     * @return the ratio of the b-bbar cross section to the hadronic cross section at sqrt_s
     */
    double getThValue();

private:
    LEP2sigmaHadron myLEP2sigmaHadron;
    
};

#endif	/* LEP2RBOTTOM_H */

