/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef WILSONCOEFFICIENTNEW_H
#define WILSONCOEFFICIENTNEW_H

#include "gslpp_vector_complex.h"
#include "WilsonTemplateNew.h"

/**
 * @class WilsonCoefficient
 * @ingroup StandardModel
 * @brief A class for the Wilson coefficients. 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class WilsonCoefficientNew : public WilsonTemplateNew<gslpp::vector<gslpp::complex> > {
public:
    WilsonCoefficientNew(unsigned int dim, schemes scheme, orders_qcd order_qcd, orders_qed order_qed = QED0);
    
    Expanded<gslpp::complex> getCoeffElement(uint i) const;
    
    void setCoeff(unsigned int i, gslpp::complex z, orders_qcd order_qcd_i, orders_qed order_qed_i = QED0);
    
    void setCoeff(const gslpp::vector<gslpp::complex>& v, orders_qcd order_qcd_i, orders_qed order_qed_i = QED0);
    
    void resetCoeff();

    gslpp::vector<gslpp::complex> getCoeff(orders_qcd order_qcd_i, orders_qed order_qed_i = QED0) const;
    
    Expanded<gslpp::vector<gslpp::complex> > getCoeff() const;
    
    void setCoeff(const Expanded<gslpp::vector<gslpp::complex> > wc);
};

#endif	/* WILSONCOEFFICIENT_H */
