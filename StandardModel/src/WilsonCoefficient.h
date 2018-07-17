/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef WILSONCOEFFICIENT_H
#define WILSONCOEFFICIENT_H

#include "gslpp_vector_complex.h"
#include "WilsonTemplate.h"

/**
 * @class WilsonCoefficient
 * @ingroup StandardModel
 * @brief A class for the Wilson coefficients. 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class WilsonCoefficient : public WilsonTemplate<gslpp::complex> {
public:
    WilsonCoefficient(unsigned int dim, schemes scheme, orders_qcd order_qcd);
    
    WilsonCoefficient(unsigned int dim, schemes scheme, orders_qcd order_qcd, orders_qed order_qed);
};

#endif	/* WILSONCOEFFICIENT_H */
