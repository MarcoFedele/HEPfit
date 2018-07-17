/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "WilsonCoefficient.h"
#include <sstream>
#include <stdexcept>

WilsonCoefficient::WilsonCoefficient(unsigned int dim, schemes scheme, orders_qcd order_qcd_i)
: WilsonTemplate<gslpp::complex>(dim, scheme, order_qcd_i)
{};

WilsonCoefficient::WilsonCoefficient(unsigned int dim, schemes scheme, orders_qcd order_qcd_i, orders_qed order_qed_i)
: WilsonTemplate<gslpp::complex>(dim, scheme, order_qcd_i, order_qed_i)
{};
