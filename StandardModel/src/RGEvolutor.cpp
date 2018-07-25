/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "RGEvolutor.h"

RGEvolutor::RGEvolutor(unsigned int dim, schemes scheme, orders_qcd order_qcd_i, orders_qed order_qed_i)
: WilsonTemplate<gslpp::matrix<double> >(dim, scheme, order_qcd_i, order_qed_i)
{}


void RGEvolutor::setEvol(unsigned int i, unsigned int  j, double x, orders_qcd order_qcd_i, orders_qed order_qed_i) 
{    
    if (i > size || j > size) {
        std::stringstream out;
        out << i << " " << j;
        throw std::runtime_error("RGEvolutor::setEvol(): matrix indices " + out.str() + " out of range"); 
    }
    if (order_qcd_i > order_qcd || order_qed_i > order_qed) {
        std::stringstream out;
        out << order_qcd_i << " and " << order_qed_i;
        throw std::runtime_error("RGEvolutor::setEvol(): order " + out.str() +" not implemented "); 
    }
    gslpp::matrix<double> tmp = wilson.getOrd(order_qcd_i, order_qed_i);
    tmp.assign(i, j, x);
    wilson.setOrd(order_qcd_i, order_qed_i, tmp);
}

void RGEvolutor::setEvol(const gslpp::matrix<double>& m, orders_qcd order_qcd_i, orders_qed order_qed_i)
{
    setWilson(m, order_qcd_i, order_qed_i);
}

Expanded<gslpp::matrix<double> >& RGEvolutor::getEvol() const
{
    return wilson;
}

double RGEvolutor::getM() const
{
    return M;
}

void RGEvolutor::setScales(double mu, double M)
{
    this->M = M;
    WilsonTemplate::setMu(mu);
    if(order_qed == NOQED)
      setEvol(gslpp::matrix<double>::Id(size), LO, NOQED);
    else
      setEvol(gslpp::matrix<double>::Id(size), LO, QED0);
}

void RGEvolutor::setM(double M)
{
    setScales(mu, M);
}

void RGEvolutor::setMu(double mu)
{
    setScales(mu, M);
}

gslpp::matrix<double> RGEvolutor::Evol(orders_qcd order_qcd_i, orders_qed order_qed_i)
{
    return getWilson(order_qcd_i, order_qed_i);
}
