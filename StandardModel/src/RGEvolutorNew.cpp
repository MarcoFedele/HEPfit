/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "RGEvolutorNew.h"

RGEvolutorNew::RGEvolutorNew(unsigned int dim, schemes scheme, orders_qcd order_qcd_i, orders_qed order_qed_i)
: WilsonTemplateNew<gslpp::matrix<double> >(dim, scheme, order_qcd_i, order_qed_i)
{}


void RGEvolutorNew::setEvol(unsigned int i, unsigned int  j, double x, orders_qcd order_qcd_i, orders_qed order_qed_i) 
{    
    if (i > size || j > size) {
        std::stringstream out;
        out << i << " " << j;
        throw std::runtime_error("RGEvolutorNew::setEvol(): matrix indices " + out.str() + " out of range"); 
    }
    if (order_qcd_i > order_qcd || order_qed_i > order_qed) {
        std::stringstream out;
        out << order_qcd_i << " and " << order_qed_i;
        throw std::runtime_error("RGEvolutorNew::setEvol(): order " + out.str() +" not implemented "); 
    }
    gslpp::matrix<double> tmp = wilson.getOrd(order_qcd_i, order_qed_i);
    tmp.assign(i, j, x);
    wilson.setOrd(order_qcd_i, order_qed_i, tmp);
}

void RGEvolutorNew::setEvol(const gslpp::matrix<double>& m, orders_qcd order_qcd_i, orders_qed order_qed_i)
{
    setWilson(m, order_qcd_i, order_qed_i);
}

const Expanded<gslpp::matrix<double> >& RGEvolutorNew::getEvol() const
{
    return wilson;
}

double RGEvolutorNew::getM() const
{
    return M;
}

void RGEvolutorNew::setScales(double mu, double M)
{
    this->M = M;
    WilsonTemplateNew::setMu(mu);
    setEvol(gslpp::matrix<double>::Id(size), LO, QED0);
}

void RGEvolutorNew::setM(double M)
{
    setScales(mu, M);
}

void RGEvolutorNew::setMu(double mu)
{
    setScales(mu, M);
}

gslpp::matrix<double> RGEvolutorNew::Evol(orders_qcd order_qcd_i, orders_qed order_qed_i)
{
    return getWilson(order_qcd_i, order_qed_i);
}
