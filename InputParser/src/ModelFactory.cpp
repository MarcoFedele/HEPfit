/* 
 * Copyright (C) 2014 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "ModelFactory.h"
#include "Axions.h"
#include "AxionsTHDM.h"
#include "FroggattNielsen.h"
#include <boost/bind.hpp>

ModelFactory::ModelFactory()
{
    modelFactory["StandardModel"] = boost::factory<StandardModel*>();
    modelFactory["Axions"] = boost::factory<Axions*>();
    modelFactory["AxionsTHDM"] = boost::factory<AxionsTHDM*>();
    modelFactory["FroggattNielsen"] = boost::factory<FroggattNielsen*>();
}

void ModelFactory::addModelToFactory(const std::string name, boost::function<StandardModel*() > funct)
{
    modelFactory[name] = funct;
}

StandardModel* ModelFactory::CreateModel(const std::string& name)
{
    if (modelFactory.find(name) == modelFactory.end())
        throw std::runtime_error("ERROR: Wrong model " + name + " passed to ModelFactory.\n");
    return (modelFactory[name]());
}
