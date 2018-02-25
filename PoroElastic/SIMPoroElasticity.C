// $Id$
//==============================================================================
//!
//! \file SIMPoroElasticity.C
//!
//! \date May 29 2016
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Simulation driver for poroelasticity problems.
//!
//==============================================================================

#include "SIMPoroElasticity.h"
#include "PoroSolutions.h"
#include "PoroMaterial.h"


template<> bool SIMPoroEl2D::parseDimSpecific (const TiXmlElement* child)
{
  if (strcasecmp(child->Value(),"anasol"))
    return false;
  else if (SIM2D::mySol)
    return true;

  std::string type;
  utl::getAttribute(child,"type",type,true);
  if (type == "terzhagi")
  {
    double height = 1.0, load = 1.0;
    utl::getAttribute(child,"height",height);
    utl::getAttribute(child,"load",load);
    IFEM::cout <<"\tAnalytical solution: Terzhagi, height = "<< height
               <<" load = "<< load << std::endl;
    Elasticity* elp = this->getIntegrand();
    PoroMaterial* pmat = static_cast<PoroMaterial*>(elp->getMaterial());
    double gacc = elp->getGravity().length();
    SIM2D::mySol = new AnaSol(new TerzhagiPressure(pmat,gacc,height,load));
  }
  else if (type == "terzhagi-stationary")
  {
    double load = 1.0;
    utl::getAttribute(child,"load",load);
    IFEM::cout <<"\tAnalytical solution: Terzhagi (stationary),"
               <<" load = "<< load << std::endl;
    Elasticity* elp = this->getIntegrand();
    PoroMaterial* pmat = static_cast<PoroMaterial*>(elp->getMaterial());
    SIM2D::mySol = new AnaSol(new ConstFunc(0.0), nullptr,
                              new StationaryTerzhagiDisplacement(pmat,load));
  }
  else
  {
    IFEM::cout <<"\tAnalytical solution: Expression"<< std::endl;
    SIM2D::mySol = new AnaSol(child);
  }

  return true;
}


template<> bool SIMPoroEl3D::parseDimSpecific (const TiXmlElement* child)
{
  if (strcasecmp(child->Value(),"anasol"))
    return false;
  else if (SIM3D::mySol)
    return true;

  IFEM::cout <<"\tAnalytical solution: Expression"<< std::endl;
  SIM3D::mySol = new AnaSol(child);

  return true;
}
