// $Id$
//==============================================================================
//!
//! \file PoroMaterial.C
//!
//! \date Apr 29 2015
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Class for poro-elastic material models.
//!
//==============================================================================


#include "PoroMaterial.h"

#include "Functions.h"
#include "IFEM.h"
#include "tinyxml.h"
#include "Utilities.h"
#include "Vec3.h"


template<>
RealFunc* PoroMaterial::FuncConstPair<RealFunc>::parse(const char* val,
                                                       const std::string& type)
{
  return utl::parseRealFunc(val, type);
}


template<>
ScalarFunc* PoroMaterial::FuncConstPair<ScalarFunc>::parse(const char* val,
                                                           const std::string& type)
{
  return utl::parseTimeFunc(val, type);
}


template<>
VecFunc* PoroMaterial::FuncConstPair<VecFunc>::parse(const char* val,
                                                     const std::string& type)
{
  return utl::parseVecFunc(val, type);
}


  template<class T>
static bool propertyParse(PoroMaterial::FuncConstPair<T>& data,
                          const TiXmlElement* elem,
                          const std::string& attr,
                          const std::string& tag)
{
  std::string constant;
  if (utl::getAttribute(elem,attr.c_str(),constant)) {
    std::stringstream str;
    str << constant;
    str >> data.constant;
    return true;
  }

  const TiXmlElement* child = elem->FirstChildElement(tag);
  if (child) {
    IFEM::cout <<" ";
    std::string type;
    utl::getAttribute(child,"type",type,true);
    const TiXmlNode* aval;
    if ((aval = child->FirstChild()))
      data.function = data.parse(aval->Value(),type);

    return data.function != nullptr;
  }

 return false;
}


void PoroMaterial::parse(const TiXmlElement* elem)
{
  propertyParse(Emod, elem, "E", "stiffness");
  propertyParse(nu, elem, "nu", "poisson");
  propertyParse(rhof, elem, "rhof", "fluiddensity");
  propertyParse(rhos, elem, "rhos", "soliddensity");

  propertyParse(thermalexpansion, elem, "alpha", "thermalexpansion");
  propertyParse(heatcapacity, elem, "cp", "heatcapacity");
  propertyParse(conductivity, elem, "kappa", "conductivity");

  propertyParse(porosity, elem, "poro", "porosity");
  propertyParse(permeability, elem, "perm", "permeability");
  propertyParse(bulkw, elem, "Kw", "waterbulk");
  propertyParse(bulks, elem, "Ks", "solidbulk");
  propertyParse(bulkm, elem, "Ko", "mediumbulk");
}


void PoroMaterial::printLog() const
{
  IFEM::cout << "\tConstitutive Properties: "
             << "\n\t\tYoung's Modulus, E = " << Emod.constant
             << "\n\t\tPoisson's Ratio, nu = " << nu.constant << std::endl;
  IFEM::cout << "\tDensities: "
             << "\n\t\tDensity of Fluid, rhof = " << rhof.constant
             << "\n\t\tDensity of Solid, rhos = " << rhos.constant << std::endl;
  IFEM::cout << "\tBulk Moduli: "
             << "\n\t\tBulk Modulus of Water, Kw = " << bulkw.constant
             << "\n\t\tBulk Modulus of Solid, Ks = " << bulks.constant
             << "\n\t\tBulk Modulus of Medium, Ko = " << bulkm.constant << std::endl;
  IFEM::cout <<"\tPorosity, n = " << porosity.constant << std::endl;
}


double PoroMaterial::getThermalExpansion (double T) const
{
  return thermalexpansion.evaluate(T);
}


double PoroMaterial::getHeatCapacity (double T) const
{
  return heatcapacity.evaluate(T);
}


double PoroMaterial::getThermalConductivity(double T) const
{
  return conductivity.evaluate(T);
}


double PoroMaterial::getPorosity(const Vec3& X) const
{
  return porosity.evaluate(X);
}


Vec3 PoroMaterial::getPermeability(const Vec3& X) const
{
  return permeability.evaluate(X);
}


double PoroMaterial::getFluidDensity(const Vec3& X) const
{
  return rhof.evaluate(X);
}


double PoroMaterial::getSolidDensity(const Vec3& X) const
{
  return rhos.evaluate(X);
}


double PoroMaterial::getBulkWater(const Vec3& X) const
{
  return bulkw.evaluate(X);
}


double PoroMaterial::getBulkSolid(const Vec3& X) const
{
  return bulks.evaluate(X);
}


double PoroMaterial::getBulkMedium(const Vec3& X) const
{
  return bulkm.evaluate(X);
}


double PoroMaterial::getStiffness(const Vec3& X) const
{
  return Emod.evaluate(X);
}


double PoroMaterial::getPoisson(const Vec3& X) const
{
  return nu.evaluate(X);
}
