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


bool PoroMaterial::formBmatrix (Matrix& Bmat, const Matrix& dNdX, size_t nsd) const
{
  const size_t nenod = dNdX.rows();
  const size_t nstrc = nsd*(nsd+1)/2;
  Bmat.resize(nstrc*nsd,nenod,true);
  if (dNdX.cols() < nsd)
  {
    std::cerr << " *** PoroElasticity::formBmatrix: Invalid dimension on dN1dX, "
              << dNdX.rows() << "x" << dNdX.cols() << "." << std::endl;
    return false;
  }

#define INDEX(i,j) i+nstrc*(j-1)

  switch (nsd) {
    case 1:
    // Strain-displacement matrix for 1D elements
    //
    //  [B] = | d/dx | * [N]

    for (size_t i = 1; i <= nenod; i++)
      Bmat(1,i) = dNdX(i,1);
    break;

    case 2:
    // Strain-displacement matrix for 2D elements
    //
    //         | d/dx   0   |
    //   [B] = |  0    d/dy | * [N]
    //         | d/dy  d/dx |

    for (size_t i = 1; i <= nenod; i++)
    {
      // Normal strain part
      Bmat(INDEX(1,1),i) = dNdX(i,1);
      Bmat(INDEX(2,2),i) = dNdX(i,2);
      // Shear strain part
      Bmat(INDEX(3,1),i) = dNdX(i,2);
      Bmat(INDEX(3,2),i) = dNdX(i,1);
    }
    break;

    case 3:
    // Strain-displacement matrix for 3D elements:
    //
    //         | d/dx   0     0   |
    //         |  0    d/dy   0   |
    //   [B] = |  0     0    d/dz | * [N]
    //         | d/dy  d/dx   0   |
    //         | d/dz   0    d/dx |
    //         |  0    d/dz  d/dy |

    for (size_t i = 1; i <= nenod; i++)
    {
      // Normal strain part
      Bmat(INDEX(1,1),i) = dNdX(i,1);
      Bmat(INDEX(2,2),i) = dNdX(i,2);
      Bmat(INDEX(3,3),i) = dNdX(i,3);
      // Shear strain part
      Bmat(INDEX(4,1),i) = dNdX(i,2);
      Bmat(INDEX(4,2),i) = dNdX(i,1);
      Bmat(INDEX(5,2),i) = dNdX(i,3);
      Bmat(INDEX(5,3),i) = dNdX(i,2);
      Bmat(INDEX(6,1),i) = dNdX(i,3);
      Bmat(INDEX(6,3),i) = dNdX(i,1);
    }
    break;

    default:
      std::cerr <<" *** PoroMaterial::formBmatrix: nsd = " << nsd << std::endl;
      return false;
  }

#undef INDEX

  Bmat.resize(nstrc,nsd*nenod);

  return true;
}


bool PoroMaterial::formElasticMatrix(Matrix& Cmat, const Vec3& X, size_t nsd) const
{
  double E = getStiffness(X);
  double nu = getPoisson(X);

  double C3 = E/(2.0+2.0*nu);
  double C2 = (E*nu)/((1.0+nu)*(1.0-2.0*nu));
  double C1 = C2 + 2.0*C3;

  const size_t nstrc = nsd*(nsd+1)/2;

  Cmat.resize(nstrc,nstrc,true);

  switch (nsd) {
    case 1:
    // Elastic modulus for 1D elements
    Cmat(1,1) = E;
    break;

    case 2: case 3:
    // Elastic matrix for 2D plane-strain and 3D elements
    for (size_t i = 1; i <= nstrc; i++) {
      for (size_t j = 1; j <= nstrc; j++) {
        if (i <= nsd && j <= nsd && i == j)
          Cmat(i,j) = C1;
        else if (i <= nsd && j <= nsd && i != j)
          Cmat(i,j) = C2;
        else if (i > nsd && i == j)
          Cmat(i,j) = C3;
      }
    }
    break;

    default:
      std::cerr <<" *** PoroMaterial::formElasticMatrix: nsd = " << nsd << std::endl;
      return false;
  }

  return true;
}
