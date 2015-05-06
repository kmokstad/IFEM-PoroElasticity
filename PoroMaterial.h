// $Id$
//==============================================================================
//!
//! \file PoroMaterial.h
//!
//! \date Apr 29 2015
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Class for poro-elastic material models.
//!
//==============================================================================

#ifndef _PORO_MATERIAL_H
#define _PORO_MATERIAL_H

#include "Function.h"
#include "Vec3.h"
#include "Vec3Oper.h"

class TiXmlElement;


/*!
  \brief Class representing a material model for a poroelastic problem.
*/

class PoroMaterial
{
public:
  /*! \brief Helper template for wrapping a constant/function pair */
    template<class Function>
  struct FuncConstPair
  {
    Function* function;                 //!< Function definition
    typename Function::Output constant; //!< Constant

    //! \brief Constructor.
    FuncConstPair() { function = nullptr; constant = 0.0; }

    //! \brief Parse an XML element. Specialized per type
    Function* parse(const char* val, const std::string& type) { return nullptr; }

    //! \brief Evaluate function.
    //! \param[in] X the value to evaluate at.
    typename Function::Output evaluate(const typename Function::Input& X) const
    {
      return function ? (*function)(X) : constant;
    }
  };

  //! \brief Empty constructor.
  PoroMaterial() {}

  //! \brief Empty destructor.
  ~PoroMaterial() {}

  //! \brief Parses material parementers from an XML element.
  void parse(const TiXmlElement*);

  //! \brief Prints out material parameters to the log stream.
  void printLog() const;

  //! \brief Evaluates the mass density of the fluid at current point.
  double getFluidDensity(const Vec3&) const;
  //! \brief Evaluates the mass density of the solid at current point.
  double getSolidDensity(const Vec3&) const;
  //! \brief Evaluates the thermal expansion coefficient for given temperature.
  double getThermalExpansion(double) const;
  //! \brief Evaluates the heat capacity at the current point
  double getHeatCapacity(double T) const;
  //! \brief Evaluates the thermal conductivity at the current point
  double getThermalConductivity(double T) const;
  //! \brief Returns porosity at the current point.
  double getPorosity(const Vec3& X) const;
  //! \brief Returns permeability at the current point.
  Vec3 getPermeability(const Vec3& X) const;
  //! \brief Returns bulk modulus of the water at the current point.
  double getBulkWater(const Vec3& X) const;
  //! \brief Returns bulk modulus of the solid at the current point.
  double getBulkSolid(const Vec3& X) const;
  //! \brief Returns bulk modulus of the medium at the current point.
  double getBulkMedium(const Vec3& X) const;
  //! \brief Returns stiffness at the current point.
  double getStiffness(const Vec3& X) const;
  //! \brief Returns Poisson's ratio at the current point.
  double getPoisson(const Vec3& X) const;
protected:
  FuncConstPair<RealFunc> Emod; //!< Young's modulus
  FuncConstPair<RealFunc> nu;   //!< Poisson's ratio
  FuncConstPair<RealFunc> rhof; //!< Fluid density
  FuncConstPair<RealFunc> rhos; //!< Solid density

  FuncConstPair<ScalarFunc> thermalexpansion; //!< Thermal expansion coefficient
  FuncConstPair<ScalarFunc> heatcapacity;     //!< Specific heat capacity
  FuncConstPair<ScalarFunc> conductivity;     //!< Thermal conductivity

  FuncConstPair<RealFunc> porosity;     //!< Porosity
  FuncConstPair<VecFunc>  permeability; //!< Permeability
  FuncConstPair<RealFunc> bulkw;        //!< Bulk modulus of water
  FuncConstPair<RealFunc> bulks;        //!< Bulk modulus of solid
  FuncConstPair<RealFunc> bulkm;        //!< Bulk modulus of medium
};

#endif
