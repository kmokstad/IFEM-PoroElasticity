// $Id$
//==============================================================================
//!
//! \file PoroElasticity.h
//!
//! \date April 16 2015
//!
//! \author Yared Bekele
//!
//! \brief Integrand implementations for time-dependent PoroElasticity problems
//!
//==============================================================================

#ifndef POROELASTICITY_H_
#define POROELASTICITY_H_

#include "BDF.h"
#include "Vec3.h"
#include "ElmMats.h"
#include "IntegrandBase.h"
#include "PoroMaterial.h"


/*!
 * \brief Class representing the integrand of the PoroElasticity problem.
*/

class PoroElasticity : public IntegrandBase
{
  /*!
   * \brief Class representing an element matrix for the PoroElasticity problem
   */
  class MixedElmMats : public ElmMats
  {
  public:
    //! \brief Default constructor
    MixedElmMats();
    //! \brief Empty destructor
    virtual ~MixedElmMats() {}
    //! \brief Returns the element level Newton matrix
    virtual const Matrix& getNewtonMatrix() const;
    //! \brief Returns the element level RHS vector
    virtual const Vector& getRHSVector() const;
  };

public:

  //! \brief The default constructor initializes all pointers to zero.
  //! \param[in] n Number of spatial dimensions
  PoroElasticity(unsigned short int n = 2, int order = 1);

  //! \brief The destructor frees the dynamically allocated data objects.
  virtual ~PoroElasticity() {}

  //! \brief Defines the traction field to use in Neumann boundary conditions
  void setTraction(TractionFunc* tf) { tracFld = tf; }

  //! \brief Defines the trction field to use in Neumann boundary conditions
  void setTraction(VecFunc* tf) { fluxFld = tf; }

  //! \brief Evaluates the boundary traction field (if any) at specified point
  //! \param[in] X Cartesian coordinate of the current integration point
  //! \param[in] n Outward-directed unit normal vector at current point
  virtual Vec3 getTraction(const Vec3& X, const Vec3& n) const;

  //! \brief Defines the gravitation vector
  //! \param[in] grav Gravity vector
  virtual void setGravity(const Vec3& gravity) { grav = gravity; }

  //! \brief Obtain current gravity vector
  const Vec3 getGravity() const { return grav; }

  //! \brief Defines the scaling factor
  //! \param[in] sc Scaling factor
  virtual void setScaling(double scaling) { sc = scaling; }

  //! \brief Obtain current scaling factor
  double getScaling() const { return sc; }

  //! \brief Defines the material properties.
  virtual void setMaterial(PoroMaterial* material) { mat = material; }

  //! Evaluates the mass density at current point
  virtual double getMassDensity(const Vec3&) const;

  //! \brief Evaluates the body force field (if any) at a specified point
  virtual Vec3 getBodyForce(const Vec3& X) const;

  using IntegrandBase::getLocalIntegral;
  //! \brief Returns a local integral container for the given element
  //! \param[in] nen1 Number of nodes on element for basis 1
  //! \param[in] nen2 Number of nodes on element for basis 2
  //! \param[in] neumann Whether or not we are assembling Neumann BCs
  virtual LocalIntegral* getLocalIntegral(size_t nen1, size_t nen2,
                                          size_t, bool neumann) const;

  //! \brief Initializes current element for numerical integration
  //! \param[in] MNPC1 Nodal point correspondence for basis 1
  //! \param[in] MNPC2 Nodal point correspondence for basis 2
  //! \param[in] n1 Number of nodes in basis 1 on this patch
  //! \param elmInt The local integral object for current element
  virtual bool initElement(const std::vector<int>& MNPC1,
                           const std::vector<int>& MNPC2,
                           size_t n1, LocalIntegral& elmInt);

  //! \brief Initializes current element for numerical boundary integration (mixed)
  //! \param[in] MNPC1 Nodal point correspondence for basis 1
  //! \param[in] MNPC2 Nodal point correspondence for basis 2
  //! \param elmInt The local integral object for current element
  virtual bool initElementBou(const std::vector<int>& MNPC1,
                              const std::vector<int>& MNPC2,
                              size_t, LocalIntegral& elmInt);

  //! \brief Calculates the strain-displacement matrix.
  //! \param[in] Bmat The strain-displacement matrix
  //! \param[in] dNdX First basis function gradients at current point
  bool formBmatrix(Matrix& Bmat, const Matrix& dNdX) const;

  //! \brief Evalutates the constitutive matrix at an integration point
  //! \param[out] Cmat Constitutive matrix at current point
  //! \param[in] E Young's modulus
  //! \param[in] nu Poisson's ratio
  bool formElasticMatrix(Matrix& Cmat, const Vec3& X) const;

  //! \brief Evaluates the integrand at an interior point
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] time Parameters for nonlinear and time-dependent simulations
  //! \param[in] X Cartesian coordinates of current integration point
  virtual bool evalIntMx(LocalIntegral& elmInt, const MxFiniteElement& fe,
                         const TimeDomain& time, const Vec3& X) const;

  //! \brief Evaluates the integrand at a boundary point
  //! \param elmInt The local interal object to receive the contributions
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] time Parameters for nonlinear and time-dependent simulations
  //! \param[in] X Cartesian coordinates of current integration point
  //! \param[in] normal Boundary normal vector at current integration point
  virtual bool evalBouMx(LocalIntegral& elmInt, const MxFiniteElement& fe,
                         const TimeDomain& time, const Vec3& X,
                         const Vec3& normal) const;

  //! \brief Evaluates the secondary solution at a result point (mixed problem).
  //! \param[out] s The solution field values at current point
  //! \param[in] fe Mixed finite element data at current point
  //! \param[in] X Cartesian coordinates of current point
  //! \param[in] MNPC1 Nodal point correspondance for the basis 1
  //! \param[in] MNPC2 Nodal point correspondance for the basis 2
  virtual bool evalSol(Vector& s, const MxFiniteElement& fe, const Vec3& X,
                       const std::vector<int>& MNPC1,
                       const std::vector<int>& MNPC2) const;

  //! \brief Finalizes the element quantities after the numerical integration
  //! \details This method is invoked once for each element, after the numerical
  //! integration loop over interior points is finished and before the resulting
  //! element quantities are assembled into their system level equivalents
  virtual bool finalizeElement(LocalIntegral&, const TimeDomain&, size_t);

  //! \brief Returns whether a mixed formulation is used
  virtual bool mixedFormulation() const { return true; }

  //! \brief Returns the number of primary/secondary solution field components
  //! \param[in] fld Which field set to consider (1=primary,2=secondary)
  virtual size_t getNoFields(int fld = 1) const;

  //! \brief Returns the name of a primary solution field component
  //! \param[in] i Field component index
  //! \param[in] prefix Name prefix for all components
  virtual const char* getField1Name(size_t i, const char* prefix = 0) const;

  //! \brief Returns the name of a secondary solution field component
  //! \param[in] i Field component index
  //! \param[in] prefix Name prefix for all components
  virtual const char* getField2Name(size_t i, const char* prefix = 0) const;

private:
  Vec3 grav;              //!< Gravitation vector
  double sc;              //!< Scaling factor

protected:
  unsigned short int nsd;                       //!< Number of space dimensions (2)
  unsigned short int eS;                        //!< Index to element load vector
  double gacc;                                  //!< Gravitational acceleration
  TractionFunc* tracFld;                        //!< Pointer to implicit boundary traction field
  VecFunc* fluxFld;                             //!< Pointer to explicit boundary traction field
  TimeIntegration::BDF bdf;                     //!< BDF time discretization parameters
  PoroMaterial* mat;                            //!< Material data
};

#endif
