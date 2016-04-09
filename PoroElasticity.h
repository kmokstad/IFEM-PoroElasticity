// $Id$
//==============================================================================
//!
//! \file PoroElasticity.h
//!
//! \date April 16 2015
//!
//! \author Yared Bekele
//!
//! \brief Integrand implementations for time-dependent PoroElasticity problems.
//!
//==============================================================================

#ifndef _PORO_ELASTICITY_H_
#define _PORO_ELASTICITY_H_

#include "Elasticity.h"
#include "ElmMats.h"


/*!
  \brief Class representing the integrand of the PoroElasticity problem.
*/

class PoroElasticity : public Elasticity
{
  /*!
   * \brief Class representing an element matrix for the mixed PoroElasticity problem
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
    //! \brief Makes the actual Newton matrix.
    //! \note Separated for reuse in finalizeElement.
    void makeNewtonMatrix(Matrix& N, size_t pp_idx) const;
  };

  /*!
   * \brief Class representing an element matrix for the non-mixed PoroElasticity problem
   */
  class NonMixedElmMats : public ElmMats
  {
  public:
    //! \brief Default constructor
    NonMixedElmMats();
    //! \brief Empty destructor
    virtual ~NonMixedElmMats() {}
    //! \brief Returns the element level Newton matrix
    virtual const Matrix& getNewtonMatrix() const;
    //! \brief Returns the element level RHS vector
    virtual const Vector& getRHSVector() const;
    //! \brief Makes the actual Newton matrix.
    //! \note Separated for reuse in finalizeElement.
    void makeNewtonMatrix(Matrix& N, size_t pp_idx) const;
  };

public:
  //! \brief The default constructor initializes all pointers to zero.
  //! \param[in] n Number of spatial dimensions
  //! \param[in] order Order of the time-integration scheme
  PoroElasticity(unsigned short int n = 3, int order = 1);
  //! \brief Empty destructor.
  virtual ~PoroElasticity() {}

  //! \brief Defines the scaling factor.
  void setScaling(double scaling) { sc = scaling; }
  //! \brief Obtain current scaling factor.
  double getScaling() const { return sc; }

  //! \brief Returns the current gravity vector.
  const Vec3 getGravity() const { return gravity; }

  //! \brief Computes the stiffness matrix for a quadrature point.
  bool evalStiffnessMatrix(Matrix& mx, const Matrix &B, const Matrix &C, double detJxW) const;

  //! \brief Computes the coupling matrix for a quadrature point.
  bool evalCouplingMatrix(Matrix &mx, const Matrix &B, const Vector &basis,
                          double scl, double alpha, const Vector &m, double detJxW) const;

  //! \brief Computes the compressibility matrix for a quadrature point.
  bool evalCompressibilityMatrix(Matrix &mx, const Vector &basis,
                                 double scl, double Minv, double detJxW) const;

  //! \brief Computes the permeability matrix for a quadrature point.
  bool evalPermeabilityMatrix(Matrix &mx, const Matrix &grad, double scl,
                              const Vec3 &permeability, double acc_dens, double detJxW) const;

  using IntegrandBase::getLocalIntegral;
  //! \brief Returns a local integral container for the given element
  //! \param[in] nen1 Number of nodes on element for basis 1
  //! \param[in] nen2 Number of nodes on element for basis 2
  //! \param[in] neumann Whether or not we are assembling Neumann BCs
  virtual LocalIntegral* getLocalIntegral(const std::vector<size_t>& nen,
                                          size_t, bool neumann) const;
  //! \brief Returns a local integral contribution object for the given element.
  //! \param[in] nen Number of nodes on element
  //! \param[in] iEl Global element number (1-based)
  //! \param[in] neumann Whether or not we are assembling Neumann BCs
  virtual LocalIntegral* getLocalIntegral(size_t nen, size_t iEl,
                                          bool neumann = false) const;

  //! \brief Initializes current element for numerical integration
  //! \param[in] MNPC1 Nodal point correspondence for basis 1
  //! \param[in] MNPC2 Nodal point correspondence for basis 2
  //! \param[in] n1 Number of nodes in basis 1 on this patch
  //! \param elmInt The local integral object for current element
  virtual bool initElement(const std::vector<int>& MNPC,
                           const std::vector<size_t>& elem_sizes,
                           const std::vector<size_t>& basis_sizes,
                           LocalIntegral& elmInt);

  //! \brief Initializes current element for numerical integration.
  //! \param[in] MNPC Matrix of nodal point correspondance for current element
  //! \param elmInt Local integral for element
  //!
  //! \details This method is invoked once before starting the numerical
  //! integration loop over the Gaussian quadrature points over an element.
  //! It is supposed to perform all the necessary internal initializations
  //! needed before the numerical integration is started for current element.
  virtual bool initElement(const std::vector<int>& MNPC, LocalIntegral& elmInt);

  //! \brief Initializes current element for numerical boundary integration (mixed)
  //! \param[in] MNPC1 Nodal point correspondence for basis 1
  //! \param[in] MNPC2 Nodal point correspondence for basis 2
  //! \param elmInt The local integral object for current element
  virtual bool initElementBou(const std::vector<int>& MNPC,
                              const std::vector<size_t>& elem_sizes,
                              const std::vector<size_t>& basis_sizes,
                              LocalIntegral& elmInt);

  //! \brief Initializes current element for boundary integration.
  //! \param[in] MNPC Matrix of nodal point correspondance for current element
  //! \param elmInt Local integral for element
  virtual bool initElementBou(const std::vector<int>& MNPC,
                              LocalIntegral& elmInt);
  //! \brief Evaluates the integrand at an interior point
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] time Parameters for nonlinear and time-dependent simulations
  //! \param[in] X Cartesian coordinates of current integration point
  virtual bool evalIntMx(LocalIntegral& elmInt, const MxFiniteElement& fe,
                         const TimeDomain& time, const Vec3& X) const;

  //! \brief Evaluates the integrand at an interior point
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] time Parameters for nonlinear and time-dependent simulations
  //! \param[in] X Cartesian coordinates of current integration point
  virtual bool evalInt(LocalIntegral& elmInt, const FiniteElement& fe,
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

  //! \brief Evaluates the integrand at a boundary point
  //! \param elmInt The local interal object to receive the contributions
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] time Parameters for nonlinear and time-dependent simulations
  //! \param[in] X Cartesian coordinates of current integration point
  //! \param[in] normal Boundary normal vector at current integration point
  virtual bool evalBou(LocalIntegral& elmInt, const FiniteElement& fe,
                       const TimeDomain& time, const Vec3& X,
                       const Vec3& normal) const;

  using IntegrandBase::evalSol;
  //! \brief Evaluates the secondary solution at a result point.
  //! \param[out] s The solution field values at current point
  //! \param[in] fe Finite element data at current point
  //! \param[in] X Cartesian coordinates of current point
  //! \param[in] MNPC Nodal point correspondance for the basis function values
  virtual bool evalSol(Vector& s, const FiniteElement& fe, const Vec3& X,
		       const std::vector<int>& MNPC) const;
  //! \brief Evaluates the secondary solution at a result point (mixed problem).
  //! \param[out] s The solution field values at current point
  //! \param[in] fe Mixed finite element data at current point
  //! \param[in] X Cartesian coordinates of current point
  //! \param[in] MNPC Nodal point correspondance for the bases
  //! \param[in] elem_sizes Size of each basis on the element
  virtual bool evalSol(Vector& s, const MxFiniteElement& fe, const Vec3& X,
                       const std::vector<int>& MNPC,
                       const std::vector<size_t>& elem_sizes) const;

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
  virtual std::string getField1Name(size_t i, const char* prefix = 0) const;

  //! \brief Returns the name of a secondary solution field component
  //! \param[in] i Field component index
  //! \param[in] prefix Name prefix for all components
  virtual std::string getField2Name(size_t i, const char* prefix = 0) const;

private:
  //! \brief Evaluates the secondary solution at a result point
  //! (shared code between mixed and non-mixed)
  //! \param[out] s The solution field values at current point
  //! \param[in] fe Finite element data at current point
  //! \param[in] X Cartesian coordinates of current point
  //! \param[in] disp The displacement coefficients
  bool evalSolCommon(Vector& s,
                     const FiniteElement& fe, const Vec3& X,
                     const Vector& disp) const;

  //! \brief Initializes the local integral
  void initLocalIntegral(ElmMats *result, size_t ndof_displ,
                         size_t ndof_press, bool neumann) const;

private:
  double sc;   //!< Scaling factor
  double gacc; //!< Gravitational acceleration
};

#endif
