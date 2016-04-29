// $Id$
//==============================================================================
//!
//! \file PoroSolutions.h
//!
//! \date Jun 22 2015
//!
//! \author Arne Morten Kvarving
//!
//! \brief Analytic solutions for PoroElasticity problems.
//!
//==============================================================================

#ifndef _PORO_SOLUTIONS_H
#define _PORO_SOLUTIONS_H

#include "Function.h"

class PoroMaterial;


class TerzhagiPressure : public RealFunc
{
public:
  //! \brief The constructor initializes the data members.
  TerzhagiPressure(const PoroMaterial* mat, double g, double h, double l)
    : myMat(mat), gacc(g), height(h), load(l) {}
  //! \brief Empty destructor.
  virtual ~TerzhagiPressure() {}

protected:
  //! \brief Evaluates the analytic pressure field at the point \a X.
  virtual double evaluate(const Vec3& X) const;

private:
  const PoroMaterial* myMat; //!< Pointer to material parameters object
  double gacc, height, load;
};


class StationaryTerzhagiDisplacement : public VecFunc
{
public:
  //! \brief The constructor initializes the data members.
  StationaryTerzhagiDisplacement(const PoroMaterial* mat, double l)
    : myMat(mat), load(l) {}
  //! \brief Empty destructor
  virtual ~StationaryTerzhagiDisplacement() {}

protected:
  //! \brief Evaluates the analytic displacement field at the point \a X.
  virtual Vec3 evaluate(const Vec3& X) const;

private:
  const PoroMaterial* myMat; //!< Pointer to material parameters object
  double load;
};

#endif
