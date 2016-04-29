// $Id$
//==============================================================================
//!
//! \file PoroSolutions.C
//!
//! \date Jun 22 2015
//!
//! \author Arne Morten Kvarving
//!
//! \brief Analytic solutions for PoroElasticity problems.
//!
//==============================================================================

#include "PoroSolutions.h"
#include "PoroMaterial.h"
#include "Vec3.h"


double TerzhagiPressure::evaluate (const Vec3& X) const
{
  double y = X.y;

  const Vec4& Xt = static_cast<const Vec4&>(X);

  double k = myMat->getPermeability(X)[0];
  double rhof = myMat->getFluidDensity(X);
  double E = myMat->getStiffness(X);
  double nu = myMat->getPoisson(X);
  double n = myMat->getPorosity(X);
  double alpha = myMat->getBiotCoeff(X);
  double Minv = myMat->getBiotModulus(X,alpha,n);

  double K = E / 3 / (1 - 2*nu);
  double G = E / 2 / (1 + nu);
  double mv = 1 / (K + 4*G/3);
  double cv = k / rhof / gacc / (alpha * alpha * mv + Minv);

  double stime = cv * Xt.t / height / height;

  double ret = 0.0;
  for (int i = 1; true; i++) {
    double trig_arg = (2*i - 1) * M_PI * y / 2 / height;
    double exp_c = exp(-(2*i-1) * (2*i-1) * M_PI * M_PI * stime / 4);
    double sign = i % 2 == 0 ? -1 : 1;
    double p_fac = 4 / (2.0*i - 1) / M_PI * load;

    double p_add = cos(trig_arg) * exp_c * p_fac * sign;
    ret += p_add;

    if (fabs(p_add) < 1e-15)
      break;
  }

  return ret;
}


Vec3 StationaryTerzhagiDisplacement::evaluate (const Vec3& X) const
{
  double y = X.y;

  double E = myMat->getStiffness(X);
  double nu = myMat->getPoisson(X);
  double c = -load/E * (1+nu) * (1-2*nu) / (1-nu);

  return Vec3(0, c*y, 0);
}
