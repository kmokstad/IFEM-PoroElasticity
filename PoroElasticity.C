// $Id$
//==============================================================================
//!
//! \file PoroElasticity.C
//!
//! \date April 16 2015
//!
//! \author Yared Bekele
//!
//! \brief Integrand implementations for time-dependent poroelasticity problems
//!
//==============================================================================

#include "PoroElasticity.h"
#include "ASMbase.h"
#include "ASMmxBase.h"
#include "FiniteElement.h"
#include "TimeDomain.h"
#include "Utilities.h"
#include "Tensor.h"
#include "ElmMats.h"
#include "ElmNorm.h"
#include "Vec3Oper.h"
#include "VTF.h"


//! \brief Enum for element level solution vectors
enum SolutionVectors
{
  U = 0, // displacement
  P = 1, // pore pressure
  NSOL = 2
};


//! \brief Enum for element level right-hand-side vectors
enum ResidualVectors
{
  Fu = 0,
  Fp = 1,
  Fres = 2,
  Fprev = 3,
  NVEC = 4
};


//! \brief Enum for element level left-hand-side matrices
enum TangentMatrices
{
  uu = 0,
  up = 1,
  pp = 2,
  Ktan = 3,
  Kprev = 4,
  NMAT = 5
};


PoroElasticity::MixedElmMats::MixedElmMats ()
{
  this->resize(NMAT,NVEC); // Number of element matrices and vectors
}


const Matrix& PoroElasticity::MixedElmMats::getNewtonMatrix () const
{
  Matrix& N = const_cast<Matrix&>(A[Ktan]);

  size_t i,j;
  size_t n = A[uu].rows();
  size_t m = A[pp].rows();

  for (i = 1; i <= n; i++)
  {
    for (j = 1; j <= n; j++)
      N(i,j) = A[uu](i,j);
    for (j = 1; j <= m; j++)
    {
      size_t k = n+j;
      N(i,k) = A[up](i,j);
      N(k,i) = A[up](i,j);
    }
  }

  for (i = 1; i <= m; i++)
    for (j = 1; j <= m; j++)
    {
      size_t ki = n+i;
      size_t kj = n+j;
      N(ki,kj) = A[pp](i,j);
    }

  return A[Ktan];
}


const Vector& PoroElasticity::MixedElmMats::getRHSVector () const
{
  Vector& F = const_cast<Vector&>(b[Fres]);
  size_t i;
  size_t n = b[Fu].size();
  size_t m = b[Fp].size();

  for (i = 1; i <= n; i++)
    F(i) = b[Fu](i);

  for (i = 1; i <= m; i++)
    F(n+i) = b[Fp](i);

  F += b[Fprev];

  return b[Fres];
}


PoroElasticity::PoroElasticity (unsigned short int n, int order) : nsd(n), gacc(9.81)
{
  primsol.resize(1+order); // Current and previous timestep solutions required
  tracFld = nullptr;
  fluxFld = nullptr;
  eS = 1;
  sc = 0.0;
}


Vec3 PoroElasticity::getTraction(const Vec3& X, const Vec3& n) const
{
  if (fluxFld)
    return (*fluxFld)(X);
  else if (tracFld)
    return (*tracFld)(X,n);
  else
    return Vec3();
}


double PoroElasticity::getMassDensity (const Vec3&) const
{
  // Density of two-phase porous medium
  double rho = rhos*(1.0 - poro) + rhof*poro;
  return rho;
}


Vec3 PoroElasticity::getBodyForce(const Vec3& X) const
{
  Vec3 f = grav;
  f *= this->getMassDensity(X);
  return f;
}


LocalIntegral* PoroElasticity::getLocalIntegral (size_t nen1, size_t nen2,
                                                 size_t, bool neumann) const
{
  const size_t nedof1 = nsd*nen1;
  const size_t nedof = nedof1 + nen2;

  ElmMats* result = new MixedElmMats();

  result->rhsOnly = neumann;
  result->withLHS = !neumann;
  result->b[Fres].resize(nedof);
  result->b[Fprev].resize(nedof);
  result->b[Fu].resize(nedof1);
  result->b[Fp].resize(nen2);

  if (!neumann)
  {
    result->A[uu].resize(nedof1,nedof1);
    result->A[up].resize(nedof1,nen2);
    result->A[pp].resize(nen2,nen2);
    result->A[Ktan].resize(nedof,nedof);
    result->A[Kprev].resize(nedof,nedof);
  }

  return result;
}


bool PoroElasticity::initElement (const std::vector<int>& MNPC1,
                                  const std::vector<int>& MNPC2, size_t n1,
                                  LocalIntegral& elmInt)
{
  if (primsol.front().empty()) return true;

  // Extract the element level solution vectors
  elmInt.vec.resize(NSOL);
  int ierr = utl::gather(MNPC1,nsd,primsol.front(),elmInt.vec[U])
           + utl::gather(MNPC2,0,1,primsol.front(),elmInt.vec[P],nsd*n1,n1);

  if (ierr == 0) return true;

  std::cerr << " *** PoroElasticity::initElement: Detected " << ierr/3
            << " node numbers out of range." << std::endl;

  return false;
}


bool PoroElasticity::initElementBou (const std::vector<int>& MNPC1,
                                     const std::vector<int>&,
                                     size_t, LocalIntegral& elmInt)
{
  return this->IntegrandBase::initElementBou(MNPC1,elmInt);
}


bool PoroElasticity::formBmatrix (Matrix& Bmat, const Matrix& dNdX) const
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

  // Strain-displacement matrix for 2D elements
  for (size_t i = 1; i <= nenod; i++)
  {
    Bmat(INDEX(1,1),i) = dNdX(i,1);
    Bmat(INDEX(2,2),i) = dNdX(i,2);
    Bmat(INDEX(3,1),i) = dNdX(i,2);
    Bmat(INDEX(3,2),i) = dNdX(i,1);
  }

#undef INDEX

  Bmat.resize(nstrc,nsd*nenod);

  return true;
}


bool PoroElasticity::formElasticMatrix(Matrix& Cmat) const
{
  double C33 = E/(2.0+2.0*nu);
  double C12 = (E*nu)/((1.0+nu)*(1.0-2.0*nu));
  double C11 = C12 + 2.0*C33;

  const size_t nstrc = nsd*(nsd+1)/2;

  Cmat.resize(nstrc,nstrc,true);

  Cmat(1,1) = C11;
  Cmat(1,2) = C12;
  Cmat(2,1) = C12;
  Cmat(2,2) = C11;
  Cmat(3,3) = C33;

  Cmat.resize(nstrc,nstrc);

  return true;
}


bool PoroElasticity::evalIntMx (LocalIntegral& elmInt,
                                const MxFiniteElement& fe,
                                const TimeDomain& time, const Vec3& X) const
{
  ElmMats& elMat = static_cast<ElmMats&>(elmInt);

  size_t i,j,k;
  Matrix Bmat, Cmat, CB;

  if (!this->formBmatrix(Bmat,fe.dN1dX))
    return false;

  if (!this->formElasticMatrix(Cmat))
    return false;

  Vec3 permeability;
  if(perm)
    permeability = (*perm)(X);
  else
    permeability = 1.0;

  double scl(sc);
  if (scl == 0.0)
    scl = sqrt(E*rhof*gacc/(permeability[0]*time.dt));

  // Biot's coefficient
  double alpha = 1.0 - (Ko/Ks);
  // Inverse of the compressibility modulus
  double Minv = ((alpha - poro)/Ks) + (poro/Kw);
  // m vector for 2D
  Vec3 m(1,1,0);

  // Integration of the stiffness matrix
  CB.multiply(Cmat,Bmat,false,false);
  CB *= -1.0 * fe.detJxW;
  elMat.A[uu].multiply(Bmat,CB,true,false,true);

  // Integration of the coupling matrix
  Matrix Kuptmp;
  const size_t nstrc = nsd*(nsd+1)/2;
  Kuptmp.resize(fe.N2.size(),nstrc);

  for (i = 1; i <= fe.N2.size(); i++)
    for (j = 1; j <= nstrc; j++)
      Kuptmp(i,j) += scl*m[j-1]*alpha*fe.N2(i)*fe.detJxW;

  elMat.A[up].multiply(Bmat,Kuptmp,true,true,true);

  // Integration of the compressibilty matrix
  Matrix Cpp;
  Cpp.resize(fe.N2.size(),fe.N2.size());
  for (i = 1; i <= fe.N2.size(); i++)
    for (j = 1; j <= fe.N2.size(); j++)
      Cpp(i,j) += scl*scl*fe.N2(i)*Minv*fe.N2(j)*fe.detJxW;

  // Integration of the permeability matrix
  Matrix Kpp;
  Kpp.resize(fe.N2.size(),fe.N2.size());
  for (i = 1; i <= fe.N2.size(); i++)
    for (j = 1; j <= fe.N2.size(); j++)
      for (k = 1; k <= nsd; k++)
        Kpp(i,j) += scl*scl*fe.dN2dX(i,k)*(permeability[k-1]/(rhof*gacc))*fe.dN2dX(j,k)*fe.detJxW;

  elMat.A[pp] += Cpp;
  elMat.A[pp].add(Kpp,time.dt);

  size_t rAuu = elMat.A[uu].rows();
  size_t rApp = elMat.A[pp].rows();

  for (i = 1; i <= rApp; i++)
    for (j = 1; j <= rApp; j++)
      elMat.A[Kprev](rAuu+i,rAuu+j) += Cpp(i,j);

  return true;
}


bool PoroElasticity::evalBouMx (LocalIntegral& elmInt,
                                const MxFiniteElement& fe,
                                const TimeDomain& time, const Vec3& X, const Vec3& normal) const
{
  if (!tracFld && !fluxFld)
  {
    std::cerr << " *** PoroElasticity::evalBouMx: No fluxes/tractions." << std::endl;
    return false;
  }
  else if (!eS)
  {
    std::cerr << " *** PoroElasticity::evalBouMx: No load vector." << std::endl;
    return false;
  }

  // Evaluate the surface traction
  Vec4 Xt = static_cast<const Vec4&>(X);
  Xt.t = time.t;
  Vec3 tr2 = this->getTraction(Xt,normal);
  Xt.t -= time.dt;
  Vec3 tr1 = this->getTraction(Xt,normal);
  Vec3 dtr;
  dtr = tr2 - tr1;

  Vec3 permeability;
  if (perm)
    permeability = (*perm)(X);

  Vec3 bf = this->getBodyForce(X);

  // Integrate the force vector fu
  ElmMats& elMat = static_cast<ElmMats&>(elmInt);
  for (size_t i = 1; i <= fe.N1.size(); i++)
    for (unsigned short int j = 1; j <= nsd; j++)
      elMat.b[Fu](nsd*(i-1)+j) += (-1.0)*(dtr[j-1]*fe.N1(i)*fe.detJxW +
          bf[j-1]*fe.N1(i)*fe.detJxW);

  // First term of fp vector in RHS, remember to add water flux term
  for (size_t i = 1; i <= fe.N2.size(); i++)
  {
    double fpvec = 0.0;
    for (size_t k = 1; k <= nsd; k++)
      fpvec += sc*fe.dN2dX(i,k)*(permeability[k-1]/gacc)*grav[k-1];
    elMat.b[Fp](i) = fpvec*time.dt*fe.detJxW;
  }

  return true;
}


bool PoroElasticity::finalizeElement (LocalIntegral& elmInt, const TimeDomain&, size_t)
{
  ElmMats& elMat = static_cast<ElmMats&>(elmInt);

  size_t i,j;

  size_t rAuu = elMat.A[uu].rows();
  size_t rApp = elMat.A[pp].rows();

  for (i = 1; i <= rAuu; i++)
  {
    for (j = 1; j <= rAuu; j++)
      elMat.A[Kprev](i,j) = elMat.A[uu](i,j);
    for (j = 1; j <= rApp; j++)
    {
      elMat.A[Kprev](i,rAuu+j) = elMat.A[up](i,j);
      elMat.A[Kprev](rAuu+j,i) = elMat.A[up](i,j);
    }
  }

  Vector prevSol;
  prevSol = elmInt.vec[U];

  prevSol.insert(prevSol.end(),elmInt.vec[P].begin(),elmInt.vec[P].end());

  elMat.b[Fprev] = elMat.A[Kprev]*prevSol;

  return true;
}


size_t PoroElasticity::getNoFields (int fld) const
{
  if (fld < 2)
    return nsd+1;
  else
    return 6;
}


const char* PoroElasticity::getField1Name (size_t i, const char* prefix) const
{
  if (i >= nsd) i = 3;

  static const char* s[4] = { "u_x", "u_y", "u_z", "p^w" };
  if(!prefix) return s[i];

  static std::string name;
  name = prefix + std::string(" ") + s[i];

  return name.c_str();
}


const char* PoroElasticity::getField2Name (size_t i, const char* prefix) const
{
  if (i >= nsd) return 0;

  static const char* s2[] = {"eps_x","eps_y","eps_xy","sig_x",
                             "sig_y","sig_z","sig_xy"}; // borked in 3D
  if (!prefix) return s2[i];

  static std::string name;
  name = prefix + std::string(" ") + s2[i];

  return name.c_str();
}


bool PoroElasticity::evalSol(Vector& s, const MxFiniteElement& fe,
                             const Vec3& X, const std::vector<int>& MNPC1,
                             const std::vector<int>& MNPC2) const
{
  Vector eV;
  utl::gather(MNPC1,nsd,s,eV);

  Vec3 V;
  for (size_t s=0;s<nsd;++s)
    V[s] = eV.dot(fe.N1,s,nsd);

  Vector eP;
  utl::gather(MNPC2,1,s,eP);

  return true;
}
