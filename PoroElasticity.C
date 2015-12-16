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


void PoroElasticity::MixedElmMats::makeNewtonMatrix(Matrix& N, bool dopp) const
{
  size_t n = A[uu].rows();

  N.fillBlock(A[uu], 1, 1);
  N.fillBlock(A[up], 1, 1+n);
  N.fillBlock(A[up], 1+n, 1, true);
  if (dopp)
    N.fillBlock(A[pp], 1+n, 1+n);
}


const Matrix& PoroElasticity::MixedElmMats::getNewtonMatrix () const
{
  Matrix& N = const_cast<Matrix&>(A[Ktan]);

  makeNewtonMatrix(N, true);

  return A[Ktan];
}


const Vector& PoroElasticity::MixedElmMats::getRHSVector () const
{
  Vector& F = const_cast<Vector&>(b[Fres]);
  size_t n = b[Fu].size();

  std::copy(b[Fu].begin(), b[Fu].end(), F.begin());
  std::copy(b[Fp].begin(), b[Fp].end(), F.begin()+n);

  F += b[Fprev];

  return b[Fres];
}


PoroElasticity::NonMixedElmMats::NonMixedElmMats ()
{
  this->resize(NMAT,NVEC); // Number of element matrices and vectors
}

void PoroElasticity::NonMixedElmMats::makeNewtonMatrix(Matrix& N, bool dopp) const
{
  size_t n = A[pp].rows();

  size_t nsd = A[uu].rows()/A[pp].rows();
  size_t nf = nsd + 1;
  if (dopp) {
    for (size_t i = 1; i <= n; ++i)
      for (size_t j = 1; j <= n; ++j)
        N(nf*i, nf*j) = A[pp](i,j);
  }

  for (size_t i = 1; i <= n; ++i) {
    for (size_t l = 1; l <= nsd; ++l) {
      for (size_t j = 1; j <= n; ++j) {
        for (size_t k = 1; k <= nsd; ++k)
          N(nf*(i-1)+l, (j-1)*nf+k) = A[uu](nsd*(i-1)+l, (j-1)*nsd+k);
        N(nf*(i-1)+l, j*nf) = N(j*nf, nf*(i-1)+l) = A[up](nsd*(i-1)+l, j);
      }
    }
  }
}


const Matrix& PoroElasticity::NonMixedElmMats::getNewtonMatrix () const
{
  Matrix& N = const_cast<Matrix&>(A[Ktan]);

  makeNewtonMatrix(N, true);

  return A[Ktan];
}


const Vector& PoroElasticity::NonMixedElmMats::getRHSVector () const
{
  Vector& F = const_cast<Vector&>(b[Fres]);

  size_t nsd = b[Fu].size()/b[Fp].size();
  utl::interleave(b[Fu], b[Fp], F, nsd, 1);

  F += b[Fprev];

  return b[Fres];
}


PoroElasticity::PoroElasticity (unsigned short int n, int order) :
  gacc(9.81), mat(nullptr)
{
  nsd = n;
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


double PoroElasticity::getMassDensity (const Vec3& X) const
{
  // Density of two-phase porous medium
  if (mat)
    return mat->getSolidDensity(X)*(1.0-mat->getPorosity(X)) +
           mat->getFluidDensity(X)*mat->getPorosity(X);

  return 0.0;
}


Vec3 PoroElasticity::getBodyForce(const Vec3& X) const
{
  Vec3 f = grav;
  f *= this->getMassDensity(X);
  return f;
}


LocalIntegral* PoroElasticity::getLocalIntegral (const std::vector<size_t>& nen,
                                                 size_t, bool neumann) const
{
  const size_t nedof1 = nsd*nen[0];
  const size_t nedof = nedof1 + nen[1];

  ElmMats* result = new MixedElmMats();

  result->rhsOnly = neumann;
  result->withLHS = !neumann;
  result->b[Fres].resize(nedof);
  result->b[Fprev].resize(nedof);
  result->b[Fu].resize(nedof1);
  result->b[Fp].resize(nen[1]);

  if (!neumann)
  {
    result->A[uu].resize(nedof1,nedof1);
    result->A[up].resize(nedof1,nen[1]);
    result->A[pp].resize(nen[1],nen[1]);
    result->A[Ktan].resize(nedof,nedof);
    result->A[Kprev].resize(nedof,nedof);
  }

  return result;
}


LocalIntegral* PoroElasticity::getLocalIntegral (size_t nen,
                                                 size_t, bool neumann) const
{
  const size_t nedof1 = nsd*nen;
  const size_t nedof = nedof1 + nen;

  ElmMats* result = new NonMixedElmMats();

  result->rhsOnly = neumann;
  result->withLHS = !neumann;
  result->b[Fres].resize(nedof);
  result->b[Fprev].resize(nedof);
  result->b[Fu].resize(nedof1);
  result->b[Fp].resize(nen);

  if (!neumann)
  {
    result->A[uu].resize(nedof1,nedof1);
    result->A[up].resize(nedof1,nen);
    result->A[pp].resize(nen,nen);
    result->A[Ktan].resize(nedof,nedof);
    result->A[Kprev].resize(nedof,nedof);
  }

  return result;
}


bool PoroElasticity::initElement (const std::vector<int>& MNPC,
                                  const std::vector<size_t>& elem_sizes,
                                  const std::vector<size_t>& basis_sizes,
                                  LocalIntegral& elmInt)
{
  if (primsol.front().empty()) return true;

  // Extract the element level solution vectors
  elmInt.vec.resize(NSOL);
  std::vector<int>::const_iterator fstart = MNPC.begin() + elem_sizes[0];
  int ierr = utl::gather(IntVec(MNPC.begin(),fstart),nsd,primsol.front(),elmInt.vec[U])
           + utl::gather(IntVec(fstart,MNPC.end()),0,1,primsol.front(),elmInt.vec[P],nsd*basis_sizes[0],basis_sizes[0]);

  if (ierr == 0) return true;

  std::cerr << " *** PoroElasticity::initElement: Detected " << ierr/3
            << " node numbers out of range." << std::endl;

  return false;
}


bool PoroElasticity::initElementBou (const std::vector<int>& MNPC,
                                     const std::vector<size_t>& elem_sizes,
                                     const std::vector<size_t>& basis_sizes,
                                     LocalIntegral& elmInt)
{
  return this->initElement(MNPC,elem_sizes,basis_sizes,elmInt);
}


bool PoroElasticity::initElement (const std::vector<int>& MNPC,
                                  LocalIntegral& elmInt)
{
  if (primsol.empty() || primsol.front().empty())
    return true;

  // Extract the element level solution vectors
  elmInt.vec.resize(NSOL);
  int ierr = 0;
  Matrix temp(nsd+1, MNPC.size());
  ierr += utl::gather(MNPC, nsd+1, primsol.front(), temp);
  Matrix temp2(nsd, MNPC.size());
  for (size_t k = 1; k <= nsd; ++k)
    temp2.fillRow(k, temp.getRow(k).ptr());
  elmInt.vec[U] = temp2;
  elmInt.vec[P] = temp.getRow(nsd+1);

  if (ierr == 0) return true;

  std::cerr << " *** PoroElasticity::initElement: Detected " << ierr/3
            << " node numbers out of range." << std::endl;

  return false;
}


bool PoroElasticity::initElementBou (const std::vector<int>& MNPC,
                                     LocalIntegral& elmInt)
{
  return this->initElement(MNPC,elmInt);
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


bool PoroElasticity::formElasticMatrix(Matrix& Cmat, const Vec3& X) const
{
  if (!mat)
    return false;

  double E = mat->getStiffness(X);
  double nu = mat->getPoisson(X);

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
  return evalInt(elmInt, fe, time, X);
}

bool PoroElasticity::evalInt (LocalIntegral& elmInt,
                              const FiniteElement& fe,
                              const TimeDomain& time, const Vec3& X) const
{
  ElmMats& elMat = static_cast<ElmMats&>(elmInt);

  if (!mat) {
    std::cerr << __FUNCTION__ << ": No material data." << std::endl;
    return false;
  }

  size_t i,j,k;
  Matrix Bmat, Cmat, CB;

  if (!this->formBmatrix(Bmat,fe.grad(1)))
    return false;

  if (!this->formElasticMatrix(Cmat,X))
    return false;

  Vec3 permeability = mat->getPermeability(X);

  double scl(sc);
  if (scl == 0.0)
    scl = sqrt(mat->getStiffness(X)*mat->getFluidDensity(X)*gacc/(permeability[0]*time.dt));

  // Biot's coefficient
  double Ko = mat->getBulkMedium(X);
  double Ks = mat->getBulkSolid(X);
  double Kw = mat->getBulkWater(X);
  double poro = mat->getPorosity(X);

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
  Kuptmp.resize(fe.basis(2).size(),nstrc);

  for (i = 1; i <= fe.basis(2).size(); i++)
    for (j = 1; j <= nstrc; j++)
      Kuptmp(i,j) += scl*m[j-1]*alpha*fe.basis(2)(i)*fe.detJxW;

  elMat.A[up].multiply(Bmat,Kuptmp,true,true,true);

  // Integration of the compressibilty matrix
  Matrix Cpp;
  Cpp.resize(fe.basis(2).size(),fe.basis(2).size());
  for (i = 1; i <= fe.basis(2).size(); i++)
    for (j = 1; j <= fe.basis(2).size(); j++)
      Cpp(i,j) += scl*scl*fe.basis(2)(i)*Minv*fe.basis(2)(j)*fe.detJxW;

  // Integration of the permeability matrix
  Matrix Kpp;
  Kpp.resize(fe.basis(2).size(),fe.basis(2).size());
  for (i = 1; i <= fe.basis(2).size(); i++)
    for (j = 1; j <= fe.basis(2).size(); j++)
      for (k = 1; k <= nsd; k++)
        Kpp(i,j) += scl*scl*fe.grad(2)(i,k)*(permeability[k-1]/(mat->getFluidDensity(X)*gacc))*fe.grad(2)(j,k)*fe.detJxW;

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
                                const TimeDomain& time, const Vec3& X,
                                const Vec3& normal) const
{
  return evalBou(elmInt, fe, time, X, normal);
}


bool PoroElasticity::evalBou (LocalIntegral& elmInt,
                              const FiniteElement& fe,
                              const TimeDomain& time, const Vec3& X,
                              const Vec3& normal) const
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

  Vec3 permeability = mat->getPermeability(X);

  Vec3 bf = this->getBodyForce(X);

  double scl(sc);
  if (scl == 0.0)
    scl = sqrt(mat->getStiffness(X)*mat->getFluidDensity(X)*gacc/(permeability[0]*time.dt));

  // Integrate the force vector fu
  ElmMats& elMat = static_cast<ElmMats&>(elmInt);
  for (size_t i = 1; i <= fe.basis(1).size(); i++)
    for (unsigned short int j = 1; j <= nsd; j++)
      elMat.b[Fu](nsd*(i-1)+j) += (-1.0)*(dtr[j-1]*fe.basis(1)(i)*fe.detJxW +
                                  bf[j-1]*fe.basis(1)(i)*fe.detJxW);

  // First term of fp vector in RHS, remember to add water flux term
  for (size_t i = 1; i <= fe.basis(2).size(); i++)
  {
    double fpvec = 0.0;
    for (size_t k = 1; k <= nsd; k++)
      fpvec += scl*fe.grad(2)(i,k)*(permeability[k-1]/gacc)*grav[k-1];
    elMat.b[Fp](i) = fpvec*time.dt*fe.detJxW;
  }

  return true;
}


bool PoroElasticity::finalizeElement (LocalIntegral& elmInt, const TimeDomain&, size_t)
{
  ElmMats& elMat = static_cast<ElmMats&>(elmInt);

  Vector prevSol;

  MixedElmMats* mMat = dynamic_cast<MixedElmMats*>(&elmInt);
  if (mMat) {
    mMat->makeNewtonMatrix(elMat.A[Kprev], false);
    prevSol = elmInt.vec[U];
    prevSol.insert(prevSol.end(),elmInt.vec[P].begin(),elmInt.vec[P].end());
  } else {
    NonMixedElmMats* mMat = dynamic_cast<NonMixedElmMats*>(&elmInt);
    mMat->makeNewtonMatrix(elMat.A[Kprev], false);
    utl::interleave(elmInt.vec[U], elmInt.vec[P], prevSol, nsd, 1);
  }

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


std::string PoroElasticity::getField1Name (size_t i, const char* prefix) const
{
  if (i == 11)
    i = 4;
  else if (i == 12)
    i = 3;
  else if (i >= nsd)
    i = 3;

  static const char* s[5] = { "u_x", "u_y", "u_z", "p^w", "u" };
  if(!prefix) return s[i];

  return prefix + std::string(" ") + s[i];
}


std::string PoroElasticity::getField2Name (size_t i, const char* prefix) const
{
  if (i >= nsd) return "";

  static const char* s2[] = {"eps_x","eps_y","eps_xy","sig_x",
                             "sig_y","sig_z","sig_xy"}; // borked in 3D
  if (!prefix) return s2[i];

  return prefix + std::string(" ") + s2[i];
}


bool PoroElasticity::evalSol(Vector& s, const MxFiniteElement& fe,
                             const Vec3& X, const std::vector<int>& MNPC,
                             const std::vector<size_t>& elem_sizes) const
{
  Vector eV;
  std::vector<int>::const_iterator fstart = MNPC.begin() + elem_sizes[0];
  utl::gather(IntVec(MNPC.begin(),fstart),nsd,s,eV);

  Vec3 V;
  for (size_t s=0;s<nsd;++s)
    V[s] = eV.dot(fe.basis(1),s,nsd);

  return true;
}


bool PoroElasticity::evalSol(Vector& s, const FiniteElement& fe,
                             const Vec3& X,
                             const std::vector<int>& MNPC) const
{
  Vector eV;
  utl::gather(MNPC,nsd+1,s,eV);

  Vec3 V;
  for (size_t s=0;s<nsd;++s)
    V[s] = eV.dot(fe.basis(1),s,nsd+1);

  return true;
}
