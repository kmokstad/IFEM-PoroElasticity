// $Id$
//==============================================================================
//!
//! \file PoroMats.C
//!
//! \date April 20 2015
//!
//! \author Eivind Fonn
//!
//! \brief Element matrix implementation for PoroElasticity.
//!
//==============================================================================

#include "PoroElasticity.h"
#include "Utilities.h"


PoroElasticity::Mats::Mats (size_t ndof_displ, size_t ndof_press,
                            bool neumann, char dynamic)
{
  this->resize(NMAT, NVEC);

  size_t ndof_total = ndof_displ + ndof_press;

  rhsOnly = neumann;
  withLHS = !neumann;
  b[Fsys].resize(ndof_total);
  b[Fu].resize(ndof_displ);
  b[Fp].resize(ndof_press);

  if (withLHS) {
    A[Nsys].resize(ndof_total, ndof_total);
    A[uu_K].resize(ndof_displ, ndof_displ);
    if (dynamic > 0)
      A[uu_M].resize(ndof_displ, ndof_displ);
    A[up_Q].resize(ndof_displ, ndof_press);
    if (dynamic == 2)
      A[up_D].resize(ndof_displ, ndof_press);
    A[pp_S].resize(ndof_press, ndof_press);
    A[pp_P].resize(ndof_press, ndof_press);
  }

  h = lambda = 0.0;
}


const Matrix& PoroElasticity::Mats::getNewtonMatrix () const
{
  Matrix& res = const_cast<Matrix&>(A[Nsys]);
  this->add_uu(A[uu_K], res);
  this->add_up(A[up_Q], res, -1.0);
  this->add_pu(A[up_Q], res);
  this->add_pp(A[pp_S], res);
  this->add_pp(A[pp_P], res, h);
#if INT_DEBUG > 2
  std::cout <<"\nPoroElasticity::Mats::getNewtonMatrix:"
            <<"\nElement stiffness matrix, K_uu"<< A[uu_K]
            <<"\nElement coupling matrix, Q_up"<< A[up_Q]
            <<"\nElement compressibility matrix, S_pp"<< A[pp_S]
            <<"\nElement permeability matrix, P_pp"<< A[pp_P]
            <<"\nElement coefficient matrix" << A[Nsys];
#endif
  return A[Nsys];
}


const Vector& PoroElasticity::Mats::getRHSVector () const
{
#if INT_DEBUG > 2
  std::cout <<"\nPoroElasticity::Mats::getRHSVector:"<< std::endl;
  if (vec.size() > Vu) std::cout <<"Element displacement, Vu"<< vec[Vu];
  if (vec.size() > Vp) std::cout <<"Element pressure, Vp"<< vec[Vp];
  std::cout <<"S_ext-S_int"<< b[Fu] <<"S_p"<< b[Fp];
#endif

  Vector temp(b[Fp]); temp *= h;                  // temp = Fp*h
  if (A.size() > up_Q && vec.size() > Vu)
    A[up_Q].multiply(vec[Vu], temp, true, true);  // temp += Q_up^t * u
  if (A.size() > pp_S && vec.size() > Vp)
    A[pp_S].multiply(vec[Vp], temp, false, true); // temp += S_pp * p

  this->form_vector(b[Fu], temp, const_cast<Vector&>(b[Fsys]));
#if INT_DEBUG > 2
  std::cout <<"\nElement right-hand-side vector"<< b[Fsys];
#endif
  return b[Fsys];
}


void PoroElasticity::MixedElmMats::add_uu (const Matrix& source, Matrix& target,
                                           double scl) const
{
  target.addBlock(source, scl, 1, 1);
}


void PoroElasticity::MixedElmMats::add_up (const Matrix& source, Matrix& target,
                                           double scl) const
{
  target.addBlock(source, scl, 1, 1 + b[Fu].size());
}


void PoroElasticity::MixedElmMats::add_pu (const Matrix& source, Matrix& target,
                                           double scl) const
{
  target.addBlock(source, scl, 1 + b[Fu].size(), 1, true);
}


void PoroElasticity::MixedElmMats::add_pp (const Matrix& source, Matrix& target,
                                           double scl) const
{
  target.addBlock(source, scl, 1 + b[Fu].size(), 1 + b[Fu].size());
}


void PoroElasticity::MixedElmMats::form_vector (const Vector& u,
                                                const Vector& p,
                                                Vector& target) const
{
  target = u;
  target.insert(target.end(), p.begin(), p.end());
}


void PoroElasticity::StdElmMats::add_uu (const Matrix& source, Matrix& target,
                                         double scl) const
{
  size_t nen = b[Fp].size();
  size_t nf  = b[Fu].size()/nen + 1;
  size_t i, j, k, l, ii, jj, kk, ll;

  for (i = kk = 1, ii = 0; i <= nen; i++, ii += nf)
    for (k = 1; k < nf; k++, kk++)
      for (j = ll = 1, jj = 0; j <= nen; j++, jj += nf)
        for (l = 1; l < nf; l++, ll++)
          target(ii+k,jj+l) += scl * source(kk,ll);
}


void PoroElasticity::StdElmMats::add_up (const Matrix& source, Matrix& target,
                                         double scl) const
{
  size_t nen = b[Fp].size();
  size_t nf  = b[Fu].size()/nen + 1;
  size_t i, j, k, ii, jj, kk;

  for (i = kk = 1, ii = 0; i <= nen; i++, ii += nf)
    for (k = 1; k < nf; k++, kk++)
      for (j = 1, jj = nf; j <= nen; j++, jj += nf)
        target(ii+k,jj) += scl * source(kk,j);
}


void PoroElasticity::StdElmMats::add_pu (const Matrix& source, Matrix& target,
                                         double scl) const
{
  size_t nen = b[Fp].size();
  size_t nf  = b[Fu].size()/nen + 1;
  size_t i, j, k, ii, jj, kk;

  for (i = kk = 1, ii = 0; i <= nen; i++, ii += nf)
    for (k = 1; k < nf; k++, kk++)
      for (j = 1, jj = nf; j <= nen; j++, jj += nf)
        target(jj,ii+k) += scl * source(kk,j);
}


void PoroElasticity::StdElmMats::add_pp (const Matrix& source, Matrix& target,
                                         double scl) const
{
  size_t nen = b[Fp].size();
  size_t nf  = b[Fu].size()/nen + 1;
  size_t i, j, ii, jj;

  for (i = 1, ii = nf; i <= nen; i++, ii += nf)
    for (j = 1, jj = nf; j <= nen; j++, jj += nf)
      target(ii,jj) += scl * source(i,j);
}


void PoroElasticity::StdElmMats::form_vector (const Vector& u,
                                              const Vector& p,
                                              Vector& target) const
{
  utl::interleave(u, p, target, b[Fu].size()/b[Fp].size(), 1);
}
