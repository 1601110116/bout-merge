/**************************************************************************
 * Laplacian solver in 2D (X-Y)
 *
 * Equation solved is: 
 *
 * d*Laplace_perp(x) +(1/c)*Grad(c) dot Grad_perp(x) + a*x = b
 *
 * Intended for use in solving n=0 component of potential
 * from inversion of vorticity equation in convective form 
 * 
 * Last Update: 2017-10-04 by H.SETO (QST)
 * 
 **************************************************************************
 * Copyright 2015 B.Dudson
 *
 * Contact: Ben Dudson, bd512@york.ac.uk
 * 
 * This file is part of BOUT++.
 *
 * BOUT++ is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * BOUT++ is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with BOUT++.  If not, see <http://www.gnu.org/licenses/>.
 *
 **************************************************************************/

#ifndef __LAPLACE_XY2_H__
#define __LAPLACE_XY2_H__

#ifndef BOUT_HAS_PETSC
// If no PETSc

#warning LaplaceXY2 requires PETSc. No LaplaceXY2 available

#include <bout/mesh.hxx>
#include <options.hxx>
#include <boutexception.hxx>

/*!
 * Create a dummy class so that code will compile
 * without PETSc, but will throw an exception if
 * LaplaceXY is used.
 */
class LaplaceXY2 {
 public:
  LaplaceXY2(Mesh *m, Options *opt = NULL) {
    throw BoutException("LaplaceXY requires PETSc. No LaplaceXY available");
  }
  void setCoefs(const Field2D &a, const Field2D &c, const Field2D &d){}
  const Field2D solve(const Field2D &rhs, const Field2D &x0) {}
};

#else // BOUT_HAS_PETSC

#include <bout/mesh.hxx>
#include <bout/petsclib.hxx>
#include <cyclic_reduction.hxx>

class LaplaceXY2 {
public:
  /*! 
   * Constructor
   */
  LaplaceXY2(Mesh *m, Options *opt = NULL);
  /*!
   * Destructor
   */
  ~LaplaceXY2();

  /*!
   * Set coefficients (a, c ,d) in equation:
   *
   * d*Laplace_perp(x) +(1/c)*Grad(c) dot Grad_perp(x) + a*x = b
   * 
   */
  void setCoefs(const Field2D &a, const Field2D &c, const Field2D &d);
  
  /*!
   * Solve Laplacian in X-Y
   * 
   * Inputs
   * ======
   *
   * rhs  - The field to be inverted. This must be allocated
   *        and contain valid data.
   * x0   - Initial guess at the solution. If this is unallocated
   *        then an initial guess of zero will be used.
   * 
   * Returns
   * =======
   * 
   * The solution as a Field2D. On failure an exception will be raised
   * 
   */
  const Field2D solve(const Field2D &rhs, const Field2D &x0);

  /*!
   * Preconditioner function
   * This is called by PETSc via a static function.
   * and should not be called by external users
   */
  int precon(Vec x, Vec y);

private:
  
  PetscLib lib;     ///< Requires PETSc library
  Mat MatA;         ///< Matrix to be inverted
  Vec xs, bs;       ///< Solution and RHS vectors
  KSP ksp;          ///< Krylov Subspace solver
  PC pc;            ///< Preconditioner

  Mesh *mesh;   ///< The mesh this operates on, provides metrics and communication
  
  // Preconditioner
  int xstart, xend;
  int nloc, nsys;
  BoutReal **acoef, **bcoef, **ccoef, **xvals, **bvals;
  CyclicReduce<BoutReal> *cr; ///< Tridiagonal solver

  // Y derivatives
  bool include_y_derivs; // Include Y derivative terms?
  
  // Boundary conditions
  bool x_inner_dirichlet; // Dirichlet on inner X boundary?
  bool x_outer_dirichlet; // Dirichlet on outer X boundary?
  bool y_bndry_dirichlet; // Dirichlet on Y boundary?
  
  /*!
   * Number of grid points on this processor
   */
  int localSize();
  
  /*!
   * Return the communicator for XY
   */
  MPI_Comm communicator();
  
  /*!
   * Return the global index of a local (x,y) coordinate
   * including guard cells.
   * Boundary cells have a global index of -1
   *
   * To do this, a Field2D (indexXY) is used to store
   * the index as a floating point number which is then rounded
   * to an integer. Guard cells are filled by communication
   * so no additional logic is needed in Mesh.
   */
  int globalIndex(int x, int y);  
  Field2D indexXY; ///< Global index (integer stored as BoutReal)


  /*!
   * parts of coefficients being independent of (a,c,d)
   */
  Field2D coefPC0, coefPC1;
  Field2D coefMC0, coefMC1;
  Field2D coefCP0, coefCP1, coefCP2;  
  Field2D coefCM0, coefCM1, coefCM2;
  Field2D coefCC0, coefCC1;
  
  /*!
   * Round a number to the nearest integer
   */
  int roundInt(BoutReal f);
  
};

#endif // BOUT_HAS_PETSC
#endif // __LAPLACE_XY2_H__