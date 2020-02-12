
#ifdef BOUT_HAS_PETSC

#include <petscksp.h>

#include <bout/invert/laplacexy2.hxx>

#include <bout/assert.hxx>

#include <boutcomm.hxx>
#include <utils.hxx>
#include <bout/sys/timer.hxx>

#include <output.hxx>

#undef __FUNCT__
#define __FUNCT__ "laplacePCapply"
static PetscErrorCode laplacePCapply(PC pc,Vec x,Vec y) {
  int ierr;
  
  // Get the context
  LaplaceXY2 *s;
  ierr = PCShellGetContext(pc,(void**)&s);CHKERRQ(ierr);
  
  PetscFunctionReturn(s->precon(x, y));
}

LaplaceXY2::LaplaceXY2(Mesh *m, Options *opt) : mesh(m) {
  Timer timer("invert");
  
  if(opt == NULL) {
    // If no options supplied, use default
    opt = Options::getRoot()->getSection("laplacexy2");
  }
  
  // Get MPI communicator
  MPI_Comm comm = BoutComm::get();
  
  // Local size
  int localN = localSize();

  // Create Vectors 
  VecCreate( comm, &xs );
  VecSetSizes( xs, localN, PETSC_DETERMINE );
  VecSetFromOptions( xs );
  VecDuplicate( xs , &bs );

  // Set size of Matrix on each processor to localN x localN
  MatCreate( comm, &MatA );                                
  MatSetSizes( MatA, localN, localN, PETSC_DETERMINE, PETSC_DETERMINE );
  MatSetFromOptions(MatA);
  
  //////////////////////////////////////////////////
  // Specify local indices. This creates a mapping
  // from local indices to index, using a Field2D object

  indexXY = -1;  // Set all points to -1, indicating out of domain

  int ind = 0;
  
  // Y boundaries
  for(RangeIterator it=mesh->iterateBndryLowerY(); !it.isDone(); it++) {
    indexXY(it.ind, mesh->ystart-1) = ind++;
  }
  for(RangeIterator it=mesh->iterateBndryUpperY(); !it.isDone(); it++) {
    indexXY(it.ind, mesh->yend+1) = ind++;
  }
  
  xstart = mesh->xstart;
  if(mesh->firstX())
    xstart -= 1; // Include X guard cells
  xend = mesh->xend;
  if(mesh->lastX())
    xend += 1;
  for(int x=xstart;x<=xend;x++)
    for(int y=mesh->ystart;y<=mesh->yend;y++) {
      indexXY(x,y) = ind++;
    }
  
  ASSERT1(ind == localN); // Reached end of range

  //////////////////////////////////////////////////
  // Allocate storage for preconditioner
  
  nloc    = xend - xstart + 1; // Number of X points on this processor
  nsys = mesh->yend - mesh->ystart + 1; // Number of separate Y slices
  
  acoef = matrix<BoutReal>(nsys, nloc);
  bcoef = matrix<BoutReal>(nsys, nloc);
  ccoef = matrix<BoutReal>(nsys, nloc);
  xvals = matrix<BoutReal>(nsys, nloc);
  bvals = matrix<BoutReal>(nsys, nloc);

  // Create a cyclic reduction object
  cr = new CyclicReduce<BoutReal>(mesh->getXcomm(), nloc);

  //////////////////////////////////////////////////
  // Pre-allocate PETSc storage
  
  PetscInt *d_nnz, *o_nnz;
  PetscMalloc( (localN)*sizeof(PetscInt), &d_nnz );
  PetscMalloc( (localN)*sizeof(PetscInt), &o_nnz );
  
  for(int i=0;i<localN;i++) {
    // Non-zero elements on this processor
    d_nnz[i] = 5; // Star pattern in 2D
    // Non-zero elements on neighboring processor
    o_nnz[i] = 0;
  }
  
  // X boundaries
  if(mesh->firstX()) {
    // Lower X boundary
    for(int y=mesh->ystart;y<=mesh->yend;y++) {
      int localIndex = indexXY(mesh->xstart-1,y);
      ASSERT1( (localIndex >= 0) && (localIndex < localN) );
    
      d_nnz[localIndex] = 2; // Diagonal sub-matrix
      o_nnz[localIndex] = 0; // Off-diagonal sub-matrix
    }
  }else {
    // On another processor
    for(int y=mesh->ystart;y<=mesh->yend;y++) {
      int localIndex = indexXY(mesh->xstart,y);
      ASSERT1( (localIndex >= 0) && (localIndex < localN) );
      d_nnz[localIndex] -= 1;
      o_nnz[localIndex] += 1;
    }
  }
  if(mesh->lastX()) {
    // Upper X boundary
    for(int y=mesh->ystart;y<=mesh->yend;y++) {
      int localIndex = indexXY(mesh->xend+1,y);
      ASSERT1( (localIndex >= 0) && (localIndex < localN) );
      d_nnz[localIndex] = 2; // Diagonal sub-matrix
      o_nnz[localIndex] = 0; // Off-diagonal sub-matrix
    }
  }else {
    // On another processor
    for(int y=mesh->ystart;y<=mesh->yend;y++) {
      int localIndex = indexXY(mesh->xend,y);
      ASSERT1( (localIndex >= 0) && (localIndex < localN) );
      d_nnz[localIndex] -= 1;
      o_nnz[localIndex] += 1;
    }
  }
  // Y boundaries
  
  for(int x=mesh->xstart; x <=mesh->xend; x++) {
    // Default to no boundary
    // NOTE: This assumes that communications in Y are to other
    //   processors. If Y is communicated with this processor (e.g. NYPE=1)
    //   then this will result in PETSc warnings about out of range allocations
    
    int localIndex = indexXY(x, mesh->ystart);
    ASSERT1( (localIndex >= 0) && (localIndex < localN) );
    //d_nnz[localIndex] -= 1;  // Note: Slightly inefficient
    o_nnz[localIndex] += 1;
    
    localIndex = indexXY(x, mesh->yend);
    ASSERT1( (localIndex >= 0) && (localIndex < localN) );
    //d_nnz[localIndex] -= 1; // Note: Slightly inefficient
    o_nnz[localIndex] += 1;
  }
  
  for(RangeIterator it=mesh->iterateBndryLowerY(); !it.isDone(); it++) {
    int localIndex = indexXY(it.ind, mesh->ystart-1);
    ASSERT1( (localIndex >= 0) && (localIndex < localN) );
    d_nnz[localIndex] = 2; // Diagonal sub-matrix
    o_nnz[localIndex] = 0; // Off-diagonal sub-matrix
    
    localIndex = indexXY(it.ind, mesh->ystart);
    ASSERT1( (localIndex >= 0) && (localIndex < localN) );
    d_nnz[localIndex] += 1;
    o_nnz[localIndex] -= 1;
  }
  for(RangeIterator it=mesh->iterateBndryUpperY(); !it.isDone(); it++) {
    int localIndex = indexXY(it.ind, mesh->yend+1);
    ASSERT1( (localIndex >= 0) && (localIndex < localN) );
    d_nnz[localIndex] = 2; // Diagonal sub-matrix
    o_nnz[localIndex] = 0; // Off-diagonal sub-matrix
    
    localIndex = indexXY(it.ind, mesh->yend);
    ASSERT1( (localIndex >= 0) && (localIndex < localN) );
    d_nnz[localIndex] += 1;
    o_nnz[localIndex] -= 1;
  }
  // Pre-allocate
  MatMPIAIJSetPreallocation( MatA, 0, d_nnz, 0, o_nnz );
  MatSetUp(MatA); 
  
  PetscFree( d_nnz );
  PetscFree( o_nnz );
  
  // Determine which row/columns of the matrix are locally owned
  int Istart, Iend;
  MatGetOwnershipRange( MatA, &Istart, &Iend );
  
  // Convert indexXY from local index to global index
  indexXY += Istart;
  
  // Now communicate to fill guard cells
  mesh->communicate(indexXY);
  
  // Calculate parts of coefficients independent of input (a,c,d)
  //
  // This part is called only once so a little less efficient implementation 
  // is employed for readability 
  //

  coefPC0.allocate();  coefPC1.allocate();
  coefMC0.allocate();  coefMC1.allocate();
  coefCP0.allocate();  coefCP1.allocate();  coefCP2.allocate();  
  coefCM0.allocate();  coefCM1.allocate();  coefCM2.allocate();  
  coefCC0.allocate();  coefCC1.allocate(); 
  for(int x=mesh->xstart; x <= mesh->xend; x++) {
    for(int y=mesh->ystart;y<=mesh->yend;y++) {

      BoutReal g11CC = mesh->g11(x,y);
      BoutReal g22CC = mesh->g22(x,y);
      BoutReal g23CC = mesh->g23(x,y);
      BoutReal g_22CC = mesh->g_22(x,y);
      BoutReal g_22CP = mesh->g_22(x,y+1);
      BoutReal g_22CM = mesh->g_22(x,y-1);
      BoutReal g_12CC = mesh->g_12(x,y);
      BoutReal g_23CC = mesh->g_23(x,y);
      BoutReal JCC = mesh->J(x,y);
      BoutReal JCP = mesh->J(x,y+1);
      BoutReal JCM = mesh->J(x,y-1);
      BoutReal G1CC = mesh->G1(x,y);
      BoutReal G2CC = mesh->G2(x,y);

      BoutReal dxCC = mesh->dx(x,y);
      BoutReal dxPC = mesh->dx(x+1,y);
      BoutReal dxMC = mesh->dx(x-1,y);
      BoutReal dyCC = mesh->dy(x,y);

      BoutReal tmp0,tmp1,tmp2,tmp3,tmp4;

      // coefPC0, coefPC1, coefMC0, coefMC1
      tmp0 = g11CC/SQ(dxCC);
      tmp1 = 0.25*(dxPC-dxMC)/dxCC; //zero if mesh is uniform 
      tmp2 = 0.5*G1CC/dxCC;
      
      coefPC0(x,y) = tmp0*(1.0 - tmp1) + tmp2;
      coefPC1(x,y) = 0.5*tmp0;
      coefMC0(x,y) = tmp0*(1.0 + tmp1) - tmp2;
      coefMC1(x,y) =-0.5*tmp0;

      // cCP0, cCP1, cCP2, cCM0, cCM1, cCM2
      tmp0 = 1.0/SQ(dyCC);
      tmp1 = g22CC - 1.0/g_22CC; //zero if mesh is uniform 
      tmp2 = (0.25/JCC)*((JCP/g_22CP) - (JCM/g_22CM)) - 0.5*dyCC*G2CC;
      tmp3 = 0.25*g11CC*g_12CC/(g_22CC*dxCC*dyCC);
      tmp4 = 0.25*g23CC*g_23CC/(g_22CC*dyCC*dyCC);
      
      coefCP0(x,y) = tmp0*(tmp1-tmp2);
      coefCP1(x,y) =-tmp3;
      coefCP2(x,y) =-tmp4;

      coefCM0(x,y) = tmp0*(tmp1+tmp2);
      coefCM1(x,y) = tmp3;
      coefCM2(x,y) = tmp4;

      // cCC0, cCC1
      coefCC0(x,y) = -2.0*g11CC/SQ(dxCC);
      coefCC1(x,y) = -2.0*(g22CC -1.0/g_22CC)/SQ(dyCC);
    }
  }

  // Now communicate to fill guard cells

  mesh->communicate(coefPC0, coefPC1);
  mesh->communicate(coefMC0, coefMC1);
  mesh->communicate(coefCP0, coefCP1, coefCP2);
  mesh->communicate(coefCM0, coefCM1, coefCM2);
  mesh->communicate(coefCC0, coefCC1);

  //////////////////////////////////////////////////
  // Set up KSP
  
  // Declare KSP Context 
  KSPCreate( comm, &ksp ); 
  
  // Configure Linear Solver
  
  bool direct;
  OPTION(opt, direct, false);
  
  if(direct) {
    KSPGetPC(ksp,&pc);
    PCSetType(pc,PCLU);
    PCFactorSetMatSolverPackage(pc,"mumps");
  }else {
    
    // Convergence Parameters. Solution is considered converged if |r_k| < max( rtol * |b| , atol )
    // where r_k = b - Ax_k. The solution is considered diverged if |r_k| > dtol * |b|.
    BoutReal rtol, atol, dtol;
    int maxits; ///< Maximum iterations
    
    OPTION(opt, rtol, 1e-5);     // Relative tolerance 
    OPTION(opt, atol, 1e-10);    // Absolute tolerance
    OPTION(opt, dtol, 1e3);      // Diverged threshold
    OPTION(opt, maxits, 100000); // Maximum iterations
    
    // Get KSP Solver Type
    string ksptype;
    opt->get("ksptype", ksptype, "gmres");
    
    // Get PC type
    string pctype;
    opt->get("pctype", pctype, "none", true);

    KSPSetType( ksp, ksptype.c_str() );
    KSPSetTolerances( ksp, rtol, atol, dtol, maxits );
    
    KSPSetInitialGuessNonzero( ksp, (PetscBool) true );
    
    KSPGetPC(ksp,&pc);
    PCSetType(pc, pctype.c_str());

    if(pctype == "shell") {
      // Using tridiagonal solver as preconditioner
      PCShellSetApply(pc,laplacePCapply);
      PCShellSetContext(pc,this);
      
      bool rightprec;
      OPTION(opt, rightprec, true);
      if(rightprec) {
        KSPSetPCSide(ksp, PC_RIGHT); // Right preconditioning
      }else
        KSPSetPCSide(ksp, PC_LEFT);  // Left preconditioning
    }
  }
  
  KSPSetFromOptions( ksp );

  ///////////////////////////////////////////////////
  // Decide boundary condititions
  if(mesh->periodicY(mesh->xstart)) {
    // Periodic in Y, so in the core
    opt->get("core_bndry_dirichlet", x_inner_dirichlet, false);
    opt->get("sol_bndry_dirichlet",  x_outer_dirichlet, false); //H.SETO
  }else {
    // Non-periodic, so in the PF region
    opt->get("pf_bndry_dirichlet", x_inner_dirichlet, true);
    opt->get("sol_bndry_dirichlet",x_outer_dirichlet, true); //H.SETO
  }
  opt->get("y_bndry_dirichlet", y_bndry_dirichlet, false);

  ///////////////////////////////////////////////////
  // Including Y derivatives?

  OPTION(opt, include_y_derivs, true);
  
  ///////////////////////////////////////////////////
  // Set the default coefficients
  setCoefs(0.0, 1.0, 1.0);
}

void LaplaceXY2::setCoefs(const Field2D &a, const Field2D &c, const Field2D &d) {
  Timer timer("invert");
  
  //////////////////////////////////////////////////
  // Set Matrix elements
  
  for(int x=mesh->xstart; x <= mesh->xend; x++) {
    for(int y=mesh->ystart;y<=mesh->yend;y++) {
      
      // stencil entries
      PetscScalar xcyc, xmyc, xpyc, xcym, xcyp;
      
      BoutReal dCC = d(x,y);
      BoutReal aCC = a(x,y);
      BoutReal cdx = (c(x+1,y)-c(x-1,y))/c(x,y);
      BoutReal cdy = (c(x,y+1)-c(x,y-1))/c(x,y);
      
      // XX component
      xmyc = coefMC0(x,y)*dCC + coefMC1(x,y)*cdx;
      xcyc = coefCC0(x,y)*dCC + aCC;
      xpyc = coefPC0(x,y)*dCC + coefPC1(x,y)*cdx;

      // Put values into the preconditioner, X derivatives only
      acoef[y - mesh->ystart][x - xstart] = xmyc;
      bcoef[y - mesh->ystart][x - xstart] = xcyc;
      ccoef[y - mesh->ystart][x - xstart] = xpyc;
      
      if( include_y_derivs ) {
        // YY component
	xcym  = coefCM0(x,y)*dCC + coefCM1(x,y)*cdx + coefCM2(x,y)*cdy;
	xcyc += coefCC1(x,y)*dCC;
	xcyp  = coefCP0(x,y)*dCC + coefCP1(x,y)*cdx + coefCP2(x,y)*cdy;

      }
      
      /////////////////////////////////////////////////
      // Now have a 5-point stencil for the Laplacian
      
      int row = globalIndex(x,y);
      
      // Set the centre (diagonal)
      MatSetValues(MatA,1,&row,1,&row,&xcyc,INSERT_VALUES);
      
      // X + 1
      int col = globalIndex(x+1, y);
      MatSetValues(MatA,1,&row,1,&col,&xpyc,INSERT_VALUES);
      
      // X - 1
      col = globalIndex(x-1, y);
      MatSetValues(MatA,1,&row,1,&col,&xmyc,INSERT_VALUES);
      
      if( include_y_derivs ) {
        // Y + 1
        col = globalIndex(x, y+1);
        MatSetValues(MatA,1,&row,1,&col,&xcyp,INSERT_VALUES);
        
        // Y - 1
        col = globalIndex(x, y-1);
        MatSetValues(MatA,1,&row,1,&col,&xcym,INSERT_VALUES);
      }
      
    }
  }
  
  // X boundaries
  if(mesh->firstX()) {
    if(x_inner_dirichlet) {

      // Dirichlet on inner X boundary
      for(int y=mesh->ystart;y<=mesh->yend;y++) {
        int row = globalIndex(mesh->xstart-1,y);
        PetscScalar val = 0.5;
        MatSetValues(MatA,1,&row,1,&row,&val,INSERT_VALUES);
        
        int col = globalIndex(mesh->xstart,y);
        MatSetValues(MatA,1,&row,1,&col,&val,INSERT_VALUES);
        
        // Preconditioner
        bcoef[y-mesh->ystart][0] = 0.5;
        ccoef[y-mesh->ystart][0] = 0.5;
      }
      
    }else {
      
      // Neumann on inner X boundary
      for(int y=mesh->ystart;y<=mesh->yend;y++) {
        int row = globalIndex(mesh->xstart-1,y);
        PetscScalar val = 1.0;
        MatSetValues(MatA,1,&row,1,&row,&val,INSERT_VALUES);
        
        int col = globalIndex(mesh->xstart,y);
        val = -1.0;
        MatSetValues(MatA,1,&row,1,&col,&val,INSERT_VALUES);
        
        // Preconditioner
        bcoef[y-mesh->ystart][0] =  1.0;
        ccoef[y-mesh->ystart][0] = -1.0;
      }
    }
  }
  if(mesh->lastX()) {//H.SETO
    if(x_outer_dirichlet) { 
      // Dirichlet on outer X boundary
    
      for(int y=mesh->ystart;y<=mesh->yend;y++) {
	int row = globalIndex(mesh->xend+1,y);
	PetscScalar val = 0.5;
	MatSetValues(MatA,1,&row,1,&row,&val,INSERT_VALUES);
      
	int col = globalIndex(mesh->xend,y);
	MatSetValues(MatA,1,&row,1,&col,&val,INSERT_VALUES);
      
	// Preconditioner
	acoef[y-mesh->ystart][mesh->xend+1 - xstart] = 0.5;
	bcoef[y-mesh->ystart][mesh->xend+1 - xstart] = 0.5;
      }
    } else {
      // Neumann on outer X boundary
      for(int y=mesh->ystart;y<=mesh->yend;y++) {
        int row = globalIndex(mesh->xend+1,y);
        PetscScalar val = 1.0;
        MatSetValues(MatA,1,&row,1,&row,&val,INSERT_VALUES);
        
        int col = globalIndex(mesh->xend,y);
        val = -1.0;
        MatSetValues(MatA,1,&row,1,&col,&val,INSERT_VALUES);
        
        // Preconditioner
        acoef[y-mesh->ystart][mesh->xend+1 - xstart] = -1.0;
        bcoef[y-mesh->ystart][mesh->xend+1 - xstart] =  1.0;
      }
    }
  }

  if(y_bndry_dirichlet) {
    // Dirichlet on Y boundaries
    for(RangeIterator it=mesh->iterateBndryLowerY(); !it.isDone(); it++) {
      int row = globalIndex(it.ind, mesh->ystart-1);
      PetscScalar val = 0.5;
      MatSetValues(MatA,1,&row,1,&row,&val,INSERT_VALUES);
      
      int col = globalIndex(it.ind, mesh->ystart);
      MatSetValues(MatA,1,&row,1,&col,&val,INSERT_VALUES);
    }
    
    for(RangeIterator it=mesh->iterateBndryUpperY(); !it.isDone(); it++) {
      int row = globalIndex(it.ind, mesh->yend+1);
      PetscScalar val = 0.5;
      MatSetValues(MatA,1,&row,1,&row,&val,INSERT_VALUES);
      
      int col = globalIndex(it.ind, mesh->yend);
      MatSetValues(MatA,1,&row,1,&col,&val,INSERT_VALUES);
    }
  }else {
    // Neumann on Y boundaries
    for(RangeIterator it=mesh->iterateBndryLowerY(); !it.isDone(); it++) {
      int row = globalIndex(it.ind, mesh->ystart-1);
      PetscScalar val = 1.0;
      MatSetValues(MatA,1,&row,1,&row,&val,INSERT_VALUES);
    
      val = -1.0;
      int col = globalIndex(it.ind, mesh->ystart);
      
      MatSetValues(MatA,1,&row,1,&col,&val,INSERT_VALUES);
    }
    
    for(RangeIterator it=mesh->iterateBndryUpperY(); !it.isDone(); it++) {
      int row = globalIndex(it.ind, mesh->yend+1);
      PetscScalar val = 1.0;
      MatSetValues(MatA,1,&row,1,&row,&val,INSERT_VALUES);
      
      val = -1.0;
      int col = globalIndex(it.ind, mesh->yend);
      MatSetValues(MatA,1,&row,1,&col,&val,INSERT_VALUES);
    }
  }
  
  // Assemble Matrix
  MatAssemblyBegin( MatA, MAT_FINAL_ASSEMBLY );
  MatAssemblyEnd( MatA, MAT_FINAL_ASSEMBLY );

  // Set the operator
#if PETSC_VERSION_GE(3,5,0)
  KSPSetOperators( ksp,MatA,MatA );
#else
  KSPSetOperators( ksp,MatA,MatA,DIFFERENT_NONZERO_PATTERN );
#endif
  
  // Set coefficients for preconditioner
  cr->setCoefs(nsys, acoef, bcoef, ccoef);
}

LaplaceXY2::~LaplaceXY2() {
  KSPDestroy( &ksp );
  VecDestroy( &xs );
  VecDestroy( &bs );
  MatDestroy( &MatA );

  // Free memory for preconditioner
  free_matrix(acoef);
  free_matrix(bcoef);
  free_matrix(ccoef);
  free_matrix(xvals);
  free_matrix(bvals);
  
  // Delete tridiagonal solver
  delete cr;
}

const Field2D LaplaceXY2::solve(const Field2D &rhs, const Field2D &x0) {
  Timer timer("invert");
  
  // Load initial guess x0 into xs and rhs into bs
  
  for(int x=mesh->xstart;x<= mesh->xend;x++) {
    for(int y=mesh->ystart;y<=mesh->yend;y++) {
      int ind = globalIndex(x,y);
      
      PetscScalar val = x0(x,y);
      VecSetValues( xs, 1, &ind, &val, INSERT_VALUES );
      
      val = rhs(x,y);
      VecSetValues( bs, 1, &ind, &val, INSERT_VALUES );
    }
  }

  if(mesh->firstX()) {
    if(x_inner_dirichlet) {
      for(int y=mesh->ystart;y<=mesh->yend;y++) {
        int ind = globalIndex(mesh->xstart-1,y);
      
        PetscScalar val = x0(mesh->xstart-1,y);
        VecSetValues( xs, 1, &ind, &val, INSERT_VALUES );
        
        val = 0.5*(x0(mesh->xstart-1,y) + x0(mesh->xstart,y));
        VecSetValues( bs, 1, &ind, &val, INSERT_VALUES );
      }
    }else {
      // Inner X boundary (Neumann)
      for(int y=mesh->ystart;y<=mesh->yend;y++) {
        int ind = globalIndex(mesh->xstart-1,y);
        
        //PetscScalar val = x0(mesh->xstart-1,y);
	PetscScalar val = x0(mesh->xstart,y);
        VecSetValues( xs, 1, &ind, &val, INSERT_VALUES );
        
        val = 0.0; //x0(mesh->xstart-1,y) - x0(mesh->xstart,y);
        VecSetValues( bs, 1, &ind, &val, INSERT_VALUES );
      }
    }
  }
  
  if(mesh->lastX()) { //H.SETO
    if(x_outer_dirichlet){  
      //Outer X boundary (Dirichlet)
      for(int y=mesh->ystart;y<=mesh->yend;y++) {
	int ind = globalIndex(mesh->xend+1,y);
	
	PetscScalar val = x0(mesh->xend+1,y);
	VecSetValues( xs, 1, &ind, &val, INSERT_VALUES );
	
	val = 0.5*(x0(mesh->xend,y) + x0(mesh->xend+1,y));
	VecSetValues( bs, 1, &ind, &val, INSERT_VALUES );
      }
    } else {
      // Outer X boundary (Neumann)
      for(int y=mesh->ystart;y<=mesh->yend;y++) {
        int ind = globalIndex(mesh->xend+1,y);
        
        //PetscScalar val = x0(mesh->xend+1,y);
	PetscScalar val = x0(mesh->xend,y);
        VecSetValues( xs, 1, &ind, &val, INSERT_VALUES );
        
        val = 0.0; //x0(mesh->xend+1,y) - x0(mesh->xend,y);
        VecSetValues( bs, 1, &ind, &val, INSERT_VALUES );
      }
    }
  }

  if(y_bndry_dirichlet) {
    for(RangeIterator it=mesh->iterateBndryLowerY(); !it.isDone(); it++) {
      int ind = globalIndex(it.ind, mesh->ystart-1);
    
      PetscScalar val = x0(it.ind,mesh->ystart-1);
      VecSetValues( xs, 1, &ind, &val, INSERT_VALUES );
    
      val = 0.5*(x0(it.ind, mesh->ystart-1) + x0(it.ind, mesh->ystart));
      VecSetValues( bs, 1, &ind, &val, INSERT_VALUES );
    }
  
    for(RangeIterator it=mesh->iterateBndryUpperY(); !it.isDone(); it++) {
      int ind = globalIndex(it.ind, mesh->yend+1);
    
      PetscScalar val = x0(it.ind,mesh->yend+1);
      VecSetValues( xs, 1, &ind, &val, INSERT_VALUES );
    
      val = 0.5*(x0(it.ind, mesh->yend+1) + x0(it.ind, mesh->yend));
      VecSetValues( bs, 1, &ind, &val, INSERT_VALUES );
    }
  } else {
    // Y boundaries Neumann
    for(RangeIterator it=mesh->iterateBndryLowerY(); !it.isDone(); it++) {
      int ind = globalIndex(it.ind, mesh->ystart-1);
    
      PetscScalar val = x0(it.ind,mesh->ystart-1);
      VecSetValues( xs, 1, &ind, &val, INSERT_VALUES );
      
      val = 0.0;
      VecSetValues( bs, 1, &ind, &val, INSERT_VALUES );
    }
  
    for(RangeIterator it=mesh->iterateBndryUpperY(); !it.isDone(); it++) {
      int ind = globalIndex(it.ind, mesh->yend+1);
      
      PetscScalar val = x0(it.ind,mesh->yend+1);
      VecSetValues( xs, 1, &ind, &val, INSERT_VALUES );
    
      val = 0.0;
      VecSetValues( bs, 1, &ind, &val, INSERT_VALUES );
    }
  }
  
  // Assemble RHS Vector
  VecAssemblyBegin(bs);
  VecAssemblyEnd(bs);

  // Assemble Trial Solution Vector
  VecAssemblyBegin(xs);
  VecAssemblyEnd(xs);
  
  // Solve the system
  KSPSolve( ksp, bs, xs );
  
  KSPConvergedReason reason;
  KSPGetConvergedReason( ksp, &reason );
  
  if(reason <= 0) {
    throw BoutException("LaplaceXY failed to converge. Reason %d", reason);
  }
  
  //////////////////////////
  // Copy data into result
  
  Field2D result;
  result.allocate();
  
  for(int x=mesh->xstart;x<= mesh->xend;x++) {
    for(int y=mesh->ystart;y<=mesh->yend;y++) {
      int ind = globalIndex(x,y);
      
      PetscScalar val;
      VecGetValues(xs, 1, &ind, &val );
      result(x,y) = val;
    }
  }
  
  // Inner X boundary
  if(mesh->firstX()) {
    for(int y=mesh->ystart;y<=mesh->yend;y++) {
      int ind = globalIndex(mesh->xstart-1,y);
      PetscScalar val;
      VecGetValues(xs, 1, &ind, &val );
      for(int x=mesh->xstart-1; x >= 0; x--)
        result(x,y) = val;
    }
  }

  // Outer X boundary
  if(mesh->lastX()) {
    for(int y=mesh->ystart;y<=mesh->yend;y++) {
      int ind = globalIndex(mesh->xend+1,y);
      PetscScalar val;
      VecGetValues(xs, 1, &ind, &val );
      for(int x=mesh->xend+1;x < mesh->ngx;x++)
        result(x,y) = val;
    }
  }
  
  // Lower Y boundary
  for(RangeIterator it=mesh->iterateBndryLowerY(); !it.isDone(); it++) {
    int ind = globalIndex(it.ind, mesh->ystart-1);
    PetscScalar val;
    VecGetValues(xs, 1, &ind, &val );
    for(int y=mesh->ystart-1;y>=0;y--)
      result(it.ind, y) = val;
  }
  
  // Upper Y boundary
  for(RangeIterator it=mesh->iterateBndryUpperY(); !it.isDone(); it++) {
    int ind = globalIndex(it.ind, mesh->yend+1);
    PetscScalar val;
    VecGetValues(xs, 1, &ind, &val );
    for(int y=mesh->yend+1;y<mesh->ngy;y++)
      result(it.ind, y) = val;
  }
  
  return result;
}

/*! Preconditioner
 * NOTE: For efficiency, this routine does not use globalIndex() 
 * in the inner loop. Instead, the indexing must be ordered in
 * exactly the same way as in the construction of indexXY
 */
int LaplaceXY2::precon(Vec input, Vec result) {
  
  // Starting index
  int ind = -1;
  
  RangeIterator itdwn=mesh->iterateBndryLowerY();
  if(!itdwn.isDone()) {
    ind = globalIndex(itdwn.ind, mesh->ystart-1);
    
    for(; !itdwn.isDone(); itdwn++) {
      PetscScalar val;
      VecGetValues(input, 1, &ind, &val ); 
      VecSetValues(result, 1, &ind, &val, INSERT_VALUES );
      ind++;
    }
  }
  RangeIterator itup=mesh->iterateBndryUpperY();
  if(!itup.isDone()) {
    if(ind == -1) {
      // No lower boundary
      ind = globalIndex(itup.ind, mesh->yend+1);
    }
    for(; !itup.isDone(); itup++) {
      PetscScalar val;
      VecGetValues(input, 1, &ind, &val ); 
      VecSetValues(result, 1, &ind, &val, INSERT_VALUES );
      ind++;
    }
  }
  if(ind == -1) {
    // No Y boundaries
    ind = globalIndex(xstart, mesh->ystart);
  }
    
  int ind0 = ind;
  // Load vector x into bvals array
  for(int x=xstart;x<=xend;x++) {
    for(int y=mesh->ystart; y<=mesh->yend;y++) {
      PetscScalar val;
      VecGetValues(input, 1, &ind, &val ); 
      bvals[y-mesh->ystart][x-xstart] = val;
      ind++;
    }
  }
  
  // Solve tridiagonal systems using CR solver
  cr->solve(nsys, bvals, xvals);

  // Save result xvals into y array
  ind = ind0;
  for(int x=xstart;x<=xend;x++) {
    for(int y=mesh->ystart; y<=mesh->yend;y++) {
      PetscScalar val = xvals[y-mesh->ystart][x-xstart];
      VecSetValues(result, 1, &ind, &val, INSERT_VALUES );
      ind++;
    }
  }
  VecAssemblyBegin(result);
  VecAssemblyEnd(result);
  return 0;
}

///////////////////////////////////////////////////////////////

int LaplaceXY2::localSize() {
  
  // Bulk of points
  int nx = mesh->xend - mesh->xstart + 1;
  int ny = mesh->yend - mesh->ystart + 1;
  
  int n = nx * ny;  
  
  // X boundaries
  if(mesh->firstX())
    n += ny;
  if(mesh->lastX())
    n += ny;
  
  // Y boundaries
  for(RangeIterator it=mesh->iterateBndryLowerY(); !it.isDone(); it++) {
    n++;
  }
  for(RangeIterator it=mesh->iterateBndryUpperY(); !it.isDone(); it++) {
    n++;
  }
  
  return n;
}

int LaplaceXY2::globalIndex(int x, int y) {
  if( (x < 0) || (x >= mesh->ngx) ||
      (y < 0) || (y >= mesh->ngy) )
    return -1; // Out of range
 
  // Get the index from a Field2D, round to integer
  return roundInt(indexXY(x,y));
}

int LaplaceXY2::roundInt(BoutReal f) {
  if(f > 0.0) {
    return (int) (f + 0.5);
  }
  return (int) (f - 0.5);
}

#endif // BOUT_HAS_PETSC
