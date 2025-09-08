// @HEADER
// *****************************************************************************
//                 Anasazi: Block Eigensolvers Package
//
// Copyright 2004 NTESS and the Anasazi contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
// This example shows how to define a custom operator and status test for Anasazi
#include "AnasaziBasicEigenproblem.hpp"
#include "AnasaziTpetraAdapter.hpp"
#include "AnasaziStatusTestDecl.hpp"
#include "AnasaziLOBPCGSolMgr.hpp"
#include "AnasaziBlockDavidsonSolMgr.hpp"
#include "AnasaziBlockKrylovSchurSolMgr.hpp"
#include "AnasaziTypes.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Tpetra_Core.hpp"
#include "MatrixMarket_Tpetra.hpp"
// These tell the compiler which namespace contains RCP, cout, etc
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::ArrayView;
using std::cout;
using std::endl;
// These define convenient aliases
typedef std::complex<double>                  Scalar;
typedef Tpetra::MultiVector<Scalar>              TMV;
typedef Tpetra::Vector<Scalar>                   Vector;
typedef Tpetra::Operator<Scalar>                 TOP;

#include <typeinfo>


int
main (int argc, char* argv[])
{
  typedef Anasazi::BasicEigenproblem<Scalar,TMV,TOP> Problem;
  typedef Anasazi::MultiVecTraits<Scalar, TMV> TMVT;
  typedef Anasazi::OperatorTraits<Scalar, TMV, TOP> TOPT;
  typedef Tpetra::CrsMatrix<Scalar> CrsMatrix;
  typedef Tpetra::MatrixMarket::Reader<CrsMatrix>  Reader;
  Tpetra::ScopeGuard tpetraScope (&argc, &argv);
  //
  // Get the default communicator
  //
  RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm ();
  const int myRank = comm->getRank();


  enum solver_types { BlockKrylovSchur, BlockDavidson, LOBPCG };
  solver_types solver_types_opt[] = { BlockKrylovSchur, BlockDavidson, LOBPCG };
  const char* solver_names[] = { "BlockKrylovSchur", "BlockDavidson", "LOBPCG" };
  // Read the command line arguments
  std::string filename ("matrix_export.mtx");
  solver_types solver_enum = solver_types::BlockKrylovSchur;
  std::string which ("SR");
  int blockSize = 60;
  double tol = 1e-15;
  bool scaled = false;
  int nev = 3;
  int max_restarts = 10000;
  int max_iterations = 1000000;

  Teuchos::CommandLineProcessor cmdp (false, true);
  cmdp.setOption ("filename", &filename, "Filename for the Matrix-Market stiffness matrix. Default: matrix_export.mtx");
  cmdp.setOption ("solver", &solver_enum, 3, solver_types_opt, solver_names, "Solver type, one of: BlockKrylovSchur (default), BlockDavidson, LOBPCG");
  cmdp.setOption ("which", &which, "Which ev are requested, one of: SR(default), LR, SM, LM");
  cmdp.setOption ("blockSize", &blockSize, "BlockSize of solver. Needs to be > nev. Default: 60");
  cmdp.setOption ("tol", &tol, "Tolerance to stop iterating. Default: 1e-16");
  cmdp.setOption ("nev", &nev, "Number of ev calculated. Default: 3");
  cmdp.setOption ("max_restarts", &max_restarts, "Maximum number of restarts. Default: 10000");
  cmdp.setOption ("max_iterations", &max_iterations, "Maximum number of iterations. Default: 1000000");

  if (cmdp.parse (argc,argv) != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL) {
    return -1;
  }
  
  // Get the matrix
  RCP<const CrsMatrix> A = Reader::readSparseFile (filename, comm);

  // set parameters for solver
  Teuchos::ParameterList MyPL;
  MyPL.set ("Which", which);
  MyPL.set ("Maximum Restarts", max_restarts);
  MyPL.set ("Maximum Iterations", max_iterations);
  MyPL.set ("Block Size", blockSize);
  MyPL.set ("Convergence Tolerance", tol );   // How small do the residuals have to be
  MyPL.set ("Relative Convergence Tolerance", scaled);
  MyPL.set ("Relative Locking Tolerance", scaled);
  MyPL.set ("Verbosity", Anasazi::TimingDetails);
  // Create a MultiVector for an initial subspace to start the solver.
  RCP<TMV> ivec = rcp (new TMV (A->getRowMap (), blockSize));
  TMVT::MvRandom (*ivec);
  // Create the eigenproblem
  RCP<Problem> MyProblem = rcp (new Problem (A, ivec));
  // Inform the eigenproblem that the matrix pencil (A,M) is symmetric
  MyProblem->setHermitian (true);
  // Set the number of eigenvalues requested
  MyProblem->setNEV (nev);
  // Tell the problem that you are finished passing it information
  MyProblem->setProblem ();
  Anasazi::ReturnType returnCode=Anasazi::Unconverged;

  // Create the eigensolver and give it your problem and parameters.
  // Tell the solver to solve the eigenproblem.
  if(solver_enum == solver_types::BlockKrylovSchur)
  {
  		Anasazi::BlockKrylovSchurSolMgr<Scalar, TMV, TOP> solver(MyProblem, MyPL);
  		returnCode = solver.solve ();
  }
  else if (solver_enum == solver_types::BlockDavidson)
  {
  		Anasazi::BlockDavidsonSolMgr<Scalar,TMV,TOP> solver (MyProblem, MyPL);
  		returnCode = solver.solve ();
  }
  else if (solver_enum == solver_types::LOBPCG)
  {
  		Anasazi::LOBPCGSolMgr<Scalar,TMV,TOP> solver (MyProblem, MyPL);
  		returnCode = solver.solve ();
  }
  else 
		std::cerr << "Unknown solver type specified. See --help \n";

  if (returnCode != Anasazi::Converged && myRank == 0) {
    cout << "The solve did NOT converge." << endl;
  } else if (myRank == 0) {
    cout << "The solve converged." << endl;
  }
  // Get the eigenvalues and eigenvectors from the eigenproblem
  Anasazi::Eigensolution<Scalar,TMV> sol = MyProblem->getSolution ();
  std::vector<Anasazi::Value<Scalar> > evals = sol.Evals;
  RCP<TMV> evecs = sol.Evecs;
  int numev = sol.numVecs;
  // Compute the residual, just as a precaution
  if (numev > 0) {
 //   std::vector<Scalar> normR (sol.numVecs);
 //   TMV Avec (A->getRowMap (), TMVT::GetNumberVecs (*evecs));
 //   TOPT::Apply (*A, *evecs, Avec);
 //   Teuchos::SerialDenseMatrix<int,Scalar> T (numev, numev);
 //   TMVT::MvTransMv (1.0, Avec, *evecs, T);
 //   TMVT::MvTimesMatAddMv (-1.0, *evecs, T, 1.0, Avec);
 //   TMVT::MvNorm (Avec, normR);
 //   if (myRank == 0) {
 //     cout.setf(std::ios_base::right, std::ios_base::adjustfield);
 //     cout<<"Actual Eigenvalues: "<<std::endl;
 //     cout<<"------------------------------------------------------"<<std::endl;
 //     cout<<std::setw(16)<<"Real Part"
 //       <<std::setw(16)<<"Error"<<std::endl;
 //     cout<<"------------------------------------------------------"<<std::endl;
 //     for (int i=0; i<numev; i++) {
 //       cout<<std::setw(16)<<T(i,i)
 //         <<std::setw(16)<<normR[i]/std::abs(T(i,i))
 //         <<std::endl;
 //     }
 //     cout<<"------------------------------------------------------"<<std::endl;
 //   }
     for (int i=0; i<numev; i++) {
      std::cout << "Eigenvalues real " << evals[i].realpart << " imaginary " << evals[i].imagpart << std::endl;
     }
  }
  return 0;
}
