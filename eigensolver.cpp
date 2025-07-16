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
//#include "AnasaziLOBPCGSolMgr.hpp"
#include "AnasaziBlockDavidsonSolMgr.hpp"
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
typedef double                                   Scalar;
typedef Tpetra::MultiVector<Scalar>              TMV;
typedef Tpetra::Vector<Scalar>                   Vector;
typedef Tpetra::Operator<Scalar>                 TOP;

int
main (int argc, char* argv[])
{
  typedef Anasazi::BasicEigenproblem<Scalar,TMV,TOP> Problem;
  typedef Anasazi::MultiVecTraits<Scalar, TMV> TMVT;
  typedef Anasazi::OperatorTraits<Scalar, TMV, TOP> TOPT;
  typedef Tpetra::CrsMatrix<> CrsMatrix;
  typedef Tpetra::MatrixMarket::Reader<CrsMatrix>  Reader;
  Tpetra::ScopeGuard tpetraScope (&argc, &argv);
  //
  // Get the default communicator
  //
  RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm ();
  const int myRank = comm->getRank();
  // Read the command line arguments
  //std::string fileA ("/u/slotnick_s2/aklinvex/matrices/anderson4.mtx");
  //Teuchos::CommandLineProcessor cmdp (false, true);
  //cmdp.setOption ("fileA", &fileA, "Filename for the Matrix-Market stiffness matrix.");
  //if (cmdp.parse (argc,argv) != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL) {
  //  return -1;
  //}
  // Get the matrix
  //RCP<const CrsMatrix> A = Reader::readSparseFile (fileA, comm);

  // The number of rows and columns in the matrix.
    const Tpetra::global_size_t numGblIndices = 100;
    // Construct a Map that puts approximately the same number of
    // equations on each processor.
    const Tpetra::Vector<>::global_ordinal_type indexBase = 0;
    RCP<const Tpetra::Map<>> map =
      rcp (new Tpetra::Map<> (numGblIndices, indexBase, comm));
    const size_t numMyElements = map->getLocalNumElements ();
    // If you like, you may get the list of global indices that the
    // calling process owns.  This is unnecessary if you don't mind
    // converting local indices to global indices.
    //
    // ArrayView<const global_ordinal_type> myGlobalElements =
    //   map->getLocalElementList ();
    if (myRank == 0) {
      cout << endl << "Creating the sparse matrix" << endl;
    }
    // Create a Tpetra sparse matrix whose rows have distribution
    // given by the Map.  We expect at most three entries per row.
    RCP<Tpetra::CrsMatrix<>> A (new Tpetra::CrsMatrix<> (map, 3));
    // Fill the sparse matrix, one row at a time.
    const Scalar two = static_cast<Scalar> (2.0);
    const Scalar negOne = static_cast<Scalar> (-1.0);
    for (Tpetra::Vector<>::local_ordinal_type lclRow = 0;
   lclRow < static_cast<Tpetra::Vector<>::local_ordinal_type> (numMyElements);
   ++lclRow) {
      const Tpetra::Vector<>::global_ordinal_type gblRow = map->getGlobalElement (lclRow);
      // A(0, 0:1) = [2, -1]
      if (gblRow == 0) {
  A->insertGlobalValues (gblRow,
             Teuchos::tuple<Tpetra::Vector<>::global_ordinal_type> (gblRow, gblRow + 1),
             Teuchos::tuple<Scalar> (two, negOne));
      }
      // A(N-1, N-2:N-1) = [-1, 2]
      else if (static_cast<Tpetra::global_size_t> (gblRow) == numGblIndices - 1) {
  A->insertGlobalValues (gblRow,
             Teuchos::tuple<Tpetra::Vector<>::global_ordinal_type> (gblRow - 1, gblRow),
             Teuchos::tuple<Scalar> (negOne, two));
      }
      // A(i, i-1:i+1) = [-1, 2, -1]
      else {
  A->insertGlobalValues (gblRow,
             Teuchos::tuple<Tpetra::Vector<>::global_ordinal_type> (gblRow - 1, gblRow, gblRow + 1),
             Teuchos::tuple<Scalar> (negOne, two, negOne));
      }
    }
    // Tell the sparse matrix that we are done adding entries to it.
    A->fillComplete ();

  // set parameters for solver
  int blockSize = 40;
  double tol = 1e-5;
  bool scaled = true;
  int nev = 25;
  Teuchos::ParameterList MyPL;
  MyPL.set ("Which", "SR");
  MyPL.set ("Maximum Restarts", 10000);
  MyPL.set ("Maximum Iterations", 1000000);
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
  // Create the eigensolver and give it your problem and parameters.
  Anasazi::BlockDavidsonSolMgr<Scalar,TMV,TOP> solver (MyProblem, MyPL);
  // Tell the solver to solve the eigenproblem.
  Anasazi::ReturnType returnCode = solver.solve ();
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
    std::vector<Scalar> normR (sol.numVecs);
    TMV Avec (A->getRowMap (), TMVT::GetNumberVecs (*evecs));
    TOPT::Apply (*A, *evecs, Avec);
    Teuchos::SerialDenseMatrix<int,Scalar> T (numev, numev);
    TMVT::MvTransMv (1.0, Avec, *evecs, T);
    TMVT::MvTimesMatAddMv (-1.0, *evecs, T, 1.0, Avec);
    TMVT::MvNorm (Avec, normR);
    if (myRank == 0) {
      cout.setf(std::ios_base::right, std::ios_base::adjustfield);
      cout<<"Actual Eigenvalues: "<<std::endl;
      cout<<"------------------------------------------------------"<<std::endl;
      cout<<std::setw(16)<<"Real Part"
        <<std::setw(16)<<"Error"<<std::endl;
      cout<<"------------------------------------------------------"<<std::endl;
      for (int i=0; i<numev; i++) {
        cout<<std::setw(16)<<T(i,i)
          <<std::setw(16)<<normR[i]/std::abs(T(i,i))
          <<std::endl;
      }
      cout<<"------------------------------------------------------"<<std::endl;
    }
  }
  return 0;
}
