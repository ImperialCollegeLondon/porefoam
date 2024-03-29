/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
	This file is part of foam-extend.

	foam-extend is free software: you can redistribute it and/or modify it
	under the terms of the GNU General Public License as published by the
	Free Software Foundation, either version 3 of the License, or (at your
	option) any later version.

	foam-extend is distributed in the hope that it will be useful, but
	WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
	General Public License for more details.

	You should have received a copy of the GNU General Public License
	along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Class
	Foam::GAMGSolver

Description
	Geometric agglomerated algebraic multigrid solver.

  Characteristics:
	  - Requires positive definite, diagonally dominant matrix.
	  - Agglomeration algorithm: selectable and optionally cached.
	  - Restriction operator: summation.
	  - Prolongation operator: injection.
	  - Smoother: Gauss-Seidel.
	  - Coarse matrix creation: central coefficient: summation of fine grid
		central coefficients with the removal of intra-cluster face;
		off-diagonal coefficient: summation of off-diagonal faces.
	  - Coarse matrix scaling: performed by correction scaling, using steepest
		descent optimisation.
	  - Type of cycle: V-cycle with optional pre-smoothing.
	  - Coarsest-level matrix solved using ICCG or BICCG.

SourceFiles
	GAMGSolver.C
	GAMGSolverCalcAgglomeration.C
	GAMGSolverMakeCoarseMatrix.C
	GAMGSolverOperations.C
	GAMGSolverSolve.C

\*---------------------------------------------------------------------------*/

#ifndef GAMGSolver_H
#define GAMGSolver_H

#include "GAMGAgglomeration.H"
#include "lduMatrix.H"
#include "labelField.H"
#include "primitiveFields.H"
#include "LUscalarMatrix.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


class GAMGSolver
:
	public lduMatrix::solver
{
	// Private data

		Switch cacheAgglomeration_;

		//- Number of pre-smoothing sweeps
		label nPreSweeps_;

		//- Number of post-smoothing sweeps
		label nPostSweeps_;

		//- Number of smoothing sweeps on finest mesh
		label nFinestSweeps_;

		//- Choose if the corrections should be scaled.
		//  By default corrections for symmetric matrices are scaled
		//  but not for asymmetric matrices.
		Switch scaleCorrection_;

		//- Direct or iteratively solve the coarsest level
		Switch directSolveCoarsest_;

		//- The agglomeration
		const GAMGAgglomeration& agglomeration_;

		//- Hierarchy of matrix levels
		PtrList<lduMatrix> matrixLevels_;

		//- Hierarchy of interfaces.
		//  Warning: Needs to be deleted explicitly.
		PtrList<lduInterfaceFieldPtrsList> interfaceLevels_;

		//- Hierarchy of interface boundary coefficients
		PtrList<FieldField<Field, scalar> > coupleLevelsBouCoeffs_;

		//- Hierarchy of interface internal coefficients
		PtrList<FieldField<Field, scalar> > coupleLevelsIntCoeffs_;

		//- LU decompsed coarsest matrix
		autoPtr<LUscalarMatrix> coarsestLUMatrixPtr_;


	// Private Member Functions

		//- Read control parameters from the control dictionary
		virtual void readControls();

		//- Make agglomeration.  Constructor helper.  HJ, 21/Oct/2008
		void makeAgglomeration();

		//- Simplified access to interface level
		const lduInterfaceFieldPtrsList& interfaceLevel
		(
			const label i
		) const;

		//- Simplified access to matrix level
		const lduMatrix& matrixLevel(const label i) const;

		//- Simplified access to interface boundary coeffs level
		const FieldField<Field, scalar>& coupleBouCoeffsLevel
		(
			const label i
		) const;

		//- Simplified access to interface internal coeffs level
		const FieldField<Field, scalar>& coupleIntCoeffsLevel
		(
			const label i
		) const;

		//- Agglomerate coarse matrix
		void agglomerateMatrix(const label fineLevelIndex);

		//- Calculate and return the scaling factor from Acf, coarseSource
		//  and coarseField.
		//  At the same time do a Jacobi iteration on the coarseField using
		//  the Acf provided after the coarseField values are used for the
		//  scaling factor.
		scalar scalingFactor
		(
			scalarField& field,
			const scalarField& source,
			const scalarField& Acf,
			const scalarField& D
		) const;

		//- Calculate Acf and calculate and return the scaling factor.
		scalar scalingFactor
		(
			scalarField& Acf,
			const lduMatrix& A,
			scalarField& field,
			const FieldField<Field, scalar>& interfaceLevelBouCoeffs,
			const lduInterfaceFieldPtrsList& interfaceLevel,
			const scalarField& source,
			const direction cmpt
		) const;


		//- Initialise the data structures for the V-cycle
		void initVcycle
		(
			PtrList<scalarField>& coarseCorrX,
			PtrList<scalarField>& coarseB,
			PtrList<lduMatrix::smoother>& smoothers
		) const;


		//- Perform a single GAMG V-cycle with pre, post and finest smoothing.
		void Vcycle
		(
			const PtrList<lduMatrix::smoother>& smoothers,
			scalarField& x,
			const scalarField& b,
			scalarField& Ax,
			scalarField& finestCorrection,
			scalarField& finestResidual,
			PtrList<scalarField>& coarseCorrX,
			PtrList<scalarField>& coarseB,
			const direction cmpt = 0
		) const;


		//- Solve the coarsest level with either an iterative or direct solver
		void solveCoarsestLevel
		(
			scalarField& coarsestCorrField,
			const scalarField& coarsestSource
		) const;


public:

	friend class GAMGPreconditioner;

	//- Runtime type information
	TypeName("GAMG");


	// Constructors

		//- Construct from lduMatrix and solver controls
		GAMGSolver
		(
			const word& fieldName,
			const lduMatrix& matrix,
			const FieldField<Field, scalar>& coupleBouCoeffs,
			const FieldField<Field, scalar>& coupleIntCoeffs,
			const lduInterfaceFieldPtrsList& interfaces,
			const dictionary& dict
		);


	// Destructor

		virtual ~GAMGSolver();


	// Member Functions

		//- Solve
		virtual lduSolverPerformance solve
		(
			scalarField& x,
			const scalarField& b,
			const direction cmpt = 0
		) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
