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
	BlockCGSolver

Description
	Preconditioned Conjugate Gradient solver.

Author
	Hrvoje Jasak, Wikki Ltd.  All rights reserved

SourceFiles
	BlockCGSolver.C

\*---------------------------------------------------------------------------*/

#ifndef BlockCGSolver_H
#define BlockCGSolver_H

#include "blockLduSolvers.H"
#include "BlockIterativeSolver.H"
#include "blockLduPrecons.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


template<class Type>
class BlockCGSolver
:
	public BlockIterativeSolver<Type>
{
	// Private data

		//- Preconditioner
		autoPtr<BlockLduPrecon<Type> > preconPtr_;


	// Private Member Functions

		//- Disallow default bitwise copy construct
		BlockCGSolver(const BlockCGSolver<Type>&);

		//- Disallow default bitwise assignment
		void operator=(const BlockCGSolver<Type>&);


public:

	//- Runtime type information
	TypeName("CG");


	// Constructors

		//- Construct from matrix components and solver data stream
		BlockCGSolver
		(
			const word& fieldName,
			const BlockLduMatrix<Type>& matrix,
			const dictionary& dict
		);


	//- Destructor
	virtual ~BlockCGSolver()
	{}


	// Member Functions

		//- Solve the matrix with this solver
		virtual BlockSolverPerformance<Type> solve
		(
			Field<Type>& x,
			const Field<Type>& b
		);

		//- Re-initialise solver after matrix coefficient update
		virtual void initMatrix()
		{
			preconPtr_->initMatrix();
		}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#	include "BlockCGSolver.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
