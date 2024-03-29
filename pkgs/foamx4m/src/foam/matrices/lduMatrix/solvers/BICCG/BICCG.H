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
	Foam::BICCG

Description
	Diagonal incomplete LU preconditioned BiCG solver derived from the general
	preconditioned BiCG solver PBiCG but with the choice of preconditioner
	pre-selected.

Deprecated
	This solver is present for backward-compatibility and the PBiCG solver
	should be used instead.

SourceFiles
	BICCG.C

\*---------------------------------------------------------------------------*/

#ifndef BICCG_H
#define BICCG_H

#include "PBiCG.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


class BICCG
:
	public PBiCG
{
	// Private Member Functions

		//- Disallow default bitwise copy construct
		BICCG(const BICCG&);

		//- Disallow default bitwise assignment
		void operator=(const BICCG&);

public:

	//- Return the dictionary constructed from the components.
	//  Needed for backward compatibility
	static dictionary solverDict
	(
		const scalar tol,
		const scalar relTol
	);

	//- Return the dictionary constructed from the old-style data-stream.
	//  Needed for backward compatibility
	static dictionary solverDict(Istream&);


	//- Runtime type information
	TypeName("BICCG");


	// Constructors

		//- Construct from matrix components and solver data stream
		BICCG
		(
			const word& fieldName,
			const lduMatrix& matrix,
			const FieldField<Field, scalar>& coupleBouCoeffs,
			const FieldField<Field, scalar>& coupleIntCoeffs,
			const lduInterfaceFieldPtrsList& interfaces,
			const dictionary& dict
		);


		//- Construct from matrix components and tolerances
		BICCG
		(
			const word& fieldName,
			const lduMatrix& matrix,
			const FieldField<Field, scalar>& coupleBouCoeffs,
			const FieldField<Field, scalar>& coupleIntCoeffs,
			const lduInterfaceFieldPtrsList& interfaces,
			const scalar tolerance,
			const scalar relTol = 0.0
		);

	// Destructor

		virtual ~BICCG()
		{}

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
