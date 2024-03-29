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
	Foam::scotchDecomp

Description
	Scotch domain decomposition

SourceFiles
	scotchDecomp.C

\*---------------------------------------------------------------------------*/

#ifndef scotchDecomp_H
#define scotchDecomp_H

#include "decompositionMethod.H"

namespace Foam
{


class scotchDecomp
:
	public decompositionMethod
{
	// Private data

		//- Mesh reference
		const polyMesh& mesh_;


	// Private Member Functions

		//- Disallow default bitwise copy construct
		scotchDecomp(const scotchDecomp&);

		//- Disallow default bitwise assignment
		void operator=(const scotchDecomp&);


		//- Check and print error message
		static void check(const int, const char*);

		label decompose
		(
			const labelList& adjncy,
			const labelList& xadj,
			const scalarField& cWeights,
			labelList& finalDecomp
		);


public:

	//- Runtime type information
	TypeName("scotch");


	// Constructors

		//- Construct given the decomposition dictionary and mesh
		scotchDecomp
		(
			const dictionary& decompositionDict,
			const polyMesh& mesh
		);


	// Destructor

		virtual ~scotchDecomp()
		{}


	// Member Functions

		//- Return mesh
		const polyMesh& mesh() const
		{
			return mesh_;
		}

		virtual bool parallelAware() const
		{
			// Scotch does not know about proc boundaries
			return false;
		}

		//- Return for every coordinate the wanted processor number. Use the
		//  mesh connectivity (if needed)
		//  Weights get normalised so the minimum value is 1 before truncation
		//  to an integer so the weights should be multiples of the minimum
		//  value. The overall sum of weights might otherwise overflow.
		virtual labelList decompose
		(
			const pointField& points,
			const scalarField& pointWeights
		);

		//- Return for every coordinate the wanted processor number. Gets
		//  passed agglomeration map (from fine to coarse cells) and coarse cell
		//  location. Can be overridden by decomposers that provide this
		//  functionality natively.
		//  See note on weights above.
		virtual labelList decompose
		(
			const labelList& agglom,
			const pointField& regionPoints,
			const scalarField& regionWeights
		);

		//- Same but with uniform weights
		virtual labelList decompose
		(
			const labelList& agglom,
			const pointField& regionPoints
		)
		{
			return decompose
			(
				agglom,
				regionPoints,
				scalarField(regionPoints.size(), 1.0)
			);
		}

		//- Return for every coordinate the wanted processor number. Explicitly
		//  provided mesh connectivity.
		//  The connectivity is equal to mesh.cellCells() except for
		//  - in parallel the cell numbers are global cell numbers (starting
		//	from 0 at processor0 and then incrementing all through the
		//	processors)
		//  - the connections are across coupled patches
		//  See note on weights above.
		virtual labelList decompose
		(
			const labelListList& globalCellCells,
			const pointField& cc,
			const scalarField& cWeights
		);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
