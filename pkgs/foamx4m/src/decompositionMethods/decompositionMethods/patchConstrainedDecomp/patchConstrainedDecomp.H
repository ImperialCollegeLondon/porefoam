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
	Foam::patchConstrainedDecomp

Description
	Base decomposition type is modified by assigning patches to given
	processors.  This is controlled by the patch entry, giving patch names
	and processor index for them.

SourceFiles
	patchConstrainedDecomp.C

\*---------------------------------------------------------------------------*/

#ifndef patchConstrainedDecomp_H
#define patchConstrainedDecomp_H

#include "decompositionMethod.H"
#include "Tuple2.H"

namespace Foam
{


class patchConstrainedDecomp
:
	public decompositionMethod
{
	// Private data types

	typedef List<Tuple2<word, label> > procWordList;


	// Private data

		//- Mesh
		const polyMesh& mesh_;

		//- Decomposition dictionary
		dictionary dict_;

		//- Base decomposition method
		autoPtr<decompositionMethod> baseDecompPtr_;

		//- Patch names to processor association
		procWordList patchConstraints_;


	// Private Member Functions

		//- Disallow default bitwise copy construct
		patchConstrainedDecomp(const patchConstrainedDecomp&);

		//- Disallow default bitwise copy construct and assignment
		void operator=(const patchConstrainedDecomp&);


public:

	//- Runtime type information
	TypeName("patchConstrained");


	// Constructors

		//- Construct given the decomposition dictionary and mesh
		patchConstrainedDecomp
		(
			const dictionary& decompositionDict,
			const polyMesh& mesh
		);


	// Destructor
		virtual ~patchConstrainedDecomp()
		{}


	// Member Functions

		//- Return parallel aware of base method
		virtual bool parallelAware() const
		{
			return baseDecompPtr_->parallelAware();
		}

		//- Decompose cells with weights
		virtual labelList decompose
		(
			const pointField& points,
			const scalarField& pointWeights
		);

		//- Decompose cells with weights with explicitly provided connectivity
		virtual labelList decompose
		(
			const labelListList& globalCellCells,
			const pointField& cc,
			const scalarField& cWeights
		)
		{
			return decompose(cc, cWeights);
		}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
