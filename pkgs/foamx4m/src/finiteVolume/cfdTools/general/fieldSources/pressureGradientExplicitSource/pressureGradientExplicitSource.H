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
	Foam::pressureGradientExplicitSource

Description
	Creates a cell set pressure gradient source

	Note: Currently only handles kinematic pressure

SourceFiles
	pressureGradientExplicitSource.C

\*---------------------------------------------------------------------------*/

#ifndef pressureGradientExplicitSource_H
#define pressureGradientExplicitSource_H

#include "autoPtr.H"
#include "topoSetSource.H"
#include "cellSet.H"
#include "fvMesh.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


class pressureGradientExplicitSource
{
	// Private data

		//- Name of the source
		const word sourceName_;

		//- Reference to the mesh
		const fvMesh& mesh_;

		//- Reference to the velocity field
		volVectorField& U_;

		//- Properties dictionary
		IOdictionary dict_;

		//- Average velocity
		vector Ubar_;

		//- Initial pressure gradient
		dimensionedScalar gradPini_;

		//- Pressure gradient
		dimensionedScalar gradP_;

		//- Flow direction
		vector flowDir_;

		//- Name of cell source
		word cellSource_;

		//- The method by which the cells will be selecetd
		autoPtr<topoSetSource> cellSelector_;

		//- The set of selected cells
		cellSet selectedCellSet_;


	// Private Member Functions

		//- Write the pressure gradient to file (for restarts etc)
		void writeGradP() const;

		//- Disallow default bitwise copy construct
		pressureGradientExplicitSource(const pressureGradientExplicitSource&);

		//- Disallow default bitwise assignment
		void operator=(const pressureGradientExplicitSource&);


public:

	// Constructors

		//- Construct from explicit source name and mesh
		pressureGradientExplicitSource
		(
			const word& sourceName,
			volVectorField& U
		);


	// Member Functions

		// Access

			//- Return a tmp field of the source
			tmp<DimensionedField<vector, volMesh> > Su() const;


		//- Correct driving force for a constant mass flow rate
		void update();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
