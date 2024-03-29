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
	Foam::reflectionVectors

Description
	Container for reflection vectors (= unit normal of nearest wall)

SourceFiles
	reflectionVectors.C

\*---------------------------------------------------------------------------*/

#ifndef reflectionVectors_H
#define reflectionVectors_H

#include "objectRegistry.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


class reflectionVectors
{

protected:

	// Protected data

		//- Unit normal vector for nearest wall
		volVectorField n_;


private:

	// Private Member Functions

		//- Disallow default bitwise copy construct
		reflectionVectors(const reflectionVectors&);

		//- Disallow default bitwise assignment
		void operator=(const reflectionVectors&);


public:

	// Constructors

		//- Construct from mesh
		reflectionVectors(const fvMesh& mesh);


	// Member Functions

		const volVectorField& n() const
		{
			return n_;
		}


		//- Correct for mesh geom/topo changes
		void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
