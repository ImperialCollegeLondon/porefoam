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
	Foam::oppositeFace

Description
	Class containing opposite face for a prismatic cell with addresing
	and a possibility of failure.

SourceFiles
	oppositeFace.C

\*---------------------------------------------------------------------------*/

#ifndef oppositeFace_H
#define oppositeFace_H

#include "face.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


class oppositeFace
:
	public face
{
	// Private data

		//- Master face index
		const label masterIndex_;

		//- Opposite face index
		const label oppositeIndex_;


public:

	// Constructors

		//- Construct from components
		oppositeFace
		(
			const face& f,
			const label masterIndex,
			const label oppositeIndex
		)
		:
			face(f),
			masterIndex_(masterIndex),
			oppositeIndex_(oppositeIndex)
		{}


	// Member Functions

		//- Master face index
		inline label masterIndex() const
		{
			return masterIndex_;
		}

		//- Slave face index
		inline label oppositeIndex() const
		{
			return oppositeIndex_;
		}

		//- Does the opposite face exist?
		inline bool found() const
		{
			return size() > 0;
		}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
