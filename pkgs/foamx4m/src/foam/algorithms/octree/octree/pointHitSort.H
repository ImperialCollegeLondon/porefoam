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
	Foam::pointHitSort

Description
	Container for sorting intersections

SourceFiles

\*---------------------------------------------------------------------------*/

#ifndef pointHitSort_H
#define pointHitSort_H

#include "pointHit.H"
#include "label.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes


class pointHitSort
{
	// Private data

		//- intersection
		pointHit inter_;

		//- original index
		label index_;

public:

	// Constructors

		//- Construct null
		pointHitSort()
		:
			inter_(false, vector::zero, GREAT),
			index_(-1)
		{}

		//- Construct from intersection, index
		pointHitSort(const pointHit& inter, const label index)
		:
			inter_(inter),
			index_(index)
		{}


	// Member Functions

		const pointHit& inter() const
		{
			return inter_;
		}

		label index() const
		{
			return index_;
		}

	// Member Operators

		bool operator==(const pointHitSort& rhs) const
		{
			return inter_.distance() == rhs.inter().distance();
		}

		bool operator>(const pointHitSort& rhs) const
		{
			return inter_.distance() > rhs.inter().distance();
		}

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
