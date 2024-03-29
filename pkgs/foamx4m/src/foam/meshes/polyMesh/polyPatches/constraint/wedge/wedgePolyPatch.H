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
	Foam::wedgePolyPatch

Description
	Wedge front and back plane patch.

SourceFiles
	wedgePolyPatch.C

\*---------------------------------------------------------------------------*/

#ifndef wedgePolyPatch_H
#define wedgePolyPatch_H

#include "polyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


class wedgePolyPatch
:
	public polyPatch
{
	// Private data

		//- Axis of the wedge
		vector axis_;

		//- Centre normal between the wedge boundaries
		vector centreNormal_;

		//- Normal to the patch
		vector patchNormal_;

		//- Face transformation tensor
		tensor faceT_;

		//- Neighbour-cell transformation tensor
		tensor cellT_;

		//- Calculate the above tensors
		void initTransforms();


public:

	//- Runtime type information
	TypeName("wedge");


	// Constructors

		//- Construct from components
		wedgePolyPatch
		(
			const word& name,
			const label size,
			const label start,
			const label index,
			const polyBoundaryMesh& bm
		);

		//- Construct from dictionary
		wedgePolyPatch
		(
			const word& name,
			const dictionary& dict,
			const label index,
			const polyBoundaryMesh& bm
		);

		//- Construct given the original patch and resetting the
		//  face list and boundary mesh information
		wedgePolyPatch
		(
			const wedgePolyPatch& pp,
			const polyBoundaryMesh& bm,
			const label index,
			const label newSize,
			const label newStart
		);

		//- Construct as copy
		wedgePolyPatch(const wedgePolyPatch&);

		//- Construct as copy, resetting the boundary mesh
		wedgePolyPatch(const wedgePolyPatch&, const polyBoundaryMesh&);

		//- Construct and return a clone, resetting the boundary mesh
		virtual autoPtr<polyPatch> clone(const polyBoundaryMesh& bm) const
		{
			return autoPtr<polyPatch>(new wedgePolyPatch(*this, bm));
		}

		//- Construct and return a clone, resetting the face list
		//  and boundary mesh
		virtual autoPtr<polyPatch> clone
		(
			const polyBoundaryMesh& bm,
			const label index,
			const label newSize,
			const label newStart
		) const
		{
			return autoPtr<polyPatch>
			(
				new wedgePolyPatch(*this, bm, index, newSize, newStart)
			);
		}


	// Member functions

		// Access

			//- Return axis of the wedge
			const vector& axis() const
			{
				return axis_;
			}

			//- Return plane normal between the wedge boundaries
			const vector& centreNormal() const
			{
				return centreNormal_;
			}

			//- Return the normal to the patch
			const vector& patchNormal() const
			{
				return patchNormal_;
			}

			//- Return face transformation tensor
			const tensor& faceT() const
			{
				return faceT_;
			}

			//- Return neighbour-cell transformation tensor
			const tensor& cellT() const
			{
				return cellT_;
			}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
