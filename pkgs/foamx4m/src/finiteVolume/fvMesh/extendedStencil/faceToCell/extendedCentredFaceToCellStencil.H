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
	Foam::extendedCentredFaceToCellStencil

Description

SourceFiles
	extendedCentredFaceToCellStencil.C

\*---------------------------------------------------------------------------*/

#ifndef extendedCentredFaceToCellStencil_H
#define extendedCentredFaceToCellStencil_H

#include "extendedFaceToCellStencil.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class faceToCellStencil;


class extendedCentredFaceToCellStencil
:
	public extendedFaceToCellStencil
{
	// Private data

		//- Swap map for getting neigbouring data
		autoPtr<mapDistribute> mapPtr_;

		//- Per face the stencil.
		labelListList stencil_;


	// Private Member Functions

		//- Disallow default bitwise copy construct
		extendedCentredFaceToCellStencil
		(
			const extendedCentredFaceToCellStencil&
		);

		//- Disallow default bitwise assignment
		void operator=(const extendedCentredFaceToCellStencil&);


public:

	// Constructors

		//- Construct from uncompacted face stencil
		explicit extendedCentredFaceToCellStencil(const faceToCellStencil&);


	// Member Functions

		//- Return reference to the parallel distribution map
		const mapDistribute& map() const
		{
			return mapPtr_();
		}

		//- Return reference to the stencil
		const labelListList& stencil() const
		{
			return stencil_;
		}

		//- After removing elements from the stencil adapt the schedule (map).
		void compact();

		//- Use map to get the data into stencil order
		template<class T>
		void collectData
		(
			const GeometricField<T, fvsPatchField, surfaceMesh>& fld,
			List<List<T> >& stencilFld
		) const
		{
			extendedFaceToCellStencil::collectData
			(
				map(),
				stencil(),
				fld,
				stencilFld
			);
		}

		//- Sum surface field contributions to create cell values
		template<class Type>
		tmp<GeometricField<Type, fvPatchField, volMesh> > weightedSum
		(
			const GeometricField<Type, fvsPatchField, surfaceMesh>& fld,
			const List<scalarList >& stencilWeights
		) const
		{
			return extendedFaceToCellStencil::weightedSum
			(
				map(),
				stencil(),
				fld,
				stencilWeights
			);
		}

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
