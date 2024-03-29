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
	Foam::extendedCentredCellToFaceStencil

Description

SourceFiles
	extendedCentredCellToFaceStencil.C

\*---------------------------------------------------------------------------*/

#ifndef extendedCentredCellToFaceStencil_H
#define extendedCentredCellToFaceStencil_H

#include "extendedCellToFaceStencil.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class cellToFaceStencil;


class extendedCentredCellToFaceStencil
:
	public extendedCellToFaceStencil
{
	// Private data

		//- Swap map for getting neigbouring data
		autoPtr<mapDistribute> mapPtr_;

		//- Per face the stencil.
		labelListList stencil_;


	// Private Member Functions

		//- Disallow default bitwise copy construct
		extendedCentredCellToFaceStencil
		(
			const extendedCentredCellToFaceStencil&
		);

		//- Disallow default bitwise assignment
		void operator=(const extendedCentredCellToFaceStencil&);


public:

	// Constructors

		//- Construct from uncompacted face stencil
		explicit extendedCentredCellToFaceStencil(const cellToFaceStencil&);


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
			const GeometricField<T, fvPatchField, volMesh>& fld,
			List<List<T> >& stencilFld
		) const
		{
			extendedCellToFaceStencil::collectData
			(
				map(),
				stencil(),
				fld,
				stencilFld
			);
		}

		//- Sum vol field contributions to create face values
		template<class Type>
		tmp<GeometricField<Type, fvsPatchField, surfaceMesh> > weightedSum
		(
			const GeometricField<Type, fvPatchField, volMesh>& fld,
			const List<scalarList >& stencilWeights
		) const
		{
			return extendedCellToFaceStencil::weightedSum
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
