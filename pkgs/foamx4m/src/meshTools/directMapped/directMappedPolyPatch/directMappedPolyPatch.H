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
	Foam::directMappedPolyPatch

Description
	Determines a mapping between patch face centres and mesh cell or face
	centres and processors they're on.

Note
	Storage is not optimal. It stores all face centres and cells on all
	processors to keep the addressing calculation simple.

SourceFiles
	directMappedPolyPatch.C

\*---------------------------------------------------------------------------*/

#ifndef directMappedPolyPatch_H
#define directMappedPolyPatch_H

#include "polyPatch.H"
#include "directMappedPatchBase.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class polyMesh;


class directMappedPolyPatch
:
	public polyPatch,
	public directMappedPatchBase
{

protected:

		//- Initialise the calculation of the patch geometry
		virtual void initGeometry();

		//- Calculate the patch geometry
		virtual void calcGeometry();

		//- Initialise the patches for moving points
		virtual void initMovePoints(const pointField&);

		//- Correct patches after moving points
		virtual void movePoints(const pointField&);

		//- Initialise the update of the patch topology
		virtual void initUpdateMesh();

		//- Update of the patch topology
		virtual void updateMesh();


public:

	//- Runtime type information
	TypeName("directMappedPatch");


	// Constructors

		//- Construct from components
		directMappedPolyPatch
		(
			const word& name,
			const label size,
			const label start,
			const label index,
			const polyBoundaryMesh& bm
		);

		//- Construct from components
		directMappedPolyPatch
		(
			const word& name,
			const label size,
			const label start,
			const label index,
			const word& sampleRegion,
			const directMappedPatchBase::sampleMode mode,
			const word& samplePatch,
			const vectorField& offset,
			const polyBoundaryMesh& bm
		);

		//- Construct from components. Uniform offset.
		directMappedPolyPatch
		(
			const word& name,
			const label size,
			const label start,
			const label index,
			const word& sampleRegion,
			const directMappedPatchBase::sampleMode mode,
			const word& samplePatch,
			const vector& offset,
			const polyBoundaryMesh& bm
		);

		//- Construct from dictionary
		directMappedPolyPatch
		(
			const word& name,
			const dictionary& dict,
			const label index,
			const polyBoundaryMesh& bm
		);

		//- Construct given the original patch and resetting the
		//  face list and boundary mesh information
		directMappedPolyPatch
		(
			const directMappedPolyPatch& pp,
			const polyBoundaryMesh& bm,
			const label index,
			const label newSize,
			const label newStart
		);

		//- Construct as copy
		directMappedPolyPatch(const directMappedPolyPatch&);

		//- Construct as copy, resetting the boundary mesh
		directMappedPolyPatch
		(
			const directMappedPolyPatch&,
			const polyBoundaryMesh&
		);

		//- Construct and return a clone, resetting the boundary mesh
		virtual autoPtr<polyPatch> clone(const polyBoundaryMesh& bm) const
		{
			return autoPtr<polyPatch>(new directMappedPolyPatch(*this, bm));
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
				new directMappedPolyPatch(*this, bm, index, newSize, newStart)
			);
		}


	//- Destructor
	virtual ~directMappedPolyPatch();


	// Member functions

		//- Write the polyPatch data as a dictionary
		virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
