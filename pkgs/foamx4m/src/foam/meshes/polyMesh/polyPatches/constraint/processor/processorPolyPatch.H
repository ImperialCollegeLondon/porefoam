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
	Foam::processorPolyPatch

Description
	Neighbour processor patch.

	Note: morph patch face ordering comes geometric or topological.
	Geometric: no cyclics allowed (assumes faces coincident)
	Topological: needs unmodified faces on both sides to correspond. Also
	needs at least one per connected patch area (so all patch faces can be
	visited from an unmodified face)

SourceFiles
	processorPolyPatch.C
	processorPolyPatchMorph.C

\*---------------------------------------------------------------------------*/

#ifndef processorPolyPatch_H
#define processorPolyPatch_H

#include "coupledPolyPatch.H"
#include "polyBoundaryMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


class processorPolyPatch
:
	public coupledPolyPatch
{
	// Private data

		int myProcNo_;
		int neighbProcNo_;

		//- Processor-neighbbour patch face centres
		vectorField neighbFaceCentres_;

		//- Processor-neighbbour patch face areas
		vectorField neighbFaceAreas_;

		//- Processor-neighbbour patch neighbour cell centres
		vectorField neighbFaceCellCentres_;

		//- Corresponding neighbouring local point label for every local point
		//  (so localPoints()[i] == neighb.localPoints()[neighbPoints_[i]])
		mutable labelList* neighbPointsPtr_;

		//- Corresponding neighbouring local edge label for every local edge
		//  (so edges()[i] == neighb.edges()[neighbEdges_[i]])
		mutable labelList* neighbEdgesPtr_;


	// Private static data

		//- Whether to use geometric or topological matching
		static bool geometricMatch_;


protected:

	// Protected Member functions

		//- Initialise the calculation of the patch addressing
		virtual void initAddressing();

		//- Calculate the patch addressing
		virtual void calcAddressing();

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
	TypeName("processor");


	// Constructors

		//- Construct from components
		processorPolyPatch
		(
			const word& name,
			const label size,
			const label start,
			const label index,
			const polyBoundaryMesh& bm,
			const int myProcNo,
			const int neighbProcNo
		);

		//- Construct from dictionary
		processorPolyPatch
		(
			const word& name,
			const dictionary& dict,
			const label index,
			const polyBoundaryMesh&
		);

		//- Construct as given the original patch and resetting the
		//  face list and boundary mesh information
		processorPolyPatch
		(
			const processorPolyPatch& pp,
			const polyBoundaryMesh& bm,
			const label index,
			const label newSize,
			const label newStart
		);

		//- Construct as copy
		processorPolyPatch(const processorPolyPatch&);

		//- Construct as copy, resetting the boundary mesh
		processorPolyPatch(const processorPolyPatch&, const polyBoundaryMesh&);

		//- Construct and return a clone, resetting the boundary mesh
		virtual autoPtr<polyPatch> clone(const polyBoundaryMesh& bm) const
		{
			return autoPtr<polyPatch>(new processorPolyPatch(*this, bm));
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
				new processorPolyPatch
				(
					refCast<const processorPolyPatch>(*this),
					bm,
					index,
					newSize,
					newStart
				)
			);
		}


	//- Destructor
	virtual ~processorPolyPatch();


	// Member functions

		//- Return processor number
		int myProcNo() const
		{
			return myProcNo_;
		}

		//- Return neigbour processor number
		int neighbProcNo() const
		{
			return neighbProcNo_;
		}

		// Communications support

			//- Return communicator used for communication
			virtual label comm() const;

			//- Return message tag to use for communication
			virtual int tag() const;


		//- Is this the master side?
		virtual bool master() const
		{
			return (myProcNo_ < neighbProcNo_);
		}

		//- Return processor-neighbbour patch face centres
		const vectorField& neighbFaceCentres() const
		{
			return neighbFaceCentres_;
		}

		//- Return processor-neighbbour patch face areas
		const vectorField& neighbFaceAreas() const
		{
			return neighbFaceAreas_;
		}

		//- Return processor-neighbbour patch neighbour cell centres
		const vectorField& neighbFaceCellCentres() const
		{
			return neighbFaceCellCentres_;
		}

		//- Return neighbour point labels. This is for my local point (-1 or)
		//  the corresponding local point on the other side. It is -1 if
		//  there are multiple corresponding points on this or the other side
		//  (can happen for cyclics being converted into proc patches)
		const labelList& neighbPoints() const;

		//- Return neighbour edge labels. This is for my local edge (-1 or) the
		//  corresponding local edge on the other side. See above for -1 cause.
		const labelList& neighbEdges() const;


		//- Initialize ordering for primitivePatch. Does not
		//  refer to *this (except for name() and type() etc.)
		virtual void initOrder(const primitivePatch&) const;

		//- Return new ordering for primitivePatch.
		//  Ordering is -faceMap: for every face
		//  index of the new face -rotation:for every new face the clockwise
		//  shift of the original face. Return false if nothing changes
		//  (faceMap is identity, rotation is 0), true otherwise.
		virtual bool order
		(
			const primitivePatch&,
			labelList& faceMap,
			labelList& rotation
		) const;

		//- Synchronise communications of ordering for primitivePatch
		//  Used in cases when no topological change happens locally,
		//  but is happening on other processors
		virtual void syncOrder() const;


		//- Write the polyPatch data as a dictionary
		virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
