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
	Foam::slidingInterface

Description
	Sliding interface mesh modifier.  Given two face zones, couple the
	master and slave side using a cutting procedure.

	The coupled faces are collected into the "coupled" zone and can become
	either internal or placed into a master and slave coupled zone.  The
	remaining faces (uncovered master or slave) are placed into the master
	and slave patch.

	The definition of the sliding interface can be either integral or partial.
	Integral interface implies that the slave side completely covers
	the master (i.e. no faces are uncovered); partial interface
	implies that the uncovered part of master/slave face zone should
	become boundary faces.

Author
	Hrvoje Jasak, Nabla Ltd. and Wikki Ltd.  All rights reserved.
	Copyright Hrvoje Jasak

SourceFiles
	slidingInterface.C
	coupleSlidingInterface.C
	decoupleSlidingInterface.C
	slidingInterfaceProjectPoints.C
	slidingInterfaceAttachedAddressing.C
	slidingInterfaceClearCouple.C
	writeSlidingInterfaceVTK.C

\*---------------------------------------------------------------------------*/

#ifndef slidingInterface_H
#define slidingInterface_H

#include "polyMeshModifier.H"
#include "primitiveFacePatch.H"
#include "polyPatchID.H"
#include "ZoneIDs.H"
#include "intersection.H"
#include "Pair.H"
#include "tolerancesSwitch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class objectHit;


class slidingInterface
:
	public polyMeshModifier
{
public:

	// Public enumerations

		//- Type of match
		enum typeOfMatch
		{
			INTEGRAL,
			PARTIAL
		};

		//- Direction names
		static const NamedEnum<typeOfMatch, 2> typeOfMatchNames_;


private:

	// Private data

		//- Master face zone ID
		faceZoneID masterFaceZoneID_;

		//- Slave face zone ID
		faceZoneID slaveFaceZoneID_;

		//- Cut point zone ID
		pointZoneID cutPointZoneID_;

		//- Cut face zone ID
		faceZoneID cutFaceZoneID_;

		//- Master patch ID
		polyPatchID masterPatchID_;

		//- Slave patch ID
		polyPatchID slavePatchID_;

		//- Type of match
		const typeOfMatch matchType_;

		//- Couple-decouple operation.
		//  If the interface is coupled, decouple it and vice versa.
		//  Used in conjuction with automatic mesh motion
		mutable Switch coupleDecouple_;

		//- State of the modifier
		mutable Switch attached_;

		//- Point projection algorithm
		intersection::algorithm projectionAlgo_;

		//- Trigger topological change
		mutable bool trigger_;


		// Private addressing data.

			//- Cut face master face.  Gives the index of face in master patch
			//  the cut face has been created from.  For a slave-only face
			//  this will be -1
			mutable labelList* cutFaceMasterPtr_;

			//- Cut face slave face.  Gives the index of face in slave patch
			//  the cut face has been created from.  For a master-only face
			//  this will be -1
			mutable labelList* cutFaceSlavePtr_;

			//- Master zone faceCells
			mutable labelList* masterFaceCellsPtr_;

			//- Slave zone faceCells
			mutable labelList* slaveFaceCellsPtr_;

			//- Master stick-out faces
			mutable labelList* masterStickOutFacesPtr_;

			//- Slave stick-out faces
			mutable labelList* slaveStickOutFacesPtr_;

			//- Retired point mapping.
			//  For every retired slave side point, gives the label of the
			//  master point that replaces it
			mutable Map<label>* retiredPointMapPtr_;

			//- Cut edge pairs
			//  For cut points created by intersection two edges,
			//  store the master-slave edge pair used in creation
			mutable Map<Pair<edge> >* cutPointEdgePairMapPtr_;

			//- Slave point hit.  The index of master point hit by the
			//  slave point in projection.  For no point hit, set to -1
			mutable labelList* slavePointPointHitsPtr_;

			//- Slave edge hit.  The index of master edge hit by the
			//  slave point in projection.
			//  For point or no edge  hit, set to -1
			mutable labelList* slavePointEdgeHitsPtr_;

			//- Slave face hit.  The index of master face hit by the
			//  slave point in projection.
			mutable List<objectHit>* slavePointFaceHitsPtr_;

			//- Master point edge hit.  The index of slave edge hit by
			//  a master point.  For no hit set to -1
			mutable labelList* masterPointEdgeHitsPtr_;

			//- Projected slave points
			mutable pointField* projectedSlavePointsPtr_;


	// Private Member Functions

		//- Disallow default bitwise copy construct
		slidingInterface(const slidingInterface&);

		//- Disallow default bitwise assignment
		void operator=(const slidingInterface&);

		//- Clear out
		void clearOut() const;


		//- Check validity of construction data
		void checkDefinition() const;

		//- Calculate attached addressing
		void calcAttachedAddressing() const;

		//- Calculate decoupled zone face-cell addressing
		void renumberAttachedAddressing(const mapPolyMesh&) const;

		//- Clear attached addressing
		void clearAttachedAddressing() const;


		// Topological changes

			//- Master faceCells
			const labelList& masterFaceCells() const;

			//- Slave faceCells
			const labelList& slaveFaceCells() const;

			//- Master stick-out faces
			const labelList& masterStickOutFaces() const;

			//- Slave stick-out faces
			const labelList& slaveStickOutFaces() const;

			//- Retired point map
			const Map<label>& retiredPointMap() const;

			//- Cut point edge pair map
			const Map<Pair<edge> >& cutPointEdgePairMap() const;

			 //- Clear addressing
			 void clearAddressing() const;

			//- Project slave points and compare with the current projection.
			//  If the projection has changed, the sliding interface
			//  changes topologically
			bool projectPoints() const;

			//- Couple sliding interface
			void coupleInterface(polyTopoChange& ref) const;

			//- Clear projection
			void clearPointProjection() const;

			//- Clear old couple
			void clearCouple(polyTopoChange& ref) const;

			//- Decouple interface (returns it to decoupled state)
			//  Note: this should not be used in normal operation of the
			//  sliding mesh, but only to return the mesh to its
			//  original state
			void decoupleInterface(polyTopoChange& ref) const;


	// Static data members

		//- Point merge tolerance
		static const debug::tolerancesSwitch pointMergeTol_;

		//- Edge merge tolerance
		static const debug::tolerancesSwitch edgeMergeTol_;

		//- Estimated number of faces an edge goes through
		static const label nFacesPerSlaveEdge_;

		//- Edge-face interaction escape limit
		static const label edgeFaceEscapeLimit_;

		//- Integral match point adjustment tolerance
		static const debug::tolerancesSwitch integralAdjTol_;

		//- Edge intersection master catch fraction
		static const debug::tolerancesSwitch edgeMasterCatchFraction_;

		//- Edge intersection co-planar tolerance
		static const debug::tolerancesSwitch edgeCoPlanarTol_;

		//- Edge end cut-off tolerance
		static const debug::tolerancesSwitch edgeEndCutoffTol_;


public:

	//- Runtime type information
	TypeName("slidingInterface");


	// Constructors

		//- Construct from components
		slidingInterface
		(
			const word& name,
			const label index,
			const polyTopoChanger& mme,
			const word& masterFaceZoneName,
			const word& slaveFaceZoneName,
			const word& cutPointZoneName,
			const word& cutFaceZoneName,
			const word& masterPatchName,
			const word& slavePatchName,
			const typeOfMatch tom,
			const bool coupleDecouple,
			const intersection::algorithm algo
		);

		//- Construct from dictionary
		slidingInterface
		(
			const word& name,
			const dictionary& dict,
			const label index,
			const polyTopoChanger& mme
		);


	// Destructor

		virtual ~slidingInterface();


	// Member Functions

		//- Return master face zone ID
		const faceZoneID& masterFaceZoneID() const;

		//- Return slave face zone ID
		const faceZoneID& slaveFaceZoneID() const;

		//- Return cut point zone ID
		const pointZoneID& cutPointZoneID() const;

		//- Return cut face zone ID
		const faceZoneID& cutFaceZoneID() const;

		//- Return true if attached
		bool attached() const
		{
			return attached_;
		}

		//- Check for topology change
		virtual bool changeTopology() const;

		//- Insert the layer addition/removal instructions
		//  into the topological change
		virtual void setRefinement(polyTopoChange&) const;

		//- Modify motion points to comply with the topological change
		virtual void modifyMotionPoints(pointField& motionPoints) const;

		//- Force recalculation of locally stored data on topological change
		virtual void updateMesh(const mapPolyMesh&);

		//- Return projected points for a slave patch
		const pointField& pointProjection() const;


		//- Write
		virtual void write(Ostream&) const;

		//- Write dictionary
		virtual void writeDict(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
