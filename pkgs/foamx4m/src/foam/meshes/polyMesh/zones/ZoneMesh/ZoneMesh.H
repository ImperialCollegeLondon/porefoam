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
	Foam::ZoneMesh

Description
	List of mesh zones

SourceFiles
	ZoneMesh.C

\*---------------------------------------------------------------------------*/

#ifndef ZoneMesh_H
#define ZoneMesh_H

#include "List.H"
#include "regIOobject.H"
#include "HashSet.H"
#include "pointFieldsFwd.H"
#include "PackedBoolList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Forward declaration of friend functions and operators

template<class ZoneType, class MeshType> class ZoneMesh;

template<class ZoneType, class MeshType>
Ostream& operator<<(Ostream&, const ZoneMesh<ZoneType, MeshType>&);


template<class ZoneType, class MeshType>
class ZoneMesh
:
	public PtrList<ZoneType>,
	public regIOobject
{
	// Private data

		//- Reference to mesh
		const MeshType& mesh_;

		//- Map of zone labels for given element
		mutable Map<label>* zoneMapPtr_;


	// Private member functions

		//- Disallow construct as copy
		ZoneMesh(const ZoneMesh&);

		//- Disallow assignment
		void operator=(const ZoneMesh<ZoneType, MeshType>&);


		//- Create zone map
		void calcZoneMap() const;

		//- Clear addressing
		void clearAddressing();


public:

	// Constructors

		//- Read constructor given IOobject and a MeshType reference
		ZoneMesh
		(
			const IOobject&,
			const MeshType&
		);

		//- Construct given size
		ZoneMesh
		(
			const IOobject&,
			const MeshType&,
			const label size
		);

		//- Construct from Istream
		ZoneMesh
		(
			const IOobject&,
			const MeshType&,
			Istream&
		);


	//- Destructor
	~ZoneMesh();


	// Member functions

		//- Return the mesh reference
		const MeshType& mesh() const
		{
			return mesh_;
		}

		//- Map of zones containing zone index for all zoned elements
		//  Return -1 if the object is not in the zone
		const Map<label>& zoneMap() const;

		//- Given a global object index, return the zone it is in.  If
		//  object does not belong to any zones, return -1
		label whichZone(const label objectIndex) const;

		//- Return a list of zone types
		wordList types() const;

		//- Return a list of zone names
		wordList names() const;

		//- Return zone indices for all matches
		labelList findIndices(const keyType&) const;

		//- Return zone index for the first match, return -1 if not found
		label findIndex(const keyType&) const;

		//- Find zone index given a name
		label findZoneID(const word& zoneName) const;

		//- Mark cells that match the zone specification
		PackedBoolList findMatching(const keyType&) const;

		//- Clear the zones
		void clear();

		//- Check zone definition. Return true if in error.
		bool checkDefinition(const bool report = false) const;

		//- Correct zone mesh after moving points
		void movePoints(const pointField&);

		//- Correct after topology update
		void updateMesh();

		//- writeData member function required by regIOobject
		bool writeData(Ostream&) const;

	// Member Operators

		//- Return const and non-const reference to ZoneType by index.
		using PtrList<ZoneType>::operator[];

		//- Return const reference to ZoneType by name.
		const ZoneType& operator[](const word&) const;

		//- Return reference to ZoneType by name.
		ZoneType& operator[](const word&);


	// Ostream operator

		friend Ostream& operator<< <ZoneType, MeshType>
		(
			Ostream&,
			const ZoneMesh<ZoneType, MeshType>&
		);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#	include "ZoneMesh.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
