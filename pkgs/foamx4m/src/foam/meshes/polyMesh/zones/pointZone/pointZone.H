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
	Foam::pointZone

Description
	A subset of mesh points.
	The labels of points in the zone can be obtained from the addressing()
	list.

	For quick check whether a point belongs to the zone use the lookup
	mechanism in pointZoneMesh, where all the zoned points are registered
	with their zone number.

SourceFiles
	pointZone.C
	newPointZone.C

\*---------------------------------------------------------------------------*/

#ifndef pointZone_H
#define pointZone_H

#include "labelList.H"
#include "typeInfo.H"
#include "dictionary.H"
#include "pointZoneMeshFwd.H"
#include "Map.H"
#include "pointField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

class pointZone;
Ostream& operator<<(Ostream&, const pointZone&);



class pointZone
:
	public labelList
{
	// Private data

		//- Name of zone
		word name_;

		//- Index of zone
		label index_;

		//- Reference to zone list
		const pointZoneMesh& zoneMesh_;


		// Demand-driven private data

			//- Map of point labels in zone for fast location lookup
			mutable Map<label>* pointLookupMapPtr_;


	// Private Member Functions

		//- Disallow default bitwise copy construct
		pointZone(const pointZone&);

		//- Return map of local point indices
		const Map<label>& pointLookupMap() const;

		//- Build map of local point indices
		void calcPointLookupMap() const;


public:

	//- Runtime type information
	TypeName("pointZone");


	// Declare run-time constructor selection tables

		declareRunTimeSelectionTable
		(
			autoPtr,
			pointZone,
			dictionary,
			(
				const word& name,
				const dictionary& dict,
				const label index,
				const pointZoneMesh& zm
			),
			(name, dict, index, zm)
		);


	// Constructors

		//- Construct from components
		pointZone
		(
			const word& name,
			const labelList& addr,
			const label index,
			const pointZoneMesh&
		);

		//- Construct from components, transferring contents
		pointZone
		(
			const word& name,
			const Xfer<labelList>& addr,
			const label index,
			const pointZoneMesh&
		);

		//- Construct from dictionary
		pointZone
		(
			const word& name,
			const dictionary&,
			const label index,
			const pointZoneMesh&
		);

		//- Construct given the original zone and resetting the
		//  point list and zone mesh information
		pointZone
		(
			const pointZone&,
			const labelList& addr,
			const label index,
			const pointZoneMesh&
		);

		//- Construct given the original zone, resetting the
		//  face list and zone mesh information
		pointZone
		(
			const pointZone&,
			const Xfer<labelList>& addr,
			const label index,
			const pointZoneMesh&
		);

		//- Construct and return a clone, resetting the zone mesh
		virtual autoPtr<pointZone> clone(const pointZoneMesh& zm) const
		{
			return autoPtr<pointZone>
			(
				new pointZone(*this, *this, index(), zm)
			);
		}

		//- Construct and return a clone, resetting the point list
		//  and zone mesh
		virtual autoPtr<pointZone> clone
		(
			const pointZoneMesh& zm,
			const label index,
			const labelList& addr
		) const
		{
			return autoPtr<pointZone>
			(
				new pointZone(*this, addr, index, zm)
			);
		}


	// Selectors

		//- Return a pointer to a new point zone
		//  created on freestore from dictionary
		static autoPtr<pointZone> New
		(
			const word& name,
			const dictionary&,
			const label index,
			const pointZoneMesh&
		);


	//- Destructor
	virtual ~pointZone();


	// Member Functions

		//- Return name
		const word& name() const
		{
			return name_;
		}

		//- Map storing the local point index for every global point
		//  index.  Used to find out the index of point in the zone from
		//  the known global point index.  If the point is not in the
		//  zone, returns -1
		label whichPoint(const label globalPointID) const;

		//- Return the index of this zone in zone list
		label index() const
		{
			return index_;
		}

		//- Return zoneMesh reference
		const pointZoneMesh& zoneMesh() const;

		//- Clear addressing
		void clearAddressing();

		//- Check zone definition. Return true if in error.
		bool checkDefinition(const bool report = false) const;

		//- Correct patch after moving points
		virtual void movePoints(const pointField&)
		{}

		//- Correct after topology update.  Note: Zone has been rebuilt
		//  and no mapping is required.  HJ, 1/Mar/2009
		void updateMesh();

		//- Write
		virtual void write(Ostream&) const;

		//- Write dictionary
		virtual void writeDict(Ostream&) const;


	// Member Operators

		//- Assign to zone clearing demand-driven data
		void operator=(const pointZone&);

		//- Assign addressing clearing demand-driven data
		void operator=(const labelList&);


	// Ostream Operator

		friend Ostream& operator<<(Ostream&, const pointZone&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
