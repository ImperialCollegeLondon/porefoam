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
	Foam::perfectInterface

Description
	Hack of attachDetach to couple patches when they perfectly align.
	Does not decouple. Used by stitchMesh app. Does geometric matching.

SourceFiles
	perfectInterface.C

\*---------------------------------------------------------------------------*/

#ifndef perfectInterface_H
#define perfectInterface_H

#include "polyMeshModifier.H"
#include "polyPatchID.H"
#include "ZoneIDs.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class Time;


class perfectInterface
:
	public polyMeshModifier
{
	// Private data

		//- Master face zone ID
		faceZoneID faceZoneID_;

		//- Master patch ID
		polyPatchID masterPatchID_;

		//- Slave patch ID
		polyPatchID slavePatchID_;

		//- Tolerance used for distance comparison (fraction of minimum edge
		//  length)
		static const scalar tol_;

	// Private Member Functions

		//- Calculate face centres on patch
		static pointField calcFaceCentres(const primitivePatch&);


		//- Disallow default bitwise copy construct
		perfectInterface(const perfectInterface&);

		//- Disallow default bitwise assignment
		void operator=(const perfectInterface&);


public:

	//- Runtime type information
	TypeName("perfectInterface");


	// Constructors

		//- Construct from components
		perfectInterface
		(
			const word& name,
			const label index,
			const polyTopoChanger& mme,
			const word& faceZoneName,
			const word& masterPatchName,
			const word& slavePatchName
		);

		//- Construct from dictionary
		perfectInterface
		(
			const word& name,
			const dictionary& dict,
			const label index,
			const polyTopoChanger& mme
		);


	// Destructor

		virtual ~perfectInterface();


	// Member Functions

		//- Check for topology change
		virtual bool changeTopology() const;

		//- Insert the layer addition/removal instructions
		//  into the topological change
		virtual void setRefinement(polyTopoChange&) const;

		//- Modify motion points to comply with the topological change
		virtual void modifyMotionPoints(pointField& motionPoints) const;

		//- Force recalculation of locally stored data on topological change
		virtual void updateMesh(const mapPolyMesh&);

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
