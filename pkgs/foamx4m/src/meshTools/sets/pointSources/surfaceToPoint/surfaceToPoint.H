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
	Foam::surfaceToPoint

Description
	A topoSetSource to selects points based on relation to surface.

	Select based on:
	- distance to surface
	- inside/outside status to surface
	  (Uses normal of nearest surface triangle so requires valid surface
	  topology!)

SourceFiles
	surfaceToPoint.C

\*---------------------------------------------------------------------------*/

#ifndef surfaceToPoint_H
#define surfaceToPoint_H

#include "topoSetSource.H"
#include "Map.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
class triSurfaceSearch;


class surfaceToPoint
:
	public topoSetSource
{

	// Private data

		//- Add usage string
		static addToUsageTable usage_;

		//- Name of surface file
		fileName surfName_;

		//- if > 0 : include points with distance to surface less than nearDist.
		scalar nearDist_;

		//- include all points on opposite sign of normal on nearest surface.
		bool includeInside_;

		//- include all points on this sign of normal on nearest surface.
		bool includeOutside_;


	// Private Member Functions

		//- Depending on surface add to or delete from pointSet.
		void combine(topoSet& set, const bool add) const;

		//- Check settings at construction time.
		void checkSettings() const;

public:

	//- Runtime type information
	TypeName("surfaceToPoint");

	// Constructors

		//- Construct from components
		surfaceToPoint
		(
			const polyMesh& mesh,
			const fileName& surfName,
			const scalar nearDist,
			const bool includeInside,
			const bool includeOutside
		);

		//- Construct from dictionary
		surfaceToPoint
		(
			const polyMesh& mesh,
			const dictionary& dict
		);

		//- Construct from Istream
		surfaceToPoint
		(
			const polyMesh& mesh,
			Istream&
		);


	// Destructor

		virtual ~surfaceToPoint();


	// Member Functions

		virtual void applyToSet
		(
			const topoSetSource::setAction action,
			topoSet&
		) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
