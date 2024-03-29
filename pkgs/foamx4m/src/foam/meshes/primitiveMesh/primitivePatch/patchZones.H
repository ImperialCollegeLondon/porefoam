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
	Foam::patchZones

Description
	Calculates zone number for every face of patch.

	Gets constructed from patch and is a labeList with zone number for
	every patch face.

SourceFiles
	patchZones.C

\*---------------------------------------------------------------------------*/

#ifndef patchZones_H
#define patchZones_H

#include "labelList.H"
#include "pointField.H"
#include "polyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


class patchZones
:
	public labelList
{
	// Private data

		//- Reference to patch
		const primitivePatch& pp_;

		//- Reference to protected edges
		const boolList& borderEdge_;

		//- Max number of zones
		label nZones_;


	// Private Member Functions

		//- Gets labels of changed faces and propagates them to the edges.
		//  Returns labels of edges changed. Does not cross edges marked in
		//  regionEdge
		labelList faceToEdge
		(
			const labelList& changedFaces,
			labelList& edgeRegion
		);

		//- Reverse of faceToEdge: gets edges and returns faces
		labelList edgeToFace(const labelList& changedEdges);

		//- Fill *this with current zone for every face reachable
		//  from faceI without crossing edge marked in borderEdge.
		void markZone(label faceI);

public:

	ClassName("patchZones");


	// Constructors

		//- Fills *this with zone of face. Zone is area
		//  reachable by edge crossing without crossing borderEdge
		//  (bool for every edge in patch).
		patchZones(const primitivePatch& pp, const boolList& borderEdge);


	// Member Functions

		//- Number of zones
		label nZones() const
		{
			return nZones_;
		}

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
