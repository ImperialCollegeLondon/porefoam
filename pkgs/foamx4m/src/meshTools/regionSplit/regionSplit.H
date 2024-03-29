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
	Foam::regionSplit

Description
	This class separates the mesh into distinct unconnected regions,
	each of which is then given a label.

SourceFiles
	regionSplit.C

\*---------------------------------------------------------------------------*/

#ifndef regionSplit_H
#define regionSplit_H

#include "polyMesh.H"
#include "demandDrivenData.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


class regionSplit
:
	public labelList
{
	// Private data

		//- Reference to mesh
		const polyMesh& mesh_;

		//- Number of regions
		label nRegions_;


	// Private Member Functions

		//- Transfer faceRegion data from one face to the other (or vice versa)
		void transferCoupledFaceRegion
		(
			const label faceI,
			const label otherFaceI,

			labelList& faceRegion,
			dynamicLabelList& newChangedFaces
		) const;

		//- Given a seed cell label, fill cellRegion/faceRegion with markValue
		//  for contiguous region around it
		void fillSeedMask
		(
			const List<labelPair>& explicitConnections,
			labelList& cellRegion,
			labelList& faceRegion,
			const label seedCellID,
			const label markValue
		) const;

		//- Calculate region split. Return number of regions.
		label calcRegionSplit
		(
			const boolList& blockedFace,
			const List<labelPair>& explicitConnections,
			labelList& cellRegion
		) const;


public:

	//- Runtime type information
	ClassName("regionSplit");


	// Constructors

		//- Construct from mesh
		regionSplit(const polyMesh&);

		//- Construct from mesh and whether face is blocked
		//  NOTE: blockedFace has to be consistent across coupled faces!
		regionSplit(const polyMesh&, const boolList& blockedFace);

		//- Construct from mesh and whether face is blocked.  Additional
		//   explicit connections between normal boundary faces.
		//  NOTE: blockedFace has to be consistent across coupled faces!
		regionSplit
		(
			const polyMesh&,
			const boolList& blockedFace,
			const List<labelPair>&
		);

	// Member Functions

		//- Return number of regions
		label nRegions() const
		{
			return nRegions_;
		}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
