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
	Foam::cellFeatures

Description
	Cell analysis class.

	Constructs feature edges and feature points, which are edges/points with
	and angle > given specification.
	Can be asked for 'superFaces' which can be used to see if a cell is a
	'splitHex'.

SourceFiles
	cellFeatures.C

\*---------------------------------------------------------------------------*/

#ifndef cellFeatures_H
#define cellFeatures_H

#include "faceList.H"
#include "labelList.H"
#include "boolList.H"
#include "HashSet.H"
#include "Map.H"
#include "DynamicList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class primitiveMesh;


class cellFeatures
{
	// Private data

		const primitiveMesh& mesh_;

		//- cos of angle between two connected faces or two connected edges on
		//  same face before edge/point is 'feature'.
		scalar minCos_;

		label cellI_;

		//- Feature edges
		labelHashSet featureEdge_;

		//- (demand driven) Faces after removing internal points&edges
		mutable faceList* facesPtr_;

		//- New to old face mapping
		mutable List<dynamicLabelList > faceMap_;


	// Private Member Functions

		bool faceAlignedEdge(const label, const label) const;

		label nextEdge
		(
			const Map<label>& toSuperFace,
			const label superFaceI,
			const label thisEdgeI,
			const label thisVertI
		) const;

		bool isCellFeatureEdge(const scalar, const label) const;

		void walkSuperFace
		(
			const label faceI,
			const label superFaceI,
			Map<label>& toSuperFace
		) const;

		void calcSuperFaces() const;


		//- Disallow default bitwise copy construct
		cellFeatures(const cellFeatures&);

		//- Disallow default bitwise assignment
		void operator=(const cellFeatures&);

public:

	// Constructors

		//- Construct from cell in mesh
		cellFeatures
		(
			const primitiveMesh&,
			const scalar minCos,    // angle to use for feature recognition.
			const label cellI
		);


	// Destructor

		~cellFeatures();


	// Member Functions

		// Access

			const labelHashSet& featureEdge() const
			{
				return featureEdge_;
			}

			const faceList& faces() const
			{
				if (!facesPtr_)
				{
					calcSuperFaces();
				}
				return *facesPtr_;
			}

			//- New to old faceMap. Guaranteed to be shrunk.
			const List<dynamicLabelList >& faceMap() const
			{
				if (!facesPtr_)
				{
					calcSuperFaces();
				}
				return faceMap_;
			}


		// Check

			//- Is edge a feature edge (uniquely determined since on cell
			//  only two faces sharing edge)
			bool isFeatureEdge(const label edgeI) const
			{
				return featureEdge().found(edgeI);
			}

			//- Are two edges connected at feature point?
			//  Is local to face since point might be seen as feature point
			//  from one face but not from another.
			bool isFeaturePoint(const label edge0, const label edge1) const;

			//- Is vertexI on faceI used by two edges that form feature
			//  point
			bool isFeatureVertex(const label faceI, const label vertI) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
