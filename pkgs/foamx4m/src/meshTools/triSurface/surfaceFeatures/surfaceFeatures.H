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
	Foam::surfaceFeatures

Description
	Holds feature edges/points of surface.

	Feature edges are stored in one list and sorted:
		0 .. externalStart_-1               : region edges
		externalStart_ .. internalStart_-1  : external edges
		internalStart_ .. size-1            : internal edges


	NOTE: angle is included angle, not feature angle and is in degrees.
	The included angle is the smallest angle between two planes. For coplanar
	faces it is 180, for straight angles it is 90. To pick up straight edges
	only use included angle of 91 degrees


SourceFiles
	surfaceFeatures.C

\*---------------------------------------------------------------------------*/

#ifndef surfaceFeatures_H
#define surfaceFeatures_H

#include "pointField.H"
#include "Map.H"
#include "HashSet.H"
#include "pointIndexHit.H"
#include "edgeList.H"
#include "typeInfo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class triSurface;


class surfaceFeatures
{
public:

		enum edgeStatus
		{
			NONE,
			REGION,
			EXTERNAL,
			INTERNAL
		};


private:

	//- label and scalar; used in path walking
	class labelScalar
	{
	public:
		label n_;
		scalar len_;

		labelScalar(const label n, const scalar len)
		:
			n_(n),
			len_(len)
		{}
	};


	// Private data

		//- Reference to surface
		const triSurface& surf_;

		//- Labels of points that are features
		labelList featurePoints_;

		//- Labels of edges that are features
		labelList featureEdges_;

		//- Start of external edges in featureEdges_
		label externalStart_;

		//- Start of internal edges in featureEdges_
		label internalStart_;


	// Private Member Functions

		//- Return nearest point on edge (start..end). Also classify nearest:
		//  index=-1: nearest on mid of edge. index=0:nearest on edge.start()
		//  index=1: nearest on edge.end().
		static pointIndexHit edgeNearest
		(
			const point& start,
			const point& end,
			const point& sample
		);


		//- Construct feature points where more than 2 feature edges meet
		void calcFeatPoints(const List<edgeStatus>&);

		//- Choose next unset feature edge.
		label nextFeatEdge
		(
			const List<edgeStatus>& edgeStat,
			const labelList& featVisited,
			const label unsetVal,
			const label prevEdgeI,
			const label vertI
		) const;

		//- Walk connected feature edges. Marks edges in featVisited.
		labelScalar walkSegment
		(
			const bool mark,
			const List<edgeStatus>& edgeStat,
			const label startEdgeI,
			const label startPointI,
			const label currentFeatI,
			labelList& featVisited
		);

public:

	ClassName("surfaceFeatures");

	// Constructors

		//- Construct from surface
		surfaceFeatures(const triSurface&);

		//- Construct from components
		surfaceFeatures
		(
			const triSurface&,
			const labelList& featurePoints,
			const labelList& featureEdges,
			const label externalStart,
			const label internalStart
		);

		//- Construct from surface, angle and min cumulative length and/or
		//  number of elements
		surfaceFeatures
		(
			const triSurface&,
			const scalar includedAngle,
			const scalar minLen = 0,
			const label minElems = 0
		);

		//- Construct from dictionary
		surfaceFeatures(const triSurface&, const dictionary& dict);

		//- Construct from file
		surfaceFeatures(const triSurface&, const fileName& fName);

		//- Construct as copy
		surfaceFeatures(const surfaceFeatures&);


	// Member Functions

		// Access

			inline const triSurface& surface() const
			{
				return surf_;
			}

			//- Return feature point list
			inline const labelList& featurePoints() const
			{
				return featurePoints_;
			}

			//- Return feature edge list
			inline const labelList& featureEdges() const
			{
				return featureEdges_;
			}

			//- start of external edges
			inline label externalStart() const
			{
				return externalStart_;
			}

			//- start of internal edges
			inline label internalStart() const
			{
				return internalStart_;
			}

			//- Return number of region edges
			inline label nRegionEdges() const
			{
				return externalStart_;
			}

			//- Return number of external edges
			inline label nExternalEdges() const
			{
				return internalStart_ - externalStart_;
			}

			//- Return number of internal edges
			inline label nInternalEdges() const
			{
				return featureEdges_.size() - internalStart_;
			}

			//- Helper function: select a subset of featureEdges_
			labelList selectFeatureEdges
			(
				const bool regionEdges,
				const bool externalEdges,
				const bool internalEdges
			) const;


		// Edit

			//- Find feature edges using provided included angle
			void findFeatures(const scalar includedAngle);

			//- Delete small sets of edges. Edges are stringed up and any
			//  string of length < minLen (or nElems < minElems) is deleted.
			void trimFeatures(const scalar minLen, const label minElems);

			//- From member feature edges to status per edge.
			List<edgeStatus> toStatus() const;

			//- Set from status per edge
			void setFromStatus(const List<edgeStatus>&);


		// Find

			//- Find nearest sample for selected surface points (usually the
			//  set of featurePoints). Return map from index in
			//  samples to surface point. Do not include points that are further
			//  than maxDist away (separate maxDist for every sample)
			Map<label> nearestSamples
			(
				const labelList& selectedPoints,
				const pointField& samples,
				const scalarField& maxDist
			) const;

			//- Find nearest sample for regularly sampled points along the
			//  selected (surface) edges. Return map from sample to edge.
			//  maxDist is distance below which gets snapped.
			//  Edge gets sampled at points sampleDist[sampleI] apart.
			//  (with a maximum of 10 samples per edge)
			Map<label> nearestSamples
			(
				const labelList& selectedEdges,
				const pointField& samples,
				const scalarField& sampleDist,
				const scalarField& maxDist,
				const scalar minSampleDist = 0.1
			) const;

			//- Like nearestSamples but now gets nearest point on
			//  sample-edge instead of nearest sample-point itself.
			//  Return map from sample edge to feature edge.
			Map<pointIndexHit> nearestEdges
			(
				const labelList& selectedEdges,
				const edgeList& sampleEdges,
				const labelList& selectedSampleEdges,
				const pointField& samplePoints,
				const scalarField& sampleDist,
				const scalarField& maxDist,
				const scalar minSampleDist = 0.1
			) const;


			//- Find nearest surface edge (out of selectedEdges) for
			//  each sample point.
			//  Sets:
			//  - edgeLabel : label of surface edge.
			//  - edgePoint : exact position of nearest point on edge.
			//  - edgeEndPoint : -1, 0, 1 depending on whether edgePoint is
			//                  on inside/start/end of edge
			void nearestSurfEdge
			(
				const labelList& selectedEdges,
				const pointField& samples,
				const vector& searchSpan,   // search span
				labelList& edgeLabel,
				labelList& edgeEndPoint,
				pointField& edgePoint
			) const;


			//- Find nearest surface edge (out of selectedEdges) for each
			// sample edge.
			//  Sets:
			//  - edgeLabel         : label of surface edge.
			//  - pointOnEdge       : exact position of nearest point on edge.
			//  - pointOnFeature    : exact position on sample edge.
			void nearestSurfEdge
			(
				const labelList& selectedEdges,
				const edgeList& sampleEdges,
				const labelList& selectedSampleEdges,
				const pointField& samplePoints,
				const vector& searchSpan,   // search span
				labelList& edgeLabel,       // label of surface edge or -1
				pointField& pointOnEdge,    // point on above edge
				pointField& pointOnFeature  // point on sample edge
			) const;


		// Write

			//- Write as dictionary
			void writeDict(Ostream&) const;

			//- Write as dictionary to file
			void write(const fileName& fName) const;

			//- Write to separate OBJ files (region, external, internal edges,
			//  feature points) for visualization
			void writeObj(const fileName& prefix) const;



	// Member Operators

		void operator=(const surfaceFeatures&);


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
