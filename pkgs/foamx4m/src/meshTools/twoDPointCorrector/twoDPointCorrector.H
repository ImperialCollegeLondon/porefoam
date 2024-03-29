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
	Foam::twoDPointCorrector

Description
	Class applies a two-dimensional correction to mesh motion point field.

	The correction guarantees that the mesh does not get twisted during motion
	and thus introduce a third dimension into a 2-D problem.

	The operation is performed by looping through all edges approximately
	normal to the plane and enforcing their orthoginality onto the plane by
	adjusting points on their ends.

SourceFiles
	twoDPointCorrector.C

\*---------------------------------------------------------------------------*/

#ifndef twoDPointCorrector_H
#define twoDPointCorrector_H

#include "pointField.H"
#include "labelList.H"
#include "vector.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward class declarations
class polyMesh;


class twoDPointCorrector
{
	// Private data

		//- Reference to moving mesh
		const polyMesh& mesh_;

		//- Is 2D correction required, i.e. is the mesh
		bool required_;

		//- 2-D plane unit normal
		mutable vector* planeNormalPtr_;

		//- Indices of edges normal to plane
		mutable labelList* normalEdgeIndicesPtr_;


	// Private Member Functions

		//- Disallow default bitwise copy construct
		twoDPointCorrector(const twoDPointCorrector&);

		//- Disallow default bitwise assignment
		void operator=(const twoDPointCorrector&);


		//- Calculate addressing
		void calcAddressing() const;

		//- Clear addressing
		void clearAddressing() const;


	// Static data members

		//- Edge orthogonality tolerance
		static const scalar edgeOrthogonalityTol;


public:

	// Constructors

		//- Construct from components
		twoDPointCorrector(const polyMesh& mesh);


	// Destructor

		~twoDPointCorrector();


	// Member Functions

		//- Is 2D correction required, i.e. is the mesh a wedge or slab
		bool required() const
		{
			return required_;
		}

		//- Return plane normal
		const vector& planeNormal() const;

		//- Return indices of normal edges.
		const labelList& normalEdgeIndices() const;

		//- Return direction normal to plane
		direction normalDir() const;

		//- Correct motion points
		void correctPoints(pointField& p) const;

		//- Update topology
		void updateMesh();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
