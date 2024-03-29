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

Description
	Class applies a two-dimensional correction to mesh motion point field.
	The correction guarantees that the mesh does not get twisted during motion
	and thus introduce a third dimension into a 2-D problem.

\*---------------------------------------------------------------------------*/

#include "twoDPointCorrector.H"
#include "polyMesh.H"
#include "wedgePolyPatch.H"
#include "emptyPolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const scalar twoDPointCorrector::edgeOrthogonalityTol = 1.0 - 1e-4;

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void twoDPointCorrector::calcAddressing() const
{
	// Find geometry normal
	planeNormalPtr_ = new vector(0, 0, 0);
	vector& pn = *planeNormalPtr_;

	bool isWedge = false;

	// Algorithm:
	// Attempt to find wedge patch and work out the normal from it.
	// If not found, find an empty patch with faces in it and use the
	// normal of the first face.  If neither is found, declare an
	// error.

	// Try and find a wedge patch
	const polyBoundaryMesh& patches = mesh_.boundaryMesh();

	forAll (patches, patchI)
	{
		if (isA<wedgePolyPatch>(patches[patchI]))
		{
			isWedge = true;

			pn = refCast<const wedgePolyPatch>(patches[patchI]).centreNormal();

			if (polyMesh::debug)
			{
				Pout << "Found normal from wedge patch " << patchI;
			}

			break;
		}
	}

	// Try to find an empty patch with faces
	if (!isWedge)
	{
		forAll (patches, patchI)
		{
			if (isA<emptyPolyPatch>(patches[patchI]) && patches[patchI].size())
			{
				pn = patches[patchI].faceAreas()[0];

				if (polyMesh::debug)
				{
					Pout << "Found normal from empty patch " << patchI;
				}

				break;
			}
		}
	}


	if (mag(pn) < VSMALL)
	{
		FatalErrorIn
		(
			"twoDPointCorrector::twoDPointCorrector(const polyMesh& mesh, "
			"const vector& n)"
		)   << "Cannot determine normal vector from patches."
			<< abort(FatalError);
	}
	else
	{
		pn /= mag(pn);
	}

	if (polyMesh::debug)
	{
		Pout << " twoDPointCorrector normal: " << pn << endl;
	}

	// Select edges to be included in check.
	normalEdgeIndicesPtr_ = new labelList(mesh_.nEdges());
	labelList& neIndices = *normalEdgeIndicesPtr_;

	const edgeList& meshEdges = mesh_.edges();

	const pointField& meshPoints = mesh_.points();

	label nNormalEdges = 0;

	forAll (meshEdges, edgeI)
	{
		vector edgeVector =
			meshEdges[edgeI].vec(meshPoints)/
			(meshEdges[edgeI].mag(meshPoints) + VSMALL);

		if (mag(edgeVector & pn) > edgeOrthogonalityTol)
		{
			// this edge is normal to the plane. Add it to the list
			neIndices[nNormalEdges++] = edgeI;
		}
	}

	neIndices.setSize(nNormalEdges);

	// Construction check: number of points in a read 2-D or wedge geometry
	// should be odd and the number of edges normal to the plane should be
	// exactly half the number of points
	if (!isWedge)
	{
		if (meshPoints.size() % 2 != 0)
		{
			WarningIn
			(
				"twoDPointCorrector::twoDPointCorrector("
				"const polyMesh& mesh, const vector& n)"
			)   << "the number of vertices in the geometry "
				<< "is odd - this should not be the case for a 2-D case. "
				<< "Please check the geometry."
				<< endl;
		}

		if (2*nNormalEdges != meshPoints.size())
		{
			WarningIn
			(
				"twoDPointCorrector::twoDPointCorrector("
				"const polyMesh& mesh, const vector& n)"
			)   << "The number of points in the mesh is not equal "
				<< "to twice the number of edges normal to the plane "
				<< "- this may be OK only for wedge geometries.\n"
				<< "    Please check the geometry or adjust "
				<< "the orthogonality tolerance.\n" << endl
				<< "Number of normal edges: " << nNormalEdges
				<< " number of points: " << meshPoints.size()
				<< endl;
		}
	}
}


void twoDPointCorrector::clearAddressing() const
{
	deleteDemandDrivenData(planeNormalPtr_);
	deleteDemandDrivenData(normalEdgeIndicesPtr_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

twoDPointCorrector::twoDPointCorrector(const polyMesh& mesh)
:
	mesh_(mesh),
	required_(mesh_.nGeometricD() == 2),
	planeNormalPtr_(nullptr),
	normalEdgeIndicesPtr_(nullptr)
{}



// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::twoDPointCorrector::~twoDPointCorrector()
{
	clearAddressing();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

direction twoDPointCorrector::normalDir() const
{
	const vector& pn = planeNormal();

	if (mag(pn.x()) >= edgeOrthogonalityTol)
	{
		return vector::X;
	}
	else if (mag(pn.y()) >= edgeOrthogonalityTol)
	{
		return vector::Y;
	}
	else if (mag(pn.z()) >= edgeOrthogonalityTol)
	{
		return vector::Z;
	}
	else
	{
		FatalErrorIn("direction twoDPointCorrector::normalDir() const")
			<< "Plane normal not aligned with the coordinate system" << nl
			<< "    pn = " << pn
			<< abort(FatalError);

		return vector::Z;
	}
}


// Return plane normal
const vector& twoDPointCorrector::planeNormal() const
{
	if (!planeNormalPtr_)
	{
		calcAddressing();
	}

	return *planeNormalPtr_;
}


// Return indices of normal edges.
const labelList& twoDPointCorrector::normalEdgeIndices() const
{
	if (!normalEdgeIndicesPtr_)
	{
		calcAddressing();
	}

	return *normalEdgeIndicesPtr_;
}


void twoDPointCorrector::correctPoints(pointField& p) const
{
	if (!required_) return;

	// Algorithm:
	// Loop through all edges. Calculate the average point position A for
	// the front and the back. Correct the position of point P (in two planes)
	// such that vectors AP and planeNormal are parallel

	// Get reference to edges
	const edgeList&  meshEdges = mesh_.edges();

	const labelList& neIndices = normalEdgeIndices();
	const vector& pn = planeNormal();

	forAll (neIndices, edgeI)
	{
		point& pStart = p[meshEdges[neIndices[edgeI]].start()];

		point& pEnd = p[meshEdges[neIndices[edgeI]].end()];

		// calculate average point position
		const point A = 0.5*(pStart + pEnd);

		// correct point locations
		pStart = A + pn*(pn & (pStart - A));
		pEnd = A + pn*(pn & (pEnd - A));
	}
}


void twoDPointCorrector::updateMesh()
{
	clearAddressing();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
