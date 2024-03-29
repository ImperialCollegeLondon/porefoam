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

\*---------------------------------------------------------------------------*/

#include "processorPolyPatch.H"
#include "addToRunTimeSelectionTable.H"
#include "dictionary.H"
#include "SubField.H"
#include "matchPoints.H"
#include "OFstream.H"
#include "polyBoundaryMesh.H"
#include "polyMesh.H"
#include "foamTime.H"
#include "transformList.H"
#include "demandDrivenData.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
	defineTypeNameAndDebug(processorPolyPatch, 0);
	addToRunTimeSelectionTable(polyPatch, processorPolyPatch, dictionary);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::processorPolyPatch::processorPolyPatch
(
	const word& name,
	const label size,
	const label start,
	const label index,
	const polyBoundaryMesh& bm,
	const int myProcNo,
	const int neighbProcNo
)
:
	coupledPolyPatch(name, size, start, index, bm),
	myProcNo_(myProcNo),
	neighbProcNo_(neighbProcNo),
	neighbFaceCentres_(),
	neighbFaceAreas_(),
	neighbFaceCellCentres_(),
	neighbPointsPtr_(nullptr),
	neighbEdgesPtr_(nullptr)
{}


Foam::processorPolyPatch::processorPolyPatch
(
	const word& name,
	const dictionary& dict,
	const label index,
	const polyBoundaryMesh& bm
)
:
	coupledPolyPatch(name, dict, index, bm),
	myProcNo_(readLabel(dict.lookup("myProcNo"))),
	neighbProcNo_(readLabel(dict.lookup("neighbProcNo"))),
	neighbFaceCentres_(),
	neighbFaceAreas_(),
	neighbFaceCellCentres_(),
	neighbPointsPtr_(nullptr),
	neighbEdgesPtr_(nullptr)
{}


Foam::processorPolyPatch::processorPolyPatch
(
	const processorPolyPatch& pp,
	const polyBoundaryMesh& bm,
	const label index,
	const label newSize,
	const label newStart
)
:
	coupledPolyPatch(pp, bm, index, newSize, newStart),
	myProcNo_(pp.myProcNo_),
	neighbProcNo_(pp.neighbProcNo_),
	neighbFaceCentres_(),
	neighbFaceAreas_(),
	neighbFaceCellCentres_(),
	neighbPointsPtr_(nullptr),
	neighbEdgesPtr_(nullptr)
{}


Foam::processorPolyPatch::processorPolyPatch
(
	const processorPolyPatch& pp
)
:
	coupledPolyPatch(pp),
	myProcNo_(pp.myProcNo_),
	neighbProcNo_(pp.neighbProcNo_),
	neighbFaceCentres_(),
	neighbFaceAreas_(),
	neighbFaceCellCentres_(),
	neighbPointsPtr_(nullptr),
	neighbEdgesPtr_(nullptr)
{}


Foam::processorPolyPatch::processorPolyPatch
(
	const processorPolyPatch& pp,
	const polyBoundaryMesh& bm
)
:
	coupledPolyPatch(pp, bm),
	myProcNo_(pp.myProcNo_),
	neighbProcNo_(pp.neighbProcNo_),
	neighbFaceCentres_(),
	neighbFaceAreas_(),
	neighbFaceCellCentres_(),
	neighbPointsPtr_(nullptr),
	neighbEdgesPtr_(nullptr)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::processorPolyPatch::~processorPolyPatch()
{
	deleteDemandDrivenData(neighbPointsPtr_);
	deleteDemandDrivenData(neighbEdgesPtr_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::processorPolyPatch::comm() const
{
	return boundaryMesh().mesh().comm();
}


int Foam::processorPolyPatch::tag() const
{
	return Pstream::msgType();
}


void Foam::processorPolyPatch::initAddressing()
{
	polyPatch::initAddressing();
}


void Foam::processorPolyPatch::calcAddressing()
{
	polyPatch::calcAddressing();
}


void Foam::processorPolyPatch::initGeometry()
{
	if (Pstream::parRun())
	{
		OPstream toNeighbProc
		(
			Pstream::blocking,
			neighbProcNo(),
			3*(sizeof(label) + size()*sizeof(vector) + sizeof(scalar))
		);

		toNeighbProc
			<< faceCentres()
			<< faceAreas()
			<< faceCellCentres();
	}
}


void Foam::processorPolyPatch::calcGeometry()
{
	if (Pstream::parRun())
	{
		{
			IPstream fromNeighbProc
			(
				Pstream::blocking,
				neighbProcNo(),
				3*(sizeof(label) + size()*sizeof(vector) + sizeof(scalar))
			);
			fromNeighbProc
				>> neighbFaceCentres_
				>> neighbFaceAreas_
				>> neighbFaceCellCentres_;
		}

		// My normals
		vectorField faceNormals(size());

		// Neighbour normals
		vectorField nbrFaceNormals(neighbFaceAreas_.size());

		// Calculate normals from areas and check

		// Cache face areas
		const vectorField::subField localFaceAreas = faceAreas();

		forAll (faceNormals, facei)
		{
			scalar magSf = mag(localFaceAreas[facei]);
			scalar nbrMagSf = mag(neighbFaceAreas_[facei]);
			scalar avSf = (magSf + nbrMagSf)/2.0;

			if (magSf < SMALL && nbrMagSf < SMALL)
			{
				// Undetermined normal. Use dummy normal to force separation
				// check.
				// (note use of sqrt(VSMALL) since that is how mag scales)
				// HR, 11/12/2013: Found a face with area = 1e-21 before
				// topo deflation. Hence must use SMALL here.
				faceNormals[facei] = point(1, 0, 0);
				nbrFaceNormals[facei] = faceNormals[facei];
			}
			else if (mag(magSf - nbrMagSf)/avSf > polyPatch::matchTol_)
			{
				FatalErrorIn
				(
					"processorPolyPatch::calcGeometry()"
				)   << "face " << facei << " area does not match neighbour by "
					<< 100*mag(magSf - nbrMagSf)/avSf
					<< "% -- possible face ordering problem." << endl
					<< "patch: " << name()
					<< " my area: " << magSf
					<< " neighbour area: " << nbrMagSf
					<< " matching tolerance: " << polyPatch::matchTol_
					<< endl
					<< "Mesh face: " << start()+facei
					<< " vertices: "
					<< UIndirectList<point>(points(), operator[](facei))()
					<< endl
					<< "Rerun with processor debug flag set for"
					<< " more information." << exit(FatalError);
			}
			else
			{
				faceNormals[facei] = localFaceAreas[facei]/magSf;
				nbrFaceNormals[facei] = neighbFaceAreas_[facei]/nbrMagSf;
			}
		}

		calcTransformTensors
		(
			faceCentres(),
			neighbFaceCentres_,
			faceNormals,
			nbrFaceNormals,
			calcFaceTol(*this, points(), faceCentres())
		);
	}
}


void Foam::processorPolyPatch::initMovePoints(const pointField& p)
{
	polyPatch::movePoints(p);
	processorPolyPatch::initGeometry();
}


void Foam::processorPolyPatch::movePoints(const pointField&)
{
	processorPolyPatch::calcGeometry();
}


void Foam::processorPolyPatch::initUpdateMesh()
{
	polyPatch::initUpdateMesh();

	deleteDemandDrivenData(neighbPointsPtr_);
	deleteDemandDrivenData(neighbEdgesPtr_);

	if (Pstream::parRun())
	{
		// Express all points as patch face and index in face.
		labelList pointFace(nPoints());
		labelList pointIndex(nPoints());

		for (label patchPointI = 0; patchPointI < nPoints(); patchPointI++)
		{
			label faceI = pointFaces()[patchPointI][0];

			pointFace[patchPointI] = faceI;

			const face& f = localFaces()[faceI];

			pointIndex[patchPointI] = findIndex(f, patchPointI);
		}

		// Express all edges as patch face and index in face.
		labelList edgeFace(nEdges());
		labelList edgeIndex(nEdges());

		for (label patchEdgeI = 0; patchEdgeI < nEdges(); patchEdgeI++)
		{
			label faceI = edgeFaces()[patchEdgeI][0];

			edgeFace[patchEdgeI] = faceI;

			const labelList& fEdges = faceEdges()[faceI];

			edgeIndex[patchEdgeI] = findIndex(fEdges, patchEdgeI);
		}

		OPstream toNeighbProc
		(
			Pstream::blocking,
			neighbProcNo(),
			8*sizeof(label)             // four headers of labelList
		  + 2*nPoints()*sizeof(label)   // two point-based labelLists
		  + 2*nEdges()*sizeof(label)    // two edge-based labelLists
		);

		toNeighbProc
			<< pointFace
			<< pointIndex
			<< edgeFace
			<< edgeIndex;
	}
}


void Foam::processorPolyPatch::updateMesh()
{
	// For completeness
	polyPatch::updateMesh();

	if (Pstream::parRun())
	{
		labelList nbrPointFace;
		labelList nbrPointIndex;
		labelList nbrEdgeFace;
		labelList nbrEdgeIndex;

		{
			// Note cannot predict exact size since opposite nPoints might
			// be different from one over here.

			// Disagree: This needs to be sorted out, so that processor patches
			// are simply internal faces and treat cyclics separately
			// HJ, 10/Mar/2011
			IPstream fromNeighbProc(Pstream::blocking, neighbProcNo());

			fromNeighbProc
				>> nbrPointFace
				>> nbrPointIndex
				>> nbrEdgeFace
				>> nbrEdgeIndex;
		}

		// Convert neighbour faces and indices into face back into
		// my edges and points.

		// Convert points.
		// ~~~~~~~~~~~~~~~

		neighbPointsPtr_ = new labelList(nPoints(), -1);
		labelList& neighbPoints = *neighbPointsPtr_;

		forAll (nbrPointFace, nbrPointI)
		{
			// Find face and index in face on this side.
			const face& f = localFaces()[nbrPointFace[nbrPointI]];
			label index = (f.size() - nbrPointIndex[nbrPointI]) % f.size();
			label patchPointI = f[index];

			if (neighbPoints[patchPointI] == -1)
			{
				// First reference of point
				neighbPoints[patchPointI] = nbrPointI;
			}
			else if (neighbPoints[patchPointI] >= 0)
			{
				// Point already visited. Mark as duplicate.
				neighbPoints[patchPointI] = -2;
			}
		}

		// Reset all duplicate entries to -1.
		forAll (neighbPoints, patchPointI)
		{
			if (neighbPoints[patchPointI] == -2)
			{
				neighbPoints[patchPointI] = -1;
			}
		}

		// Convert edges.
		// ~~~~~~~~~~~~~~

		neighbEdgesPtr_ = new labelList(nEdges(), -1);
		labelList& neighbEdges = *neighbEdgesPtr_;

		forAll (nbrEdgeFace, nbrEdgeI)
		{
			// Find face and index in face on this side.
			const labelList& f = faceEdges()[nbrEdgeFace[nbrEdgeI]];
			label index = (f.size() - nbrEdgeIndex[nbrEdgeI] - 1) % f.size();
			label patchEdgeI = f[index];

			if (neighbEdges[patchEdgeI] == -1)
			{
				// First reference of edge
				neighbEdges[patchEdgeI] = nbrEdgeI;
			}
			else if (neighbEdges[patchEdgeI] >= 0)
			{
				// Edge already visited. Mark as duplicate.
				neighbEdges[patchEdgeI] = -2;
			}
		}

		// Reset all duplicate entries to -1.
		forAll (neighbEdges, patchEdgeI)
		{
			if (neighbEdges[patchEdgeI] == -2)
			{
				neighbEdges[patchEdgeI] = -1;
			}
		}

		// Remove any addressing used for shared points/edges calculation
		primitivePatch::clearOut();
	}
}


const Foam::labelList& Foam::processorPolyPatch::neighbPoints() const
{
	if (!neighbPointsPtr_)
	{
		FatalErrorIn("processorPolyPatch::neighbPoints() const")
			<< "No extended addressing calculated for patch " << name()
			<< abort(FatalError);
	}

	return *neighbPointsPtr_;
}


const Foam::labelList& Foam::processorPolyPatch::neighbEdges() const
{
	if (!neighbEdgesPtr_)
	{
		FatalErrorIn("processorPolyPatch::neighbEdges() const")
			<< "No extended addressing calculated for patch " << name()
			<< abort(FatalError);
	}

	return *neighbEdgesPtr_;
}


void Foam::processorPolyPatch::initOrder(const primitivePatch& pp) const
{
	if (!Pstream::parRun())
	{
		return;
	}

	// Master side sends the data and slave side rotates the faces
	// HJ, 10/Mar/2011
	if (master())
	{
		pointField ctrs(calcFaceCentres(pp, pp.points()));

		// Anchor point is the start point of all faces.  HJ, 10/Mar/2011
		pointField anchors(getAnchorPoints(pp, pp.points()));

		// Now send all info over to the neighbour
		OPstream toNeighbour(Pstream::blocking, neighbProcNo());
		toNeighbour << ctrs << anchors;
	}
}


// Return new ordering. Ordering is -faceMap: for every face index
// the new face -rotation:for every new face the clockwise shift
// of the original face. Return false if nothing changes (faceMap
// is identity, rotation is 0)
bool Foam::processorPolyPatch::order
(
	const primitivePatch& pp,
	labelList& faceMap,
	labelList& rotation
) const
{
	if (!Pstream::parRun())
	{
		return false;
	}

	// Moved from initOrder to simplify communications
	// HJ, 27/Nov/2009
	if (debug)
	{
		fileName nm
		(
			boundaryMesh().mesh().time().path()
		   /name() + "_faces.obj"
		);
		Pout<< "processorPolyPatch::order : Writing my " << pp.size()
			<< " faces to OBJ file " << nm << endl;
		writeOBJ(nm, pp, pp.points());

		// Calculate my face centres
		pointField ctrs(calcFaceCentres(pp, pp.points()));

		OFstream localStr
		(
			boundaryMesh().mesh().time().path()
		   /name() + "_localFaceCentres.obj"
		);
		Pout<< "processorPolyPatch::order : "
			<< "Dumping " << ctrs.size()
			<< " local faceCentres to " << localStr.name() << endl;

		forAll (ctrs, faceI)
		{
			writeOBJ(localStr, ctrs[faceI]);
		}
	}

	faceMap.setSize(pp.size());
	faceMap = -1;

	rotation.setSize(pp.size());
	rotation = 0;

	if (master())
	{
		// Do nothing (i.e. identical mapping, zero rotation).
		// See comment at top.
		forAll (faceMap, patchFaceI)
		{
			faceMap[patchFaceI] = patchFaceI;
		}

		return false;
	}
	else
	{
		// Slave side: receive points
		vectorField masterCtrs;
		vectorField masterAnchors;

		// Receive data from neighbour
		{
			IPstream fromNeighbour(Pstream::blocking, neighbProcNo());
			fromNeighbour >> masterCtrs >> masterAnchors;
		}

		// Calculate my face centres
		pointField ctrs(calcFaceCentres(pp, pp.points()));

		// Calculate typical distance from face centre
		scalarField tols(calcFaceTol(pp, pp.points(), ctrs));

		if (debug || masterCtrs.size() != pp.size())
		{
			{
				OFstream nbrStr
				(
					boundaryMesh().mesh().time().path()
				   /name() + "_nbrFaceCentres.obj"
				);

				Pout<< "processorPolyPatch::order : "
					<< "Dumping neighbour faceCentres to " << nbrStr.name()
					<< endl;

				forAll (masterCtrs, faceI)
				{
					writeOBJ(nbrStr, masterCtrs[faceI]);
				}
			}

			if (masterCtrs.size() != pp.size())
			{
				FatalErrorIn
				(
					"processorPolyPatch::order(const primitivePatch&"
					", labelList&, labelList&) const"
				)   << "in patch:" << name() << " : "
					<< "Local size of patch is " << pp.size() << " (faces)."
					<< endl
					<< "Received from neighbour " << masterCtrs.size()
					<< " faceCentres!"
					<< abort(FatalError);
			}
		}

		// Geometric match of face centre vectors
		// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

		// 1. Try existing ordering and transformation
		bool matchedAll = false;

		if
		(
			separated()
		 && (separation().size() == 1 || separation().size() == pp.size())
		)
		{
			vectorField transformedCtrs;

			const vectorField& v = separation();

			if (v.size() == 1)
			{
				transformedCtrs = masterCtrs - v[0];
			}
			else
			{
				transformedCtrs = masterCtrs - v;
			}

			matchedAll = matchPoints
			(
				ctrs,
				transformedCtrs,
				tols,
				true,
				faceMap
			);

			if (matchedAll)
			{
				// Use transformed centers from now on
				masterCtrs = transformedCtrs;

				// Transform anchors
				if (v.size() == 1)
				{
					masterAnchors -= v[0];
				}
				else
				{
					masterAnchors -= v;
				}
			}
		}
		else if
		(
		   !parallel()
		 && (forwardT().size() == 1 || forwardT().size() == pp.size())
		)
		{
			vectorField transformedCtrs = masterCtrs;
			transformList(forwardT(), transformedCtrs);

			matchedAll = matchPoints
			(
				ctrs,
				transformedCtrs,
				tols,
				true,
				faceMap
			);

			if (matchedAll)
			{
				// Use transformed centers from now on
				masterCtrs = transformedCtrs;

				// Transform anchors
				transformList(forwardT(), masterAnchors);
			}
		}


		// 2. Try zero separation automatic matching
		if (!matchedAll)
		{
			matchedAll = matchPoints(ctrs, masterCtrs, tols, true, faceMap);
		}

		if (!matchedAll || debug)
		{
			// Dump faces
			fileName str
			(
				boundaryMesh().mesh().time().path()
			   /name()/name() + "_faces.obj"
			);
			Pout<< "processorPolyPatch::order :"
				<< " Writing faces to OBJ file " << str.name() << endl;
			writeOBJ(str, pp, pp.points());

			OFstream ccStr
			(
				boundaryMesh().mesh().time().path()
			   /name() + "_faceCentresConnections.obj"
			);

			Pout<< "processorPolyPatch::order :"
				<< " Dumping newly found match as lines between"
				<< " corresponding face centres to OBJ file " << ccStr.name()
				<< endl;

			label vertI = 0;

			forAll (ctrs, faceI)
			{
				label masterFaceI = faceMap[faceI];

				if (masterFaceI != -1)
				{
					const point& c0 = masterCtrs[masterFaceI];
					const point& c1 = ctrs[faceI];
					writeOBJ(ccStr, c0, c1, vertI);
				}
			}
		}

		if (!matchedAll)
		{
			FatalErrorIn
//             SeriousErrorIn
			(
				"processorPolyPatch::order(const primitivePatch&"
				", labelList&, labelList&) const"
			)   << "in patch:" << name() << " : "
				<< "Cannot match vectors to faces on both sides of patch"
				<< endl
				<< "    masterCtrs[0]:" << masterCtrs[0] << endl
				<< "    ctrs[0]:" << ctrs[0] << endl
				<< "    Please check your topology changes or maybe you have"
				<< " multiple separated (from cyclics) processor patches"
				<< endl
				<< "    Continuing with incorrect face ordering from now on!"
				<< abort(FatalError);

			return false;
		}

		// Set rotation
		forAll (faceMap, oldFaceI)
		{
			// The face f will be at newFaceI (after morphing) and we want its
			// anchorPoint (= f[0]) to align with the anchorpoint for the
			// corresponding face on the other side.

			label newFaceI = faceMap[oldFaceI];

			const point& wantedAnchor = masterAnchors[newFaceI];

			rotation[newFaceI] = getRotation
			(
				pp.points(),
				pp[oldFaceI],
				wantedAnchor,
				tols[oldFaceI]
			);

			if (rotation[newFaceI] == -1)
			{
				FatalErrorIn
//                 SeriousErrorIn
				(
					"processorPolyPatch::order(const primitivePatch&"
					", labelList&, labelList&) const"
				)   << "in patch " << name()
					<< " : "
					<< "Cannot find point on face " << pp[oldFaceI]
					<< " with vertices "
					<< UIndirectList<point>(pp.points(), pp[oldFaceI])()
					<< " that matches point " << wantedAnchor
					<< " when matching the halves of processor patch "
					<< name()
					<< "Continuing with incorrect face ordering from now on!"
					<< abort(FatalError);

				return false;
			}
		}

		forAll (faceMap, faceI)
		{
			if (faceMap[faceI] != faceI || rotation[faceI] != 0)
			{
				return true;
			}
		}

		return false;
	}
}


void Foam::processorPolyPatch::syncOrder() const
{
	if (!Pstream::parRun())
	{
		return;
	}

	// Read and discard info from owner side from the neighbour
	if (slave())
	{
		vectorField masterCtrs;
		vectorField masterAnchors;

		// Receive data from neighbour
		{
			IPstream fromNeighbour(Pstream::blocking, neighbProcNo());
			fromNeighbour >> masterCtrs >> masterAnchors;
		}
	}
}


void Foam::processorPolyPatch::write(Ostream& os) const
{
	polyPatch::write(os);
	os.writeKeyword("myProcNo") << myProcNo_
		<< token::END_STATEMENT << nl;
	os.writeKeyword("neighbProcNo") << neighbProcNo_
		<< token::END_STATEMENT << nl;
}


// ************************************************************************* //
