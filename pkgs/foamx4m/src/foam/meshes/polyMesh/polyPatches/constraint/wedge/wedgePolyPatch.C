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

#include "wedgePolyPatch.H"
#include "addToRunTimeSelectionTable.H"
#include "SubField.H"
#include "transform.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
	defineTypeNameAndDebug(wedgePolyPatch, 0);

	addToRunTimeSelectionTable(polyPatch, wedgePolyPatch, word);
	addToRunTimeSelectionTable(polyPatch, wedgePolyPatch, dictionary);
}

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

void Foam::wedgePolyPatch::initTransforms()
{
	if (empty())
	{
		return;
	}

	const pointField& points = this->points();

	patchNormal_ = operator[](0).normal(points);
	patchNormal_ /= mag(patchNormal_);

	centreNormal_ =
		vector
		(
			sign(patchNormal_.x())*(max(mag(patchNormal_.x()), 0.5) - 0.5),
			sign(patchNormal_.y())*(max(mag(patchNormal_.y()), 0.5) - 0.5),
			sign(patchNormal_.z())*(max(mag(patchNormal_.z()), 0.5) - 0.5)
		);
	centreNormal_ /= mag(centreNormal_);

	if
	(
		mag(centreNormal_.x() + centreNormal_.y() + centreNormal_.z())
		< (1 - SMALL)
	)
	{
		FatalErrorIn
		(
			"wedgePolyPatch::wedgePolyPatch(const polyPatch&, "
			"const fvBoundaryMesh&)"
		)   << "wedge " << name()
			<< " centre plane does not align with a coordinate plane by "
			<< 1
			 - mag(centreNormal_.x() + centreNormal_.y() + centreNormal_.z())
			<< exit(FatalError);
	}

	axis_ = centreNormal_ ^ patchNormal_;
	scalar magAxis = mag(axis_);
	axis_ /= magAxis;

	if (magAxis < SMALL)
	{
		FatalErrorIn
		(
			"wedgePolyPatch::initTransforms()"
		)   << "wedge " << name()
			<< " plane aligns with a coordinate plane." << nl
			<< "    The wedge plane should make a small angle (~2.5deg)"
			   " with the coordinate plane" << nl
			<< "    and the the pair of wedge planes should be symmetric"
			<< " about the coordinate plane." << nl
			<< "    Normal of face " << 0 << " is " << patchNormal_
			<< " , implied coordinate plane direction is " << centreNormal_
			<< exit(FatalError);
	}

	faceT_ = rotationTensor(centreNormal_, patchNormal_);
	cellT_ = faceT_ & faceT_;
}


// * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * * * * //

Foam::wedgePolyPatch::wedgePolyPatch
(
	const word& name,
	const label size,
	const label start,
	const label index,
	const polyBoundaryMesh& bm
)
:
	polyPatch(name, size, start, index, bm)
{
	initTransforms();
}


Foam::wedgePolyPatch::wedgePolyPatch
(
	const word& name,
	const dictionary& dict,
	const label index,
	const polyBoundaryMesh& bm
)
:
	polyPatch(name, dict, index, bm)
{
	initTransforms();
}


Foam::wedgePolyPatch::wedgePolyPatch
(
	const wedgePolyPatch& pp,
	const polyBoundaryMesh& bm,
	const label index,
	const label newSize,
	const label newStart
)
:
	polyPatch(pp, bm, index, newSize, newStart)
{
	initTransforms();
}


Foam::wedgePolyPatch::wedgePolyPatch
(
	const wedgePolyPatch& pp
)
:
	polyPatch(pp)
{
	initTransforms();
}


Foam::wedgePolyPatch::wedgePolyPatch
(
	const wedgePolyPatch& pp,
	const polyBoundaryMesh& bm
)
:
	polyPatch(pp, bm)
{
	initTransforms();
}


// ************************************************************************* //
