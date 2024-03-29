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

#include "dynamicInkJetFvMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
	defineTypeNameAndDebug(dynamicInkJetFvMesh, 0);
	addToRunTimeSelectionTable(dynamicFvMesh, dynamicInkJetFvMesh, IOobject);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dynamicInkJetFvMesh::dynamicInkJetFvMesh(const IOobject& io)
:
	dynamicFvMesh(io),
	dynamicMeshCoeffs_
	(
		IOdictionary
		(
			IOobject
			(
				"dynamicMeshDict",
				io.time().constant(),
				*this,
				IOobject::MUST_READ_IF_MODIFIED,
				IOobject::NO_WRITE
			)
		).subDict(typeName + "Coeffs")
	),
	amplitude_(readScalar(dynamicMeshCoeffs_.lookup("amplitude"))),
	frequency_(readScalar(dynamicMeshCoeffs_.lookup("frequency"))),
	refPlaneX_(readScalar(dynamicMeshCoeffs_.lookup("refPlaneX"))),
	stationaryPoints_
	(
		IOobject
		(
			"points",
			io.time().constant(),
			meshSubDir,
			*this,
			IOobject::MUST_READ,
			IOobject::NO_WRITE
		)
	)
{
	Info<< "Performing a dynamic mesh calculation: " << endl
		<< "amplitude: " << amplitude_
		<< " frequency: " << frequency_
		<< " refPlaneX: " << refPlaneX_ << endl;
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::dynamicInkJetFvMesh::~dynamicInkJetFvMesh()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::dynamicInkJetFvMesh::update()
{
	scalar scalingFunction =
		0.5*(::cos(2*mathematicalConstant::pi*frequency_*time().value()) - 1.0);

	Info<< "Mesh scaling. Time = " << time().value() << " scaling: "
		<< scalingFunction << endl;

	pointField newPoints = stationaryPoints_;

	newPoints.replace
	(
		vector::X,
		stationaryPoints_.component(vector::X)*
		(
			1.0
		  + pos
			(
			  - (stationaryPoints_.component(vector::X))
			  - refPlaneX_
			)*amplitude_*scalingFunction
		)
	);

	fvMesh::movePoints(newPoints);

	// Mesh motion only - return false
	return false;
}


// ************************************************************************* //
