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

#include "cloud.H"
#include "foamTime.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(Foam::cloud, 0);

const Foam::word Foam::cloud::prefix("lagrangian");
Foam::word Foam::cloud::defaultName("defaultCloud");

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cloud::cloud(const objectRegistry& obr, const word& cloudName)
:
	objectRegistry
	(
		IOobject
		(
			( cloudName.size() ? cloudName : defaultName ),
			obr.time().timeName(),
			prefix,
			obr,
			IOobject::NO_READ,
			IOobject::AUTO_WRITE
		)
	)
{}


Foam::autoPtr<Foam::cloudDistribute> Foam::cloud::cloudDist
(
	const labelList& cellToProc,
	const labelListList& procCellAddressing,
	const labelListList& procFaceAddressing
)
{
	NotImplemented;
	return autoPtr<cloudDistribute>(nullptr);
}


Foam::labelList Foam::cloud::nParticlesPerCell() const
{
	NotImplemented;
	return labelList(0);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::cloud::~cloud()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::cloud::size() const
{
	NotImplemented;
	return 0;
}


void Foam::cloud::autoMap(const mapPolyMesh&)
{
	NotImplemented;
}


// ************************************************************************* //
