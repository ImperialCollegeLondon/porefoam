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

#include "viscosityModel.H"
#include "volFields.H"
#include "fvcGrad.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
	defineTypeNameAndDebug(viscosityModel, 0);
	defineRunTimeSelectionTable(viscosityModel, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::viscosityModel::viscosityModel
(
	const word& name,
	const dictionary& viscosityProperties,
	const volVectorField& U,
	const surfaceScalarField& phi
)
:
	name_(name),
	viscosityProperties_(viscosityProperties),
	U_(U),
	phi_(phi)
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::viscosityModel::strainRate() const
{
	// Bug fix: sqrt(2) inconsistency.  HJ, 8/Dec/2009
	return sqrt(2.0)*mag(symm(fvc::grad(U_)));
}


bool Foam::viscosityModel::read(const dictionary& viscosityProperties)
{
	viscosityProperties_ = viscosityProperties;

	return true;
}


// ************************************************************************* //
