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

#include "multiMode.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
	defineTypeNameAndDebug(multiMode, 0);
	addToRunTimeSelectionTable(viscoelasticLaw, multiMode, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::multiMode::multiMode
(
	const word& name,
	const volVectorField& U,
	const surfaceScalarField& phi,
	const dictionary& dict
)
:
	viscoelasticLaw(name, U, phi),
	tau_
	(
		IOobject
		(
			"tau" + name,
			U.time().timeName(),
			U.mesh(),
			IOobject::NO_READ,
			IOobject::AUTO_WRITE
		),
		U.mesh(),
		dimensionedSymmTensor
		(
			"zero",
			dimensionSet(1, -1, -2, 0, 0, 0, 0),
			symmTensor::zero
		)
	),
	models_()
{
	PtrList<entry> modelEntries(dict.lookup("models"));
	models_.setSize(modelEntries.size());

	forAll (models_, modelI)
	{
		models_.set
		(
			modelI,
			viscoelasticLaw::New
			(
				modelEntries[modelI].keyword(),
				U,
				phi,
				modelEntries[modelI].dict()
			)
		);
	}
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::fvVectorMatrix>
Foam::multiMode::divTau(volVectorField& U) const
{
	tmp<fvVectorMatrix> divMatrix = models_[0].divTau(U);

	for (label i = 1; i < models_.size(); i++)
	{
		divMatrix() += models_[i].divTau(U);
	}

	return divMatrix;
}


Foam::tmp<Foam::volSymmTensorField> Foam::multiMode::tau() const
{
	tau_ *= dimensionedScalar(0);

	for (label i = 0; i < models_.size(); i++)
	{
		tau_ += models_[i].tau();
	}

	return tau_;
}


void Foam::multiMode::correct()
{
	forAll (models_, i)
	{
		Info<< "Model mode "  << i+1 << endl;
		models_[i].correct();
	}

	tau();
}


// ************************************************************************* //
