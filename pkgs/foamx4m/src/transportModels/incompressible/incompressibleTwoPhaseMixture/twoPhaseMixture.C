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

#include "twoPhaseMixture.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"
#include "fvc.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

//- Calculate and return the laminar viscosity
word twoPhaseMixture::getPhaseName(const word& key) const
{
	if (isDict(key))
	{
		return key;
	}
	else
	{
		return word(lookup(key));
	}
}

void twoPhaseMixture::calcNu()
{
	nuModel1_->correct();
	nuModel2_->correct();

	volScalarField limitedAlpha1
	(
		"limitedAlpha1",
		min(max(alpha1_, scalar(0)), scalar(1))
	);

	// Average kinematic viscosity calculated from dynamic viscosity
	nu_ = mu()/(limitedAlpha1*rho1_ + (scalar(1) - limitedAlpha1)*rho2_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

twoPhaseMixture::twoPhaseMixture
(
	const volVectorField& U,
	const surfaceScalarField& phi,
	const word& alpha1Name
)
:
	transportModel(U, phi),

	phase1Name_(getPhaseName("phase1")),
	phase2Name_(getPhaseName("phase2")),

	nuModel1_
	(
		viscosityModel::New
		(
			"nu1",
			subDict(phase1Name_),
			U,
			phi
		)
	),
	nuModel2_
	(
		viscosityModel::New
		(
			"nu2",
			subDict(phase2Name_),
			U,
			phi
		)
	),

	rho1_(nuModel1_->viscosityProperties().lookup("rho")),
	rho2_(nuModel2_->viscosityProperties().lookup("rho")),

	U_(U),
	phi_(phi),

	alpha1_(U_.db().lookupObject<const volScalarField> (alpha1Name)),

	nu_
	(
		IOobject
		(
			"nu",
			U_.time().timeName(),
			U_.db()
		),
		U_.mesh(),
		dimensionedScalar("nu", dimensionSet(0, 2, -1, 0, 0), 0),
		calculatedFvPatchScalarField::typeName
	)
{
	calcNu();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

tmp<volScalarField> twoPhaseMixture::rho() const
{
	volScalarField limitedAlpha1 = min(max(alpha1_, scalar(0)), scalar(1));

	return tmp<volScalarField>
	(
		new volScalarField
		(
			"rho_twoPhaseMixture",
			limitedAlpha1*rho1_
		  + (scalar(1) - limitedAlpha1)*rho2_
		)
	);
}


tmp<volScalarField> twoPhaseMixture::mu() const
{
	volScalarField limitedAlpha1 = min(max(alpha1_, scalar(0)), scalar(1));

	return tmp<volScalarField>
	(
		new volScalarField
		(
			"mu_twoPhaseMixture",
			limitedAlpha1*rho1_*nuModel1_->nu()
		  + (scalar(1) - limitedAlpha1)*rho2_*nuModel2_->nu()
		)
	);
}


tmp<surfaceScalarField> twoPhaseMixture::muf() const
{
	surfaceScalarField alpha1f =
		min(max(fvc::interpolate(alpha1_), scalar(0)), scalar(1));

	return tmp<surfaceScalarField>
	(
		new surfaceScalarField
		(
			"muf_twoPhaseMixture",
			alpha1f*rho1_*fvc::interpolate(nuModel1_->nu())
		  + (scalar(1) - alpha1f)*rho2_*fvc::interpolate(nuModel2_->nu())
		)
	);
}


tmp<surfaceScalarField> twoPhaseMixture::nuf() const
{
	surfaceScalarField alpha1f =
		min(max(fvc::interpolate(alpha1_), scalar(0)), scalar(1));

	return tmp<surfaceScalarField>
	(
		new surfaceScalarField
		(
			"nuf_twoPhaseMixture",
			(
				alpha1f*rho1_*fvc::interpolate(nuModel1_->nu())
			  + (scalar(1) - alpha1f)*rho2_*fvc::interpolate(nuModel2_->nu())
			)/(alpha1f*rho1_ + (scalar(1) - alpha1f)*rho2_)
		)
	);
}


bool twoPhaseMixture::read()
{
	if (transportModel::read())
	{
		if
		(
			nuModel1_().read(subDict(phase1Name_))
		 && nuModel2_().read(subDict(phase2Name_))
		)
		{
			nuModel1_->viscosityProperties().lookup("rho") >> rho1_;
			nuModel2_->viscosityProperties().lookup("rho") >> rho2_;

			return true;
		}
		else
		{
			return false;
		}
	}
	else
	{
		return false;
	}
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
