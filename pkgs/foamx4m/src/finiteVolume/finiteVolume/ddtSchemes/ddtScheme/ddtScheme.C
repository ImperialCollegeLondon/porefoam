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

#include "fv.H"
#include "HashTable.H"
#include "surfaceInterpolate.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

template<class Type>
tmp<ddtScheme<Type> > ddtScheme<Type>::New
(
	const fvMesh& mesh,
	Istream& schemeData
)
{
	if (fv::debug)
	{
		Info<< "ddtScheme<Type>::New(const fvMesh&, Istream&) : "
			   "constructing ddtScheme<Type>"
			<< endl;
	}

	if (schemeData.eof())
	{
		FatalIOErrorIn
		(
			"ddtScheme<Type>::New(const fvMesh&, Istream&)",
			schemeData
		)   << "Ddt scheme not specified" << nl << nl
			<< "Valid ddt schemes are :" << endl
			<< IstreamConstructorTablePtr_->sortedToc()
			<< exit(FatalIOError);
	}

	const word schemeName(schemeData);

	typename IstreamConstructorTable::iterator cstrIter =
		IstreamConstructorTablePtr_->find(schemeName);

	if (cstrIter == IstreamConstructorTablePtr_->end())
	{
		FatalIOErrorIn
		(
			"ddtScheme<Type>::New(const fvMesh&, Istream&)",
			schemeData
		)   << "Unknown ddt scheme " << schemeName << nl << nl
			<< "Valid ddt schemes are :" << endl
			<< IstreamConstructorTablePtr_->sortedToc()
			<< exit(FatalIOError);
	}

	return cstrIter()(mesh, schemeData);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
ddtScheme<Type>::~ddtScheme()
{}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
tmp<surfaceScalarField> ddtScheme<Type>::fvcDdtPhiCoeff
(
	const GeometricField<Type, fvPatchField, volMesh>& U,
	const fluxFieldType& phi,
	const fluxFieldType& phiCorr
)
{
	tmp<surfaceScalarField> tddtCouplingCoeff = scalar(1)
	  - min
		(
			mag(phiCorr)
		   /(mag(phi) + dimensionedScalar("small", phi.dimensions(), SMALL)),
			scalar(1)
		);

	surfaceScalarField& ddtCouplingCoeff = tddtCouplingCoeff();

	forAll (U.boundaryField(), patchi)
	{
		if (U.boundaryField()[patchi].fixesValue())
		{
			ddtCouplingCoeff.boundaryField()[patchi] = 0.0;
		}
	}

	if (debug > 1)
	{
		Info<< "ddtCouplingCoeff mean max min = "
			<< gAverage(ddtCouplingCoeff.internalField())
			<< " " << gMax(ddtCouplingCoeff.internalField())
			<< " " << gMin(ddtCouplingCoeff.internalField())
			<< endl;
	}

	return tddtCouplingCoeff;
}


template<class Type>
tmp<surfaceScalarField> ddtScheme<Type>::fvcDdtPhiCoeff
(
	const GeometricField<Type, fvPatchField, volMesh>& U,
	const fluxFieldType& phi
)
{
	dimensionedScalar rDeltaT = 1.0/mesh().time().deltaT();

	tmp<surfaceScalarField> tddtCouplingCoeff = scalar(1)
	  - min
		(
			mag(phi - (mesh().Sf() & fvc::interpolate(U)))
		   /(mag(phi) + dimensionedScalar("small", phi.dimensions(), VSMALL)),
		   //(rDeltaT*mesh().magSf()/mesh().deltaCoeffs()),
			scalar(1)
		);

	surfaceScalarField& ddtCouplingCoeff = tddtCouplingCoeff();

	forAll (U.boundaryField(), patchi)
	{
		if (U.boundaryField()[patchi].fixesValue())
		{
			ddtCouplingCoeff.boundaryField()[patchi] = 0.0;
		}
	}

	if (debug > 1)
	{
		Info<< "ddtCouplingCoeff mean max min = "
			<< gAverage(ddtCouplingCoeff.internalField())
			<< " " << gMax(ddtCouplingCoeff.internalField())
			<< " " << gMin(ddtCouplingCoeff.internalField())
			<< endl;
	}

	return tddtCouplingCoeff;
}


template<class Type>
tmp<surfaceScalarField> ddtScheme<Type>::fvcDdtPhiCoeff
(
	const volScalarField& rho,
	const GeometricField<Type, fvPatchField, volMesh>& rhoU,
	const fluxFieldType& phi
)
{
	dimensionedScalar rDeltaT = 1.0/mesh().time().deltaT();

	tmp<surfaceScalarField> tddtCouplingCoeff = scalar(1)
	  - min
		(
			mag(phi - (mesh().Sf() & fvc::interpolate(rhoU)))
		   /(
				mag(phi) + dimensionedScalar("small", phi.dimensions(), VSMALL)
				//fvc::interpolate(rho)*rDeltaT
				//*mesh().magSf()/mesh().deltaCoeffs()
			),
			scalar(1)
		);

	surfaceScalarField& ddtCouplingCoeff = tddtCouplingCoeff();

	forAll (rhoU.boundaryField(), patchi)
	{
		if (rhoU.boundaryField()[patchi].fixesValue())
		{
			ddtCouplingCoeff.boundaryField()[patchi] = 0.0;
		}
	}

	if (debug > 1)
	{
		Info<< "ddtCouplingCoeff mean max min = "
			<< gAverage(ddtCouplingCoeff.internalField())
			<< " " << gMax(ddtCouplingCoeff.internalField())
			<< " " << gMin(ddtCouplingCoeff.internalField())
			<< endl;
	}

	return tddtCouplingCoeff;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
