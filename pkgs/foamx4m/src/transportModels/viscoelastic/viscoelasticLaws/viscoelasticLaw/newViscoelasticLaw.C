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

#include "error.H"
#include "viscoelasticModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::viscoelasticLaw> Foam::viscoelasticLaw::New
(
	const word& name,
	const volVectorField& U,
	const surfaceScalarField& phi,
	const dictionary& dict
)
{
	word typeName = dict.lookup("type");

	Info<< "Selecting viscoelastic model " << typeName << endl;

	dictionaryConstructorTable::iterator cstrIter =
		dictionaryConstructorTablePtr_->find(typeName);

	if (cstrIter == dictionaryConstructorTablePtr_->end())
	{
		FatalErrorIn
		(
			"viscoelasticLaw::New(const word& name, const volVectorField&, "
			"const surfaceScalarField&)"
		)   << "Unknown viscoelasticLaw type " << typeName
			<< endl << endl
			<< "Valid viscoelasticLaw types are :" << endl
			<< dictionaryConstructorTablePtr_->sortedToc()
			<< exit(FatalError);
	}

	return autoPtr<viscoelasticLaw>(cstrIter()(name, U, phi, dict));
}


// ************************************************************************* //
