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

#include "token.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const char* const Foam::token::typeName = "token";
Foam::token Foam::token::undefinedToken;

defineTypeNameAndDebug(Foam::token::compound, 0);
defineRunTimeSelectionTable(Foam::token::compound, Istream);


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

void Foam::token::parseError(const char* expected) const
{
	FatalIOError
		<< "Parse error, expected a " << expected
		<< ", found \n    " << info() << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::token::compound::~compound()
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::token::compound> Foam::token::compound::New
(
	const word& compoundType,
	Istream& is
)
{
	IstreamConstructorTable::iterator cstrIter =
		IstreamConstructorTablePtr_->find(compoundType);

	if (cstrIter == IstreamConstructorTablePtr_->end())
	{
		FatalIOErrorInFunction(is)
			<< "Unknown compound type " << compoundType << nl << nl
			<< "Valid compound types:" << endl
			<< IstreamConstructorTablePtr_->sortedToc()
			<< abort(FatalIOError);
	}

	return autoPtr<Foam::token::compound>(cstrIter()(is));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::token::compound::isCompound(const word& name)
{
	return
	(
		IstreamConstructorTablePtr_
	 && IstreamConstructorTablePtr_->found(name)
	);
}


Foam::token::compound& Foam::token::transferCompoundToken(const Istream& is)
{
	if (type_ == COMPOUND)
	{
		if (compoundTokenPtr_->empty())
		{
			FatalIOErrorInFunction(is)
				<< "compound has already been transfered from token\n    "
				<< info() << abort(FatalIOError);
		}
		else
		{
			compoundTokenPtr_->empty() = true;
		}

		return *compoundTokenPtr_;
	}
	else
	{
		parseError("compound");
		return *compoundTokenPtr_;
	}
}


// ************************************************************************* //
