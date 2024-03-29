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

#include "Function1.H"
#include "foamTime.H"

// * * * * * * * * * * * * * * * * Constructor * * * * * * * * * * * * * * * //

template<class Type>
Foam::Function1<Type>::Function1(const word& entryName)
:
	name_(entryName)
{}


template<class Type>
Foam::Function1<Type>::Function1(const Function1<Type>& de)
:
	refCount(),
	name_(de.name_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::Function1<Type>::~Function1()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
const Foam::word& Foam::Function1<Type>::name() const
{
	return name_;
}


template<class Type>
void Foam::Function1<Type>::convertTimeBase(const Time&)
{}


template<class Type>
Type Foam::Function1<Type>::value(const scalar x) const
{
	NotImplemented;

	return pTraits<Type>::zero;
}


template<class Type>
Type Foam::Function1<Type>::integrate(const scalar x1, const scalar x2) const
{
	NotImplemented;

	return pTraits<Type>::zero;
}


template<class Type>
Foam::tmp<Foam::Field<Type> > Foam::Function1<Type>::value
(
	const scalarField& x
) const
{
	tmp<Field<Type> > tfld(new Field<Type>(x.size()));
	Field<Type>& fld = tfld();

	forAll(x, i)
	{
		fld[i] = this->value(x[i]);
	}
	return tfld;
}


template<class Type>
Foam::tmp<Foam::Field<Type> > Foam::Function1<Type>::integrate
(
	const scalarField& x1,
	const scalarField& x2
) const
{
	tmp<Field<Type> > tfld(new Field<Type>(x1.size()));
	Field<Type>& fld = tfld();

	forAll(x1, i)
	{
		fld[i] = this->integrate(x1[i], x2[i]);
	}
	return tfld;
}


template<class Type>
void Foam::Function1<Type>::writeData(Ostream& os) const
{
	os.writeKeyword(name_) << type();
}


// * * * * * * * * * * * * * *  IOStream operators * * * * * * * * * * * * * //

template<class Type>
Foam::Ostream& Foam::operator<<
(
	Ostream& os,
	const Function1<Type>& f1
)
{
	// Check state of Ostream
	os.check
	(
		"Ostream& operator<<(Ostream&, const Function1<Type>&)"
	);

	os  << f1.name_;
	f1.writeData(os);

	return os;
}


// ************************************************************************* //
