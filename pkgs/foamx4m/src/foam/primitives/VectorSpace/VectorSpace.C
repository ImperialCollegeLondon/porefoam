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

#include "VectorSpace.H"
#include "IOstreams.H"

#include <sstream>

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from Istream
template<class Form, class Cmpt, int nCmpt>
Foam::VectorSpace<Form, Cmpt, nCmpt>::VectorSpace
(
	Istream& is
)
{
	// Read beginning of VectorSpace<Cmpt>
	is.readBegin("VectorSpace<Form, Cmpt, nCmpt>");

	for (int i=0; i<nCmpt; i++)
	{
		is >> v_[i];
	}

	// Read end of VectorSpace<Cmpt>
	is.readEnd("VectorSpace<Form, Cmpt, nCmpt>");

	// Check state of Istream
	is.check("VectorSpace<Form, Cmpt, nCmpt>::VectorSpace(Istream&)");
}


// Return a string representation
template<class Form, class Cmpt, int nCmpt>
Foam::word Foam::name
(
	const VectorSpace<Form, Cmpt, nCmpt>& vs
)
{
	std::ostringstream buf;

	buf << '(';

	for (int i=0; i<nCmpt-1; i++)
	{
		buf << vs.v_[i] << ',';
	}

	buf << vs.v_[nCmpt-1] << ')';

	return buf.str();
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class Form, class Cmpt, int nCmpt>
Foam::Istream& Foam::operator>>
(
	Istream& is,
	VectorSpace<Form, Cmpt, nCmpt>& vs
)
{
	// Read beginning of VectorSpace<Cmpt, nCmpt>
	is.readBegin("VectorSpace<Form, Cmpt, nCmpt>");

	for (int i=0; i<nCmpt; i++)
	{
		is >> vs.v_[i];
	}

	// Read end of VectorSpace<Cmpt, nCmpt>
	is.readEnd("VectorSpace<Form, Cmpt, nCmpt>");

	// Check state of Istream
	is.check("operator>>(Istream&, VectorSpace<Form, Cmpt, nCmpt>&)");

	return is;
}


template<class Form, class Cmpt, int nCmpt>
Foam::Ostream& Foam::operator<<
(
	Ostream& os,
	const VectorSpace<Form, Cmpt, nCmpt>& vs
)
{
	os << token::BEGIN_LIST;

	for (int i=0; i<nCmpt-1; i++)
	{
		os << vs.v_[i] << token::SPACE;
	}

	os << vs.v_[nCmpt-1] << token::END_LIST;

	// Check state of Ostream
	os.check("operator<<(Ostream&, const VectorSpace<Form, Cmpt, nCmpt>&)");

	return os;
}


// ************************************************************************* //
