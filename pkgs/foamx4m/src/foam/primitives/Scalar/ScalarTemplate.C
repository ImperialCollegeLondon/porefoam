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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

const char* const pTraits<Scalar>::typeName = "scalar";
const Scalar pTraits<Scalar>::zero = 0.0;
const Scalar pTraits<Scalar>::one = 1.0;
const Scalar pTraits<Scalar>::min = -ScalarVGREAT;
const Scalar pTraits<Scalar>::max = ScalarVGREAT;
const Scalar pTraits<Scalar>::rootMin = -ScalarROOTVGREAT;
const Scalar pTraits<Scalar>::rootMax = ScalarROOTVGREAT;

const char* pTraits<Scalar>::componentNames[] = { "x" };

pTraits<Scalar>::pTraits(const Scalar& p)
:
	p_(p)
{}


pTraits<Scalar>::pTraits(Istream& is)
{
	is >> p_;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

word name(const Scalar val)
{
	std::ostringstream buf;
	buf << val;
	return buf.str();
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Scalar readScalar(Istream& is)
{
	Scalar rs;
	is >> rs;

	return rs;
}


Istream& operator>>(Istream& is, Scalar& s)
{
	token t(is);

	if (!t.good())
	{
		is.setBad();
		return is;
	}

	if (t.isNumber())
	{
		s = t.number();
	}
	else
	{
		is.setBad();
		FatalIOErrorIn("operator>>(Istream&, Scalar&)", is)
			<< "wrong token type - expected Scalar, found " << t.info()
			<< exit(FatalIOError);

		return is;
	}

	// Check state of Istream
	is.check("Istream& operator>>(Istream&, Scalar&)");

	return is;
}


Ostream& operator<<(Ostream& os, const Scalar s)
{
	os.write(s);
	os.check("Ostream& operator<<(Ostream&, const Scalar&)");
	return os;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
