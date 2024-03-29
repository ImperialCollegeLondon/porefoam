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

Description
	Stream operators for token

\*---------------------------------------------------------------------------*/

#include "error.H"
#include "token.H"

#include "IOstreams.H"
#include "scalar.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::token::token(Istream& is)
:
	type_(UNDEFINED)
{
	is.read(*this);
}


// * * * * * * * * * * * * IOstream operators  * * * * * * * * * * * * * * * //

Foam::Istream& Foam::operator>>(Istream& is, token& t)
{
	t.clear();
	return is.read(t);
}


Foam::Ostream& Foam::operator<<(Ostream& os, const token& t)
{
	switch (t.type_)
	{
		case token::UNDEFINED:
			os << "UNDEFINED";
			WarningIn("Ostream& operator<<(Ostream&, const token&)")
				<< "Undefined token" << endl;
		break;

		case token::PUNCTUATION:
			os << t.punctuationToken_;
		break;

		case token::WORD:
			os << *t.wordTokenPtr_;
		break;

		case token::STRING:
			os << *t.stringTokenPtr_;
		break;

		case token::LABEL:
			os << t.labelToken_;
		break;

		case token::FLOAT_SCALAR:
			os << t.floatScalarToken_;
		break;

		case token::DOUBLE_SCALAR:
			os << t.doubleScalarToken_;
		break;

		case token::LONG_DOUBLE_SCALAR:
			os << t.longDoubleScalarToken_;
		break;

		case token::COMPOUND:
			os << *t.compoundTokenPtr_;
		break;

		case token::ERROR:
			os << "FATALERROR";
			WarningIn("Ostream& operator<<(Ostream&, const token&)")
				<< "Error token" << endl;
		break;

		default:
			os << "UNKNOWN";
			SeriousErrorIn("Ostream& operator<<(Ostream&, const token&)")
				<< "Unknown token"
				<< endl;
	}

	// Check state of stream
	os.check("Ostream& operator<<(Ostream&, const token&)");

	return os;
}


ostream& Foam::operator<<(ostream& os, const token::punctuationToken& pt)
{
	return os << char(pt);
}


Foam::Ostream& Foam::operator<<(Ostream& os, const token::punctuationToken& pt)
{
	return os << char(pt);
}


Foam::Ostream& Foam::operator<<(Ostream& os, const token::compound& ct)
{
	os << ct.type() << token::SPACE;
	ct.write(os);

	return os;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

ostream& Foam::operator<<(ostream& os, const InfoProxy<token>& ip)
{
	const token& t = ip.t_;

	os  << "on line " << t.lineNumber();

	switch (t.type())
	{
		case token::UNDEFINED:
			os  << " an undefined token";
		break;

		case token::PUNCTUATION:
			os  << " the punctuation token " << '\'' << t.pToken() << '\'';
		break;

		case token::WORD:
			os  << " the word " << '\'' << t.wordToken() << '\'';
		break;

		case token::STRING:
			os  << " the string " << t.stringToken();
		break;

		case token::LABEL:
			os  << " the label " << t.labelToken();
		break;

		case token::FLOAT_SCALAR:
			os  << " the floatScalar " << t.floatScalarToken();
		break;

		case token::DOUBLE_SCALAR:
			os  << " the doubleScalar " << t.doubleScalarToken();
		break;

		case token::LONG_DOUBLE_SCALAR:
			os  << " the longDoubleScalar " << t.doubleScalarToken();
		break;

		case token::COMPOUND:
		{
			if (t.compoundToken().empty())
			{
				os  << " the empty compound of type "
					<< t.compoundToken().type();
			}
			else
			{
				os  << " the compound of type "
					<< t.compoundToken().type();
			}
		}
		break;

		case token::ERROR:
			os  << " an error";
		break;

		default:
			os  << " an unknown token type " << '\'' << int(t.type()) << '\'';
	}

	return os;
}


// template specialization
namespace Foam
{

#if defined (__GNUC__)
template<>
#endif
Ostream& operator<<(Ostream& os, const InfoProxy<token>& ip)
{
	const token& t = ip.t_;

	os  << "on line " << t.lineNumber();

	switch (t.type())
	{
		case token::UNDEFINED:
			os  << " an undefined token";
		break;

		case token::PUNCTUATION:
			os  << " the punctuation token " << '\'' << t.pToken() << '\'';
		break;

		case token::WORD:
			os  << " the word " << '\'' << t.wordToken() << '\'';
		break;

		case token::STRING:
			os  << " the string " << t.stringToken();
		break;

		case token::LABEL:
			os  << " the label " << t.labelToken();
		break;

		case token::FLOAT_SCALAR:
			os  << " the floatScalar " << t.floatScalarToken();
		break;

		case token::DOUBLE_SCALAR:
			os  << " the doubleScalar " << t.doubleScalarToken();
		break;

		case token::LONG_DOUBLE_SCALAR:
			os  << " the longDoubleScalar " << t.longDoubleScalarToken();
		break;

		case token::COMPOUND:
		{
			if (t.compoundToken().empty())
			{
				os  << " the empty compound of type "
					<< t.compoundToken().type();
			}
			else
			{
				os  << " the compound of type "
					<< t.compoundToken().type();
			}
		}
		break;

		case token::ERROR:
			os  << " an error";
		break;

		default:
			os  << " an unknown token type "  << '\'' << int(t.type()) << '\'';
	}

	return os;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
