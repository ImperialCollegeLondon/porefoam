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

#include "LPtrList.H"
#include "Istream.H"
#include "Ostream.H"
#include "INew.H"

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

template<class LListBase, class T>
template<class INew>
void Foam::LPtrList<LListBase, T>::read(Istream& is, const INew& iNew)
{
	is.fatalCheck
	(
		"LPtrList<LListBase, T>::read(Istream&, const INew&)"
	);

	token firstToken(is);

	is.fatalCheck
	(
		"LPtrList<LListBase, T>::read(Istream&, const INew&) : "
		"reading first token"
	);

	if (firstToken.isLabel())
	{
		label s = firstToken.labelToken();

		// Read beginning of contents
		char delimiter = is.readBeginList("LPtrList<LListBase, T>");

		if (s)
		{
			if (delimiter == token::BEGIN_LIST)
			{
				for (label i=0; i<s; ++i)
				{
					this->append(iNew(is).ptr());

					is.fatalCheck
					(
					    "LPtrList<LListBase, T>::read(Istream&, const INew&) : "
					    "reading entry"
					);
				}
			}
			else
			{
				T* tPtr = iNew(is).ptr();
				this->append(tPtr);

				is.fatalCheck
				(
					"LPtrList<LListBase, T>::read(Istream&, const INew&) : "
					"reading entry"
				);

				for (label i=1; i<s; ++i)
				{
					this->append(tPtr->clone().ptr());
				}
			}
		}

		// Read end of contents
		is.readEndList("LPtrList<LListBase, T>");
	}
	else if (firstToken.isPunctuation())
	{
		if (firstToken.pToken() != token::BEGIN_LIST)
		{
			FatalIOErrorIn
			(
				"LPtrList<LListBase, T>::read(Istream&, const INew&)",
				is
			)   << "incorrect first token, '(', found " << firstToken.info()
				<< exit(FatalIOError);
		}

		token lastToken(is);
		is.fatalCheck("LPtrList<LListBase, T>::read(Istream&, const INew&)");

		while
		(
		   !(
				lastToken.isPunctuation()
			 && lastToken.pToken() == token::END_LIST
			)
		)
		{
			is.putBack(lastToken);
			this->append(iNew(is).ptr());

			is >> lastToken;
			is.fatalCheck
			(
				"LPtrList<LListBase, T>::read(Istream&, const INew&)"
			);
		}
	}
	else
	{
		FatalIOErrorIn
		(
			"LPtrList<LListBase, T>::read(Istream&, const INew&)",
			is
		)   << "incorrect first token, expected <int> or '(', found "
			<< firstToken.info()
			<< exit(FatalIOError);
	}

	is.fatalCheck("LPtrList<LListBase, T>::read(Istream&, const INew&)");
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class LListBase, class T>
template<class INew>
Foam::LPtrList<LListBase, T>::LPtrList(Istream& is, const INew& iNew)
{
	this->read(is, iNew);
}


template<class LListBase, class T>
Foam::LPtrList<LListBase, T>::LPtrList(Istream& is)
{
	this->read(is, INew<T>());
}


// * * * * * * * * * * * * * * * Istream Operator  * * * * * * * * * * * * * //

template<class LListBase, class T>
Foam::Istream& Foam::operator>>(Istream& is, LPtrList<LListBase, T>& L)
{
	L.clear();
	L.read(is, INew<T>());

	return is;
}


// * * * * * * * * * * * * * * * Ostream Operators * * * * * * * * * * * * * //

template<class LListBase, class T>
Foam::Ostream& Foam::operator<<(Ostream& os, const LPtrList<LListBase, T>& lst)
{
	// Write size
	os << nl << lst.size();

	// Write beginning of contents
	os << nl << token::BEGIN_LIST << nl;

	// Write contents
	for
	(
		typename LPtrList<LListBase, T>::const_iterator iter = lst.begin();
		iter != lst.end();
		++iter
	)
	{
		os << iter() << nl;
	}

	// Write end of contents
	os << token::END_LIST;

	// Check state of IOstream
	os.check("Ostream& operator<<(Ostream&, const LPtrList<LListBase, T>&)");

	return os;
}

// ************************************************************************* //
