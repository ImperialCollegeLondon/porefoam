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

\*---------------------------------------------------------------------------*/

#include "IOstreams.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline objectMap::objectMap()
:
	index_(-1),
	masterObjects_(0)
{}


inline objectMap::objectMap(const label index, const labelList& master)
:
	index_(index),
	masterObjects_(master)
{}


inline objectMap::objectMap(Istream& is)
{
	// Read beginning of objectMap
	is.readBegin("objectMap");

	is >> index_ >> static_cast<labelList&>(masterObjects_);

	// Read master of objectMap
	is.readEnd("objectMap");

	// Check state of Istream
	is.check("objectMap::objectMap(Istream&)");
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

label& objectMap::index()
{
	return index_;
}


inline label objectMap::index() const
{
	return index_;
}


inline labelList& objectMap::masterObjects()
{
	return masterObjects_;
}


inline const labelList& objectMap::masterObjects() const
{
	return masterObjects_;
}


// * * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * //

inline bool operator==(const objectMap& a, const objectMap& b)
{
	return
	(
		(a.index_ == b.index_) && (a.masterObjects_ == b.masterObjects_)
	);
}


inline bool operator!=(const objectMap& a, const objectMap& b)
{
	return (!(a == b));
}


// * * * * * * * * * * * * * * * Ostream Operator *  * * * * * * * * * * * * //

inline Ostream& operator<<(Ostream& os, const objectMap& a)
{
	os  << token::BEGIN_LIST
		<< a.index_ << token::SPACE
		<< a.masterObjects_
		<< token::END_LIST;

	// Check state of Ostream
	os.check("Ostream& operator<<(Ostream&, const objectMap&)");

	return os;
}


inline Istream& operator>>(Istream& is, objectMap& a)
{
	is.readBegin("objectMap");
	is  >> a.index_ >> a.masterObjects_;
	is.readEnd("objectMap");

	// Check state of Istream
	is.check("Istream& operator>>(Istream&, objectMap&)");

	return is;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // Master namespace Foam

// ************************************************************************* //
