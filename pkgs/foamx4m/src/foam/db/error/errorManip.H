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

Class
	Foam::errorManip

Description
	Error stream manipulators for exit and abort which may terminate the
	program or throw an exception depending if the exception
	handling has been switched on (off by default).

Usage
	\code
		error << "message1" << "message2" << FoamDataType << exit(error, errNo);
		error << "message1" << "message2" << FoamDataType << abort(error);
	\endcode

\*---------------------------------------------------------------------------*/

#ifndef errorManip_H
#define errorManip_H

#include "error.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators
template<class Err> class errorManip;
template<class Err> Ostream& operator<<(Ostream&, errorManip<Err>);

template<class Err, class T> class errorManipArg;
template<class Err, class T>
Ostream& operator<<(Ostream&, errorManipArg<Err, T>);



template<class Err>
class errorManip
{
	void (Err::*fPtr_)();
	Err& err_;

public:

	errorManip(void (Err::*fPtr)(), Err& t)
	:
		fPtr_(fPtr),
		err_(t)
	{}

	friend Ostream& operator<< <Err>(Ostream& os, errorManip<Err> m);
};


template<class Err>
inline Ostream& operator<<(Ostream& os, errorManip<Err> m)
{
	(m.err_.*m.fPtr_)();
	return os;
}



//- errorManipArg
template<class Err, class T>
class errorManipArg
{
	void (Err::*fPtr_)(const T);
	Err& err_;
	T arg_;

public:

	errorManipArg(void (Err::*fPtr)(const T), Err& t, const T i)
	:
		fPtr_(fPtr),
		err_(t),
		arg_(i)
	{}

	friend Ostream& operator<< <Err, T>(Ostream& os, errorManipArg<Err, T> m);
};


template<class Err, class T>
inline Ostream& operator<<(Ostream& os, errorManipArg<Err, T> m)
{
	(m.err_.*m.fPtr_)(m.arg_);
	return os;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

inline errorManipArg<error, int> exit(error& err, const int errNo = 1)
{
	return errorManipArg<error, int>(&error::exit, err, errNo);
}

inline errorManip<error> abort(error& err)
{
	return errorManip<error>(&error::abort, err);
}


inline errorManipArg<IOerror, int> exit(IOerror& err, const int errNo = 1)
{
	return errorManipArg<IOerror, int>(&IOerror::exit, err, errNo);
}

inline errorManip<IOerror> abort(IOerror& err)
{
	return errorManip<IOerror>(&IOerror::abort, err);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
