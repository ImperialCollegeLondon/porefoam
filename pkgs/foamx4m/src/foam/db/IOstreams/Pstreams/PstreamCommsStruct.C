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

#include "Pstream.H"
#include "boolList.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::Pstream::commsStruct::commsStruct()
:
	above_(-1),
	below_(0),
	allBelow_(0),
	allNotBelow_(0)
{}


Foam::Pstream::commsStruct::commsStruct
(
	const label above,
	const labelList& below,
	const labelList& allBelow,
	const labelList& allNotBelow
)
:
	above_(above),
	below_(below),
	allBelow_(allBelow),
	allNotBelow_(allNotBelow)
{}


Foam::Pstream::commsStruct::commsStruct
(
	const label nProcs,
	const label myProcID,
	const label above,
	const labelList& below,
	const labelList& allBelow
)
:
	above_(above),
	below_(below),
	allBelow_(allBelow),
	allNotBelow_(nProcs - allBelow.size() - 1)
{
	boolList inBelow(nProcs, false);

	forAll (allBelow, belowI)
	{
		inBelow[allBelow[belowI]] = true;
	}

	label notI = 0;
	forAll (inBelow, procI)
	{
		if ((procI != myProcID) && !inBelow[procI])
		{
			allNotBelow_[notI++] = procI;
		}
	}
	if (notI != allNotBelow_.size())
	{
		FatalErrorIn("commsStruct") << "problem!" << Foam::abort(FatalError);
	}
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

bool Foam::Pstream::commsStruct::operator==(const commsStruct& comm) const
{
	return
	(
		(above_ == comm.above())
	 && (below_ == comm.below())
	 && (allBelow_ == allBelow())
	 && (allNotBelow_ == allNotBelow())
	);
}


bool Foam::Pstream::commsStruct::operator!=(const commsStruct& comm) const
{
	return !operator==(comm);
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const Pstream::commsStruct& comm)
{
	os  << comm.above_ << token::SPACE
		<< comm.below_ << token::SPACE
		<< comm.allBelow_ << token::SPACE
		<< comm.allNotBelow_;

	os.check
	(
		"Ostream& operator<<(Ostream&, const commsStruct&)"
	);

	return os;
}


// ************************************************************************* //
