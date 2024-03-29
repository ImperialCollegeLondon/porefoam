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

#include "IOOutputFilter.H"
#include "foamTime.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class OutputFilter>
Foam::IOOutputFilter<OutputFilter>::IOOutputFilter
(
	const word& outputFilterName,
	const IOobject& ioDict,
	const bool readFromFiles
)
:
	IOdictionary(ioDict),
	OutputFilter(outputFilterName, ioDict.db(), *this, readFromFiles)
{}


template<class OutputFilter>
Foam::IOOutputFilter<OutputFilter>::IOOutputFilter
(
	const word& outputFilterName,
	const objectRegistry& obr,
	const word& dictName,
	const IOobject::readOption rOpt,
	const bool readFromFiles
)
:
	IOdictionary
	(
		IOobject
		(
			dictName,
			obr.time().system(),
			obr,
			rOpt,
			IOobject::NO_WRITE
		)
	),
	OutputFilter(outputFilterName, obr, *this, readFromFiles)
{}


template<class OutputFilter>
Foam::IOOutputFilter<OutputFilter>::IOOutputFilter
(
	const word& outputFilterName,
	const objectRegistry& obr,
	const fileName& dictName,
	const IOobject::readOption rOpt,
	const bool readFromFiles
)
:
	IOdictionary
	(
		IOobject
		(
			dictName,
			obr,
			rOpt,
			IOobject::NO_WRITE
		)
	),
	OutputFilter(outputFilterName, obr, *this, readFromFiles)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class OutputFilter>
Foam::IOOutputFilter<OutputFilter>::~IOOutputFilter()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class OutputFilter>
bool Foam::IOOutputFilter<OutputFilter>::read()
{
	if (regIOobject::read())
	{
		OutputFilter::read(*this);
		return true;
	}
	else
	{
		return false;
	}
}


template<class OutputFilter>
void Foam::IOOutputFilter<OutputFilter>::write()
{
	OutputFilter::write();
}


// ************************************************************************* //
