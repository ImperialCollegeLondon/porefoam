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
	Foam::IOdictionary

Description
	IOdictionary is derived from dictionary and IOobject to give the dictionary
	automatic IO functionality via the objectRegistry.  To facilitate IO,
	IOdictioanry is provided with a constructor from IOobject and writeData and
	write functions.

SourceFiles
	IOdictionary.C
	IOdictionaryIO.C

\*---------------------------------------------------------------------------*/

#ifndef IOdictionary_H
#define IOdictionary_H

#include "dictionary.H"
#include "regIOobject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


class IOdictionary
:
	public regIOobject,
	public dictionary
{
	// Private data

		static bool writeDictionaries;


	// Private Member Functions

		//- Read dictionary from file
		void readFile(const bool);


public:

	TypeName("dictionary");


	// Constructors

		//- Construct given an IOobject
		IOdictionary(const IOobject&);

		//- Construct given an IOobject and dictionary
		IOdictionary(const IOobject&, const dictionary&);

		//- Construct given an IOobject and Istream
		IOdictionary(const IOobject&, Istream&);


	//- Destructor
	virtual ~IOdictionary();


	// Member functions

		//- Name function is needed to disambiguate those inherited
		//  from regIOobject and dictionary
		const word& name() const;

		//- ReadData function required for regIOobject read operation
		bool readData(Istream&);

		//- WriteData function required for regIOobject write operation
		bool writeData(Ostream&) const;


	// Member operators

		//- Assignment of other IOdictionary's entries to this IOdictionary
		void operator=(const IOdictionary&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
