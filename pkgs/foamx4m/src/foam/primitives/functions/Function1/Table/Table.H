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
	Foam::Function1Types::Table

Description
	Templated table container data entry. Items are stored in a list of
	Tuple2's. First column is always stored as scalar entries. Data is read
	in Tuple2 form, e.g. for an entry \<entryName\> that is (scalar, vector):

	\verbatim
		<entryName>   table
		(
			(0.0 (1 2 3))
			(1.0 (4 5 6))
		);
	\endverbatim

SourceFiles
	Table.C

\*---------------------------------------------------------------------------*/

#ifndef Table_H
#define Table_H

#include "Function1.H"
#include "Tuple2.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Function1Types
{


template<class Type>
class Table
:
	public TableBase<Type>
{
	// Private Member Functions

		//- Disallow default bitwise assignment
		void operator=(const Table<Type>&);


public:

	//- Runtime type information
	TypeName("table");


	// Constructors

		//- Construct from entry name and Istream
		Table(const word& entryName, const dictionary& dict);

		//- Copy constructor
		Table(const Table<Type>& tbl);

		//- Construct and return a clone
		virtual tmp<Function1<Type> > clone() const
		{
			return tmp<Function1<Type> >(new Table<Type>(*this));
		}


	//- Destructor
	virtual ~Table();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Function1Types
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
	#include "Table.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
