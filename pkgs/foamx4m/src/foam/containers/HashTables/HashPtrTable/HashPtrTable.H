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
	Foam::HashPtrTable

Description
	A HashTable specialization for hashing pointers.

SourceFiles
	HashPtrTable.C
	HashPtrTableIO.C

\*---------------------------------------------------------------------------*/

#ifndef HashPtrTable_H
#define HashPtrTable_H

#include "HashTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class Istream;
class Ostream;

// Forward declaration of friend functions and operators

template<class T, class Key, class Hash> class HashPtrTable;

template<class T, class Key, class Hash>
Istream& operator>>(Istream&, HashPtrTable<T, Key, Hash>&);

template<class T, class Key, class Hash>
Ostream& operator<<(Ostream&, const HashPtrTable<T, Key, Hash>&);



template<class T, class Key=word, class Hash=string::hash>
class HashPtrTable
:
	public HashTable<T*, Key, Hash>
{
	// Private Member Functions

		//- Read from Istream using given Istream constructor class
		template<class INew>
		void read(Istream&, const INew& inewt);

		//- Read from dictionary using given dictionary constructor class
		template<class INew>
		void read(const dictionary& dict, const INew& inewt);



public:

	typedef typename HashTable<T*, Key, Hash>::iterator iterator;
	typedef typename HashTable<T*, Key, Hash>::const_iterator const_iterator;


	// Constructors

		//- Construct given initial table size
		HashPtrTable(const label size = 128);

		//- Construct from Istream using given Istream constructor class
		template<class INew>
		HashPtrTable(Istream&, const INew&);

		//- Construct from Istream using default Istream constructor class
		HashPtrTable(Istream&);

		//- Construct from dictionary using default dictionary constructor
		//  class
		HashPtrTable(const dictionary&);

		//- Construct as copy
		HashPtrTable(const HashPtrTable<T, Key, Hash>&);


	//- Destructor
	~HashPtrTable();


	// Member Functions

		// Edit

			//- Remove and return the pointer specified by given iterator
			T* remove(iterator&);

			//- Erase an hashedEntry specified by given iterator
			bool erase(iterator&);

			//- Clear all entries from table
			void clear();

			//- Write
			void write(Ostream& os) const;


	// Member Operators

		void operator=(const HashPtrTable<T, Key, Hash>&);


	// IOstream Operators

		friend Istream& operator>> <T, Key, Hash>
		(
			Istream&,
			HashPtrTable<T, Key, Hash>&
		);

		friend Ostream& operator<< <T, Key, Hash>
		(
			Ostream&,
			const HashPtrTable<T, Key, Hash>&
		);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#	include "HashPtrTable.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
