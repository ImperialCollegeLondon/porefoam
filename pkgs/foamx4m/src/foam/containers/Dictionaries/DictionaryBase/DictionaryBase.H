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
	Foam::DictionaryBase

Description
	Base dictionary class templated on both the form of doubly-linked list
	it uses as well as the type it holds.

	The double templating allows for the instantiation of forms with or
	without storage management.

Note
	The IDLListType parameter should itself be a template but this confused
	gcc 2.95.2 so it has to be instantiated for T when an instantiation of
	DictionaryBase is requested

See Also
	Dictionary and UDictionary

SourceFiles
	DictionaryBase.C
	DictionaryBaseIO.C

\*---------------------------------------------------------------------------*/

#ifndef DictionaryBase_H
#define DictionaryBase_H

#include "HashTable.H"
#include "wordList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class IDLListType, class T>
class DictionaryBase;

template<class IDLListType, class T>
Ostream& operator<<(Ostream&, const DictionaryBase<IDLListType, T>&);



template<class IDLListType, class T>
class DictionaryBase
:
	public IDLListType
{
	// Private data

		//- HashTable of the entries held on the IDLListType for quick lookup
		HashTable<T*> hashedTs_;


	// Private Member functions

		// Add the IDLListType entries into the HashTable
		void addEntries();


public:

	// Constructors

		//- Null constructor
		DictionaryBase();

		//- Copy construct
		DictionaryBase(const DictionaryBase&);

		//- Construct from Istream using given Istream constructor class
		template<class INew>
		DictionaryBase(Istream&, const INew&);

		//- Construct from Istream using default Istream constructor class
		DictionaryBase(Istream&);


	// Member functions

		// Search and lookup

			//- Search DictionaryBase for given keyword
			bool found(const word&) const;

			//- Find and return an entry if present, otherwise return nullptr
			const T* lookupPtr(const word&) const;

			//- Find and return an entry if present, otherwise return nullptr
			T* lookupPtr(const word&);

			//- Find and return entry
			const T* lookup(const word&) const;

			//- Find and return entry
			T* lookup(const word&);

			//- Return the table of contents
			wordList toc() const;


		// Editing

			//- Add at head of dictionary
			void insert(const word&, T*);

			//- Add at tail of dictionary
			void append(const word&, T*);

			//- Remove and return entry specified by keyword.
			//  Return nullptr if the keyword was not found.
			T* remove(const word&);

			//- Clear the dictionary
			void clear();

			//- Transfer the contents of the argument into this DictionaryBase
			//  and annull the argument.
			void transfer(DictionaryBase<IDLListType, T>&);

	// Member operators

		void operator=(const DictionaryBase&);

		//- Find and return entry
		const T* operator[](const word& key) const
		{
			return lookup(key);
		}

		//- Find and return entry
		T* operator[](const word& key)
		{
			return lookup(key);
		}


	// Ostream operator

		friend Ostream& operator<< <IDLListType, T>
		(
			Ostream&,
			const DictionaryBase<IDLListType, T>&
		);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#	include "DictionaryBase.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
