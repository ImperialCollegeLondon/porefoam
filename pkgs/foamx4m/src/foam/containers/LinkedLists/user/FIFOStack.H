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
	Foam::FIFOStack

Description
	A FIFO stack based on a singly-linked list.
	Operations are push(), pop(), top(), bottom() and empty().

SourceFiles
	FIFOStack.C

\*---------------------------------------------------------------------------*/

#ifndef FIFOStack_H
#define FIFOStack_H

#include "SLList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


template<class T>
class FIFOStack
:
	public SLList<T>
{

public:

	// Constructors

		//- Construct null
		FIFOStack()
		{}

		//- Construct given initial T
		FIFOStack(T a)
		:
			SLList<T>(a)
		{}

		//- Construct from Istream
		FIFOStack(Istream& is)
		:
			SLList<T>(is)
		{}


	// Member Functions

		// Access

			//- Return a copy of the top element
			T top() const
			{
				return this->last();
			}

			//- Return a copy of the bottom element
			T bottom() const
			{
				return this->first();
			}


		// Edit

			//- Push an element onto the stack
			void push(const T& a)
			{
				this->append(a);
			}

			//- Pop the bottom element off the stack
			T pop()
			{
				return this->removeHead();
			}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
