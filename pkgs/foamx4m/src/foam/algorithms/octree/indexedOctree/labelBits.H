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
	Foam::labelBits

Description
	A 29bits label and 3bits direction packed into single label

SourceFiles

\*---------------------------------------------------------------------------*/

#ifndef labelBits_H
#define labelBits_H

#include "label.H"
#include "uLabel.H"
#include "direction.H"
#include "error.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{



class labelBits
{
	// Private data

		label data_;

		inline static label pack(const label val, const direction bits)
		{
#			ifdef FULLDEBUG
			if (bits > 7 || (((val<<3)>>3) != val))
			{
				FatalErrorIn
				(
					"labelBits::pack(const label, const direction)"
				)   << "Direction " << bits << " outside range 0..7"
					<< abort(FatalError);
			}
#			endif

			return (val<<3) | bits;
		}

public:

	// Constructors

		//- Construct null
		inline labelBits()
		{}

		//- Construct from components
		inline labelBits(const label val, const direction bits)
		:
			 data_(pack(val, bits))
		{}

		//- Construct from Istream
		inline labelBits(Istream& is)
		{
			is >> data_;
		}



	// Member Functions

		inline label val() const
		{
			return data_ >> 3;
		}

		inline direction bits() const
		{
			return data_ & 7;
		}

		inline void setVal(const label val)
		{
			data_ = pack(val, bits());
		}

		inline void setBits(const direction bits)
		{
			data_ = pack(val(), bits);
		}


	// Member Operators

		friend inline bool operator==(const labelBits& a, const labelBits& b)
		{
			return a.data_ == b.data_;
		}

		friend inline bool operator!=(const labelBits& a, const labelBits& b)
		{
			return !(a == b);
		}

	// IOstream Operators

		friend inline Istream& operator>>(Istream& is, labelBits& lb)
		{
			return is >> lb.data_;
		}

		friend inline Ostream& operator<<(Ostream& os, const labelBits& lb)
		{
			return os << lb.data_;
		}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
