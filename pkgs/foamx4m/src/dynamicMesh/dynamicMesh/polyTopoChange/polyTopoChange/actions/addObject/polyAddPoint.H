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
	Foam::polyAddPoint

Description
	Class containing data for point addition.

\*---------------------------------------------------------------------------*/

#ifndef polyAddPoint_H
#define polyAddPoint_H

#include "label.H"
#include "point.H"
#include "topoAction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


class polyAddPoint
:
	public topoAction
{
	// Private data

		//- Point to add
		point p_;

		//- Master point
		label masterPointID_;

		//- Point zone ID
		label zoneID_;

		//- Does the point support a cell
		bool inCell_;


	// Private Member Functions

		//- Check
		void check() const
		{
			if (zoneID_ < 0 && !inCell_)
			{
				FatalErrorIn
				(
					"polyAddPoint\n"
					"(\n"
					"    const point& p,\n"
					"    const label masterPointID,\n"
					"    const label zoneID,\n"
					"    const bool inCell\n"
					")"
				)   << "Point is not in a cell and not in a zone.  "
					<< "This is not allowed.\n"
					<< "point: " << p_
					<< " master: " << masterPointID_
					<< " zone: " << zoneID_
					<< abort(FatalError);
			}
		}


public:

	// Static data members

		//- Runtime type information
		TypeName("addPoint");


	// Constructors

		//- Construct null.  Used only for list construction
		polyAddPoint()
		:
			p_(vector::zero),
			masterPointID_(-1),
			zoneID_(-1),
			inCell_(false)
		{}

		//- Construct from components
		polyAddPoint
		(
			const point& p,
			const label masterPointID,
			const label zoneID,
			const bool inCell
		)
		:
			p_(p),
			masterPointID_(masterPointID),
			zoneID_(zoneID),
			inCell_(inCell)
		{
			check();
		}

		//- Construct from Istream
		polyAddPoint(Istream& is)
		:
			p_(is),
			masterPointID_(readLabel(is)),
			zoneID_(readLabel(is)),
			inCell_(readBool(is))
		{
			check();
		}

		//- Construct and return a clone
		virtual autoPtr<topoAction> clone() const
		{
			return autoPtr<topoAction>(new polyAddPoint(*this));
		}


	// Default Destructor


	// Member Functions

		//- Point location
		const point& newPoint() const
		{
			return p_;
		}

		//- Master point label
		label masterPointID() const
		{
			return masterPointID_;
		}

		//- Is the point appended with no master
		bool appended() const
		{
			return masterPointID_ < 0;
		}

		//- Does the point belong to a zone?
		bool isInZone() const
		{
			return zoneID_ >= 0;
		}

		//- Point zone ID
		label zoneID() const
		{
			return zoneID_;
		}

		//- Does the point support a cell
		bool inCell() const
		{
			return inCell_;
		}


	// IOstream Operators

		friend Ostream& operator<<(Ostream& os , const polyAddPoint& pap)
		{
			os << pap.p_ << pap.masterPointID_ << pap.zoneID_ << pap.inCell_;

			return os;
		}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
