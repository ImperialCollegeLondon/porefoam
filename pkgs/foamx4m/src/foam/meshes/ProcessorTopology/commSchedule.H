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
	Foam::commSchedule

Description
	Determines the order in which a set of processors should communicate
	with one another.

	The communication order should
	  - have maximum overlap
	  - allow blocking communication without deadlock

	Does a very simple scheduling which assumes same time for all operations.

	After construction:
	  - schedule() gives the order in which the input communication should occur
	  - procSchedule()[procI] gives per procI

	Does not care whether 'talking' is first send, second receive or maybe
	full swap. This is all responsability of caller. See ProcessorTopology
	class for use in scheduling processor boundary swaps.

SourceFiles
	commSchedule.C

\*---------------------------------------------------------------------------*/

#ifndef commSchedule_H
#define commSchedule_H

#include "dynamicLabelList.H"
#include "labelPair.H"
#include "labelList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


class commSchedule
{
	// Private data

		//- Order in which input communication has been scheduled
		labelList schedule_;

		//- Per processor the order in which communication has been scheduled
		labelListList procSchedule_;

	// Private Member Functions

		//- Count the number of outstanding communications for a single
		//  processor
		label outstandingComms(const labelList&, dynamicLabelList&) const;


public:

	ClassName("commSchedule");

	// Constructors

		//- Construct from wanted communication. Wanted communication is between
		//  two processors. Can be a one-way communication or
		//  two-way communication, that is up to the caller. This class just
		//  determines an order for it such that any processor is only talking
		//  to one other at a time. After construction:
		//  - schedule is the order in which comms is done.
		//  - procSchedule[procI] is for procI the order in which comms is done.
		commSchedule(const label nProcs, const List<labelPair>& comms);


	// Member Functions

		//- order in which comms is scheduled
		const labelList& schedule() const
		{
			return schedule_;
		}

		//- Per processor the order in which communication has been scheduled
		const labelListList& procSchedule() const
		{
			return procSchedule_;
		}


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
