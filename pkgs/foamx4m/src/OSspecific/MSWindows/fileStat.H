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
	Foam::fileStat

Description
	Wrapper for stat() system call.

Warning
	on Linux (an maybe on others) a stat() of an nfs mounted (remote)
	file does never timeout and cannot be interrupted!
	So e.g. Foam::ping first and hope nfs is running.

SourceFiles
	fileStat.C

\*---------------------------------------------------------------------------*/

#ifndef fileStat_H
#define fileStat_H

#include <sys/stat.h>
#include <sys/types.h>

#include "label.H"
#include "fileName.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


class fileStat
{
	// Private data

		struct stat status_;

		bool isValid_;


public:

	// Constructors

		//- Empty constructor
		fileStat();

		//- Construct from components
		fileStat(const fileName& fName, const unsigned int maxTime=0);

		//- Construct from Istream
		fileStat(Istream&);


	// Member Functions

		// Access

			//- Raw status
			const struct stat& status() const
			{
				return status_;
			}

			//- Did constructor fail
			bool isValid() const
			{
				return isValid_;
			}


		// Check

			//- compare two fileStats for same device
			bool sameDevice(const fileStat& stat2) const;

			//- compare two fileStats for same Inode
			bool sameINode(const fileStat& stat2) const;

			//- compare state against inode
			bool sameINode(const label iNode) const;


	// IOstream Operators

		friend Istream& operator>>(Istream&, fileStat&);
		friend Ostream& operator<<(Ostream&, const fileStat&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
