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
	Foam::cloudDistribute

Description
	Lagrangian field decomposer.

SourceFiles
	cloudDistribute.C

\*---------------------------------------------------------------------------*/

#ifndef cloudDistribute_H
#define cloudDistribute_H

#include "labelIOList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


class cloudDistribute
{
	// Private Member Functions

		//- Disallow default bitwise copy construct
		cloudDistribute(const cloudDistribute&);

		//- Disallow default bitwise assignment
		void operator=(const cloudDistribute&);


public:

	// Constructors

		//- Null constructor
		cloudDistribute();


	//- Destructor
	virtual ~cloudDistribute();


	// Member Functions

		// Send particles
		virtual void send(Ostream& toProc, const label procIndex);

		// Receive particles
		virtual void receive(Istream& fromProc, const label procIndex);

		// Receive particles
		virtual void rebuild
		(
			const PtrList<labelIOList>& cellProcAddressing,
			const PtrList<labelIOList>& faceProcAddressing
		);

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
