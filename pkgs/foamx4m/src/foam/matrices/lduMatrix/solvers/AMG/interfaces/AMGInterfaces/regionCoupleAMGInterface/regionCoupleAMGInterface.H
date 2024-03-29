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
	Foam::regionCoupleAMGInterface

Description
	AMG agglomerated regionCouple interface.

Author
	Hrvoje Jasak, Wikki Ltd.  All rights reserved.

SourceFiles
	regionCoupleAMGInterface.C

\*---------------------------------------------------------------------------*/

#ifndef regionCoupleAMGInterface_H
#define regionCoupleAMGInterface_H

#include "ggiAMGInterface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


class regionCoupleAMGInterface
:
	public ggiAMGInterface
{

	// Private data

		//- Is the interface coupled?
		const bool coupled_;


public:

	//- Runtime type information
	TypeName("regionCouple");


	// Constructors

		//- Construct from fine level interface,
		//  local and neighbour restrict addressing
		regionCoupleAMGInterface
		(
			const lduPrimitiveMesh& lduMesh,
			const lduInterfacePtrsList& coarseInterfaces,
			const lduInterface& fineInterface,
			const labelField& localRestrictAddressing,
			const labelField& neighbourRestrictAddressing
		);


	//- Destructor
	virtual ~regionCoupleAMGInterface();


	// Member Functions

		// Access

			//- Return true if interface is coupled
			virtual bool coupled() const
			{
				return coupled_;
			}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
