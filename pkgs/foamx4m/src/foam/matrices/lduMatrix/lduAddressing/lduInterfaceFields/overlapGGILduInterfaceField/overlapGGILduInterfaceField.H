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
	overlapGGILduInterfaceField

Description
	Abstract base class for overlapGgi coupled interface fields.

Author:
	Hrvoje Jasak, Wikki Ltd.  All rights reserved.

SourceFiles
	overlapGGILduInterfaceField.C

\*---------------------------------------------------------------------------*/

#ifndef overlapGGILduInterfaceField_H
#define overlapGGILduInterfaceField_H

#include "primitiveFieldsFwd.H"
#include "typeInfo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


class overlapGGILduInterfaceField
{

public:

	//- Runtime type information
	TypeName("overlapGGILduInterfaceField");


	// Constructors

		//- Construct given coupled patch
		overlapGGILduInterfaceField()
		{}


	//- Destructor
	virtual ~overlapGGILduInterfaceField();


	// Member Functions

		// Access

			//- Is the transform required
			virtual bool doTransform() const = 0;

			//- Return face transformation tensor
			virtual const tensorField& forwardT() const = 0;

			//- Return neighbour-cell transformation tensor
			virtual const tensorField& reverseT() const = 0;

			//- Return rank of component for transform
			virtual int rank() const = 0;


		//- Transform given patch internal field
		void transformCoupleField
		(
			scalarField& psiInternal,
			const direction cmpt
		) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
