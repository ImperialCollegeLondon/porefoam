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
	Foam::ggiSAMGInterfaceField

Description
	AMG selected ggi interface field.

Author
	Hrvoje Jasak, Wikki Ltd.  All rights reserved.

SourceFiles
	ggiSAMGInterfaceField.C

\*---------------------------------------------------------------------------*/

#ifndef ggiSAMGInterfaceField_H
#define ggiSAMGInterfaceField_H

#include "SAMGInterfaceField.H"
#include "ggiSAMGInterface.H"
#include "ggiLduInterfaceField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


class ggiSAMGInterfaceField
:
	public SAMGInterfaceField,
	virtual public ggiLduInterfaceField
{
	// Private data

		//- Local reference cast into the ggi interface
		const ggiSAMGInterface& ggiInterface_;

		//- Is the transform required
		bool doTransform_;

		//- Rank of component for transformation
		int rank_;


	// Private Member Functions

		//- Disallow default bitwise copy construct
		ggiSAMGInterfaceField(const ggiSAMGInterfaceField&);

		//- Disallow default bitwise assignment
		void operator=(const ggiSAMGInterfaceField&);


public:

	//- Runtime type information
	TypeName("ggi");


	// Constructors

		//- Construct from AMG interface and fine level interface field
		ggiSAMGInterfaceField
		(
			const SAMGInterface& SAMGCp,
			const lduInterfaceField& fineInterfaceField
		);


	//- Destructor
	virtual ~ggiSAMGInterfaceField();


	// Member Functions

		// Access

			//- Return size
			label size() const
			{
				return ggiInterface_.size();
			}


		// Coupled interface matrix update

			//- Transform neighbour field
			virtual void transformCoupleField
			(
				scalarField& pnf,
				const direction cmpt
			) const
			{
				ggiLduInterfaceField::transformCoupleField(pnf, cmpt);
			}

			//- Initialise neighbour matrix update
			virtual void initInterfaceMatrixUpdate
			(
				const scalarField& psiInternal,
				scalarField& result,
				const lduMatrix& m,
				const scalarField& coeffs,
				const direction cmpt,
				const Pstream::commsTypes commsType,
				const bool switchToLhs
			) const;

			//- Update result field based on interface functionality
			virtual void updateInterfaceMatrix
			(
				const scalarField& psiInternal,
				scalarField& result,
				const lduMatrix&,
				const scalarField& coeffs,
				const direction cmpt,
				const Pstream::commsTypes commsType,
				const bool switchToLhs
			) const;


		//- Ggi interface functions

			//- Does the interface field perform the transfromation
			virtual bool doTransform() const
			{
				return doTransform_;
			}

			//- Return face transformation tensor
			virtual const tensorField& forwardT() const
			{
				return ggiInterface_.forwardT();
			}

			//- Return neighbour-cell transformation tensor
			virtual const tensorField& reverseT() const
			{
				return ggiInterface_.reverseT();
			}

			//- Return rank of component for transform
			virtual int rank() const
			{
				return rank_;
			}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
