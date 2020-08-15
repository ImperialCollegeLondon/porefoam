/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2010 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
	This file is partly based on OpenFOAM code.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::hysteresisContactAngleFvPatchScalarField

Description
    A dynamic alphaContactAngle scalar boundary condition
    (alphaContactAngleFvPatchScalarField)

SourceFiles
    hysteresisContactAngleFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef hysteresisContactAngleFvPatchScalarField_H
#define hysteresisContactAngleFvPatchScalarField_H

#include "alphaContactAngleFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class hysteresisContactAngleFvPatch Declaration
\*---------------------------------------------------------------------------*/

class hysteresisContactAngleFvPatchScalarField
:
    public alphaContactAngleFvPatchScalarField
{
	// Private data

		//- Equilibrium contact angle
		//scalar theta0_;

		//- Dynamic contact angle velocity scale
		//scalar uTheta_;

		//- Limiting advancing contact angle
		scalar thetaA_;

		//- Limiting receeding contact angle
		scalar thetaR_;

		//- Limiting advancing contact angle
		scalar thetaAOW_;

		//- Limiting receeding contact angle
		scalar thetaROW_;


		//scalarField minAlphaHist;
		scalarField maxAlphaHist;

public:

	//- Runtime type information
	TypeName("hysteresisContactAngle");


	// Constructors

		//- Construct from patch and internal field
		hysteresisContactAngleFvPatchScalarField
		(
			const fvPatch&,
			const DimensionedField<scalar, volMesh>&
		);

		//- Construct from patch, internal field and dictionary
		hysteresisContactAngleFvPatchScalarField
		(
			const fvPatch&,
			const DimensionedField<scalar, volMesh>&,
			const dictionary&
		);

		//- Construct by mapping given
		//  hysteresisContactAngleFvPatchScalarField
		//  onto a new patch
		hysteresisContactAngleFvPatchScalarField
		(
			const hysteresisContactAngleFvPatchScalarField&,
			const fvPatch&,
			const DimensionedField<scalar, volMesh>&,
			const fvPatchFieldMapper&
		);

		//- Construct as copy
		hysteresisContactAngleFvPatchScalarField
		(
			const hysteresisContactAngleFvPatchScalarField&
		);

		//- Construct and return a clone
		virtual tmp<fvPatchScalarField> clone() const
		{
			return tmp<fvPatchScalarField>
			(
				 new hysteresisContactAngleFvPatchScalarField(*this)
			);
		}

		//- Construct as copy setting internal field reference
		hysteresisContactAngleFvPatchScalarField
		(
			const hysteresisContactAngleFvPatchScalarField&,
			const DimensionedField<scalar, volMesh>&
		);

		//- Construct and return a clone setting internal field reference
		virtual tmp<fvPatchScalarField> clone
		(
			const DimensionedField<scalar, volMesh>& iF
		) const
		{
			return tmp<fvPatchScalarField>
			(
				 new hysteresisContactAngleFvPatchScalarField(*this, iF)
			);
		}


  // Member functions
	  virtual void evaluate
	  (
			const Pstream::commsTypes commsType=Pstream::blocking
	  );

		//- Evaluate and return dynamic contact-angle
		virtual tmp<scalarField> theta
		(
			const vectorField & Up,
			const vectorField & nHat,
			const vectorField & nWall
		) const;

		//- Write
		virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
