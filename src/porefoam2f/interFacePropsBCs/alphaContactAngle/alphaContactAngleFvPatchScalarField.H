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
    Foam::alphaContactAngleFvPatchScalarField

Description
    Abstract base class for alphaContactAngle boundary conditions.

    Derived classes must implement the theta() fuction which returns the
    wall contact angle field.

    The essential entry "limit" controls the gradient of alpha1 on the wall:

    limit none;         // Calculate the gradient from the contact-angle without
                        // limiter

    limit gradient;     // Limit the wall-gradient such that alpha1 remains
                        // bounded on the wall

    limit alpha;        // Bound the calculated alpha1 on the wall

    limit zeroGradient; // Set the gradient of alpha1 to 0 on the wall
                        // i.e. reproduce previous behaviour

    Note that if any of the first three options are used the boundary condition
    on p_rgh must set to guarantee that the flux is corrected to be zero at the
    wall e.g.

        walls
        {
            type            fixedFluxPressure;
            adjoint         no;
        }

    If "limit zeroGradient;" is used the pressure BCs can be left as before.

SourceFiles
    alphaContactAngleFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef alphaContactAngleFvPatchScalarField_H
#define alphaContactAngleFvPatchScalarField_H

#include "fixedGradientFvPatchFields.H"
#include "fvsPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#ifndef FOAMX // used in official openfoam
 #define BLOCKING Pstream::commsTypes::blocking
#else
 #define BLOCKING Pstream::blocking
#endif
namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class alphaContactAngleFvPatch Declaration
\*---------------------------------------------------------------------------*/

class alphaContactAngleFvPatchScalarField
:
    public fixedGradientFvPatchScalarField
{

public:

		//- Runtime type information
		TypeName("alphaContactAngle");

		//- Alpha limit options
		enum limitControls
		{
		  lcNone,
		  lcGradient,
		  lcZeroGradient,
		  lcAlpha
		};
    
		static bool reset;
		//static bool reset2; 

		static const NamedEnum<limitControls, 4> limitControlNames_;
		limitControls limit_;

		scalar transitionFactor_;
    
    // Constructors

        //- Construct from patch and internal field
        alphaContactAngleFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        alphaContactAngleFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given alphaContactAngleFvPatchScalarField
        //  onto a new patch
        alphaContactAngleFvPatchScalarField
        (
            const alphaContactAngleFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        alphaContactAngleFvPatchScalarField
        (
            const alphaContactAngleFvPatchScalarField&
        );
        //- Construct as copy setting internal field reference
        alphaContactAngleFvPatchScalarField
        (
            const alphaContactAngleFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );


    // Member functions

        //- Return the contact angle
        virtual tmp<scalarField> theta
        (
            const vectorField& Up,
            const vectorField& nHat,
				const vectorField & nWall
        ) const = 0;

        //- Evaluate the patch field
        virtual void evaluate
        (
            const Pstream::commsTypes commsType=BLOCKING
        );

        //- Write
        virtual void write(Ostream&) const;
        
		void noLimit() {limit_ = lcNone;}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
