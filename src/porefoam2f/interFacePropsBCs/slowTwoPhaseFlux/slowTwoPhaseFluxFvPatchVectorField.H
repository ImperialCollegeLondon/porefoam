/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    Foam::slowTwoPhaseFluxFvPatchVectorField

Description
    Velocity inlet/outlet boundary condition for patches where the pressure is
    specified.  zero-gradient is applied for outflow (as defined by the flux)
    and for inflow the velocity is obtained from the flux with a direction
    normal to the patch faces.

SourceFiles
    slowTwoPhaseFluxFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef slowTwoPhaseFluxFvPatchVectorField_H
#define slowTwoPhaseFluxFvPatchVectorField_H

#include "fvPatchFields.H"
#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
          Class slowTwoPhaseFluxFvPatch Declaration
\*---------------------------------------------------------------------------*/

class slowTwoPhaseFluxFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
    // Private data

	//word phiName_;
	//word rhoName_;


	//- Mean value
	scalar flowRate0_;
	scalar flowRate1_;
	scalar gradientFactor0_;
	scalar gradientFactor1_;
	scalar pdFactor_;
	scalar pcFactor_;

	label curTimeIndex_;

	  //- Relaxation factor 0 -> fixedValue
	  //scalar relaxationFactor_;
	  
public:

    //- Runtime type information
    TypeName("slowTwoPhaseFlux");


    // Constructors

        //- Construct from patch and internal field
        slowTwoPhaseFluxFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        slowTwoPhaseFluxFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given slowTwoPhaseFluxFvPatchVectorField
        //  onto a new patch
        slowTwoPhaseFluxFvPatchVectorField
        (
            const slowTwoPhaseFluxFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        slowTwoPhaseFluxFvPatchVectorField
        (
            const slowTwoPhaseFluxFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new slowTwoPhaseFluxFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        slowTwoPhaseFluxFvPatchVectorField
        (
            const slowTwoPhaseFluxFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField>
            (
                new slowTwoPhaseFluxFvPatchVectorField(*this, iF)
            );
        }


    // Member functions

        // Access



        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;


    // Member operators

        //virtual void operator=(const fvPatchField<vector>& pvf);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
