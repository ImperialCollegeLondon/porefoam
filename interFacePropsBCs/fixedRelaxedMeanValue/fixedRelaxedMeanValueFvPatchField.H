/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is partly based on OpenFOAM code.

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
    fixedRelaxedMeanValueFvPatchField

Description

SourceFiles
    fixedRelaxedMeanValueFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef fixedRelaxedMeanValueFvPatchField_H
#define fixedRelaxedMeanValueFvPatchField_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
               Class fixedRelaxedMeanValueFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class fixedRelaxedMeanValueFvPatchField
:
    public fixedValueFvPatchField<Type>
{
    // Private data

        //- Mean value
        Type meanValue1_;
        Type meanValue2_;

        //- Current time index
        label curTimeIndex_;
        
        //- Relaxation factor 0 -> fixedValue
        scalar relaxationFactor_;
        scalar thicknessFactor_;


public:

    //- Runtime type information
    TypeName("fixedRelaxedMeanValue");


    // Constructors

        //- Construct from patch and internal field
        fixedRelaxedMeanValueFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        fixedRelaxedMeanValueFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given fixedRelaxedMeanValueFvPatchField
        //  onto a new patch
        fixedRelaxedMeanValueFvPatchField
        (
            const fixedRelaxedMeanValueFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<Type> > clone() const
        {
            return tmp<fvPatchField<Type> >
            (
                new fixedRelaxedMeanValueFvPatchField<Type>(*this)
            );
        }

        //- Construct as copy setting internal field reference
        fixedRelaxedMeanValueFvPatchField
        (
            const fixedRelaxedMeanValueFvPatchField<Type>&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<Type> > clone
        (
            const DimensionedField<Type, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<Type> >
            (
                new fixedRelaxedMeanValueFvPatchField<Type>(*this, iF)
            );
        }


    // Member functions

        // Access

            //- Return the ref value
            //const Type& meanValue() const
            //{
                //return meanValue_;
            //}

            //- Return reference to the ref value to allow adjustment
            //Type& meanValue()
            //{
                //return meanValue_;
            //}


        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchField<Type>&,
                const labelList&
            );


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "fixedRelaxedMeanValueFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
