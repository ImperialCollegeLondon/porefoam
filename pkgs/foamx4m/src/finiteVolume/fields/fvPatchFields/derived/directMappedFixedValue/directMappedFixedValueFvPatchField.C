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

\*---------------------------------------------------------------------------*/

#include "directMappedFixedValueFvPatchField.H"
#include "directMappedPatchBase.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
void directMappedFixedValueFvPatchField<Type>::mapField()
{
	typedef GeometricField<Type, fvPatchField, volMesh> fieldType;

	// Get the scheduling information from the directMappedPatchBase
	const directMappedPatchBase& mpp = refCast<const directMappedPatchBase>
	(
		directMappedFixedValueFvPatchField<Type>::patch().patch()
	);
	const mapDistribute& distMap = mpp.map();

	// Force recalculation of schedule
	distMap.schedule();

	const fvMesh& nbrMesh = refCast<const fvMesh>(mpp.sampleMesh());
	const word& fldName = this->dimensionedInternalField().name();

	// Result of obtaining remote values
	if (debug)
	{
		Info<< "direct mapping field "
			<< this->dimensionedInternalField().name() << endl;
	}

	switch (mpp.mode())
	{
		case directMappedPatchBase::NEARESTCELL:
		{
			if (mpp.sameRegion())
			{
				newValues_ = this->internalField();
			}
			else
			{
				newValues_ = nbrMesh.lookupObject<fieldType>
				(
					fldName
				).internalField();
			}
			mapDistribute::distribute
			(
				Pstream::defaultComms(),
				distMap.schedule(),
				distMap.constructSize(),
				distMap.subMap(),
				distMap.constructMap(),
				newValues_
			);

			break;
		}
		case directMappedPatchBase::NEARESTPATCHFACE:
		{
			const label nbrPatchID = nbrMesh.boundaryMesh().findPatchID
			(
				mpp.samplePatch()
			);
			if (nbrPatchID < 0)
			{
				FatalErrorIn
				(
					"void directMappedFixedValueFvPatchField<Type>::"
					"updateCoeffs()"
				)<< "Unable to find sample patch " << mpp.samplePatch()
				 << " in region " << mpp.sampleRegion()
				 << " for patch " << this->patch().name() << nl
				 << abort(FatalError);
			}

			const fieldType& nbrField = nbrMesh.lookupObject<fieldType>
			(
				fldName
			);
			newValues_ = nbrField.boundaryField()[nbrPatchID];
			mapDistribute::distribute
			(
				Pstream::defaultComms(),
				distMap.schedule(),
				distMap.constructSize(),
				distMap.subMap(),
				distMap.constructMap(),
				newValues_
			);

			break;
		}
		case directMappedPatchBase::NEARESTFACE:
		{
			Field<Type> allValues(nbrMesh.nFaces(), pTraits<Type>::zero);

			const fieldType& nbrField = nbrMesh.lookupObject<fieldType>
			(
				fldName
			);
			forAll(nbrField.boundaryField(), patchI)
			{
				const fvPatchField<Type>& pf =
					nbrField.boundaryField()[patchI];
				label faceStart = pf.patch().patch().start();

				forAll(pf, faceI)
				{
					allValues[faceStart++] = pf[faceI];
				}
			}

			mapDistribute::distribute
			(
				Pstream::defaultComms(),
				distMap.schedule(),
				distMap.constructSize(),
				distMap.subMap(),
				distMap.constructMap(),
				allValues
			);

			newValues_ = this->patch().patchSlice(allValues);

			break;
		}
		default:
		{
			FatalErrorIn
			(
				"directMappedFixedValueFvPatchField<Type>::updateCoeffs()"
			)<< "Unknown sampling mode: " << mpp.mode()
			 << nl << abort(FatalError);
		}
	}

	if (setAverage_)
	{
		Type averagePsi =
			gSum(this->patch().magSf()*newValues_)
		   /gSum(this->patch().magSf());

		if (mag(averagePsi)/mag(average_) > 0.5)
		{
			newValues_ *= mag(average_)/mag(averagePsi);
		}
		else
		{
			newValues_ += (average_ - averagePsi);
		}
	}

	if (debug)
	{
		Info<< "directMapped on field:" << fldName
			<< " patch:" << this->patch().name()
			<< "  avg:" << gAverage(newValues_)
			<< "  min:" << gMin(newValues_)
			<< "  max:" << gMax(newValues_)
			<< endl;
	}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
directMappedFixedValueFvPatchField<Type>::directMappedFixedValueFvPatchField
(
	const fvPatch& p,
	const DimensionedField<Type, volMesh>& iF
)
:
	fixedValueFvPatchField<Type>(p, iF),
	setAverage_(false),
	average_(pTraits<Type>::zero),
	curTimeIndex_(-1)
{}


template<class Type>
directMappedFixedValueFvPatchField<Type>::directMappedFixedValueFvPatchField
(
	const fvPatch& p,
	const DimensionedField<Type, volMesh>& iF,
	const dictionary& dict
)
:
	fixedValueFvPatchField<Type>(p, iF, dict),
	setAverage_(readBool(dict.lookup("setAverage"))),
	average_(pTraits<Type>(dict.lookup("average"))),
	curTimeIndex_(-1)
{
	if (!isA<directMappedPatchBase>(this->patch().patch()))
	{
		FatalErrorIn
		(
			"directMappedFixedValueFvPatchField<Type>::"
			"directMappedFixedValueFvPatchField\n"
			"(\n"
			"    const fvPatch& p,\n"
			"    const DimensionedField<Type, volMesh>& iF,\n"
			"    const dictionary& dict\n"
			")\n"
		)   << "\n	patch type '" << p.type()
			<< "' not type '" << directMappedPatchBase::typeName << "'"
			<< "\n	for patch " << p.name()
			<< " of field " << this->dimensionedInternalField().name()
			<< " in file " << this->dimensionedInternalField().objectPath()
			<< exit(FatalError);
	}
}


template<class Type>
directMappedFixedValueFvPatchField<Type>::directMappedFixedValueFvPatchField
(
	const directMappedFixedValueFvPatchField<Type>& ptf,
	const fvPatch& p,
	const DimensionedField<Type, volMesh>& iF,
	const fvPatchFieldMapper& mapper
)
:
	fixedValueFvPatchField<Type>(ptf, p, iF, mapper),
	setAverage_(ptf.setAverage_),
	average_(ptf.average_),
	curTimeIndex_(-1)
{
	if (!isA<directMappedPatchBase>(this->patch().patch()))
	{
		FatalErrorIn
		(
			"directMappedFixedValueFvPatchField<Type>::"
			"directMappedFixedValueFvPatchField\n"
			"(\n"
			"    const directMappedFixedValueFvPatchField<Type>&,\n"
			"    const fvPatch&,\n"
			"    const Field<Type>&,\n"
			"    const fvPatchFieldMapper&\n"
			")\n"
		)   << "\n	patch type '" << p.type()
			<< "' not type '" << directMappedPatchBase::typeName << "'"
			<< "\n	for patch " << p.name()
			<< " of field " << this->dimensionedInternalField().name()
			<< " in file " << this->dimensionedInternalField().objectPath()
			<< exit(FatalError);
	}
}


template<class Type>
directMappedFixedValueFvPatchField<Type>::directMappedFixedValueFvPatchField
(
	const directMappedFixedValueFvPatchField<Type>& ptf
)
:
	fixedValueFvPatchField<Type>(ptf),
	setAverage_(ptf.setAverage_),
	average_(ptf.average_),
	curTimeIndex_(-1)
{}


template<class Type>
directMappedFixedValueFvPatchField<Type>::directMappedFixedValueFvPatchField
(
	const directMappedFixedValueFvPatchField<Type>& ptf,
	const DimensionedField<Type, volMesh>& iF
)
:
	fixedValueFvPatchField<Type>(ptf, iF),
	setAverage_(ptf.setAverage_),
	average_(ptf.average_),
	curTimeIndex_(-1)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void directMappedFixedValueFvPatchField<Type>::initEvaluate
(
	const Pstream::commsTypes commsType
)
{
	// Note: All parallel comms need to happen at this stage
	// HJ, 13/Mar/2012

	// Only map
	if (curTimeIndex_ != this->db().time().timeIndex())
	{
		mapField();

		curTimeIndex_ = this->db().time().timeIndex();
	}
}


template<class Type>
void directMappedFixedValueFvPatchField<Type>::evaluate
(
	const Pstream::commsTypes commsType
)
{
	this->operator==(newValues_);

	fixedValueFvPatchField<Type>::evaluate();
}


template<class Type>
void directMappedFixedValueFvPatchField<Type>::write(Ostream& os) const
{
	fvPatchField<Type>::write(os);
	os.writeKeyword("setAverage") << setAverage_ << token::END_STATEMENT << nl;
	os.writeKeyword("average") << average_ << token::END_STATEMENT << nl;
	this->writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
