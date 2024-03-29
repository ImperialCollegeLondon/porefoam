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

#include "angularOscillatingDisplacementPointPatchVectorField.H"
#include "pointPatchFields.H"
#include "addToRunTimeSelectionTable.H"
#include "foamTime.H"
#include "polyMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

angularOscillatingDisplacementPointPatchVectorField::
angularOscillatingDisplacementPointPatchVectorField
(
	const pointPatch& p,
	const DimensionedField<vector, pointMesh>& iF
)
:
	fixedValuePointPatchVectorField(p, iF),
	axis_(vector::zero),
	origin_(vector::zero),
	angle0_(0.0),
	amplitude_(0.0),
	omega_(0.0),
	p0_(p.localPoints())
{}


angularOscillatingDisplacementPointPatchVectorField::
angularOscillatingDisplacementPointPatchVectorField
(
	const pointPatch& p,
	const DimensionedField<vector, pointMesh>& iF,
	const dictionary& dict
)
:
	fixedValuePointPatchVectorField(p, iF, dict),
	axis_(dict.lookup("axis")),
	origin_(dict.lookup("origin")),
	angle0_(readScalar(dict.lookup("angle0"))),
	amplitude_(readScalar(dict.lookup("amplitude"))),
	omega_(readScalar(dict.lookup("omega")))
{
	if (!dict.found("value"))
	{
		updateCoeffs();
	}

	if (dict.found("p0"))
	{
		p0_ = vectorField("p0", dict , p.size());
	}
	else
	{
		p0_ = p.localPoints();
	}
}


angularOscillatingDisplacementPointPatchVectorField::
angularOscillatingDisplacementPointPatchVectorField
(
	const angularOscillatingDisplacementPointPatchVectorField& ptf,
	const pointPatch& p,
	const DimensionedField<vector, pointMesh>& iF,
	const PointPatchFieldMapper& mapper
)
:
	fixedValuePointPatchVectorField(ptf, p, iF, mapper),
	axis_(ptf.axis_),
	origin_(ptf.origin_),
	angle0_(ptf.angle0_),
	amplitude_(ptf.amplitude_),
	omega_(ptf.omega_),
	p0_(ptf.p0_, mapper)
{}


angularOscillatingDisplacementPointPatchVectorField::
angularOscillatingDisplacementPointPatchVectorField
(
	const angularOscillatingDisplacementPointPatchVectorField& ptf,
	const DimensionedField<vector, pointMesh>& iF
)
:
	fixedValuePointPatchVectorField(ptf, iF),
	axis_(ptf.axis_),
	origin_(ptf.origin_),
	angle0_(ptf.angle0_),
	amplitude_(ptf.amplitude_),
	omega_(ptf.omega_),
	p0_(ptf.p0_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void angularOscillatingDisplacementPointPatchVectorField::autoMap
(
	const PointPatchFieldMapper& m
)
{
	fixedValuePointPatchVectorField::autoMap(m);
	p0_.autoMap(m);
}

void angularOscillatingDisplacementPointPatchVectorField::rmap
(
	const PointPatchField
		<pointPatchField, pointMesh, pointPatch, DummyMatrix, vector>& ptf,
	const labelList& addr
)
{
	const angularOscillatingDisplacementPointPatchVectorField& aODptf =
			refCast<const angularOscillatingDisplacementPointPatchVectorField>(ptf);

	fixedValuePointPatchVectorField::rmap(aODptf, addr);
	p0_.rmap(aODptf.p0_, addr);
}

void angularOscillatingDisplacementPointPatchVectorField::updateCoeffs()
{
	if (this->updated())
	{
		return;
	}

	const polyMesh& mesh = this->dimensionedInternalField().mesh()();
	const Time& t = mesh.time();

	scalar angle = angle0_ + amplitude_*sin(omega_*t.value());
	vector axisHat = axis_/mag(axis_);
	vectorField p0Rel = p0_ - origin_;

	vectorField::operator=
	(
		p0Rel*(cos(angle) - 1)
	  + (axisHat ^ p0Rel*sin(angle))
	  + (axisHat & p0Rel)*(1 - cos(angle))*axisHat
	);

	fixedValuePointPatchVectorField::updateCoeffs();
}


void angularOscillatingDisplacementPointPatchVectorField::write
(
	Ostream& os
) const
{
	pointPatchVectorField::write(os);

	os.writeKeyword("axis")
		<< axis_ << token::END_STATEMENT << nl;
	os.writeKeyword("origin")
		<< origin_ << token::END_STATEMENT << nl;
	os.writeKeyword("angle0")
		<< angle0_ << token::END_STATEMENT << nl;
	os.writeKeyword("amplitude")
		<< amplitude_ << token::END_STATEMENT << nl;
	os.writeKeyword("omega")
		<< omega_ << token::END_STATEMENT << nl;
	p0_.writeEntry("p0", os);
	writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePointPatchTypeField
(
	pointPatchVectorField,
	angularOscillatingDisplacementPointPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
