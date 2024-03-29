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

#include "timeVaryingMappedFixedValueFvPatchField.H"
#include "foamTime.H"
#include "triSurfaceTools.H"
#include "triSurface.H"
#include "vector2D.H"
#include "OFstream.H"
#include "AverageIOField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
timeVaryingMappedFixedValueFvPatchField<Type>::
timeVaryingMappedFixedValueFvPatchField
(
	const fvPatch& p,
	const DimensionedField<Type, volMesh>& iF
)
:
	fixedValueFvPatchField<Type>(p, iF),
	fieldTableName_(iF.name()),
	setAverage_(false),
	referenceCS_(nullptr),
	nearestVertex_(0),
	nearestVertexWeight_(0),
	sampleTimes_(0),
	startSampleTime_(-1),
	startSampledValues_(0),
	startAverage_(pTraits<Type>::zero),
	endSampleTime_(-1),
	endSampledValues_(0),
	endAverage_(pTraits<Type>::zero)
{
	if (debug)
	{
		Pout<< "timeVaryingMappedFixedValue :"
			<< " construct from fvPatch and internalField" << endl;
	}
}


template<class Type>
timeVaryingMappedFixedValueFvPatchField<Type>::
timeVaryingMappedFixedValueFvPatchField
(
	const timeVaryingMappedFixedValueFvPatchField<Type>& ptf,
	const fvPatch& p,
	const DimensionedField<Type, volMesh>& iF,
	const fvPatchFieldMapper& mapper
)
:
	fixedValueFvPatchField<Type>(ptf, p, iF, mapper),
	fieldTableName_(ptf.fieldTableName_),
	setAverage_(ptf.setAverage_),
	referenceCS_(nullptr),
	nearestVertex_(0),
	nearestVertexWeight_(0),
	sampleTimes_(0),
	startSampleTime_(-1),
	startSampledValues_(0),
	startAverage_(pTraits<Type>::zero),
	endSampleTime_(-1),
	endSampledValues_(0),
	endAverage_(pTraits<Type>::zero)
{
	if (debug)
	{
		Pout<< "timeVaryingMappedFixedValue"
			<< " : construct from mappedFixedValue and mapper" << endl;
	}
}


template<class Type>
timeVaryingMappedFixedValueFvPatchField<Type>::
timeVaryingMappedFixedValueFvPatchField
(
	const fvPatch& p,
	const DimensionedField<Type, volMesh>& iF,
	const dictionary& dict
)
:
	fixedValueFvPatchField<Type>(p, iF),
	fieldTableName_(iF.name()),
	setAverage_(readBool(dict.lookup("setAverage"))),
	referenceCS_(nullptr),
	nearestVertex_(0),
	nearestVertexWeight_(0),
	sampleTimes_(0),
	startSampleTime_(-1),
	startSampledValues_(0),
	startAverage_(pTraits<Type>::zero),
	endSampleTime_(-1),
	endSampledValues_(0),
	endAverage_(pTraits<Type>::zero)
{
	if (debug)
	{
		Pout<< "timeVaryingMappedFixedValue : construct from dictionary"
			<< endl;
	}

	if (dict.found("fieldTableName"))
	{
		dict.lookup("fieldTableName") >> fieldTableName_;
	}

	if (dict.found("value"))
	{
		fvPatchField<Type>::operator==(Field<Type>("value", dict, p.size()));
	}
	else
	{
		updateCoeffs();
	}
}


template<class Type>
timeVaryingMappedFixedValueFvPatchField<Type>::
timeVaryingMappedFixedValueFvPatchField
(
	const timeVaryingMappedFixedValueFvPatchField<Type>& ptf
)
:
	fixedValueFvPatchField<Type>(ptf),
	fieldTableName_(ptf.fieldTableName_),
	setAverage_(ptf.setAverage_),
	referenceCS_(ptf.referenceCS_),
	nearestVertex_(ptf.nearestVertex_),
	nearestVertexWeight_(ptf.nearestVertexWeight_),
	sampleTimes_(ptf.sampleTimes_),
	startSampleTime_(ptf.startSampleTime_),
	startSampledValues_(ptf.startSampledValues_),
	startAverage_(ptf.startAverage_),
	endSampleTime_(ptf.endSampleTime_),
	endSampledValues_(ptf.endSampledValues_),
	endAverage_(ptf.endAverage_)
{
	if (debug)
	{
		Pout<< "timeVaryingMappedFixedValue : copy construct"
			<< endl;
	}
}


template<class Type>
timeVaryingMappedFixedValueFvPatchField<Type>::
timeVaryingMappedFixedValueFvPatchField
(
	const timeVaryingMappedFixedValueFvPatchField<Type>& ptf,
	const DimensionedField<Type, volMesh>& iF
)
:
	fixedValueFvPatchField<Type>(ptf, iF),
	fieldTableName_(ptf.fieldTableName_),
	setAverage_(ptf.setAverage_),
	referenceCS_(ptf.referenceCS_),
	nearestVertex_(ptf.nearestVertex_),
	nearestVertexWeight_(ptf.nearestVertexWeight_),
	sampleTimes_(ptf.sampleTimes_),
	startSampleTime_(ptf.startSampleTime_),
	startSampledValues_(ptf.startSampledValues_),
	startAverage_(ptf.startAverage_),
	endSampleTime_(ptf.endSampleTime_),
	endSampledValues_(ptf.endSampledValues_),
	endAverage_(ptf.endAverage_)
{
	if (debug)
	{
		Pout<< "timeVaryingMappedFixedValue :"
			<< " copy construct, resetting internal field" << endl;
	}
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void timeVaryingMappedFixedValueFvPatchField<Type>::autoMap
(
	const fvPatchFieldMapper& m
)
{
	fixedValueFvPatchField<Type>::autoMap(m);
	if (startSampledValues_.size() > 0)
	{
		startSampledValues_.autoMap(m);
		endSampledValues_.autoMap(m);
	}
}


template<class Type>
void timeVaryingMappedFixedValueFvPatchField<Type>::rmap
(
	const fvPatchField<Type>& ptf,
	const labelList& addr
)
{
	fixedValueFvPatchField<Type>::rmap(ptf, addr);

	const timeVaryingMappedFixedValueFvPatchField<Type>& tiptf =
		refCast<const timeVaryingMappedFixedValueFvPatchField<Type> >(ptf);

	startSampledValues_.rmap(tiptf.startSampledValues_, addr);
	endSampledValues_.rmap(tiptf.endSampledValues_, addr);
}


template<class Type>
void timeVaryingMappedFixedValueFvPatchField<Type>::readSamplePoints()
{
	// Read the sample points

	pointIOField samplePoints
	(
		IOobject
		(
			"points",
			this->db().time().constant(),
			"boundaryData"/this->patch().name(),
			this->db(),
			IOobject::MUST_READ,
			IOobject::AUTO_WRITE,
			false
		)
	);

	const fileName samplePointsFile = samplePoints.filePath();

	if (debug)
	{
		Info<< "timeVaryingMappedFixedValueFvPatchField :"
			<< " Read " << samplePoints.size() << " sample points from "
			<< samplePointsFile << endl;
	}

	// Determine coordinate system from samplePoints

	if (samplePoints.size() < 3)
	{
		FatalErrorIn
		(
			"timeVaryingMappedFixedValueFvPatchField<Type>::readSamplePoints()"
		)   << "Only " << samplePoints.size() << " points read from file "
			<< samplePoints.objectPath() << nl
			<< "Need at least three non-colinear samplePoints"
			<< " to be able to interpolate."
			<< "\n    on patch " << this->patch().name()
			<< " of field " << this->dimensionedInternalField().name()
			<< " in file " << this->dimensionedInternalField().objectPath()
			<< exit(FatalError);
	}

	const point& p0 = samplePoints[0];

	// Find point separate from p0
	vector e1;
	label index1 = -1;

	for (label i = 1; i < samplePoints.size(); i++)
	{
		e1 = samplePoints[i] - p0;

		scalar magE1 = mag(e1);

		if (magE1 > SMALL)
		{
			e1 /= magE1;
			index1 = i;
			break;
		}
	}

	// Find point that makes angle with n1
	label index2 = -1;
	vector e2;
	vector n;

	for (label i = index1+1; i < samplePoints.size(); i++)
	{
		e2 = samplePoints[i] - p0;

		scalar magE2 = mag(e2);

		if (magE2 > SMALL)
		{
			e2 /= magE2;

			n = e1^e2;

			scalar magN = mag(n);

			if (magN > SMALL)
			{
				index2 = i;
				n /= magN;
				break;
			}
		}
	}

	if (index2 == -1)
	{
		FatalErrorIn
		(
			"timeVaryingMappedFixedValueFvPatchField<Type>::readSamplePoints()"
		)   << "Cannot find points that make valid normal." << nl
			<< "Need at least three sample points which are not in a line."
			<< "\n    on patch " << this->patch().name()
			<< " of points " << samplePoints.name()
			<< " in file " << samplePoints.objectPath()
			<< exit(FatalError);
	}

	if (debug)
	{
		Info<< "timeVaryingMappedFixedValueFvPatchField :"
			<< " Used points " << p0 << ' ' << samplePoints[index1]
			<< ' ' << samplePoints[index2]
			<< " to define coordinate system with normal " << n << endl;
	}

	referenceCS_.reset
	(
		new coordinateSystem
		(
			"reference",
			p0,  // origin
			n,   // normal
			e1   // 0-axis
		)
	);

	tmp<vectorField> tlocalVertices
	(
		referenceCS().localPosition(samplePoints)
	);
	const vectorField& localVertices = tlocalVertices();

	// Determine triangulation
	List<vector2D> localVertices2D(localVertices.size());
	forAll(localVertices, i)
	{
		localVertices2D[i][0] = localVertices[i][0];
		localVertices2D[i][1] = localVertices[i][1];
	}

	triSurface s(triSurfaceTools::delaunay2D(localVertices2D));

	tmp<pointField> localFaceCentres
	(
		referenceCS().localPosition
		(
			this->patch().patch().faceCentres()
		)
	);

	if (debug)
	{
		Pout<< "readSamplePoints :"
			<<" Dumping triangulated surface to triangulation.stl" << endl;
		s.write(this->db().time().path()/"triangulation.stl");

		OFstream str(this->db().time().path()/"localFaceCentres.obj");
		Pout<< "readSamplePoints :"
			<< " Dumping face centres to " << str.name() << endl;

		forAll(localFaceCentres(), i)
		{
			const point& p = localFaceCentres()[i];
			str<< "v " << p.x() << ' ' << p.y() << ' ' << p.z() << nl;
		}
	}

	// Determine interpolation onto face centres.
	triSurfaceTools::calcInterpolationWeights
	(
		s,
		localFaceCentres,   // points to interpolate to
		nearestVertex_,
		nearestVertexWeight_
	);



	// Read the times for which data is available

	const fileName samplePointsDir = samplePointsFile.path();

	sampleTimes_ = Time::findTimes(samplePointsDir);

	if (debug)
	{
		Info<< "timeVaryingMappedFixedValueFvPatchField : In directory "
			<< samplePointsDir << " found times " << timeNames(sampleTimes_)
			<< endl;
	}
}


template<class Type>
wordList timeVaryingMappedFixedValueFvPatchField<Type>::timeNames
(
	const instantList& times
)
{
	wordList names(times.size());

	forAll(times, i)
	{
		names[i] = times[i].name();
	}
	return names;
}


template<class Type>
void timeVaryingMappedFixedValueFvPatchField<Type>::findTime
(
	const fileName& instance,
	const fileName& local,
	const scalar timeVal,
	label& lo,
	label& hi
) const
{
	lo = startSampleTime_;
	hi = -1;

	for (label i = startSampleTime_+1; i < sampleTimes_.size(); i++)
	{
		if (sampleTimes_[i].value() > timeVal)
		{
			break;
		}
		else
		{
			lo = i;
		}
	}

	if (lo == -1)
	{
		FatalErrorIn("findTime")
			<< "Cannot find starting sampling values for current time "
			<< timeVal << nl
			<< "Have sampling values for times "
			<< timeNames(sampleTimes_) << nl
			<< "In directory "
			<<  this->db().time().constant()/"boundaryData"/this->patch().name()
			<< "\n    on patch " << this->patch().name()
			<< " of field " << fieldTableName_
			<< exit(FatalError);
	}

	if (lo < sampleTimes_.size()-1)
	{
		hi = lo+1;
	}


	if (debug)
	{
		if (hi == -1)
		{
			Pout<< "findTime : Found time " << timeVal << " after"
				<< " index:" << lo << " time:" << sampleTimes_[lo].value()
				<< endl;
		}
		else
		{
			Pout<< "findTime : Found time " << timeVal << " inbetween"
				<< " index:" << lo << " time:" << sampleTimes_[lo].value()
				<< " and index:" << hi << " time:" << sampleTimes_[hi].value()
				<< endl;
		}
	}
}


template<class Type>
void timeVaryingMappedFixedValueFvPatchField<Type>::checkTable()
{
	// Initialise
	if (startSampleTime_ == -1 && endSampleTime_ == -1)
	{
		readSamplePoints();
	}

	// Find current time in sampleTimes
	label lo = -1;
	label hi = -1;

	findTime
	(
		this->db().time().constant(),
		"boundaryData"/this->patch().name(),
		this->db().time().value(),
		lo,
		hi
	);

	// Update sampled data fields.

	if (lo != startSampleTime_)
	{
		startSampleTime_ = lo;

		if (startSampleTime_ == endSampleTime_)
		{
			// No need to reread since are end values
			if (debug)
			{
				Pout<< "checkTable : Setting startValues to (already read) "
					<<   "boundaryData"
					    /this->patch().name()
					    /sampleTimes_[startSampleTime_].name()
					<< endl;
			}
			startSampledValues_ = endSampledValues_;
			startAverage_ = endAverage_;
		}
		else
		{
			if (debug)
			{
				Pout<< "checkTable : Reading startValues from "
					<<   "boundaryData"
					    /this->patch().name()
					    /sampleTimes_[lo].name()
					<< endl;
			}


			// Reread values and interpolate
			AverageIOField<Type> vals
			(
				IOobject
				(
					fieldTableName_,
					this->db().time().constant(),
					"boundaryData"
				   /this->patch().name()
				   /sampleTimes_[startSampleTime_].name(),
					this->db(),
					IOobject::MUST_READ,
					IOobject::AUTO_WRITE,
					false
				)
			);

			startAverage_ = vals.average();
			startSampledValues_ = interpolate(vals);
		}
	}

	if (hi != endSampleTime_)
	{
		endSampleTime_ = hi;

		if (endSampleTime_ == -1)
		{
			// endTime no longer valid. Might as well clear endValues.
			if (debug)
			{
				Pout<< "checkTable : Clearing endValues" << endl;
			}
			endSampledValues_.clear();
		}
		else
		{
			if (debug)
			{
				Pout<< "checkTable : Reading endValues from "
					<<   "boundaryData"
					    /this->patch().name()
					    /sampleTimes_[endSampleTime_].name()
					<< endl;
			}
			// Reread values and interpolate
			AverageIOField<Type> vals
			(
				IOobject
				(
					fieldTableName_,
					this->db().time().constant(),
					"boundaryData"
				   /this->patch().name()
				   /sampleTimes_[endSampleTime_].name(),
					this->db(),
					IOobject::MUST_READ,
					IOobject::AUTO_WRITE,
					false
				)
			);
			endAverage_ = vals.average();
			endSampledValues_ = interpolate(vals);
		}
	}
}


template<class Type>
tmp<Field<Type> > timeVaryingMappedFixedValueFvPatchField<Type>::interpolate
(
	const Field<Type>& sourceFld
) const
{
	tmp<Field<Type> > tfld(new Field<Type>(nearestVertex_.size()));
	Field<Type>& fld = tfld();

	forAll(fld, i)
	{
		const FixedList<label, 3>& verts = nearestVertex_[i];
		const FixedList<scalar, 3>& w = nearestVertexWeight_[i];

		if (verts[2] == -1)
		{
			if (verts[1] == -1)
			{
				// Use vertex0 only
				fld[i] = sourceFld[verts[0]];
			}
			else
			{
				// Use vertex 0,1
				fld[i] =
					w[0]*sourceFld[verts[0]]
				  + w[1]*sourceFld[verts[1]];
			}
		}
		else
		{
			fld[i] =
				w[0]*sourceFld[verts[0]]
			  + w[1]*sourceFld[verts[1]]
			  + w[2]*sourceFld[verts[2]];
		}
	}
	return tfld;
}


template<class Type>
void timeVaryingMappedFixedValueFvPatchField<Type>::updateCoeffs()
{
	if (this->updated())
	{
		return;
	}

	checkTable();

	// Interpolate between the sampled data

	Type wantedAverage;

	if (endSampleTime_ == -1)
	{
		// only start value
		if (debug)
		{
			Pout<< "updateCoeffs : Sampled, non-interpolated values"
				<< " from start time:"
				<< sampleTimes_[startSampleTime_].name() << nl;
		}

		this->operator==(startSampledValues_);
		wantedAverage = startAverage_;
	}
	else
	{
		scalar start = sampleTimes_[startSampleTime_].value();
		scalar end = sampleTimes_[endSampleTime_].value();

		scalar s = (this->db().time().value()-start)/(end-start);

		if (debug)
		{
			Pout<< "updateCoeffs : Sampled, interpolated values"
				<< " between start time:"
				<< sampleTimes_[startSampleTime_].name()
				<< " and end time:" << sampleTimes_[endSampleTime_].name()
				<< " with weight:" << s << endl;
		}

		this->operator==((1-s)*startSampledValues_ + s*endSampledValues_);
		wantedAverage = (1-s)*startAverage_ + s*endAverage_;
	}

	// Enforce average. Either by scaling (if scaling factor > 0.5) or by
	// offsetting.
	if (setAverage_)
	{
		const Field<Type>& fld = *this;

		Type averagePsi =
			gSum(this->patch().magSf()*fld)
		   /gSum(this->patch().magSf());

		if (debug)
		{
			Pout<< "updateCoeffs :"
				<< " actual average:" << averagePsi
				<< " wanted average:" << wantedAverage
				<< endl;
		}

		if (mag(averagePsi) < VSMALL)
		{
			// Field too small to scale. Offset instead.
			const Type offset = wantedAverage - averagePsi;
			if (debug)
			{
				Pout<< "updateCoeffs :"
					<< " offsetting with:" << offset << endl;
			}
			this->operator==(fld+offset);
		}
		else
		{
			const scalar scale = mag(wantedAverage)/mag(averagePsi);

			if (debug)
			{
				Pout<< "updateCoeffs :"
					<< " scaling with:" << scale << endl;
			}
			this->operator==(scale*fld);
		}
	}

	if (debug)
	{
		Pout<< "updateCoeffs : set fixedValue to min:" << gMin(*this)
			<< " max:" << gMax(*this) << endl;
	}

	fixedValueFvPatchField<Type>::updateCoeffs();
}


template<class Type>
void timeVaryingMappedFixedValueFvPatchField<Type>::write(Ostream& os) const
{
	fvPatchField<Type>::write(os);
	os.writeKeyword("setAverage") << setAverage_ << token::END_STATEMENT << nl;

	if (fieldTableName_ != this->dimensionedInternalField().name())
	{
		os.writeKeyword("fieldTableName") << fieldTableName_ << token::END_STATEMENT << nl;
	}

	this->writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
