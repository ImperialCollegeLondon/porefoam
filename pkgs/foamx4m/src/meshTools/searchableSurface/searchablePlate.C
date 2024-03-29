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

#include "searchablePlate.H"
#include "addToRunTimeSelectionTable.H"
#include "SortableList.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(searchablePlate, 0);
addToRunTimeSelectionTable(searchableSurface, searchablePlate, dict);

}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::direction Foam::searchablePlate::calcNormal(const point& span)
{
	direction normalDir = 3;

	for (direction dir = 0; dir < vector::nComponents; dir++)
	{
		if (span[dir] < 0)
		{
			FatalErrorIn("searchablePlate::calcNormal()")
				<< "Span should have two positive and one zero entry. Now:"
				<< span << exit(FatalError);
		}
		else if (span[dir] < VSMALL)
		{
			if (normalDir == 3)
			{
				normalDir = dir;
			}
			else
			{
				// Multiple zero entries. Flag and exit.
				normalDir = 3;
				break;
			}
		}
	}

	if (normalDir == 3)
	{
		FatalErrorIn("searchablePlate::calcNormal()")
			<< "Span should have one and only zero entry. Now:" << span
			<< exit(FatalError);
	}

	return normalDir;
}


// Returns miss or hit with face (always 0)
Foam::pointIndexHit Foam::searchablePlate::findNearest
(
	const point& sample,
	const scalar nearestDistSqr
) const
{
	// For every component direction can be
	// left of min, right of max or inbetween.
	// - outside points: project first one x plane (either min().x()
	// or max().x()), then onto y plane and finally z. You should be left
	// with intersection point
	// - inside point: find nearest side (compare to mid point). Project onto
	//   that.

	// Project point on plane.
	pointIndexHit info(true, sample, 0);
	info.rawPoint()[normalDir_] = origin_[normalDir_];

	// Clip to edges if outside
	for (direction dir = 0; dir < vector::nComponents; dir++)
	{
		if (dir != normalDir_)
		{
			if (info.rawPoint()[dir] < origin_[dir])
			{
				info.rawPoint()[dir] = origin_[dir];
			}
			else if (info.rawPoint()[dir] > origin_[dir]+span_[dir])
			{
				info.rawPoint()[dir] = origin_[dir]+span_[dir];
			}
		}
	}

	// Check if outside. Optimisation: could do some checks on distance already
	// on components above
	if (magSqr(info.rawPoint() - sample) > nearestDistSqr)
	{
		info.setMiss();
		info.setIndex(-1);
	}

	return info;
}


Foam::pointIndexHit Foam::searchablePlate::findLine
(
	const point& start,
	const point& end
) const
{
	pointIndexHit info
	(
		true,
		vector::zero,
		0
	);

	const vector dir(end-start);

	if (mag(dir[normalDir_]) < VSMALL)
	{
		info.setMiss();
		info.setIndex(-1);
	}
	else
	{
		scalar t = (origin_[normalDir_]-start[normalDir_]) / dir[normalDir_];

		if (t < 0 || t > 1)
		{
			info.setMiss();
			info.setIndex(-1);
		}
		else
		{
			info.rawPoint() = start+t*dir;
			info.rawPoint()[normalDir_] = origin_[normalDir_];

			// Clip to edges
			for (direction dir = 0; dir < vector::nComponents; dir++)
			{
				if (dir != normalDir_)
				{
					if (info.rawPoint()[dir] < origin_[dir])
					{
					    info.setMiss();
					    info.setIndex(-1);
					    break;
					}
					else if (info.rawPoint()[dir] > origin_[dir]+span_[dir])
					{
					    info.setMiss();
					    info.setIndex(-1);
					    break;
					}
				}
			}
		}
	}

	// Debug
	if (info.hit())
	{
		treeBoundBox bb(origin_, origin_+span_);
		bb.min()[normalDir_] -= 1E-6;
		bb.max()[normalDir_] += 1E-6;

		if (!bb.contains(info.hitPoint()))
		{
			FatalErrorIn("searchablePlate::findLine(..)")
				<< "bb:" << bb << endl
				<< "origin_:" << origin_ << endl
				<< "span_:" << span_ << endl
				<< "normalDir_:" << normalDir_ << endl
				<< "hitPoint:" << info.hitPoint()
				<< abort(FatalError);
		}
	}

	return info;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::searchablePlate::searchablePlate
(
	const IOobject& io,
	const point& origin,
	const vector& span
)
:
	searchableSurface(io),
	origin_(origin),
	span_(span),
	normalDir_(calcNormal(span_))
{
	if (debug)
	{
		Info<< "searchablePlate::searchablePlate :"
			<< " origin:" << origin_
			<< " origin+span:" << origin_+span_
			<< " normal:" << vector::componentNames[normalDir_]
			<< endl;
	}
}


Foam::searchablePlate::searchablePlate
(
	const IOobject& io,
	const dictionary& dict
)
:
	searchableSurface(io),
	origin_(dict.lookup("origin")),
	span_(dict.lookup("span")),
	normalDir_(calcNormal(span_))
{
	if (debug)
	{
		Info<< "searchablePlate::searchablePlate :"
			<< " origin:" << origin_
			<< " origin+span:" << origin_+span_
			<< " normal:" << vector::componentNames[normalDir_]
			<< endl;
	}
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::searchablePlate::~searchablePlate()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::wordList& Foam::searchablePlate::regions() const
{
	if (regions_.empty())
	{
		regions_.setSize(1);
		regions_[0] = "region0";
	}
	return regions_;
}


void Foam::searchablePlate::findNearest
(
	const pointField& samples,
	const scalarField& nearestDistSqr,
	List<pointIndexHit>& info
) const
{
	info.setSize(samples.size());

	forAll(samples, i)
	{
		info[i] = findNearest(samples[i], nearestDistSqr[i]);
	}
}


void Foam::searchablePlate::findLine
(
	const pointField& start,
	const pointField& end,
	List<pointIndexHit>& info
) const
{
	info.setSize(start.size());

	forAll(start, i)
	{
		info[i] = findLine(start[i], end[i]);
	}
}


void Foam::searchablePlate::findLineAny
(
	const pointField& start,
	const pointField& end,
	List<pointIndexHit>& info
) const
{
	findLine(start, end, info);
}


void Foam::searchablePlate::findLineAll
(
	const pointField& start,
	const pointField& end,
	List<List<pointIndexHit> >& info
) const
{
	List<pointIndexHit> nearestInfo;
	findLine(start, end, nearestInfo);

	info.setSize(start.size());
	forAll(info, pointI)
	{
		if (nearestInfo[pointI].hit())
		{
			info[pointI].setSize(1);
			info[pointI][0] = nearestInfo[pointI];
		}
		else
		{
			info[pointI].clear();
		}
	}
}


void Foam::searchablePlate::getRegion
(
	const List<pointIndexHit>& info,
	labelList& region
) const
{
	region.setSize(info.size());
	region = 0;
}


void Foam::searchablePlate::getNormal
(
	const List<pointIndexHit>& info,
	vectorField& normal
) const
{
	normal.setSize(info.size());
	normal = vector::zero;
	forAll(normal, i)
	{
		normal[i][normalDir_] = 1.0;
	}
}


void Foam::searchablePlate::getVolumeType
(
	const pointField& points,
	List<volumeType>& volType
) const
{
	FatalErrorIn
	(
		"searchableCollection::getVolumeType(const pointField&"
		", List<volumeType>&) const"
	)   << "Volume type not supported for plate."
		<< exit(FatalError);
}


// ************************************************************************* //
