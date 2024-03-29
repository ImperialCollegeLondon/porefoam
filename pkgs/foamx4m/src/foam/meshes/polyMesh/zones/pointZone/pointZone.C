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

Description
	A subset of mesh points.

\*---------------------------------------------------------------------------*/

#include "pointZone.H"
#include "addToRunTimeSelectionTable.H"
#include "pointZoneMesh.H"
#include "polyMesh.H"
#include "primitiveMesh.H"
#include "demandDrivenData.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
	defineTypeNameAndDebug(pointZone, 0);
	defineRunTimeSelectionTable(pointZone, dictionary);
	addToRunTimeSelectionTable(pointZone, pointZone, dictionary);
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

const Foam::Map<Foam::label>& Foam::pointZone::pointLookupMap() const
{
	if (!pointLookupMapPtr_)
	{
		calcPointLookupMap();
	}

	return *pointLookupMapPtr_;
}


void Foam::pointZone::calcPointLookupMap() const
{
	if (debug)
	{
		Info<< "void pointZone::calcPointLookupMap() const : "
			<< "Calculating point lookup map"
			<< endl;
	}

	if (pointLookupMapPtr_)
	{
		FatalErrorIn
		(
			"void pointZone::calcPointLookupMap() const"
		)   << "point lookup map already calculated"
			<< abort(FatalError);
	}

	const labelList& addr = *this;

	pointLookupMapPtr_ = new Map<label>(2*addr.size());
	Map<label>& plm = *pointLookupMapPtr_;

	forAll (addr, pointI)
	{
		plm.insert(addr[pointI], pointI);
	}

	if (debug)
	{
		Info<< "void pointZone::calcPointLookupMap() const : "
			<< "Finished calculating point lookup map"
			<< endl;
	}
}


void Foam::pointZone::clearAddressing()
{
	deleteDemandDrivenData(pointLookupMapPtr_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::pointZone::pointZone
(
	const word& name,
	const labelList& addr,
	const label index,
	const pointZoneMesh& zm
)
:
	labelList(addr),
	name_(name),
	index_(index),
	zoneMesh_(zm),
	pointLookupMapPtr_(nullptr)
{}


Foam::pointZone::pointZone
(
	const word& name,
	const Xfer<labelList>& addr,
	const label index,
	const pointZoneMesh& zm
)
:
	labelList(addr),
	name_(name),
	index_(index),
	zoneMesh_(zm),
	pointLookupMapPtr_(nullptr)
{}


// Construct from dictionary
Foam::pointZone::pointZone
(
	const word& name,
	const dictionary& dict,
	const label index,
	const pointZoneMesh& zm
)
:
	labelList(dict.lookup("pointLabels")),
	name_(name),
	index_(index),
	zoneMesh_(zm),
	pointLookupMapPtr_(nullptr)
{}


// Construct given the original zone and resetting the
// point list and zone mesh information
Foam::pointZone::pointZone
(
	const pointZone& pz,
	const labelList& addr,
	const label index,
	const pointZoneMesh& zm
)
:
	labelList(addr),
	name_(pz.name()),
	index_(index),
	zoneMesh_(zm),
	pointLookupMapPtr_(nullptr)
{}


Foam::pointZone::pointZone
(
	const pointZone& pz,
	const Xfer<labelList>& addr,
	const label index,
	const pointZoneMesh& zm
)
:
	labelList(addr),
	name_(pz.name()),
	index_(index),
	zoneMesh_(zm),
	pointLookupMapPtr_(nullptr)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::pointZone::~pointZone()
{
	clearAddressing();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::pointZone::whichPoint(const label globalPointID) const
{
	const Map<label>& plm = pointLookupMap();

	Map<label>::const_iterator plmIter = plm.find(globalPointID);

	if (plmIter == plm.end())
	{
		return -1;
	}
	else
	{
		return plmIter();
	}
}


const Foam::pointZoneMesh& Foam::pointZone::zoneMesh() const
{
	return zoneMesh_;
}


void Foam::pointZone::updateMesh()
{
	clearAddressing();
}


bool Foam::pointZone::checkDefinition(const bool report) const
{
	const labelList& addr = *this;

	bool boundaryError = false;

	forAll(addr, i)
	{
		if (addr[i] < 0 || addr[i] >= zoneMesh_.mesh().allPoints().size())
		{
			boundaryError = true;

			if (report)
			{
				SeriousErrorIn
				(
					"bool pointZone::checkDefinition("
					"const bool report) const"
				)   << "Zone " << name()
					<< " contains invalid point label " << addr[i] << nl
					<< "Valid point labels are 0.."
					<< zoneMesh_.mesh().allPoints().size() - 1 << endl;
			}
		}
	}
	return boundaryError;
}


void Foam::pointZone::write(Ostream& os) const
{
	os  << nl << name()
		<< nl << static_cast<const labelList&>(*this);
}


void Foam::pointZone::writeDict(Ostream& os) const
{
	os  << nl << name() << nl << token::BEGIN_BLOCK << incrIndent << nl
		<< indent << "type " << type() << token::END_STATEMENT << nl;

	writeEntry("pointLabels", os);

	os  << decrIndent << token::END_BLOCK << endl;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void Foam::pointZone::operator=(const pointZone& cz)
{
	clearAddressing();
	labelList::operator=(cz);
}


void Foam::pointZone::operator=(const labelList& addr)
{
	clearAddressing();
	labelList::operator=(addr);
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const pointZone& p)
{
	p.write(os);
	os.check("Ostream& operator<<(Ostream& f, const pointZone& p");
	return os;
}


// ************************************************************************* //
