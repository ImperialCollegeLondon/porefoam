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

\*---------------------------------------------------------------------------*/

#include "graph.H"
#include "OFstream.H"
#include "IOmanip.H"
#include "Pair.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(graph::writer, 0);
defineRunTimeSelectionTable(graph::writer, word);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void graph::readCurves(Istream& is)
{
	List<xy> xyData(is);

	x_.setSize(xyData.size());
	scalarField y(xyData.size());

	forAll (xyData, i)
	{
		x_[i] = xyData[i].x_;
		y[i] = xyData[i].y_;
	}

	insert(yName_, new curve(yName_, curve::curveStyle::CONTINUOUS, y));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

graph::graph
(
	const string& title,
	const string& xName,
	const string& yName,
	const scalarField& x
)
:
	title_(title),
	xName_(xName),
	yName_(yName),
	x_(x)
{}


graph::graph
(
	const string& title,
	const string& xName,
	const string& yName,
	const scalarField& x,
	const scalarField& y
)
:
	title_(title),
	xName_(xName),
	yName_(yName),
	x_(x)
{
	insert(yName, new curve(yName, curve::curveStyle::CONTINUOUS, y));
}


graph::graph
(
	const string& title,
	const string& xName,
	const string& yName,
	Istream& is
)
:
	title_(title),
	xName_(xName),
	yName_(yName)
{
	readCurves(is);
}


graph::graph(Istream& is)
:
	title_(is),
	xName_(is),
	yName_(is)
{
	readCurves(is);
}


const scalarField& graph::y() const
{
	if (size() != 1)
	{
		FatalErrorIn("const scalarField& graph::y() const")
			<< "y field requested for graph containing " << size()
			<< "ys" << exit(FatalError);
	}

	return *begin()();
}

scalarField& graph::y()
{
	if (size() != 1)
	{
		FatalErrorIn("scalarField& graph::y()")
			<< "y field requested for graph containing " << size()
			<< "ys" << exit(FatalError);
	}

	return *begin()();
}


autoPtr<graph::writer> graph::writer::New(const word& graphFormat)
{
	if (!wordConstructorTablePtr_)
	{
		FatalErrorIn
		(
			"graph::writer::New(const word&)"
		)   << "Graph writer table is empty"
			<< exit(FatalError);
	}

	wordConstructorTable::iterator cstrIter =
		wordConstructorTablePtr_->find(graphFormat);

	if (cstrIter == wordConstructorTablePtr_->end())
	{
		FatalErrorIn
		(
			"graph::writer::New(const word&)"
		)   << "Unknown graph format " << graphFormat
			<< endl << endl
			<< "Valid graph formats are : " << endl
			<< wordConstructorTablePtr_->sortedToc()
			<< exit(FatalError);
	}

	return autoPtr<graph::writer>(cstrIter()());
}


void graph::writer::writeXY
(
	const scalarField& x,
	const scalarField& y,
	Ostream& os
) const
{
	forAll(x, xi)
	{
		os << setw(10) << x[xi] << token::SPACE << setw(10) << y[xi]<< endl;
	}
}


void graph::writeTable(Ostream& os) const
{
	forAll(x_, xi)
	{
		os << setw(10) << x_[xi];

		for
		(
			graph::const_iterator iter = begin();
			iter != end();
			++iter
		)
		{
			os << token::SPACE << setw(10) << (*iter())[xi];
		}
		os << endl;
	}
}


void graph::write(Ostream& os, const word& format) const
{
	writer::New(format)().write(*this, os);
}


void graph::write(const fileName& fName, const word& format) const
{
	autoPtr<writer> graphWriter(writer::New(format));

	OFstream graphFile(fName + '.' + graphWriter().ext());

	if (graphFile.good())
	{
		write(graphFile, format);
	}
	else
	{
		WarningIn("graph::write(const word& format, const fileName& dir)")
			<< "Could not open graph file " << graphFile.name()
			<< endl;
	}
}


Ostream& operator<<(Ostream& os, const graph& g)
{
	// Write size of list and start contents delimiter
	os << nl << g.x().size() << nl << token::BEGIN_LIST << nl;

	g.writeTable(os);

	// Write end of contents delimiter
	os << token::END_LIST;

	os.check("Ostream& operator<<(Ostream&, const graph&)");
	return os;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
