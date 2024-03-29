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

#include "gnuplotGraph.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(Foam::gnuplotGraph, 0);
const Foam::word Foam::gnuplotGraph::ext_("gplt");

namespace Foam
{
	typedef graph::writer graphWriter;
	addToRunTimeSelectionTable(graphWriter, gnuplotGraph, word);
};


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::gnuplotGraph::write(const graph& g, Ostream& os) const
{
	os  << "#set term postscript color" << endl
		<< "set output \"" << word(g.title()) << ".ps\"" << endl
		<< "set title " << g.title() << " 0,0" << endl << "show title" << endl
		<< "set xlabel " << g.xName() << " 0,0" << endl << "show xlabel" << endl
		<< "set ylabel " << g.yName() << " 0,0" << endl << "show ylabel" << endl
		<< "plot";

	bool firstField = true;

	for (graph::const_iterator iter = g.begin(); iter != g.end(); ++iter)
	{
		if (!firstField)
		{
			os << ',';
		}
		firstField = false;

		os  << "'-' title " << iter()->name() << " with lines";
	}
	os << "; pause -1" << endl;


	for (graph::const_iterator iter = g.begin(); iter != g.end(); ++iter)
	{
		os << endl;
		writeXY(g.x(), *iter(), os);
	}
}


// ************************************************************************* //
