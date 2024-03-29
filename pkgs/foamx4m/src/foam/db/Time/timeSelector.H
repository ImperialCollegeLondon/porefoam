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

Class
	Foam::timeSelector

Description
	A List of scalarRange for selecting times.

	The timeSelector provides a convenient means of selecting multiple
	times. A typical use would be the following:

	\verbatim
	timeSelector::addOptions();
	// add other options
	#include "setRootCase.H"
	#include "createTime.H"
	instantList timeDirs = timeSelector::select0(runTime, args);
	...
	forAll(timeDirs, timeI)
	{
		...
	}
	\endverbatim

	The result program would receive \b -time, @b -latestTime, @b -constant
	and \b -noZero options. The @b -constant option explicitly includes the
	\c constant/ directory in the time list and the \b -noZero option
	explicitly excludes the \c 0/ directory from the time list.

	There may however also be many cases in which neither the \c constant/
	directory nor the \c 0/ directory contain particularly relevant
	information. This might occur, for example, when post-processing
	results. In this case, addOptions is called with optional boolean
	arguments.

	\verbatim
	timeSelector::addOptions(false, true);
	\endverbatim

	The first argument avoids adding the \b -constant option. The second
	argument adds an additional \b -withZero option and also prevents the
	\c 0/ directory from being included in the default time range and in the
	\b -latestTime selection.

SourceFiles
	timeSelector.C

\*---------------------------------------------------------------------------*/

#ifndef timeSelector_H
#define timeSelector_H

#include "scalarRanges.H"
#include "instantList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class argList;
class Time;


class timeSelector
:
	public scalarRanges
{
public:

	// Constructors

		//- Construct Null
		timeSelector();

		//- Construct from Istream
		timeSelector(Istream&);


	// Member Functions

		//- Return true if the given instant is within the ranges
		bool selected(const instant&) const;

		//- Return the set of selected instants in the given list that are
		//  within the ranges
		boolList selected(const instantList&) const;

		//- Select a list of Time values that are within the ranges
		instantList select(const instantList&) const;

		//- Select a list of Time values that are within the ranges
		void inplaceSelect(instantList&) const;

		//- Add the options handled by timeSelector to argList::validOptions
		//
		// \param constant
		//   Add the \b -constant option to include the \c constant/ directory
		//
		// \param withZero
		//   Enable the \b -withZero option and alter the normal time selection
		//   behaviour (and \b -latestTime behaviour) to exclude the \c 0/
		//   directory. The \c 0/ directory will only be included when
		//   \b -withZero is specified.
		//   The \b -noZero option has precedence over the @b -withZero option.
		static void addOptions
		(
			const bool constant=true,
			const bool withZero=false
		);

		//- Return the set of times selected based on the argList options
		static instantList select
		(
			const instantList&,
			const argList& args,
			const word& constantName = "constant"
		);

		//- Return the set of times selected based on the argList options
		//  also set the runTime to the first instance or the
		//  \c constant/ directory if no instances are specified or available
		static instantList select0
		(
			Time& runTime,
			const argList& args
		);

		//- If any time option provided return the set of times (as select0)
		//  otherwise return just the current time.
		//  Also set the runTime to the first instance
		static instantList selectIfPresent
		(
			Time& runTime,
			const argList& args
		);

		//- Return the set of times selected based on the argList options
		//  including support for \b -newTimes in which times are selected
		//  if the file <fName> does not exist in the time directory.
		//  Also set the runTime to the first instance or the
		//  \c constant/ directory if no instances are specified or available
		static instantList select
		(
			Time& runTime,
			const argList& args,
			const word& fName
		);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
