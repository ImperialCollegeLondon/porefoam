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

#include "timeSelector.H"
#include "ListOps.H"
#include "argList.H"
#include "foamTime.H"
#include "IStringStream.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::timeSelector::timeSelector()
:
	scalarRanges()
{}


Foam::timeSelector::timeSelector(Istream& is)
:
	scalarRanges(is)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::timeSelector::selected(const instant& value) const
{
	return scalarRanges::selected(value.value());
}


Foam::List<bool> Foam::timeSelector::selected(const instantList& Times) const
{
	boolList lst(Times.size(), false);

	// Check ranges, avoid false positive on constant/
	forAll(Times, timeI)
	{
		if (Times[timeI].name() != "constant" && selected(Times[timeI]))
		{
			lst[timeI] = true;
		}
	}

	// Check specific values
	forAll(*this, rangeI)
	{
		if (operator[](rangeI).isExact())
		{
			scalar target = operator[](rangeI).value();

			int nearestIndex = -1;
			scalar nearestDiff = Foam::GREAT;

			forAll(Times, timeI)
			{
				if (Times[timeI].name() == "constant") continue;

				scalar diff = fabs(Times[timeI].value() - target);
				if (diff < nearestDiff)
				{
					nearestDiff = diff;
					nearestIndex = timeI;
				}
			}

			if (nearestIndex >= 0)
			{
				lst[nearestIndex] = true;
			}
		}
	}

	return lst;
}


Foam::instantList Foam::timeSelector::select(const instantList& Times)
const
{
	return subset(selected(Times), Times);
}


void Foam::timeSelector::inplaceSelect(instantList& Times) const
{
	inplaceSubset(selected(Times), Times);
}


void Foam::timeSelector::addOptions
(
	const bool constant,
	const bool withZero
)
{
	if (constant)
	{
		argList::addBoolOption
		(
			"constant",
			"include the 'constant/' dir in the times list"
		);
	}
	if (withZero)
	{
		argList::addBoolOption
		(
			"withZero",
			"include the '0/' dir in the times list"
		);
	}
	argList::addBoolOption
	(
		"noZero",
		"exclude the '0/' dir from the times list, "
		"has precedence over the -withZero option"
	);
	argList::addBoolOption
	(
		"latestTime",
		"select the latest time"
	);
	argList::addBoolOption
	(
		"newTimes",
		"select the new times"
	);
	argList::addOption
	(
		"time",
		"ranges",
		"comma-separated time ranges - eg, ':10,20,40:70,1000:'"
	);
}


Foam::instantList Foam::timeSelector::select
(
	const instantList& timeDirs,
	const argList& args,
	const word& constantName
)
{
	if (timeDirs.size())
	{
		boolList selectTimes(timeDirs.size(), true);

		// Determine locations of constant/ and 0/ directories
		label constantIdx = -1;
		label zeroIdx = -1;

		forAll(timeDirs, timeI)
		{
			if (timeDirs[timeI].name() == constantName)
			{
				constantIdx = timeI;
			}
			else if (timeDirs[timeI].value() == 0)
			{
				zeroIdx = timeI;
			}

			if (constantIdx >= 0 && zeroIdx >= 0)
			{
				break;
			}
		}

		// Determine latestTime selection (if any)
		// This must appear before the -time option processing
		label latestIdx = -1;
		if (args.optionFound("latestTime"))
		{
			selectTimes = false;
			latestIdx = timeDirs.size() - 1;

			// Avoid false match on constant/
			if (latestIdx == constantIdx)
			{
				latestIdx = -1;
			}
		}

		if (args.optionFound("time"))
		{
			// Can match 0/, but can never match constant/
			selectTimes = timeSelector
			(
				args.optionLookup("time")()
			).selected(timeDirs);
		}

		// Add in latestTime (if selected)
		if (latestIdx >= 0)
		{
			selectTimes[latestIdx] = true;
		}

		if (constantIdx >= 0)
		{
			// Only add constant/ if specifically requested
			selectTimes[constantIdx] = args.optionFound("constant");
		}

		// Special treatment for 0/
		if (zeroIdx >= 0)
		{
			if (args.optionFound("noZero"))
			{
				// Exclude 0/ if specifically requested
				selectTimes[zeroIdx] = false;
			}
			else if (argList::validOptions.found("withZero"))
			{
				// With -withZero enabled, drop 0/ unless specifically requested
				selectTimes[zeroIdx] = args.optionFound("withZero");
			}
		}

		return subset(selectTimes, timeDirs);
	}
	else
	{
		return timeDirs;
	}
}


Foam::instantList Foam::timeSelector::select0
(
	Time& runTime,
	const argList& args
)
{
	instantList timeDirs
	(
		timeSelector::select
		(
			runTime.times(),
			args,
			runTime.constant()
		)
	);

	if (timeDirs.empty())
	{
		WarningIn(args.executable())
			<< "No time specified or available, selecting 'constant'"
			<< endl;

		timeDirs.append(instant(0, runTime.constant()));
	}

	runTime.setTime(timeDirs[0], 0);

	return timeDirs;
}


Foam::instantList Foam::timeSelector::selectIfPresent
(
	Time& runTime,
	const argList& args
)
{
	if
	(
		args.optionFound("latestTime")
	 || args.optionFound("time")
	 || args.optionFound("constant")
	 || args.optionFound("noZero")
	 || args.optionFound("withZero")
	)
	{
		return select0(runTime, args);
	}
	else
	{
		// No timeSelector option specified. Do not change runTime.
		return instantList(1, instant(runTime.value(), runTime.timeName()));
	}
}


Foam::instantList Foam::timeSelector::select
(
	Time& runTime,
	const argList& args,
	const word& fName
)
{
	instantList timeDirs(timeSelector::select0(runTime, args));

	if (timeDirs.size() && args.optionFound("newTimes"))
	{
		boolList selectTimes(timeDirs.size(), true);

		forAll(timeDirs, timeI)
		{
			selectTimes[timeI] =
				!exists(runTime.path()/timeDirs[timeI].name()/fName);
		}

		return subset(selectTimes, timeDirs);
	}
	else
	{
		return timeDirs;
	}
}


// ************************************************************************* //
