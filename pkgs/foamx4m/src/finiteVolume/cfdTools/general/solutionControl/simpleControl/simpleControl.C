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

#include "simpleControl.H"
#include "foamTime.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
	defineTypeNameAndDebug(simpleControl, 0);
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::simpleControl::read()
{
	solutionControl::read(true);
}


bool Foam::simpleControl::criteriaSatisfied()
{
	if (residualControl_.empty())
	{
		return false;
	}

	bool achieved = true;
	bool checked = false;	// safety that some checks were indeed performed

	const dictionary& solverDict = mesh_.solutionDict().solverPerformanceDict();
	forAllConstIter(dictionary, solverDict, iter)
	{
		const word& variableName = iter().keyword();
		const label fieldI = applyToField(variableName);

		if (fieldI != -1)
		{
			scalar lastResidual = 0;
			const scalar residual =
				maxResidual(variableName, iter().stream(), lastResidual);

			checked = true;

			bool absCheck = residual < residualControl_[fieldI].absTol;
			achieved = achieved && absCheck;

			if (debug)
			{
				Info<< algorithmName_ << " solution statistics:" << endl;

				Info<< "    " << variableName << ": tolerance = " << residual
					<< " (" << residualControl_[fieldI].absTol << ")"
					<< endl;
			}
		}
	}

	return checked && achieved;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::simpleControl::simpleControl(fvMesh& mesh)
:
	solutionControl(mesh, "SIMPLE"),
	initialised_(false)
{
	read();

	Info<< nl;

	if (residualControl_.empty())
	{
		Info<< algorithmName_ << ": no convergence criteria found. "
			<< "Calculations will run for " << mesh_.time().endTime().value()
			<< " steps." << nl << endl;
	}
	else
	{
		Info<< algorithmName_ << ": convergence criteria" << nl;
		forAll(residualControl_, i)
		{
			Info<< "    field " << residualControl_[i].name << token::TAB
				<< " tolerance " << residualControl_[i].absTol
				<< nl;
		}
		Info<< endl;
	}
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::simpleControl::~simpleControl()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::simpleControl::loop()
{
	read();

	Time& time = const_cast<Time&>(mesh_.time());

	if (initialised_)
	{
		if (criteriaSatisfied())
		{
			Info<< nl << algorithmName_ << " solution converged in "
				<< time.value() << " iterations" << nl << endl;

			// Set to finalise calculation
			time.writeAndEnd();
		}
		else
		{
			storePrevIterFields();
		}
	}
	else
	{
		initialised_ = true;
		storePrevIterFields();
	}


	return time.loop();
}


// ************************************************************************* //
