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

#include "SDA.H"
#include "addToRunTimeSelectionTable.H"
#include "mathematicalConstants.H"

using namespace Foam::mathematicalConstant;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solidBodyMotionFunctions
{
	defineTypeNameAndDebug(SDA, 0);
	addToRunTimeSelectionTable(solidBodyMotionFunction, SDA, dictionary);
};
};


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::septernion
Foam::solidBodyMotionFunctions::SDA::calcTransformation(const scalar t) const
{
	scalar Tpi = Tp_ + dTp_*(t/dTi_);   // Current roll period [sec]
	scalar wr = 2*pi/Tpi; // Current freq [/sec]

	// Current phase for roll [rad]
	scalar r = dTp_/dTi_;
	scalar u = Tp_ + r*t;
	scalar phr = 2*pi*((Tp_/u - 1) + log(mag(u)) - log(Tp_))/r;

	// Current phase for Sway [rad]
	scalar phs = phr + pi;

	// Current phase for Heave [rad]
	scalar phh = phr + pi/2;

	scalar rollA = max(rollAmax_*exp(-sqr(Tpi - Tpn_)/(2*Q_)), rollAmin_);

	vector T
	(
		0,
		swayA_*(sin(wr*t + phs) - sin(phs)),
		heaveA_*(sin(wr*t + phh) - sin(phh))
	);
	quaternion R(rollA*sin(wr*t + phr), 0, 0);
	septernion TR(septernion(CofG_ + T)*R*septernion(-CofG_));

	return TR;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidBodyMotionFunctions::SDA::SDA
(
	const dictionary& SBMFCoeffs,
	const Time& runTime
)
:
	solidBodyMotionFunction(SBMFCoeffs, runTime),
	CofG_(SBMFCoeffs_.lookup("CofG"))
{
	read(SBMFCoeffs);
}


// * * * * * * * * * * * * * * * * Destructors * * * * * * * * * * * * * * * //

Foam::solidBodyMotionFunctions::SDA::~SDA()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::septernion Foam::solidBodyMotionFunctions::SDA::transformation() const
{
	scalar t = time_.value();

	septernion TR = calcTransformation(t);

	Info<< "solidBodyMotionFunctions::SDA::transformation(): "
		<< "Time = " << t << " transformation: " << TR << endl;

	return TR;
}


Foam::septernion Foam::solidBodyMotionFunctions::SDA::velocity() const
{
	// Velocity is calculated as a difference
	scalar t = time_.value();
	scalar dt = time_.deltaT().value();

	const septernion velocity
	(
		(calcTransformation(t).t() - calcTransformation(t - dt).t())/dt,
		(calcTransformation(t).r()/calcTransformation(t - dt).r())/dt
	);

	return velocity;
}


bool Foam::solidBodyMotionFunctions::SDA::read(const dictionary& SBMFCoeffs)
{
	solidBodyMotionFunction::read(SBMFCoeffs);

	SBMFCoeffs_.lookup("CofG") >> CofG_;
	SBMFCoeffs_.lookup("lamda") >> lamda_;
	SBMFCoeffs_.lookup("rollAmax") >> rollAmax_;
	SBMFCoeffs_.lookup("rollAmin") >> rollAmin_;
	SBMFCoeffs_.lookup("heaveA") >> heaveA_;
	SBMFCoeffs_.lookup("swayA") >> swayA_;
	SBMFCoeffs_.lookup("Q") >> Q_;
	SBMFCoeffs_.lookup("Tp") >> Tp_;
	SBMFCoeffs_.lookup("Tpn") >> Tpn_;
	SBMFCoeffs_.lookup("dTi") >> dTi_;
	SBMFCoeffs_.lookup("dTp") >> dTp_;

	// Rescale parameters according to the given scale parameter
	if (lamda_ > 1 + SMALL)
	{
		heaveA_ /= lamda_;
		swayA_ /= lamda_;
		Tp_ /= sqrt(lamda_);
		Tpn_ /= sqrt(lamda_);
		dTi_ /= sqrt(lamda_);
		dTp_ /= sqrt(lamda_);
	}

	return true;
}


// ************************************************************************* //
