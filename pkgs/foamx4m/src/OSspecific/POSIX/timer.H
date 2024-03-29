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
	Foam::timer

Description
	Implements a timeout mechanism via sigalarm.

	Example usage:
	\code
		timer myTimer(5);     // 5 sec
		..
		if (timedOut(myTimer))
		{
			// timed out
		}
		else
		{
			// do something possible blocking
		}
	\endcode

	Constructor set signal handler on sigalarm and alarm(). Destructor
	clears these.

	timedOut is macro because setjmp can't be in member function of timer.
	?something to do with stack frames.

Warning
	The setjmp restores complete register state so including local vars
	held in regs. So if in blocking part something gets calced in a stack
	based variable make sure it is declared 'volatile'.

SourceFiles
	timer.C

\*---------------------------------------------------------------------------*/

#ifndef timer_H
#define timer_H

#include "className.H"

#include <signal.h>
#include <setjmp.h>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Check it a timeout has occured
//  keep setjmp in same stack frame so no function calls
#define timedOut(x) \
	(((x).newTimeOut_ > 0) ? setjmp(Foam::timer::envAlarm) : false)

namespace Foam
{


class timer
{
	// Private data

		//- old signal masks
		static struct sigaction oldAction_;

		//- old alarm() value
		static unsigned int oldTimeOut_;


	// Private Member Functions

		//- alarm handler
		static void signalHandler(int);


public:

	// Public data

		//- Declare name of the class and its debug switch
		ClassName("timer");

		//- current time out value. Needed by macro timedOut
		unsigned int newTimeOut_;

		//- state for setjmp. Needed by macro timedOut
		static jmp_buf envAlarm;


	// Constructors

		//- Construct from components.
		//  newTimeOut=0 makes it do nothing.
		timer(const unsigned int newTimeOut);


	// Destructor

		~timer();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
