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
	Foam::TimePaths

Description
	A class for addressing time paths without using the Time class.

SourceFiles
	TimePaths.C

\*---------------------------------------------------------------------------*/

#ifndef TimePaths_H
#define TimePaths_H

#include "fileName.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{



class TimePaths
{
	// Private data

		//- Processor case
		bool processorCase_;

		//- Root path
		fileName rootPath_;

		//- Global case name
		fileName globalCaseName_;

		//- Case name
		fileName case_;

		//- System directory name
		word system_;

		//- Constant directory name
		word constant_;


public:

	// Constructors

		//- Construct given database name, rootPath and casePath
		TimePaths
		(
			const fileName& rootPath,
			const fileName& caseName,
			const word& systemName,
			const word& constantName
		);


		//- Construct given database name, rootPath and casePath
		TimePaths
		(
			const bool processorCase,
			const fileName& rootPath,
			const fileName& globalCaseName,
			const fileName& caseName,
			const word& systemName,
			const word& constantName
		);


	// Member functions

		//- Return true if this is a processor case
		bool processorCase() const
		{
			return processorCase_;
		}

		//- Return root path
		const fileName& rootPath() const
		{
			return rootPath_;
		}

		//- Return global case name
		const fileName& globalCaseName() const
		{
			return globalCaseName_;
		}

		//- Return case name
		const fileName& caseName() const
		{
			return case_;
		}

		//- Return system name
		const word& system() const
		{
			return system_;
		}

		//- Return system name for the case
		//  which for parallel runs returns ../system()
		fileName caseSystem() const;

		//- Return constant name
		const word& constant() const
		{
			return constant_;
		}

		//- Return constant name for the case
		//  which for parallel runs returns ../constant()
		fileName caseConstant() const;

		//- Return path
		fileName path() const
		{
			return rootPath()/caseName();
		}

		//- Return system path
		fileName systemPath() const
		{
			return path()/system();
		}

		//- Return constant path
		fileName constantPath() const
		{
			return path()/constant();
		}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
