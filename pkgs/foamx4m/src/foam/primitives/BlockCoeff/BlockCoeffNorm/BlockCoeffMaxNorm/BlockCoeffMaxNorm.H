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
	BlockCoeffMaxNorm

Description
	Class for BlockCoeff maximum element magnitude norm.

Author
	Klas Jareteg, 2013-03-07

SourceFiles
	BlockCoeffMaxNorm.C

\*---------------------------------------------------------------------------*/

#ifndef BlockCoeffMaxNorm_H
#define BlockCoeffMaxNorm_H

#include "blockCoeffs.H"
#include "blockCoeffNorms.H"
#include "BlockCoeffNorm.H"
#include "coeffFields.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * Forward declaration of template friend fuctions * * * * * * * //



template<class Type>
class BlockCoeffMaxNorm
:
	public BlockCoeffNorm<Type>
{
	// Private data

		//- Control data dictionary
		dictionary dict_;


protected:

	// Protected Member Functions

		//- Return dictionary
		const dictionary& dict() const
		{
			return dict_;
		}


public:

	//- Runtime type information
	TypeName("maxNorm");


	// Constructors

		//- Construct from dictionary
		BlockCoeffMaxNorm
		(
			const dictionary& dict
		);

	//- Destructor
	virtual ~BlockCoeffMaxNorm()
	{}


	// Member functions

		//- Solve
		virtual scalar normalize
		(
			const BlockCoeff<Type>& a
		);

		virtual void normalize
		(
			Field<scalar>& b,
			const CoeffField<Type>& a
		);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#	include "BlockCoeffMaxNorm.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
