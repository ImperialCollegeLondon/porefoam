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
	Foam::SphericalTensor2D

Description
	Templated 2D sphericalTensor derived from VectorSpace adding construction
	from 1 component, element access using ii() member function and the
	inner-product (dot-product) and outer-product operators.

SourceFiles
	SphericalTensor2DTemplateI.H

\*---------------------------------------------------------------------------*/

#ifndef SphericalTensor2DTemplate_H
#define SphericalTensor2DTemplate_H

#include "VectorSpace.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


template<class Cmpt>
class SphericalTensor2D
:
	public VectorSpace<SphericalTensor2D<Cmpt>, Cmpt, 1>
{

public:

	// Member constants

		enum
		{
			rank = 2 // Rank of SphericalTensor2D is 2
		};


	// Static data members

		static const char* const typeName;
		static const char* componentNames[];
		static const SphericalTensor2D zero;
		static const SphericalTensor2D one;
		static const SphericalTensor2D max;
		static const SphericalTensor2D min;
		static const SphericalTensor2D I;
		static const SphericalTensor2D oneThirdI;
		static const SphericalTensor2D twoThirdsI;


	//- Component labeling enumeration
	enum components { II };


	// Constructors

		//- Construct null
		inline SphericalTensor2D();

		//- Construct given VectorSpace
		inline SphericalTensor2D
		(
			const VectorSpace<SphericalTensor2D<Cmpt>, Cmpt, 1>&
		);

		//- Construct given the component
		inline SphericalTensor2D(const Cmpt& tii);

		//- Construct from Istream
		inline SphericalTensor2D(Istream&);


	// Member Functions

		// Access

			inline const Cmpt& ii() const;
			inline Cmpt& ii();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Include inline implementations
#include "SphericalTensor2DTemplateI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
