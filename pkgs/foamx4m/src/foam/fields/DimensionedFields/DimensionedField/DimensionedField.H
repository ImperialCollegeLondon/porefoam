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
	Foam::DimensionedField

Description
	Field with dimensions and associated with geometry type GeoMesh which is
	used to size the field and a reference to it is maintained.

SourceFiles
	DimensionedFieldI.H
	DimensionedField.C
	DimensionedFieldIO.C

\*---------------------------------------------------------------------------*/

#ifndef DimensionedField_H
#define DimensionedField_H

#include "foamTime.H"
#include "regIOobject.H"
#include "Field.H"
#include "dimensionedType.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Type, class GeoMesh> class DimensionedField;

template<class Type, class GeoMesh> Ostream& operator<<
(
	Ostream&,
	const DimensionedField<Type, GeoMesh>&
);

template<class Type, class GeoMesh> Ostream& operator<<
(
	Ostream&,
	const tmp<DimensionedField<Type, GeoMesh> >&
);



template<class Type, class GeoMesh>
class DimensionedField
:
	public regIOobject,
	public Field<Type>
{

public:

	// Public typedefs

		//- Field mesh type (eg. volMesh, surfaceMesh, areaMesh)
		typedef GeoMesh GeoMeshType;

		//- Actual mesh type on which the DimensionedField is instantiated
		//  (eg fvMesh, faMesh)
		typedef typename GeoMesh::Mesh Mesh;

		//- Field element type (eg. scalar, vector)
		typedef Type PrimitiveType;

		//- Field element component type
		typedef typename Field<Type>::cmptType cmptType;


private:

	// Private data

		//- Reference to mesh
		const Mesh& mesh_;

		//- Dimension set for this field
		dimensionSet dimensions_;


public:

	//- Runtime type information
	TypeName("DimensionedField");


	// Static Member Functions

		//- Return a null DimensionedField
		//  Dangerous, but no null constructor.  Reconsider.  HJ, 10/Apr/2013
		inline static const DimensionedField<Type, GeoMesh>& null();


	// Constructors

		//- Construct from components
		DimensionedField
		(
			const IOobject&,
			const Mesh& mesh,
			const dimensionSet&,
			const Field<Type>&
		);

		//- Construct from components
		//  Used for temporary fields which are initialised after construction
		DimensionedField
		(
			const IOobject&,
			const Mesh& mesh,
			const dimensionSet&
		);

		//- Construct from components
		DimensionedField
		(
			const IOobject&,
			const Mesh& mesh,
			const dimensioned<Type>&
		);

		//- Construct from Istream
		DimensionedField
		(
			const IOobject&,
			const Mesh& mesh,
			const word& fieldDictEntry="value"
		);

		void readField
		(
			const dictionary& fieldDict,
			const word& fieldDictEntry="value"
		);

		//- Construct as copy
		DimensionedField
		(
			const DimensionedField<Type, GeoMesh>&
		);

		//- Construct as copy or re-use as specified.
		DimensionedField
		(
			DimensionedField<Type, GeoMesh>&,
			bool reUse
		);

		//- Construct by transferring the DimensionedField
		DimensionedField
		(
			const Xfer<DimensionedField<Type, GeoMesh> >&
		);

		//- Construct as copy of tmp<DimensionedField> deleting argument
#		ifdef ConstructFromTmp
		DimensionedField
		(
			const tmp<DimensionedField<Type, GeoMesh> >&
		);
#		endif

		//- Construct as copy resetting IO parameters
		DimensionedField
		(
			const IOobject&,
			const DimensionedField<Type, GeoMesh>&
		);

		//- Construct as copy resetting name
		DimensionedField
		(
			const word& newName,
			const DimensionedField<Type, GeoMesh>&
		);

		//- Construct as copy resetting name and re-use as specified.
		DimensionedField
		(
			const word& newName,
			DimensionedField<Type, GeoMesh>&,
			bool reUse
		);

		//- Construct by transferring the DimensionedField with a new name
		DimensionedField
		(
			const word& newName,
			const Xfer<DimensionedField<Type, GeoMesh> >&
		);

		//- Construct as copy resetting name
#		ifdef ConstructFromTmp
		DimensionedField
		(
			const word& newName,
			const tmp<DimensionedField<Type, GeoMesh> >&
		);
#		endif

		//- Clone
		tmp<DimensionedField<Type, GeoMesh> > clone() const;


	//- Destructor
	virtual ~DimensionedField();


	// Member Functions

		//- Return mesh
		inline const Mesh& mesh() const;

		//- Return dimensions
		inline const dimensionSet& dimensions() const;

		//- Return non-const access to dimensions
		inline dimensionSet& dimensions();

		inline const Field<Type>& field() const;

		inline Field<Type>& field();

		//- Return a component field of the field
		tmp<DimensionedField<cmptType, GeoMesh> > component
		(
			const direction
		) const;

		//- Replace a component field of the field
		void replace
		(
			const direction,
			const DimensionedField<cmptType, GeoMesh>&
		);

		//- Replace a component field of the field
		void replace
		(
			const direction,
			const tmp<DimensionedField<cmptType, GeoMesh> >&
		);

		//- Return the field transpose (only defined for second rank tensors)
		tmp<DimensionedField<Type, GeoMesh> > T() const;

		//- Calculate and return arithmetic average
		dimensioned<Type> average() const;

		//- Calculate and return weighted average
		dimensioned<Type> weightedAverage
		(
			const DimensionedField<scalar, GeoMesh>&
		) const;

		//- Calculate and return weighted average
		dimensioned<Type> weightedAverage
		(
			const tmp<DimensionedField<scalar, GeoMesh> >&
		) const;


		// Write

			bool writeData(Ostream&, const word& fieldDictEntry) const;

			bool writeData(Ostream&) const;


	// Member Operators

		void operator=(const DimensionedField<Type, GeoMesh>&);
		void operator=(const tmp<DimensionedField<Type, GeoMesh> >&);
		void operator=(const dimensioned<Type>&);

		void operator+=(const DimensionedField<Type, GeoMesh>&);
		void operator+=(const tmp<DimensionedField<Type, GeoMesh> >&);

		void operator-=(const DimensionedField<Type, GeoMesh>&);
		void operator-=(const tmp<DimensionedField<Type, GeoMesh> >&);

		void operator*=(const DimensionedField<scalar, GeoMesh>&);
		void operator*=(const tmp<DimensionedField<scalar, GeoMesh> >&);

		void operator/=(const DimensionedField<scalar, GeoMesh>&);
		void operator/=(const tmp<DimensionedField<scalar, GeoMesh> >&);

		void operator+=(const dimensioned<Type>&);
		void operator-=(const dimensioned<Type>&);

		void operator*=(const dimensioned<scalar>&);
		void operator/=(const dimensioned<scalar>&);


	// Ostream Operators

		friend Ostream& operator<< <Type, GeoMesh>
		(
			Ostream&,
			const DimensionedField<Type, GeoMesh>&
		);

		friend Ostream& operator<< <Type, GeoMesh>
		(
			Ostream&,
			const tmp<DimensionedField<Type, GeoMesh> >&
		);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "DimensionedFieldI.H"
#include "DimensionedFieldFunctions.H"

#ifdef NoRepository
#	include "DimensionedField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
