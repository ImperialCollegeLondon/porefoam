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

#include "GAMGAgglomeration.H"
#include "lduMesh.H"
#include "lduMatrix.H"
#include "foamTime.H"
#include "dlLibraryTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
	defineTypeNameAndDebug(GAMGAgglomeration, 0);
	defineRunTimeSelectionTable(GAMGAgglomeration, lduMesh);
	defineRunTimeSelectionTable(GAMGAgglomeration, lduMatrix);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::GAMGAgglomeration::compactLevels(const label nCreatedLevels)
{
	nCells_.setSize(nCreatedLevels);
	restrictAddressing_.setSize(nCreatedLevels);
	meshLevels_.setSize(nCreatedLevels);
	interfaceLevels_.setSize(nCreatedLevels + 1);
}


bool Foam::GAMGAgglomeration::continueAgglomerating
(
	const label nCoarseCells
) const
{
	// Check the need for further agglomeration on all processors
	bool contAgg = nCoarseCells >= nCellsInCoarsestLevel_;
	reduce(contAgg, andOp<bool>());
	return contAgg;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::GAMGAgglomeration::GAMGAgglomeration
(
	const lduMesh& mesh,
	const dictionary& dict
)
:
	MeshObject<lduMesh, GAMGAgglomeration>(mesh),

	maxLevels_(50),

	nCellsInCoarsestLevel_
	(
		readLabel(dict.lookup("nCellsInCoarsestLevel"))
	),

	nCells_(maxLevels_),
	restrictAddressing_(maxLevels_),
	faceRestrictAddressing_(maxLevels_),

	meshLevels_(maxLevels_),
	interfaceLevels_(maxLevels_ + 1)
{}


const Foam::GAMGAgglomeration& Foam::GAMGAgglomeration::New
(
	const lduMesh& mesh,
	const dictionary& dict
)
{
	if
	(
		!mesh.thisDb().foundObject<GAMGAgglomeration>
		(
			GAMGAgglomeration::typeName
		)
	)
	{
		word agglomeratorType(dict.lookup("agglomerator"));

		const_cast<Time&>(mesh.thisDb().time()).libs().open
		(
			dict,
			"geometricGAMGAgglomerationLibs",
			lduMeshConstructorTablePtr_
		);

		lduMeshConstructorTable::iterator cstrIter =
			lduMeshConstructorTablePtr_->find(agglomeratorType);

		if (cstrIter == lduMeshConstructorTablePtr_->end())
		{
			FatalErrorIn
			(
				"GAMGAgglomeration::New"
				"(const lduMesh& mesh, const dictionary& dict)"
			)   << "Unknown GAMGAgglomeration type "
				<< agglomeratorType << ".\n"
				<< "Valid algebraic GAMGAgglomeration types are :"
				<< lduMatrixConstructorTablePtr_->sortedToc() << endl
				<< "Valid algebraic GAMGAgglomeration types are :"
				<< lduMeshConstructorTablePtr_->sortedToc()
				<< exit(FatalError);
		}

		return store(cstrIter()(mesh, dict).ptr());
	}
	else
	{
		return mesh.thisDb().lookupObject<GAMGAgglomeration>
		(
			GAMGAgglomeration::typeName
		);
	}
}


const Foam::GAMGAgglomeration& Foam::GAMGAgglomeration::New
(
	const lduMatrix& matrix,
	const dictionary& dict
)
{
	const lduMesh& mesh = matrix.mesh();

	if
	(
		!mesh.thisDb().foundObject<GAMGAgglomeration>
		(
			GAMGAgglomeration::typeName
		)
	)
	{
		word agglomeratorType(dict.lookup("agglomerator"));

		const_cast<Time&>(mesh.thisDb().time()).libs().open
		(
			dict,
			"algebraicGAMGAgglomerationLibs",
			lduMatrixConstructorTablePtr_
		);

		if
		(
			!lduMatrixConstructorTablePtr_
		 || !lduMatrixConstructorTablePtr_->found(agglomeratorType)
		)
		{
			return New(mesh, dict);
		}
		else
		{
			lduMatrixConstructorTable::iterator cstrIter =
				lduMatrixConstructorTablePtr_->find(agglomeratorType);

			return store(cstrIter()(matrix, dict).ptr());
		}
	}
	else
	{
		return mesh.thisDb().lookupObject<GAMGAgglomeration>
		(
			GAMGAgglomeration::typeName
		);
	}
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::GAMGAgglomeration::~GAMGAgglomeration()
{
	// Clear the interface storage by hand.
	// It is a list of ptrs not a PtrList for consistency of the interface
	for (label leveli=1; leveli<interfaceLevels_.size(); leveli++)
	{
		lduInterfacePtrsList& curLevel = interfaceLevels_[leveli];

		forAll (curLevel, i)
		{
			if (curLevel.set(i))
			{
				delete curLevel(i);
			}
		}
	}
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::lduMesh& Foam::GAMGAgglomeration::meshLevel
(
	const label i
) const
{
	if (i == 0)
	{
		return mesh();
	}
	else
	{
		return meshLevels_[i - 1];
	}
}


const Foam::lduInterfacePtrsList& Foam::GAMGAgglomeration::interfaceLevel
(
	const label i
) const
{
	return interfaceLevels_[i];
}


// ************************************************************************* //
