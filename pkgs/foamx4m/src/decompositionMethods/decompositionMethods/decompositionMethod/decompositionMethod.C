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

InClass
	decompositionMethod

\*---------------------------------------------------------------------------*/

#include "decompositionMethod.H"
#include "cyclicPolyPatch.H"
#include "syncTools.H"
#include "globalIndex.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
	defineTypeNameAndDebug(decompositionMethod, 0);
	defineRunTimeSelectionTable(decompositionMethod, dictionary);
	defineRunTimeSelectionTable(decompositionMethod, dictionaryMesh);
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::decompositionMethod::calcCSR
(
	const labelListList& cellCells,
	labelList& adjncy,
	labelList& xadj
)
{
	// Count number of internal faces
	label nConnections = 0;

	forAll(cellCells, coarseI)
	{
		nConnections += cellCells[coarseI].size();
	}

	// Create the adjncy array as twice the size of the total number of
	// internal faces
	adjncy.setSize(nConnections);

	xadj.setSize(cellCells.size() + 1);


	// Fill in xadj
	// ~~~~~~~~~~~~
	label freeAdj = 0;

	forAll(cellCells, coarseI)
	{
		xadj[coarseI] = freeAdj;

		const labelList& cCells = cellCells[coarseI];

		forAll(cCells, i)
		{
			adjncy[freeAdj++] = cCells[i];
		}
	}
	xadj[cellCells.size()] = freeAdj;
}



void Foam::decompositionMethod::calcCSR
(
	const polyMesh& mesh,
	labelList& adjncy,
	labelList& xadj
)
{
	// Make Metis CSR (Compressed Storage Format) storage
	//   adjncy      : contains neighbours (= edges in graph)
	//   xadj(celli) : start of information in adjncy for celli

	xadj.setSize(mesh.nCells() + 1);

	// Initialise the number of internal faces of the cells to twice the
	// number of internal faces
	label nInternalFaces = 2*mesh.nInternalFaces();

	// Check the boundary for coupled patches and add to the number of
	// internal faces
	const polyBoundaryMesh& pbm = mesh.boundaryMesh();

	forAll(pbm, patchi)
	{
		if (isA<cyclicPolyPatch>(pbm[patchi]))
		{
			nInternalFaces += pbm[patchi].size();
		}
	}

	// Create the adjncy array the size of the total number of internal and
	// coupled faces
	adjncy.setSize(nInternalFaces);

	// Fill in xadj
	// ~~~~~~~~~~~~
	label freeAdj = 0;

	for (label cellI = 0; cellI < mesh.nCells(); cellI++)
	{
		xadj[cellI] = freeAdj;

		const labelList& cFaces = mesh.cells()[cellI];

		forAll(cFaces, i)
		{
			label faceI = cFaces[i];

			if
			(
				mesh.isInternalFace(faceI)
			 || isA<cyclicPolyPatch>(pbm[pbm.whichPatch(faceI)])
			)
			{
				freeAdj++;
			}
		}
	}
	xadj[mesh.nCells()] = freeAdj;


	// Fill in adjncy
	// ~~~~~~~~~~~~~~

	labelList nFacesPerCell(mesh.nCells(), 0);

	// Internal faces
	for (label faceI = 0; faceI < mesh.nInternalFaces(); faceI++)
	{
		label own = mesh.faceOwner()[faceI];
		label nei = mesh.faceNeighbour()[faceI];

		adjncy[xadj[own] + nFacesPerCell[own]++] = nei;
		adjncy[xadj[nei] + nFacesPerCell[nei]++] = own;
	}

	// Coupled faces. Only cyclics done.
	forAll(pbm, patchi)
	{
		if (isA<cyclicPolyPatch>(pbm[patchi]))
		{
			const unallocLabelList& faceCells = pbm[patchi].faceCells();

			label sizeby2 = faceCells.size()/2;

			for (label facei=0; facei<sizeby2; facei++)
			{
				label own = faceCells[facei];
				label nei = faceCells[facei + sizeby2];

				adjncy[xadj[own] + nFacesPerCell[own]++] = nei;
				adjncy[xadj[nei] + nFacesPerCell[nei]++] = own;
			}
		}
	}
}


void Foam::decompositionMethod::calcDistributedCSR
(
	const polyMesh& mesh,
	labelList& adjncy,
	labelList& xadj
)
{
	// Create global cell numbers
	// ~~~~~~~~~~~~~~~~~~~~~~~~~~

	globalIndex globalCells(mesh.nCells());


	//
	// Make Metis Distributed CSR (Compressed Storage Format) storage
	//   adjncy      : contains cellCells (= edges in graph)
	//   xadj(celli) : start of information in adjncy for celli
	//


	const labelList& faceOwner = mesh.faceOwner();
	const labelList& faceNeighbour = mesh.faceNeighbour();
	const polyBoundaryMesh& patches = mesh.boundaryMesh();


	// Get renumbered owner on other side of coupled faces
	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	labelList globalNeighbour(mesh.nFaces()-mesh.nInternalFaces());

	forAll(patches, patchI)
	{
		const polyPatch& pp = patches[patchI];

		if (pp.coupled())
		{
			label faceI = pp.start();
			label bFaceI = pp.start() - mesh.nInternalFaces();

			forAll(pp, i)
			{
				globalNeighbour[bFaceI++] = globalCells.toGlobal
				(
					faceOwner[faceI++]
				);
			}
		}
	}

	// Get the cell on the other side of coupled patches
	syncTools::swapBoundaryFaceList(mesh, globalNeighbour, false);


	// Count number of faces (internal + coupled)
	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	// Number of faces per cell
	labelList nFacesPerCell(mesh.nCells(), 0);

	// Number of coupled faces
	label nCoupledFaces = 0;

	for (label faceI = 0; faceI < mesh.nInternalFaces(); faceI++)
	{
		nFacesPerCell[faceOwner[faceI]]++;
		nFacesPerCell[faceNeighbour[faceI]]++;
	}
	// Handle coupled faces
	forAll(patches, patchI)
	{
		const polyPatch& pp = patches[patchI];

		if (pp.coupled())
		{
			label faceI = pp.start();

			forAll(pp, i)
			{
				nCoupledFaces++;
				nFacesPerCell[faceOwner[faceI++]]++;
			}
		}
	}


	// Fill in xadj
	// ~~~~~~~~~~~~

	xadj.setSize(mesh.nCells() + 1);

	label freeAdj = 0;

	for (label cellI = 0; cellI < mesh.nCells(); cellI++)
	{
		xadj[cellI] = freeAdj;

		freeAdj += nFacesPerCell[cellI];
	}
	xadj[mesh.nCells()] = freeAdj;



	// Fill in adjncy
	// ~~~~~~~~~~~~~~

	adjncy.setSize(2*mesh.nInternalFaces() + nCoupledFaces);

	nFacesPerCell = 0;

	// For internal faces is just offsetted owner and neighbour
	for (label faceI = 0; faceI < mesh.nInternalFaces(); faceI++)
	{
		label own = faceOwner[faceI];
		label nei = faceNeighbour[faceI];

		adjncy[xadj[own] + nFacesPerCell[own]++] = globalCells.toGlobal(nei);
		adjncy[xadj[nei] + nFacesPerCell[nei]++] = globalCells.toGlobal(own);
	}
	// For boundary faces is offsetted coupled neighbour
	forAll(patches, patchI)
	{
		const polyPatch& pp = patches[patchI];

		if (pp.coupled())
		{
			label faceI = pp.start();
			label bFaceI = pp.start()-mesh.nInternalFaces();

			forAll(pp, i)
			{
				label own = faceOwner[faceI];
				adjncy[xadj[own] + nFacesPerCell[own]++] =
					globalNeighbour[bFaceI];

				faceI++;
				bFaceI++;
			}
		}
	}
}


void Foam::decompositionMethod::calcCellCells
(
	const polyMesh& mesh,
	const labelList& fineToCoarse,
	const label nCoarse,
	labelListList& cellCells
)
{
	if (fineToCoarse.size() != mesh.nCells())
	{
		FatalErrorIn
		(
			"decompositionMethod::calcCellCells"
			"(const labelList&, labelListList&) const"
		)   << "Only valid for mesh agglomeration." << exit(FatalError);
	}

	List<dynamicLabelList > dynCellCells(nCoarse);

	forAll(mesh.faceNeighbour(), faceI)
	{
		label own = fineToCoarse[mesh.faceOwner()[faceI]];
		label nei = fineToCoarse[mesh.faceNeighbour()[faceI]];

		if (own != nei)
		{
			if (findIndex(dynCellCells[own], nei) == -1)
			{
				dynCellCells[own].append(nei);
			}
			if (findIndex(dynCellCells[nei], own) == -1)
			{
				dynCellCells[nei].append(own);
			}
		}
	}

	cellCells.setSize(dynCellCells.size());
	forAll(dynCellCells, coarseI)
	{
		cellCells[coarseI].transfer(dynCellCells[coarseI]);
	}
}


void Foam::decompositionMethod::fixCyclics
(
	const polyMesh& mesh,
	labelList& decomp
)
{
	const polyBoundaryMesh& pbm = mesh.boundaryMesh();

	label nFixedCyclics = 0;

	// Coupled faces. Only cyclics done.
	do
	{
		nFixedCyclics = 0;

		forAll(pbm, patchi)
		{
			if (isA<cyclicPolyPatch>(pbm[patchi]))
			{
				const unallocLabelList& faceCells = pbm[patchi].faceCells();

				const label sizeby2 = faceCells.size()/2;

				for (label facei=0; facei<sizeby2; facei++)
				{
					const label own = faceCells[facei];
					const label nei = faceCells[facei + sizeby2];

					if(decomp[own] < decomp[nei])
					{
						decomp[own] = decomp[nei];
						nFixedCyclics++;
					}
					else if(decomp[own] > decomp[nei])
					{
						decomp[nei] = decomp[own];
						nFixedCyclics++;
					}
				}
			}
		}

		if(nFixedCyclics > 0)
		{
			WarningIn
			(
				"decompositionMethod::fixCyclics"
				"(const polyMesh& mesh, labelList& decomp)"
			)   << "Fixed " << nFixedCyclics << " disconnected cyclic faces";
		}
	}
	while (nFixedCyclics > 0);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::decompositionMethod> Foam::decompositionMethod::New
(
	const dictionary& decompositionDict
)
{
	// HR 11.02.18: Library table is member of runtime (like in
	// vanilla). Therefore need runtime here and we get it by
	// traversing to the top of the dict in the hope that it is
	// an IOdictionary.
	const dictionary& topDict = decompositionDict.topDict();
	if (isA<IOdictionary>(topDict))
	{
		Time& time = const_cast<Time&>
		(
			reinterpret_cast<const Time&>(topDict)
		);

		loadExternalLibraries(time);
	}

	word methodName(decompositionDict.lookup("method"));

	Info<< "Selecting decompositionMethod "<< methodName << endl;

	dictionaryConstructorTable::iterator cstrIter =
		dictionaryConstructorTablePtr_->find(methodName);

	if (cstrIter == dictionaryConstructorTablePtr_->end())
	{
		FatalErrorIn
		(
			"decompositionMethod::New"
			"(const dictionary& decompositionDict)"
		)   << "Unknown decompositionMethod "
			<< methodName << endl << endl
			<< "Valid decompositionMethods are : " << endl
			<< dictionaryConstructorTablePtr_->sortedToc()
			<< exit(FatalError);
	}

	return autoPtr<decompositionMethod>(cstrIter()(decompositionDict));
}


Foam::autoPtr<Foam::decompositionMethod> Foam::decompositionMethod::New
(
	const dictionary& decompositionDict,
	const polyMesh& mesh
)
{
	loadExternalLibraries(const_cast<Time&>(mesh.time()));

	word methodName(decompositionDict.lookup("method"));

	Info<< "Selecting decompositionMethod "
		<< methodName << endl;

	dictionaryMeshConstructorTable::iterator cstrIter =
		dictionaryMeshConstructorTablePtr_->find(methodName);

	if (cstrIter == dictionaryMeshConstructorTablePtr_->end())
	{
		FatalErrorIn
		(
			"decompositionMethod::New"
			"(const dictionary& decompositionDict, "
			"const polyMesh& mesh)"
		)   << "Unknown decompositionMethod "
			<< methodName << endl << endl
			<< "Valid decompositionMethods are : " << endl
			<< dictionaryMeshConstructorTablePtr_->sortedToc()
			<< exit(FatalError);
	}

	return autoPtr<decompositionMethod>(cstrIter()(decompositionDict, mesh));
}


void Foam::decompositionMethod::loadExternalLibraries(Time& time)
{
	wordList libNames(3);
	libNames[0]=word("scotchDecomp");
	libNames[1]=word("metisDecomp");
	libNames[2]=word("parMetisDecomp");

	forAll(libNames,i)
	{
		const word libName("lib"+libNames[i]+DOTSO);

		if(!time.libs().open(libName))
		{
			WarningIn("decompositionMethod::loadExternalLibraries()")
				<< "Loading of decomposition library " << libName
					<< " unsuccesful. Some decomposition methods may not be "
					<< " available"
					<< endl;
		}
	}
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::labelList Foam::decompositionMethod::decompose
(
	const pointField& points
)
{
	scalarField weights(points.size(), 1);

	return decompose(points, weights);
}


Foam::labelList Foam::decompositionMethod::decompose
(
	const labelList& fineToCoarse,
	const pointField& coarsePoints
)
{
	// Decompose based on agglomerated points
	labelList coarseDistribution(decompose(coarsePoints));

	// Rework back into decomposition for original mesh_
	labelList fineDistribution(fineToCoarse.size());

	forAll(fineDistribution, i)
	{
		fineDistribution[i] = coarseDistribution[fineToCoarse[i]];
	}

	return fineDistribution;
}


Foam::labelList Foam::decompositionMethod::decompose
(
	const labelList& fineToCoarse,
	const pointField& coarsePoints,
	const scalarField& coarseWeights
)
{
	// Decompose based on agglomerated points
	labelList coarseDistribution(decompose(coarsePoints, coarseWeights));

	// Rework back into decomposition for original mesh_
	labelList fineDistribution(fineToCoarse.size());

	forAll(fineDistribution, i)
	{
		fineDistribution[i] = coarseDistribution[fineToCoarse[i]];
	}

	return fineDistribution;
}


// ************************************************************************* //
