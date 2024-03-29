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

	From scotch forum:

	By: Francois PELLEGRINI RE: Graph mapping 'strategy' string [ reply ]
	2008-08-22 10:09 Strategy handling in Scotch is a bit tricky. In order
	not to be confused, you must have a clear view of how they are built.
	Here are some rules:

	1- Strategies are made up of "methods" which are combined by means of
	"operators".

	2- A method is of the form "m{param=value,param=value,...}", where "m"
	is a single character (this is your first error: "f" is a method name,
	not a parameter name).

	3- There exist different sort of strategies : bipartitioning strategies,
	mapping strategies, ordering strategies, which cannot be mixed. For
	instance, you cannot build a bipartitioning strategy and feed it to a
	mapping method (this is your second error).

	To use the "mapCompute" routine, you must create a mapping strategy, not
	a bipartitioning one, and so use stratGraphMap() and not
	stratGraphBipart(). Your mapping strategy should however be based on the
	"recursive bipartitioning" method ("b"). For instance, a simple (and
	hence not very efficient) mapping strategy can be :

	"b{sep=f}"

	which computes mappings with the recursive bipartitioning method "b",
	this latter using the Fiduccia-Mattheyses method "f" to compute its
	separators.

	If you want an exact partition (see your previous post), try
	"b{sep=fx}".

	However, these strategies are not the most efficient, as they do not
	make use of the multi-level framework.

	To use the multi-level framework, try for instance:

	"b{sep=m{vert=100,low=h,asc=f}x}"

	The current default mapping strategy in Scotch can be seen by using the
	"-vs" option of program gmap. It is, to date:

	b
	{
		job=t,
		map=t,
		poli=S,
		sep=
		(
			m
			{
				asc=b
				{
					bnd=d{pass=40,dif=1,rem=1}f{move=80,pass=-1,bal=0.005},
					org=f{move=80,pass=-1,bal=0.005},
					width=3
				},
				low=h{pass=10}f{move=80,pass=-1,bal=0.0005},
				type=h,
				vert=80,
				rat=0.8
			}
		  | m
			{
				asc=b
				{
					bnd=d{pass=40,dif=1,rem=1}f{move=80,pass=-1,bal=0.005},
					org=f{move=80,pass=-1,bal=0.005},
					width=3
				},
				low=h{pass=10}f{move=80,pass=-1,bal=0.0005},
				type=h,
				vert=80,
				rat=0.8
			}
		)
	}


\*---------------------------------------------------------------------------*/

#include "scotchDecomp.H"
#include "addToRunTimeSelectionTable.H"
#include "floatScalar.H"
#include "foamTime.H"
#include "OFstream.H"

extern "C"
{
#include "scotch.h"
}


// Hack: scotch generates floating point errors so need to switch of error
//	   trapping!
#if defined(linux) || defined(linuxAMD64) || defined(linuxIA64)
#	define LINUX
#endif

#if defined(LINUX) && defined(__GNUC__)
#	define LINUX_GNUC
#endif

#ifdef LINUX_GNUC
#	ifndef __USE_GNU
#	   define __USE_GNU
#	endif
#	include <fenv.h>
#endif


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
	defineTypeNameAndDebug(scotchDecomp, 0);

	addToRunTimeSelectionTable
	(
		decompositionMethod,
		scotchDecomp,
		dictionaryMesh
	);
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::scotchDecomp::check(const int retVal, const char* str)
{
	if (retVal)
	{
		FatalErrorIn("scotchDecomp::decompose(..)")
			<< "Call to scotch routine " << str << " failed."
			<< exit(FatalError);
	}
}


// Call scotch with options from dictionary.
Foam::label Foam::scotchDecomp::decompose
(
	const labelList& adjncy,
	const labelList& xadj,
	const scalarField& cWeights,

	labelList& finalDecomp
)
{
	// Dump graph
	if (decompositionDict_.found("scotchCoeffs"))
	{
		const dictionary& scotchCoeffs =
			decompositionDict_.subDict("scotchCoeffs");

		if (scotchCoeffs.found("writeGraph"))
		{
			Switch writeGraph(scotchCoeffs.lookup("writeGraph"));

			if (writeGraph)
			{
				OFstream str(mesh_.time().path() / mesh_.name() + ".grf");

				Info<< "Dumping Scotch graph file to " << str.name() << endl
					<< "Use this in combination with gpart." << endl;

				label version = 0;
				str << version << nl;
				// Numer of vertices
				str << xadj.size()-1 << ' ' << adjncy.size() << nl;
				// Numbering starts from 0
				label baseval = 0;
				// Has weights?
				label hasEdgeWeights = 0;
				label hasVertexWeights = 0;
				label numericflag = 10*hasEdgeWeights+hasVertexWeights;
				str << baseval << ' ' << numericflag << nl;
				for (label cellI = 0; cellI < xadj.size()-1; cellI++)
				{
					label start = xadj[cellI];
					label end = xadj[cellI+1];
					str << end-start;

					for (label i = start; i < end; i++)
					{
						str << ' ' << adjncy[i];
					}
					str << nl;
				}
			}
		}
	}


	// Strategy
	// ~~~~~~~~

	// Default.
	SCOTCH_Strat stradat;
	check(SCOTCH_stratInit(&stradat), "SCOTCH_stratInit");

	if (decompositionDict_.found("scotchCoeffs"))
	{
		const dictionary& scotchCoeffs =
			decompositionDict_.subDict("scotchCoeffs");


		string strategy;
		if (scotchCoeffs.readIfPresent("strategy", strategy))
		{
			if (debug)
			{
				Info<< "scotchDecomp : Using strategy " << strategy << endl;
			}
			SCOTCH_stratGraphMap(&stradat, strategy.c_str());
			//fprintf(stdout, "S\tStrat=");
			//SCOTCH_stratSave(&stradat, stdout);
			//fprintf(stdout, "\n");
		}
	}


	// Graph
	// ~~~~~

	labelList velotab;


	// Check for externally provided cellweights and if so initialise weights
	scalar minWeights = gMin(cWeights);
	if (cWeights.size() > 0)
	{
		if (minWeights <= 0)
		{
			WarningIn
			(
				"scotchDecomp::decompose"
				"(const pointField&, const scalarField&)"
			)   << "Illegal minimum weight " << minWeights
				<< endl;
		}

		if (cWeights.size() != xadj.size()-1)
		{
			FatalErrorIn
			(
				"scotchDecomp::decompose"
				"(const pointField&, const scalarField&)"
			)   << "Number of cell weights " << cWeights.size()
				<< " does not equal number of cells " << xadj.size()-1
				<< exit(FatalError);
		}

		// Convert to integers.
		velotab.setSize(cWeights.size());
		forAll(velotab, i)
		{
			velotab[i] = int(cWeights[i]/minWeights);
		}
	}



	SCOTCH_Graph grafdat;
	check(SCOTCH_graphInit(&grafdat), "SCOTCH_graphInit");
	check
	(
		SCOTCH_graphBuild
		(
			&grafdat,
			0,					  // baseval, c-style numbering
			xadj.size()-1,		  // vertnbr, nCells
			xadj.begin(),		   // verttab, start index per cell into adjncy
			&xadj[1],			   // vendtab, end index  ,,
			velotab.begin(),		// velotab, vertex weights
			nullptr,				   // vlbltab
			adjncy.size(),		  // edgenbr, number of arcs
			adjncy.begin(),		 // edgetab
			nullptr					// edlotab, edge weights
		),
		"SCOTCH_graphBuild"
	);
	check(SCOTCH_graphCheck(&grafdat), "SCOTCH_graphCheck");


	// Architecture
	// ~~~~~~~~~~~~
	// (fully connected network topology since using switch)

	SCOTCH_Arch archdat;
	check(SCOTCH_archInit(&archdat), "SCOTCH_archInit");

	labelList processorWeights;
	if (decompositionDict_.found("scotchCoeffs"))
	{
		const dictionary& scotchCoeffs =
			decompositionDict_.subDict("scotchCoeffs");

		scotchCoeffs.readIfPresent("processorWeights", processorWeights);
	}
	if (processorWeights.size())
	{
		if (debug)
		{
			Info<< "scotchDecomp : Using procesor weights " << processorWeights
				<< endl;
		}
		check
		(
			SCOTCH_archCmpltw(&archdat, nProcessors_, processorWeights.begin()),
			"SCOTCH_archCmpltw"
		);
	}
	else
	{
		check
		(
			SCOTCH_archCmplt(&archdat, nProcessors_),
			"SCOTCH_archCmplt"
		);
	}


	//SCOTCH_Mapping mapdat;
	//SCOTCH_graphMapInit(&grafdat, &mapdat, &archdat, nullptr);
	//SCOTCH_graphMapCompute(&grafdat, &mapdat, &stradat); /* Perform mapping */
	//SCOTCH_graphMapExit(&grafdat, &mapdat);


	// Hack:switch off fpu error trapping
#	ifdef LINUX_GNUC
	int oldExcepts = fedisableexcept
	(
		FE_DIVBYZERO
	  | FE_INVALID
	  | FE_OVERFLOW
	);
#	endif

	finalDecomp.setSize(xadj.size()-1);
	finalDecomp = 0;
	check
	(
		SCOTCH_graphMap
		(
			&grafdat,
			&archdat,
			&stradat,		   // const SCOTCH_Strat *
			finalDecomp.begin() // parttab
		),
		"SCOTCH_graphMap"
	);

#	ifdef LINUX_GNUC
	feenableexcept(oldExcepts);
#	endif



	//finalDecomp.setSize(xadj.size()-1);
	//check
	//(
	//	SCOTCH_graphPart
	//	(
	//		&grafdat,
	//		nProcessors_,	   // partnbr
	//		&stradat,		   // const SCOTCH_Strat *
	//		finalDecomp.begin() // parttab
	//	),
	//	"SCOTCH_graphPart"
	//);

	// Release storage for graph
	SCOTCH_graphExit(&grafdat);
	// Release storage for strategy
	SCOTCH_stratExit(&stradat);
	// Release storage for network topology
	SCOTCH_archExit(&archdat);

	return 0;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::scotchDecomp::scotchDecomp
(
	const dictionary& decompositionDict,
	const polyMesh& mesh
)
:
	decompositionMethod(decompositionDict),
	mesh_(mesh)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::labelList Foam::scotchDecomp::decompose
(
	const pointField& points,
	const scalarField& pointWeights
)
{
	if (points.size() != mesh_.nCells())
	{
		FatalErrorIn
		(
			"scotchDecomp::decompose(const pointField&, const scalarField&)"
		)
			<< "Can use this decomposition method only for the whole mesh"
			<< endl
			<< "and supply one coordinate (cellCentre) for every cell." << endl
			<< "The number of coordinates " << points.size() << endl
			<< "The number of cells in the mesh " << mesh_.nCells()
			<< exit(FatalError);
	}

	// Make Metis CSR (Compressed Storage Format) storage
	//   adjncy	  : contains neighbours (= edges in graph)
	//   xadj(celli) : start of information in adjncy for celli
	labelList adjncy;
	labelList xadj;
	calcCSR(mesh_, adjncy, xadj);

	// Decompose using default weights
	labelList finalDecomp;
	decompose(adjncy, xadj, pointWeights, finalDecomp);

	// Copy back to labelList
	labelList decomp(finalDecomp.size());

	forAll(decomp, i)
	{
		decomp[i] = finalDecomp[i];
	}

	fixCyclics(mesh_, decomp);

	return decomp;
}


Foam::labelList Foam::scotchDecomp::decompose
(
	const labelList& fineToCoarse,
	const pointField& coarsePoints,
	const scalarField& coarseWeights
)
{
	if (fineToCoarse.size() != mesh_.nCells())
	{
		FatalErrorIn
		(
			"parMetisDecomp::decompose(const labelList&, const pointField&)"
		)   << "Size of cell-to-coarse map " << fineToCoarse.size()
			<< " differs from number of cells in mesh " << mesh_.nCells()
			<< exit(FatalError);
	}

	// Make Metis CSR (Compressed Storage Format) storage
	//   adjncy	  : contains neighbours (= edges in graph)
	//   xadj(celli) : start of information in adjncy for celli
	labelList adjncy;
	labelList xadj;
	{
		// Get cellCells on coarse mesh.
		labelListList cellCells;

		calcCellCells
		(
			mesh_,
			fineToCoarse,
			coarsePoints.size(),
			cellCells
		);

		calcCSR(cellCells, adjncy, xadj);
	}

	// Decompose using weights
	labelList finalDecomp;
	decompose(adjncy, xadj, coarseWeights, finalDecomp);

	// Rework back into decomposition for original mesh_
	labelList fineDistribution(fineToCoarse.size());

	forAll(fineDistribution, i)
	{
		fineDistribution[i] = finalDecomp[fineToCoarse[i]];
	}

	fixCyclics(mesh_, fineDistribution);

	return fineDistribution;
}


Foam::labelList Foam::scotchDecomp::decompose
(
	const labelListList& globalCellCells,
	const pointField& cc,
	const scalarField& cWeights
)
{
	if (cc.size() != globalCellCells.size())
	{
		FatalErrorIn
		(
			"scotchDecomp::decompose"
			"(const labelListList&, const pointField&, const scalarField&)"
		)   << "Inconsistent number of cells (" << globalCellCells.size()
			<< ") and number of cell centres (" << cc.size()
			<< ")." << exit(FatalError);
	}


	// Make Metis CSR (Compressed Storage Format) storage
	//   adjncy	  : contains neighbours (= edges in graph)
	//   xadj(celli) : start of information in adjncy for celli

	labelList adjncy;
	labelList xadj;
	calcCSR(globalCellCells, adjncy, xadj);

	// Decompose using weights
	labelList finalDecomp;
	decompose(adjncy, xadj, cWeights, finalDecomp);

	// Copy back to labelList
	labelList decomp(finalDecomp.size());

	forAll(decomp, i)
	{
		decomp[i] = finalDecomp[i];
	}

	fixCyclics(mesh_, decomp);

	return decomp;
}


// ************************************************************************* //
