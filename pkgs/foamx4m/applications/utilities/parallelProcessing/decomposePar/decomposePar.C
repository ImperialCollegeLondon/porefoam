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

Application
	decomposePar

Description
	Automatically decomposes a mesh and fields of a case for parallel
	execution of FOAM.

Usage

	- decomposePar [OPTION]

	@param -cellDist \n
	Write the cell distribution as a labelList for use with 'manual'
	decomposition method and as a volScalarField for post-processing.

	@param -region regionName \n
	Decompose named region. Does not check for existence of processor*.

	@param -copyUniform \n
	Copy any @a uniform directories too.

	@param -fields \n
	Use existing geometry decomposition and convert fields only.

	@param -filterPatches \n
	Remove empty patches when decomposing the geometry.

	@param -force \n
	Remove any existing @a processor subdirectories before decomposing the
	geometry.

	@param -ifRequired \n
	Only decompose the geometry if the number of domains has changed from a
	previous decomposition. No @a processor subdirectories will be removed
	unless the @a -force option is also specified. This option can be used
	to avoid redundant geometry decomposition (eg, in scripts), but should
	be used with caution when the underlying (serial) geometry or the
	decomposition method etc. have been changed between decompositions.

\*---------------------------------------------------------------------------*/

#include "OSspecific.H"
#include "fvCFD.H"
#include "IOobjectList.H"
#include "processorFvPatchFields.H"
#include "domainDecomposition.H"
#include "labelIOField.H"
#include "scalarIOField.H"
#include "vectorIOField.H"
#include "sphericalTensorIOField.H"
#include "symmTensorIOField.H"
#include "tensorIOField.H"

//#include "tetPointFields.H"
#include "pointFields.H"

#include "readFields.H"
#include "fvFieldDecomposer.H"
#include "pointFieldDecomposer.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
	argList::noParallel();
#	include "addRegionOption.H"
	argList::validOptions.insert("cellDist", "");
	argList::validOptions.insert("copyUniform", "");
	argList::validOptions.insert("fields", "");
	argList::validOptions.insert("filterPatches", "");
	argList::validOptions.insert("force", "");
	argList::validOptions.insert("ifRequired", "");

#	include "setRootCase.H"

	word regionName = fvMesh::defaultRegion;
	word regionDir = word::null;

	if (args.optionFound("region"))
	{
		regionName = args.option("region");
		regionDir = regionName;
		Info<< "Decomposing mesh " << regionName << nl << endl;
	}


	bool writeCellDist = args.optionFound("cellDist");
	bool copyUniform = args.optionFound("copyUniform");
	bool decomposeFieldsOnly = args.optionFound("fields");
	bool filterPatches = args.optionFound("filterPatches");
	bool forceOverwrite = args.optionFound("force");
	bool ifRequiredDecomposition = args.optionFound("ifRequired");

#	include "createTime.H"

	Info<< "Time = " << runTime.timeName() << endl;

	// Determine the existing processor count directly
	label nProcs = 0;
	while
	(
		isDir
		(
			runTime.path()
		   /(word("processor") + name(nProcs))
		   /runTime.constant()
		   /regionDir
		   /polyMesh::meshSubDir
		)
	)
	{
		++nProcs;
	}

	// get requested numberOfSubdomains
	label nDomains = 0;
	{
		IOdictionary decompDict
		(
			IOobject
			(
				"decomposeParDict",
				runTime.time().system(),
				regionDir,		  // use region if non-standard
				runTime,
				IOobject::MUST_READ,
				IOobject::NO_WRITE,
				false
			)
		);

		decompDict.lookup("numberOfSubdomains") >> nDomains;
	}

	if (decomposeFieldsOnly)
	{
		// Sanity check on previously decomposed case
		if (nProcs != nDomains)
		{
			FatalErrorIn(args.executable())
				<< "Specified -fields, but the case was decomposed with "
				<< nProcs << " domains"
				<< nl
				<< "instead of " << nDomains
				<< " domains as specified in decomposeParDict"
				<< nl
				<< exit(FatalError);
		}
	}
	else if (nProcs)
	{
		bool procDirsProblem = true;

		if (regionName != fvMesh::defaultRegion)
		{
			decomposeFieldsOnly = false;
			procDirsProblem = false;
		}


		if (ifRequiredDecomposition && nProcs == nDomains)
		{
			// we can reuse the decomposition
			decomposeFieldsOnly = true;
			procDirsProblem = false;
			forceOverwrite = false;

			Info<< "Using existing processor directories" << nl;
		}

		if (forceOverwrite)
		{
			Info<< "Removing " << nProcs
				<< " existing processor directories" << endl;

			// remove existing processor dirs
			// reverse order to avoid gaps if someone interrupts the process
			for (label procI = nProcs-1; procI >= 0; --procI)
			{
				fileName procDir
				(
					runTime.path()/(word("processor") + name(procI))
				);

				rmDir(procDir);
			}

			procDirsProblem = false;
		}

		if (procDirsProblem)
		{
			FatalErrorIn(args.executable())
				<< "Case is already decomposed with " << nProcs
				<< " domains, use the -force option or manually" << nl
				<< "remove processor directories before decomposing. e.g.,"
				<< nl
				<< "    rm -rf " << runTime.path().c_str() << "/processor*"
				<< nl
				<< exit(FatalError);
		}
	}

	Info<< "Create mesh for region " << regionName << endl;
    fvMesh mesh
	(
		IOobject
		(
			regionName,
			runTime.timeName(),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
		)
	);

    domainDecomposition meshDecomp
    (
        mesh,
        IOdictionary
        (
            IOobject
            (
                "decomposeParDict",
                runTime.system(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        )
    );


	// Decompose the mesh
	if (!decomposeFieldsOnly)
	{
        meshDecomp.decomposeMesh(filterPatches);

        meshDecomp.writeDecomposition();

		if (writeCellDist)
		{
            const labelList& procIds = meshDecomp.cellToProc();

			// Write the decomposition as labelList for use with 'manual'
			// decomposition method.
			labelIOList cellDecomposition
			(
				IOobject
				(
					"cellDecomposition",
					mesh.facesInstance(),
					mesh,
					IOobject::NO_READ,
					IOobject::NO_WRITE,
					false
				),
				procIds
			);
			cellDecomposition.write();

			Info<< nl << "Wrote decomposition to "
				<< cellDecomposition.objectPath()
				<< " for use in manual decomposition." << endl;

			// Write as volScalarField for post-processing
			volScalarField cellDist
			(
				IOobject
				(
					"cellDist",
					runTime.timeName(),
                    mesh.dbDir(),
					mesh,
					IOobject::NO_READ,
					IOobject::NO_WRITE
				),
				mesh,
				dimensionedScalar("cellDist", dimless, 0),
				zeroGradientFvPatchScalarField::typeName
			);

			forAll(procIds, celli)
			{
			   cellDist[celli] = procIds[celli];
			}

			cellDist.write();

			Info<< nl << "Wrote decomposition as volScalarField to "
				<< cellDist.name() << " for use in post-processing."
				<< endl;
		}
	}


	// Search for list of objects for this time
	IOobjectList objects(mesh, runTime.timeName());

	// Construct the vol fields
	// ~~~~~~~~~~~~~~~~~~~~~~~~
	PtrList<volScalarField> volScalarFields;
	readFields(mesh, objects, volScalarFields);

	PtrList<volVectorField> volVectorFields;
	readFields(mesh, objects, volVectorFields);

	PtrList<volSphericalTensorField> volSphericalTensorFields;
	readFields(mesh, objects, volSphericalTensorFields);

	PtrList<volSymmTensorField> volSymmTensorFields;
	readFields(mesh, objects, volSymmTensorFields);

	PtrList<volTensorField> volTensorFields;
	readFields(mesh, objects, volTensorFields);


	// Construct the surface fields
	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	PtrList<surfaceScalarField> surfaceScalarFields;
	readFields(mesh, objects, surfaceScalarFields);
	PtrList<surfaceVectorField> surfaceVectorFields;
	readFields(mesh, objects, surfaceVectorFields);
	PtrList<surfaceSphericalTensorField> surfaceSphericalTensorFields;
	readFields(mesh, objects, surfaceSphericalTensorFields);
	PtrList<surfaceSymmTensorField> surfaceSymmTensorFields;
	readFields(mesh, objects, surfaceSymmTensorFields);
	PtrList<surfaceTensorField> surfaceTensorFields;
	readFields(mesh, objects, surfaceTensorFields);


	// Construct the point fields
	// ~~~~~~~~~~~~~~~~~~~~~~~~~~
    const pointMesh& pMesh = pointMesh::New(mesh);

	PtrList<pointScalarField> pointScalarFields;
	readFields(pMesh, objects, pointScalarFields);

	PtrList<pointVectorField> pointVectorFields;
	readFields(pMesh, objects, pointVectorFields);

	PtrList<pointSphericalTensorField> pointSphericalTensorFields;
	readFields(pMesh, objects, pointSphericalTensorFields);

	PtrList<pointSymmTensorField> pointSymmTensorFields;
	readFields(pMesh, objects, pointSymmTensorFields);

	PtrList<pointTensorField> pointTensorFields;
	readFields(pMesh, objects, pointTensorFields);





	// Any uniform data to copy/link?
	fileName uniformDir("uniform");

	if (isDir(runTime.timePath()/uniformDir))
	{
		Info<< "Detected additional non-decomposed files in "
			<< runTime.timePath()/uniformDir
			<< endl;
	}
	else
	{
		uniformDir.clear();
	}

	Info<< endl;

	// Split the fields over processors
    for (label procI = 0; procI < meshDecomp.nProcs(); procI++)
	{
		Info<< "Processor " << procI << ": field transfer" << endl;

		// open the database
		Time processorDb
		(
			Time::controlDictName,
			args.rootPath(),
			args.caseName()/fileName(word("processor") + name(procI))
		);

		processorDb.setTime(runTime);

		// Remove files remnants that can cause horrible problems
		// - mut and nut are used to mark the new turbulence models,
		//   their existence prevents old models from being upgraded
		// 1.6.x merge.  HJ, 25/Aug/2010
		{
			fileName timeDir(processorDb.path()/processorDb.timeName());

			rm(timeDir/"mut");
			rm(timeDir/"nut");
			rm(timeDir/"mut.gz");
			rm(timeDir/"nut.gz");
		}

		// read the mesh
		fvMesh procMesh
		(
			IOobject
			(
				regionName,
				processorDb.timeName(),
				processorDb
            )
        );

        labelIOList pointProcAddressing
        (
            IOobject
            (
                "pointProcAddressing",
                procMesh.facesInstance(),
                procMesh.meshSubDir,
                procMesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );

        labelIOList faceProcAddressing
        (
            IOobject
            (
                "faceProcAddressing",
                procMesh.facesInstance(),
                procMesh.meshSubDir,
                procMesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
			)
		);

		labelIOList cellProcAddressing
		(
			IOobject
			(
				"cellProcAddressing",
				procMesh.facesInstance(),
				procMesh.meshSubDir,
				procMesh,
				IOobject::MUST_READ,
				IOobject::NO_WRITE
			)
		);

		labelIOList boundaryProcAddressing
		(
			IOobject
			(
				"boundaryProcAddressing",
				procMesh.facesInstance(),
				procMesh.meshSubDir,
				procMesh,
				IOobject::MUST_READ,
				IOobject::NO_WRITE
			)
		);

		// FV fields
		if
		(
			volScalarFields.size()
		 || volVectorFields.size()
		 || volSphericalTensorFields.size()
		 || volSymmTensorFields.size()
		 || volTensorFields.size()
		 || surfaceScalarFields.size()
		 || surfaceVectorFields.size()
		 || surfaceSphericalTensorFields.size()
		 || surfaceSymmTensorFields.size()
		 || surfaceTensorFields.size()
		)
		{
			fvFieldDecomposer fieldDecomposer
			(
				mesh,
				procMesh,
				faceProcAddressing,
				cellProcAddressing,
				boundaryProcAddressing
			);

			fieldDecomposer.decomposeFields(volScalarFields);
			fieldDecomposer.decomposeFields(volVectorFields);
			fieldDecomposer.decomposeFields(volSphericalTensorFields);
			fieldDecomposer.decomposeFields(volSymmTensorFields);
			fieldDecomposer.decomposeFields(volTensorFields);

			fieldDecomposer.decomposeFields(surfaceScalarFields);
			fieldDecomposer.decomposeFields(surfaceVectorFields);
			fieldDecomposer.decomposeFields(surfaceSphericalTensorFields);
			fieldDecomposer.decomposeFields(surfaceSymmTensorFields);
			fieldDecomposer.decomposeFields(surfaceTensorFields);
		}


		// Point fields
		if
		(
			pointScalarFields.size()
		 || pointVectorFields.size()
		 || pointSphericalTensorFields.size()
		 || pointSymmTensorFields.size()
		 || pointTensorFields.size()
		)
		{
            const pointMesh& procPMesh = pointMesh::New(procMesh, true);

			pointFieldDecomposer fieldDecomposer
			(
				pMesh,
				procPMesh,
				pointProcAddressing,
				boundaryProcAddressing
			);

			fieldDecomposer.decomposeFields(pointScalarFields);
			fieldDecomposer.decomposeFields(pointVectorFields);
			fieldDecomposer.decomposeFields(pointSphericalTensorFields);
			fieldDecomposer.decomposeFields(pointSymmTensorFields);
			fieldDecomposer.decomposeFields(pointTensorFields);
		}




		// Any non-decomposed data to copy?
		if (uniformDir.size())
		{
			const fileName timePath = processorDb.timePath();

            if (copyUniform || meshDecomp.distributed())
			{
				cp
				(
					runTime.timePath()/uniformDir,
					timePath/uniformDir
				);
			}
			else
			{
				// Link with relative paths
				const string parentPath = string("..")/"..";

				fileName currentDir(cwd());
				chDir(timePath);
				if (!exists(uniformDir))
				{
					ln
					(
						parentPath/runTime.timeName()/uniformDir,
						uniformDir
					);
					chDir(currentDir);
				}
			}
		}
	}




    Info<< "\nEnd.\n" << endl;

    return 0;
}


// ************************************************************************* //
