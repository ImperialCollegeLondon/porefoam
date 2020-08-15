namespace Foam
{
	Foam::pointScalarField pointInterpolate( const surfaceScalarField& alphaf, const pointMesh &  pMesh )
	{
        const fvMesh& mesh = alphaf.mesh();
		//const labelListList& pointCells = mesh.pointCells();
		const pointField& points = mesh.points();

		pointScalarField alphaP
		(
		  IOobject   ( "alphaP",  mesh.time().timeName(),  mesh ),
		 pMesh,  dimensionedScalar("alphaP",dimless,0.0)
		);
		pointScalarField alphaPSumWeights
		(
		  IOobject  (  "alphaPSumWeights",  alphaf.time().timeName(),  alphaf.mesh()  ),
		 pMesh, dimensionedScalar("alphaPSumWeights",dimless,0.0)
		);


        //surfaceScalarField alphaf = fvc::interpolate(alpha1);
		const faceList &  faces=mesh.faces();
		const surfaceVectorField & Cf=mesh.Cf();




		forAll(alphaf, fI)
		{
			vector Cfi = Cf[fI];
			scalar alphafi = alphaf[fI];
			const face& ff = faces[fI];
			forAll(ff, ip)
				{
					label pI = ff[ip];
					scalar weight=mag((Cfi - points[pI]));
					alphaPSumWeights[pI] += weight;
					alphaP[pI] += weight * alphafi;
				}
		}


		const fvBoundaryMesh& patches = mesh.boundary();
		forAll(patches, patchI)
		{
          label pStart=patches[patchI].patch().start();
          const vectorField & Cfp = Cf.boundaryField()[patchI];
          if (!patches[patchI].coupled())
          {
            forAll(Cfp, pfI)  //if (pmagdS[pfI]>1.0e-8)
            {
					const face& ff = faces[pStart+pfI];
					forAll(ff, ip)
					{
						label pI = ff[ip];
						alphaPSumWeights[pI] *= 1.0e-32;
						alphaP[pI] *= 1.0e-32;
					}
            }
          }
		}

		forAll(patches, patchI)
		{
			label pStart=patches[patchI].patch().start();
			const vectorField & Cfp = Cf.boundaryField()[patchI];
			const scalarField & alphafp = alphaf.boundaryField()[patchI];
			if (patches[patchI].coupled())
			{
				forAll(Cfp, pfI)  //if (pmagdS[pfI]>1.0e-8)
				{
					vector Cfi = Cfp[pfI];
					scalar alphafi = alphafp[pfI];
					const face& ff = faces[pStart+pfI];
					forAll(ff, ip)
					{
						label pI = ff[ip];
						scalar weight=mag((Cfi - points[pI]));
						alphaPSumWeights[pI] += weight;
						alphaP[pI] += weight * alphafi;
					}
				}
			}
			else
			{
				forAll(Cfp, pfI)  //if (pmagdS[pfI]>1.0e-8)
				{
					vector Cfi = Cfp[pfI];
					scalar alphafi = alphafp[pfI];
					const face& ff = faces[pStart+pfI];
					forAll(ff, ip)
					{
						label pI = ff[ip];
						scalar weight=0.0001*mag((Cfi - points[pI]));
						alphaPSumWeights[pI] += weight;
						alphaP[pI] += weight * alphafi;
					}
				}
			}
		}

		syncTools::syncPointList(mesh, alphaPSumWeights, plusEqOp<scalar>(), 0.0, false);
		syncTools::syncPointList(mesh, alphaP, plusEqOp<scalar>(), 0.0, false);


		return alphaP/(alphaPSumWeights+1.0e-37);
	}
}
namespace Foam
{
    Foam::pointVectorField pointInterpolate
    (
        const volVectorField& U,
        const pointMesh &  pMesh
    )
    {
        const fvMesh& mesh = U.mesh();
        //const labelListList& pointCells = mesh.pointCells();
         //= UP.mesh();
        const pointField& points = mesh.points();

        pointVectorField UP
                        (
                          IOobject ( "UP",  U.time().timeName(),  U.mesh() ),
                         pMesh,  dimensionedVector("UP",U.dimensions(),vector::zero)
                       );
        pointScalarField  UPSumWeights
                        (
                          IOobject  ( "UPSumWeights", U.time().timeName(),  U.mesh() ),
                         pMesh,  dimensionedScalar("UPSumWeights",dimless,0.0)
                       );


        surfaceVectorField Uf = fvc::interpolate(U);
        const faceList &  faces=mesh.faces();
        const surfaceVectorField & Cf=mesh.Cf();




        forAll(Uf, fI)
        {
            vector Cfi = Cf[fI];
            vector Ufi = Uf[fI];
            const face& ff = faces[fI];
            forAll(ff, ip)
				{
					label pI = ff[ip];
					scalar weight=mag((Cfi - points[pI]));
					UPSumWeights[pI] += weight;
					UP[pI] += weight * Ufi;
				}
        }


        const fvBoundaryMesh& patches = mesh.boundary();
        forAll(patches, patchI)
        {
          label pStart=patches[patchI].patch().start();
          const vectorField & Cfp = Cf.boundaryField()[patchI];
          const vectorField & Ufp = U.boundaryField()[patchI];
          if (patches[patchI].coupled())
          {
            forAll(Cfp, pfI)  //if (pmagdS[pfI]>1.0e-8)
            {
					vector Cfi = Cfp[pfI];
					vector Ufi = Ufp[pfI];
					const face& ff = faces[pStart+pfI];
					forAll(ff, ip)
					{
						label pI = ff[ip];
						scalar weight=mag((Cfi - points[pI]));
						UPSumWeights[pI] += weight;
						UP[pI] += weight * Ufi;
					}
            }
          }
          else
          {
            forAll(Cfp, pfI)  //if (pmagdS[pfI]>1.0e-8)
            {
					vector Cfi = Cfp[pfI];
					vector Ufi = Ufp[pfI];
					const face& ff = faces[pStart+pfI];
					forAll(ff, ip)
					{
						label pI = ff[ip];
						scalar weight=0.0001*mag((Cfi - points[pI]));
						UPSumWeights[pI] += weight;
						UP[pI] += weight * Ufi;
					}
            }
          }
        }

        syncTools::syncPointList(mesh, UPSumWeights, plusEqOp<scalar>(), 0.0, false);
        syncTools::syncPointList(mesh, UP, plusEqOp<vector>(), vector::zero, false);


        return UP/(UPSumWeights+1.0e-37);
    }
}
