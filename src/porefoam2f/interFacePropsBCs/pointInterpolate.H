namespace Foam
{
	Foam::pointScalarField pointInterpolate( const surfaceScalarField& alphaf, const pointMesh &  pMesh )
	{
        const fvMesh& msh = alphaf.mesh();
		//const labelListList& pointCells = msh.pointCells();
		const pointField& points = msh.points();

		pointScalarField alphaP
		(
		  IOobject   ( "alphaP",  msh.time().timeName(),  msh ),
		 pMesh,  dimensionedScalar("alphaP",dimless,0.0)
		);
		pointScalarField alphaPSumWeights
		(
		  IOobject  (  "alphaPSumWeights",  alphaf.time().timeName(),  msh  ),
		 pMesh, dimensionedScalar("alphaPSumWeights",dimless,0.0)
		);


        //surfaceScalarField alphaf = fvc::interpolate(alpha1);
		const faceList &  faces=msh.faces();
		const surfaceVectorField & Cf=msh.Cf();




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


		const fvBoundaryMesh& patches = msh.boundary();
		forAll(patches, bI)
		{
          label pStart=patches[bI].patch().start();
          const vectorField & Cfp = Cf.boundaryField()[bI];
          if (!patches[bI].coupled())
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

		forAll(patches, bI)
		{
			label pStart=patches[bI].patch().start();
			const vectorField & Cfp = Cf.boundaryField()[bI];
			const scalarField & alphafp = alphaf.boundaryField()[bI];
			if (patches[bI].coupled())
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

		syncTools::syncPointList(msh, alphaPSumWeights, plusEqOp<scalar>(), 0.0, false);
		syncTools::syncPointList(msh, alphaP, plusEqOp<scalar>(), 0.0, false);


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
        const fvMesh& msh = U.mesh();
        //const labelListList& pointCells = msh.pointCells();
        const pointField& points = msh.points();

        pointVectorField UP
                        (
                          IOobject ( "UP",  U.time().timeName(), msh ),
                         pMesh,  dimensionedVector("UP",U.dimensions(),vector::zero)
                       );
        pointScalarField  UPSumWeights
                        (
                          IOobject  ( "UPSumWeights", U.time().timeName(), msh ),
                         pMesh,  dimensionedScalar("UPSumWeights",dimless,0.0)
                       );


        surfaceVectorField Uf = fvc::interpolate(U);
        const faceList &  faces=msh.faces();
        const surfaceVectorField & Cf=msh.Cf();




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


        const fvBoundaryMesh& patches = msh.boundary();
        forAll(patches, bI)
        {
          label pStart=patches[bI].patch().start();
          const vectorField & Cfp = Cf.boundaryField()[bI];
          const vectorField & Ufp = U.boundaryField()[bI];
          if (patches[bI].coupled())
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

        syncTools::syncPointList(msh, UPSumWeights, plusEqOp<scalar>(), 0.0, false);
        syncTools::syncPointList(msh, UP, plusEqOp<vector>(), vector::zero, false);


        return UP/(UPSumWeights+1.0e-37);
    }
}
