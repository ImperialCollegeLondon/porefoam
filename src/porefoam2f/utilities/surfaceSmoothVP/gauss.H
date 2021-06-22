

    const pointField & points=surf1.points();
    const pointField& Cf=surf1.faceCentres();

// 1.  to re-implement, to avoid crash with multiply-connected surfaces
    vectorField fNorms(faces.size(),vector(0,0,0));
    forAll(faces,faceI)
    {
        const face & f=faces[faceI];
        forAll(f, pI)  fNorms[faceI] += (newPoints[f.prevLabel(pI)]-newPoints[f[pI]])^(newPoints[f.nextLabel(pI)]-newPoints[f[pI]]);
	}
    fNorms=fNorms/(mag(fNorms)+1e-12);
///. point-normal vectors
    vectorField pNorms(pointPoints.size(),vector(0,0,0));
    forAll(faces,faceI)
    {
        const face & f=faces[faceI];
        forAll(f, pI)
        {
			edge e=f.faceEdge(pI);
			vector Ce=0.5*(newPoints[e[0]] + newPoints[e[1]]);

			pNorms[e[0]] -= (newPoints[e[0]]-Cf[faceI]) ^ (Ce-Cf[faceI]);
			pNorms[e[1]] += (newPoints[e[1]]-Cf[faceI]) ^ (Ce-Cf[faceI]);

        }
    }

	scalarField pAreas(mag(pNorms)+1e-18);

	pNorms=pNorms/(mag(pNorms)+1e-12);



	//vectorField pointsOld(newPoints);


//----------------------------------------------------------------------------------------------------------------------------------------------------------------------

	for(int iter3=0;iter3<kernelRadius;iter3++)
	{
		vectorField newPointsTmp(newPoints);
		forAll(pointPoints, vertI)
		{
			vector avgPos(vector::zero);
			scalar sumWeights(1e-14);

			const labelLoop& neiPoints = pointPoints[vertI];
			forAll(neiPoints, neipI)
			{
				scalar weight=0.05+(pNorms[neiPoints[neipI]] & pNorms[vertI]);
				weight*=weight; weight*=weight;
				avgPos += weight * newPointsTmp[neiPoints[neipI]];
				sumWeights += weight;
			}
			avgPos /= sumWeights;//myEdges.size();

			newPoints[vertI] = (1-relax)*newPointsTmp[vertI] + relax*avgPos;
		}
	}
    

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------
//===========================================================================
Info<< "VP, avgDisp: "; cout.flush();



	vectorField  displacements(newPoints-points);
	Info<< " -> "<<max(mag(displacements)); cout.flush();

      for(int iter3=0;iter3<kernelRadius;iter3++)
      {

        vectorField dispTmp(displacements);

        forAll(pointPoints, vertI)
        {
			vector avgDispV(displacements[vertI]);
            //scalar avgDispS((pNorms[vertI] & dispTmp[vertI]));
            //scalar sumWeightsDisp(1e-14);
            scalar sumWeights(1.0000000001);

            const labelLoop& neiPoints = pointPoints[vertI];
            forAll(neiPoints, neipI)
            {
                scalar weight=0.05+(pNorms[neiPoints[neipI]] & pNorms[vertI]);
                weight*=weight; weight*=weight;
                avgDispV += weight * dispTmp[neiPoints[neipI]];
                //avgDispS += weight * (pNorms[neiPoints[neipI]] & dispTmp[neiPoints[neipI]]);
                sumWeights += weight;
            }
			displacements[vertI]=avgDispV/sumWeights;
			//displacements[vertI]+=0.2*avgDispS/sumWeights*pNorms[vertI];
        }
		Info<< " -> "<<max(mag(displacements)); cout.flush();
    }
    
    newPoints -= 1.0625*displacements;

	Info<< " -> "<<max(mag((newPoints-points)))<< "  "; cout.flush();



	(Info<< " -  ").flush();



  if(iter==iters-1  +1000 )
  {
	
	 vectorField Disp2 = newPoints-points;
	 
	std::ofstream  offnorms("norms.txt");
	offnorms<<"\nPOINT_DATA "<<pointPoints.size()<<std::endl;
	offnorms<<"FIELD attributes 3"<<std::endl;

	offnorms<<"\n\n norms 3 "<<pointPoints.size()<<" float"<<std::endl<<std::endl;
	forAll(pointPoints, vertI)
	{
		if ((vertI+1)%20==0)	offnorms<<"\n";
		offnorms<<pNorms[vertI][0]<<" "<<pNorms[vertI][1]<<" "<<pNorms[vertI][2]<<"   ";    
	}
	
	offnorms<<"\n\n dispVP 3 "<<pointPoints.size()<<" float"<<std::endl<<std::endl;
	forAll(pointPoints, vertI)
	{
		if ((vertI+1)%20==0)	offnorms<<"\n";
		offnorms<<displacements[vertI][0]<<" "<<displacements[vertI][1]<<" "<<displacements[vertI][2]<<"   ";    
	}


	offnorms<<"\n\n DispGauss 3 "<<pointPoints.size()<<" float"<<std::endl<<std::endl;
	forAll(pointPoints, vertI)
	{
		if ((vertI+1)%20==0)	offnorms<<"\n";
		offnorms<<Disp2[vertI][0]<<" "<<Disp2[vertI][1]<<" "<<Disp2[vertI][2]<<"   ";    
	}

	offnorms.close();
	
	
  }

	//surf1.clearGeom();
	//surf1.movePoints(newPoints);

