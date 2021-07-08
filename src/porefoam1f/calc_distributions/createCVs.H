
	volScalarField VC
	(	IOobject
		(	"VC",  runTime.timeName(),  mesh
		),  mesh,  dimensionedScalar("VC",dimVolume,0.),  "fixedValue"
	);	
	VC.internalField()=mesh.V();




	label iBack = mesh.boundaryMesh().findPatchID("inlet");
	if (iBack < 0)
	{	iBack = mesh.boundaryMesh().findPatchID("Left");
		if (iBack < 0) FatalError	<< "Unable to find patch inlet nor Left" << nl	<< exit(FatalError);
	}

	label iFront = mesh.boundaryMesh().findPatchID("outlet");
	if (iFront < 0)
	{	iFront = mesh.boundaryMesh().findPatchID("Right");
		if (iFront < 0) FatalError	<< "Unable to find patch outlet nor Right" << nl << exit(FatalError);
	}





	IOdictionary postProcDict 
	(	IOobject
		(	"postProcessDict",	"system",	mesh,
			IOobject::READ_IF_PRESENT
	)	);





	scalarField CVBounds1(postProcDict.lookupOrDefault("CVxBounds1", scalarField(1,-1e31)));
	scalarField CVBounds2(postProcDict.lookupOrDefault("CVxBounds2", scalarField(1,1e31)));


	scalarField CVyBounds1(postProcDict.lookupOrDefault("CVyBounds1", scalarField(CVBounds1.size(),-1e31)));
	scalarField CVyBounds2(postProcDict.lookupOrDefault("CVyBounds2", scalarField(CVBounds1.size(),1e31)));
	scalarField CVzBounds1(postProcDict.lookupOrDefault("CVzBounds1", scalarField(CVBounds1.size(),-1e31)));
	scalarField CVzBounds2(postProcDict.lookupOrDefault("CVzBounds2", scalarField(CVBounds1.size(),1e31)));
	wordList cvTags(postProcDict.lookupOrDefault("cvTags", wordList(CVBounds1.size(),"")));

	Info <<"postProcessing x bounds  "<<CVBounds1<<"  -  "<<CVBounds2<<endl;


	labelList CVxlVals=( postProcDict.lookupOrDefault("CVoxelValues", labelList(CVBounds1.size(),1)));

	labelField PPRegions(mesh.V().size(),-1);

  {
	boundBox bbox(mesh.points());

	volScalarField PPRegionsScalar
	(
		IOobject
		(
			"PPRegions",
			runTime.timeName(),
			mesh,
			IOobject::NO_READ,
			IOobject::AUTO_WRITE
		),
		mesh,
		dimensionedScalar("PPRegions",dimless,0.),
		"zeroGradient"
	);



	if(word(postProcDict.lookupOrDefault("relativeBounds", word("true")))==word("true"))
	{
		CVBounds1=bbox.min()[0]+CVBounds1*(bbox.max()[0]-bbox.min()[0]);
		CVBounds2=bbox.min()[0]+CVBounds2*(bbox.max()[0]-bbox.min()[0]);
		CVyBounds1=bbox.min()[1]+CVyBounds1*(bbox.max()[1]-bbox.min()[1]);
		CVyBounds2=bbox.min()[1]+CVyBounds2*(bbox.max()[1]-bbox.min()[1]);
		CVzBounds1=bbox.min()[2]+CVzBounds1*(bbox.max()[2]-bbox.min()[2]);
		CVzBounds2=bbox.min()[2]+CVzBounds2*(bbox.max()[2]-bbox.min()[2]);
	}

	Info<<"\nboundingBox: "<<bbox.min()<<"   -   "<<bbox.max()<<endl;
	Info<<"\nCVxBounds1: "<<CVBounds1<<"   :   "<<CVBounds2<<endl;
	Info<<"CVyBounds1: "<<CVyBounds1<<"   :   "<<CVyBounds2<<endl;
	Info<<"CVzBounds1: "<<CVzBounds1<<"   :   "<<CVzBounds2<<endl<<endl;

	labelList nCVxlVals(CVxlVals.size(),0);
	label nGrowCVs( postProcDict.lookupOrDefault("nGrowCVs", 0) );
	Info <<"nGrowCVs:"<<nGrowCVs<<endl;
	Info <<"CVoxelValues: "<<CVxlVals<<endl;

	if (CVBounds1.size()!=CVxlVals.size())
		Info <<"\n\nWarning in postProcessDict:\n"
		<<" lengths of arrays assigned to CVxBounds1 and CVoxelValues are not equal:  \n"
		<<CVBounds1.size()<<" != "<<CVxlVals.size()
		<<"extra values to CVxBounds1 will be ignored"<<endl<<endl;
	if (CVBounds1.size() < CVxlVals.size())   {Info <<"\n\Error in postProcessDict: CVBounds1.size() < CVxlVals.size()\n"<<endl; exit(-1);}

	fileName imageName( postProcDict.lookupOrDefault("labelImage", fileName("")) );
	Info << "labelImage:" << imageName<<endl;
	voxelImageT<int> CV_img(1,1,1,1); CV_img.X0Ch()-=dbl3(100000000.,100000000.,100000000.);CV_img.dxCh()*=1000000000.;


	if(!imageName.empty()) { readConvertFromHeader(CV_img,runTime.system()+"/../"+imageName); }
	else Info<<"NO labelImage given"<<endl;;

	bool  skipImage = CV_img.nx()<=1;

	//CV_img.write(runTime.system()+"/../voxelI mageT_CV_img.dat");
	int3 nnn = CV_img.size3();
	dbl3 xmin=CV_img.X0();  if(postProcDict.lookupOrDefault("resetX0", word("true"))==word("true")) xmin*=0.;//1e-6 ;      ///.   WARNING IGNORING ORIGIN
	dbl3 dx=CV_img.dx();   //dx*=1e-6;
	std::cout<<"  dx: "<<dx<<";  xmin: "<<xmin<<std::endl;
	int nRejecteds(0);
	int nBndErrs=0;
	const vectorField & C =	mesh.C().internalField();
	forAll(C,c)
	{
		int i=(C[c][0]-xmin[0]+0.*dx[0])/dx[0];
		int j=(C[c][1]-xmin[1]+0.*dx[1])/dx[1];
		int k=(C[c][2]-xmin[2]+0.*dx[2])/dx[2];
		if (k>=0 && k<nnn[2] && j>=0 && j<nnn[1] && i>=0 && i<nnn[0] )
		{
			label CV_ind=CV_img(i,j,k);

			if (nGrowCVs && k>nGrowCVs && k<nnn[2]-nGrowCVs && j>nGrowCVs && j<nnn[1]-nGrowCVs && i>nGrowCVs && i<nnn[0]-nGrowCVs)
			{
				CV_ind = maxNei(CV_img,i,j,k,-nGrowCVs,nGrowCVs);
			}
			if(CV_ind>=0)
			{
				forAll(CVxlVals,iSam)
				{
					if ( (skipImage || CV_ind == CVxlVals[iSam]) && 
						(C[c][0]>=CVBounds1[iSam] && C[c][0]<=CVBounds2[iSam] && 
						 C[c][1]>=CVyBounds1[iSam] && C[c][1]<=CVyBounds2[iSam] && 
						 C[c][2]>=CVzBounds1[iSam] && C[c][2]<=CVzBounds2[iSam] )
					   )
					 {
						PPRegions[c]=CV_ind;
						PPRegionsScalar[c]=CV_ind;
						++(nCVxlVals[iSam]);
				     }
					 else
					 {  //Info<<CV_ind<<"  "; 
					    ++nRejecteds;
					 }
				}
			}
			else				Info<<" ! ";

		}
		else { if (!nBndErrs) Info<<endl<<endl<<"\n\n *** Error out of bounds.... ***   @:"<<C[c]<<endl<<endl<<endl; 
			++nBndErrs; }

	}
	Info<<" \n\nPostProcessing cells: nCollected cells: "<<nCVxlVals<<endl<<",  nRejecteds :"<<nRejecteds<<",  nBndErrs :"<<nBndErrs<<endl<<endl;
	PPRegionsScalar.write();

  }






	std::ofstream data_out_for_plot;
	if (Pstream::master())
	{	data_out_for_plot.open("grads_SP.csv",std::fstream::app);
		writePostProcHeader("grads_SP_hdr.csv", CVBounds1.size());
	}

