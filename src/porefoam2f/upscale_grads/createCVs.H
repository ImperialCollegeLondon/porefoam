
	volScalarField VC
	(	IOobject("VC",  runTime.timeName(),  mesh)
		,  mesh,  dimensionedScalar("VC",dimVolume,0.0),  "fixedValue"
	);	
	VC.internalField()=mesh.V();




	label iBack = mesh.boundaryMesh().findPatchID("inlet");
	if (iBack < 0)
	{
		iBack = mesh.boundaryMesh().findPatchID("Left");
		ensure(iBack >=0,  "Unable to find patch inlet or Left", 2);
	}

	label iFront = mesh.boundaryMesh().findPatchID("outlet");
	if (iFront < 0)
	{
		iFront = mesh.boundaryMesh().findPatchID("Right");
		ensure(iFront >=0,  "Unable to find patch outlet or Right", 2);
	}





	IOdictionary postProcessDict 
	(	IOobject
		(	"postProcessDict",	"system",	mesh,
			IOobject::READ_IF_PRESENT
	)	);





	scalarField CVBounds1(postProcessDict.lookupOrDefault("CVxBounds1", scalarField(1,-1.0e31)));
	scalarField CVBounds2(postProcessDict.lookupOrDefault("CVxBounds2", scalarField(1,1.0e31)));

	Info <<"postProcessing x bounds  "<<CVBounds1<<"  -  "<<CVBounds2<<endl;



	
	volScalarField Cx=mesh.C().component(0);



	labelList regionVoxelValues=( postProcessDict.lookupOrDefault("CVoxelValues", labelList(CVBounds1.size(),1)));

	labelField PPRegions(mesh.V().size(),0);


	scalarField CVyBounds1(postProcessDict.lookupOrDefault("CVyBounds1", scalarField(CVBounds1.size(),-1.0e31)));
	scalarField CVyBounds2(postProcessDict.lookupOrDefault("CVyBounds2", scalarField(CVBounds1.size(),1.0e31)));
	scalarField CVzBounds1(postProcessDict.lookupOrDefault("CVzBounds1", scalarField(CVBounds1.size(),-1.0e31)));
	scalarField CVzBounds2(postProcessDict.lookupOrDefault("CVzBounds2", scalarField(CVBounds1.size(),1.0e31)));

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
		dimensionedScalar("PPRegions",dimless,0.0),
		"zeroGradient"
	);



	if(word(postProcessDict.lookupOrDefault("relativeBounds", word("true")))==word("true"))
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

	labelList nregionVoxelValues(regionVoxelValues.size(),0);
	label nGrowCVs( postProcessDict.lookupOrDefault("nGrowCVs", 0) );
	Info <<"nGrowCVs:"<<nGrowCVs<<endl;
	Info <<"regionVoxelValues: "<<regionVoxelValues<<endl;

	if (CVBounds1.size()!=regionVoxelValues.size())
		Info <<"\n\nWarning in postProcessDict:\n"
		<<" lengths of arrays assigned to CVxBounds1 and CVoxelValues are not equal:  \n"
		<<CVBounds1.size()<<" != "<<regionVoxelValues.size()
		<<"extra values to CVxBounds1 will be ignored"<<endl<<endl;
	if (CVBounds1.size() < regionVoxelValues.size())   {Info <<"\n\Error in postProcessDict: CVBounds1.size() < regionVoxelValues.size()\n"<<endl; exit(-1);}

	fileName imageName( postProcessDict.lookupOrDefault("labelImage", fileName("")) );
	Info << "labelImage:" << imageName<<endl;
	voxelImageT<int> CV_img(1,1,1,1); CV_img.X0Ch()-=dbl3(100000000.0,100000000.0,100000000.0);CV_img.dxCh()*=1000000000.0;



#ifdef _USE_MPI_
	MPI_Comm temp_comm;  	MPI_Group temp_group;
	MPI_Comm_group( MPI_COMM_WORLD, &temp_group );
	MPI_Comm_create( MPI_COMM_WORLD, temp_group, &temp_comm);
	int  myrank;	MPI_Comm_rank(temp_comm, &myrank );

	if(myrank==0)
	{
#endif
		if(!imageName.empty()) { readConvertFromHeader(CV_img,runTime.system()+"/../"+imageName); }
		else Info<<"NO labelImage given"<<endl;
#ifdef _USE_MPI_
	}

	if(Pstream::parRun())
	{
		if(MPI_Bcast(&CV_img.nij_, 1, MPI_UNSIGNED_LONG_LONG, 0, temp_comm )) cout<<" pid:"<<myrank<<" Error nij_"<<std::endl;
		if(MPI_Bcast(&CV_img.nnn_, 3, MPI_INT, 0, temp_comm )) cout<<" pid:"<<myrank<<" Error nnn_"<<std::endl;
		if(MPI_Bcast(&CV_img.X0Ch(), 3, MPI_DOUBLE, 0, temp_comm )) cout<<" pid:"<<myrank<<" Error X0Ch"<<std::endl;
		if(MPI_Bcast(&CV_img.dxCh(), 3, MPI_DOUBLE, 0, temp_comm )) cout<<" pid:"<<myrank<<" Error dxCh"<<std::endl;
		//if(MPI_Barrier(temp_comm )) cout<<" pid:"<<myrank<<" Error MPI_Barrier"<<std::endl;
		if(myrank!=0) CV_img.data_.resize(CV_img.nij_*size_t(CV_img.nnn_[2]));

		cout<<" pid:"<<myrank<<" CV_img.data_.size(): "<<CV_img.data_.size()<<std::endl;;
		MPI_Bcast(&CV_img.data_[0], CV_img.data_.size(), MPI_INT, 0, temp_comm );
	}

	MPI_Group_free(&temp_group);
	MPI_Comm_free(&temp_comm);
#endif


	bool  skipImage = CV_img.nx()<=1;

	//CV_img.write(runTime.system()+"/../voxelImag eT_CV_img.dat");
	int3 nnn = CV_img.size3();
	dbl3 xmin=CV_img.X0(); //xmin*=1.0e-6 ;      ///.   WARNING IGNORING ORIGIN
	dbl3 dx=CV_img.dx();   //dx*=1.0e-6;
	int nRejecteds(0);
	const vectorField & C =	mesh.C().internalField();
	forAll(C,c)
	{
		int i=(C[c][0]-xmin[0])/dx[0];//+0.0
		int j=(C[c][1]-xmin[1])/dx[1];//+0.0
		int k=(C[c][2]-xmin[2])/dx[2];//+0.0
		if (k>=0 && k<nnn[2] && j>=0 && j<nnn[1] && i>=0 && i<nnn[0] )
		{
			label CV_ind=CV_img(i,j,k);

			if (nGrowCVs && k>nGrowCVs && k<nnn[2]-nGrowCVs && j>nGrowCVs && j<nnn[1]-nGrowCVs && i>nGrowCVs && i<nnn[0]-nGrowCVs)
			{
				CV_ind = maxNei(CV_img,i,j,k,-nGrowCVs,nGrowCVs);
			}
			if(CV_ind>=0)
			{
				forAll(regionVoxelValues,iSam)
				{
					if ( (skipImage || CV_ind == regionVoxelValues[iSam]) && 
						(C[c][0]>=CVBounds1[iSam] && C[c][0]<=CVBounds2[iSam] && 
						 C[c][1]>=CVyBounds1[iSam] && C[c][1]<=CVyBounds2[iSam] && 
						 C[c][2]>=CVzBounds1[iSam] && C[c][2]<=CVzBounds2[iSam] )
					   )
					 {
						PPRegions[c]=CV_ind;
						PPRegionsScalar[c]=CV_ind;
						++(nregionVoxelValues[iSam]);
				     }
					 else
					 {  //Info<<CV_ind<<"  "; 
					    ++nRejecteds;
					 }
				}
			}
			else				Info<<" ! ";

		}
		else Info<<" Error ";

		
	}
	Info<<" \n\nPostProcessing cells: nCollected cells: "<<nregionVoxelValues<<endl<<",  nRejecteds :"<<nRejecteds<<endl<<endl;
	PPRegionsScalar.write();

  }





	volScalarField clip
	(	IOobject
		(	"clip",	runTime.timeName(),	mesh,
			IOobject::NO_READ,	IOobject::NO_WRITE
		),	mesh,	dimensionedScalar("clip",dimless,0.0),	"fixedValue"
	);

	std::ofstream data_out_for_plot;
	if (Pstream::master())
	{
		data_out_for_plot.open("grads.csv",std::fstream::app);
		writePostProcHeader("grads_hdr.csv", CVBounds1.size());
	}

	{
		const vectorField & C =	mesh.C().internalField();

		label iSam=0;
		clip==0.0;
		forAll(C,c) 
		{
				if (PPRegions[c]==regionVoxelValues[iSam]&& 
					(C[c][0]>=CVBounds1[iSam] && C[c][0]<=CVBounds2[iSam] && 
					 C[c][1]>=CVyBounds1[iSam] && C[c][1]<=CVyBounds2[iSam] && 
					 C[c][2]>=CVzBounds1[iSam] && C[c][2]<=CVzBounds2[iSam] ))
					clip[c]=1.0;
		}
	}
