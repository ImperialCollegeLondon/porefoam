

const double minZfrac=0.00;
const double maxZfrac=1.;



class snapShot_zt
{
	scalarField data_;
 public:


    scalar& alpha(){return data_[0];};
    scalar& q(){return data_[1];};
    scalar& vol(){return data_[2];};
    scalar& f_1(){return data_[3];};
    scalar& dpddz(){return data_[4];};
    scalar& dpcdz(){return data_[5];};
    scalar& dpcdz_1(){return data_[6];};
    scalar& dpddz_1(){return data_[7];};
    scalar& viscz(){return data_[8];};
    scalar& viscz_1(){return data_[9];};
    scalar& phiu(){return data_[10];};
    scalar& phiu_1(){return data_[11];};

    scalar& delPdelZ(){return data_[12];};
    scalar& delPcelZ(){return data_[13];};

    scalar& viscInterf_1(){return data_[14];};



    scalar& viscE(){return data_[15];};
    scalar& viscE_1(){return data_[16];};
    scalar& dpEc(){return data_[17];};
    scalar& dpEc_1(){return data_[18];};
    scalar& dpEd(){return data_[19];};
    scalar& dpEd_1(){return data_[20];};
    scalar& phiE(){return data_[21];};
    scalar& phiE_1(){return data_[22];};




    scalar& Pc(){return data_[24];};
    scalar& zDropAvg(){return data_[25];};
    scalar& zDrop1(){return data_[26];};
    scalar& zDrop2(){return data_[27];};
    scalar& x1(){return data_[28];};
    scalar& x2(){return data_[29];};
    scalar& mu1(){return data_[30];};
    scalar& mu2(){return data_[31];};

    //scalar f_2(){return 1-f_1();};
    //scalar dpcdz_2(){return dpddz()-dpcdz_1();};
    //scalar dpddz_2(){return dpcdz()-dpddz_1();};
    //scalar viscz_2(){return viscz()-dpddz_1();};


	scalarField& data(){return data_;};

	snapShot_zt();
};

void writePostProcHeader(std::string fnam, int nSlices)  {
	std::string Names[32]={
	"alpha",
	"U",
	"vol",
	"f_1",
	"dpddz",
	"dpcdz",
	"dpcdz_1",
	"dpddz_1",
	"viscz",
	"viscz_1",
	"phiu",
	"phiu_1",
	"delPdelZ",
	"delPcelZ",
	"viscInterf_1",
	"viscE",
	"viscE_1",
	"dpEc",
	"dpEc_1",
	"dpEd",
	"dpEd_1",
	"phiE",
	"phiE_1",
	"ZERO",
	"Pc",
	"xDropAvg",
	"xDrop1",
	"xDrop2",
	"x1",
	"x2",
	"mu1",
	"mu2"
	}	;


	std::ofstream out (fnam);

	out << "t ";
	out << "maxMagU ";
	out << "aAvg ";
	out << "aAvgL ";
	out << "aAvgR ";
	for (int i=0; i<3; ++i)    out << "avgUAlpha1_"<<i<<" ";
	for (int i=0; i<3; ++i)    out << "avgUAlpha2_"<<i<<" ";
	out << "QIn ";
	out << "QOut ";
	out << "Dp ";
	out << "Dpc ";
	out << "pcAvg ";
	out << "ADarcy ";
	for (int i=0; i<nSlices; ++i)
	{
	 for(int j=0; j<32; ++j)   out << "S"<<i+1<<"-"<<Names[j]<<" ";
	}
	out << std::endl ;
	out.close();
	Info<<"wrote "<< fnam <<endl;

}

snapShot_zt::snapShot_zt() : data_(32,0.) {	};

class snapShot_t
{
 public:
    double t;
    scalar maxMagU;
    scalar aAvg;
    scalar aAvgL;
    scalar aAvgR;
    vector avgUAlpha1;
    vector avgUAlpha2;
    scalar QIn;
    scalar QOut;
    scalar Dp;
    scalar Dpc;
    scalar pcAvg;
    scalar ADarcy;

    List<snapShot_zt> slices;

	snapShot_t()	{ };
	snapShot_t(unsigned int numSlices)	{	reset(numSlices);	};

	void reset(unsigned int numSlices)
	{
		t=0.;
		maxMagU=0.;
		aAvg=0.;
		aAvgL=0.;
		aAvgR=0.;
		avgUAlpha1=vector(0.,0.,0.);
		avgUAlpha2=vector(0.,0.,0.);
		QIn=0.;
		QOut=0.;
		Dp=0.;
		Dpc=0.;
		pcAvg=0.;
		ADarcy=0.;
		slices=List<snapShot_zt>(numSlices);
	};

	void read(std::ifstream& in)
	{
		in >> t;
		in >> maxMagU;
		in >> aAvg;
		in >> aAvgL;
		in >> aAvgR;
		for (int i=0; i<3; ++i)    in >> avgUAlpha1[i];
		for (int i=0; i<3; ++i)    in >> avgUAlpha2[i];
		in >> QIn;
		in >> QOut;
		in >> Dp;
		in >> Dpc;
		in >> pcAvg;
		in >> ADarcy;
		for (int i=0; i<slices.size(); ++i)
		 forAll( slices[i].data(), j )    in >> slices[i].data()[j];
	}
	void write(std::ofstream& out)
	{
		out << t<<",";
		out << maxMagU<<",";
		out << aAvg<<",";
		out << aAvgL<<",";
		out << aAvgR<<",";
		for (int i=0; i<3; ++i)    out << avgUAlpha1[i]<<",";
		for (int i=0; i<3; ++i)    out << avgUAlpha2[i]<<",";
		out << QIn<<",";
		out << QOut<<",";
		out << Dp<<",";
		out << Dpc<<",";
		out << pcAvg<<",";
		out << ADarcy<<",";
		for (int i=0; i<slices.size(); ++i)
		{
		 forAll( slices[i].data(), j )   out << slices[i].data()[j]<<",";
		}
		out << std::endl ;

	}


};

// not being used anymore, but kept in case
class tzData
{
 public:
    List<snapShot_t> timeData;


    tzData(label size, label numSlices);
    void smoothData()
    {
		Info<<"Smoothing over Time and z using Gaussian kernel";
		for(label i=1; i<timeData.size()-1; ++i)
		for(label j=1; j<timeData[0].slices.size()-1; ++j)
		{
			timeData[i].slices[j].data()=0.125*
				(
					4*timeData[i].slices[j].data()+
					timeData[i+1].slices[j].data()+
					timeData[i-1].slices[j].data()+
					timeData[i].slices[j+1].data()+
					timeData[i].slices[j-1].data()
				);
		};
	};
    void smoothData_z()
    {
		Info<<"Smoothing over Time and z using Gaussian kernel";
		for(label i=1; i<timeData.size()-1; ++i)
		for(label j=1; j<timeData[0].slices.size()-1; ++j)
		{
			timeData[i].slices[j].data()=0.25*
				(
					2*timeData[i].slices[j].data()+
					timeData[i].slices[j+1].data()+
					timeData[i].slices[j-1].data()
				);
		};
	};
    void smoothData_t()
    {
		Info<<"Smoothing over Time and z using Gaussian kernel";
		for(label i=1; i<timeData.size()-1; ++i)
		for(label j=1; j<timeData[0].slices.size()-1; ++j)
		{
			timeData[i].slices[j].data()=0.25*
				(
					2*timeData[i].slices[j].data()+
					timeData[i+1].slices[j].data()+
					timeData[i-1].slices[j].data()
				);
		};
	};

	void read()
	{
		std::ifstream myFile ("data_out_for_plot");
		forAll(timeData,i)
		{
			Info<<i<<endl;
			timeData[i].read(myFile);
		}
	};


	void write(std::string name)
	{
		writePostProcHeader(name+"_header",timeData[0].slices.size());

		std::ofstream myFile (name);
		forAll(timeData,i)
		{
			timeData[i].write(myFile);
			Info<<"wrote "<< name <<endl;
		}
	};



};
tzData::tzData(label size, label numSlices):
	timeData(size)
{
	forAll(timeData, i) timeData[i].reset(numSlices);
};
