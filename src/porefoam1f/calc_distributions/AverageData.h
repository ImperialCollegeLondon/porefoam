

const double minZfrac=0.00;
const double maxZfrac=1.;


	
class snapShot_zt 
{
	scalarField data_;
 public:
	//static std::string Names[12];
	
 
    scalar& q(){return data_[0];};
    scalar& vol(){return data_[1];};

    scalar& delPdelX(){return data_[2];};

    scalar& dpdE(){return data_[3];};
    scalar& viscE(){return data_[4];};
    scalar& phiE(){return data_[5];};

    scalar& dpddz(){return data_[6];};
    scalar& viscz(){return data_[7];};
    scalar& phiu(){return data_[8];};
    
    scalar& xDropAvg(){return data_[9];};
    scalar& x1(){return data_[10];};
    scalar& x2(){return data_[11];};       
    


	scalarField& data(){return data_;};	
	
	snapShot_zt();   
};

void writePostProcHeader(std::string fnam, int nSlices)
{
	std::string Names[12]={
	"U", 
	"vol", 
	"delPdelX", 
	"dpdE", 
	"viscE", 
	"phiE", 
	"dpddz", 
	"viscz", 
	"phiu",
	"xDropAvg", 
	"x1", 
	"x2"
	}	;


	std::ofstream out (fnam);

	out << "t ";
	out << "maxMagU ";
	out << "QIn ";
	out << "QOut ";
	out << "Dp ";
	out << "ADarcy ";
	for (int i=0; i<nSlices; ++i) 
	{
	 for(int j=0; j<12; ++j)   out << "S"<<i+1<<"-"<<Names[j]<<" ";
	}
	out << std::endl ;
	out.close();
	Info<<"wrote "<< fnam <<endl;

}

snapShot_zt::snapShot_zt() : data_(12,0.) {	};

class snapShot_t 
{
 public:
    double t;
    scalar maxMagU;
    scalar QIn;
    scalar QOut;
    scalar Dp;
    scalar ADarcy;
    
    List<snapShot_zt> slices;

	snapShot_t(unsigned int numSlices)	{	reset(numSlices);	};

	void reset(unsigned int numSlices)
	{
		t=0.;
		maxMagU=0.;
		QIn=0.;
		QOut=0.;
		Dp=0.;
		ADarcy=0.;
		slices=List<snapShot_zt>(numSlices);
	};
	

	void write(std::ofstream& out)
	{
		out << t<<",";
		out << maxMagU<<",";
		out << QIn<<",";
		out << QOut<<",";
		out << Dp<<",";
		out << ADarcy<<",";
		for (int i=0; i<slices.size(); ++i) 
		{
		 forAll(slices[i].data(), j)   out << slices[i].data()[j]<<",";
		}
		out << std::endl ;
		
	}


};








std::string removesuffix(const std::string& str)
{
  return str.substr(0,str.find_last_of("."));
}

std::string dirname(const std::string& str)
{
  return str.substr(0,str.find_last_of("/\\"));
}

std::string basename(const std::string& str)
{
  return str.substr(str.find_last_of("/\\")+1);
}





