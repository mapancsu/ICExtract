#ifdef ICEXTRACT_STATIC_DEFINE
#  define ICEXTRACT_EXPORT
#  define ICEXTRACT_NO_EXPORT
#else
#  ifndef ICEXTRACT_EXPORT
#    ifdef ICEXTRACT_EXPORTS
/* We are building this library */
#      ifdef _MSC_VER
#      		define ICEXTRACT_EXPORT __declspec(dllexport)
#      else
#			define ICEXTRACT_EXPORT
#      endif
#    else
/* We are using this library */
#      ifdef _MSC_VER
#      define ICEXTRACT_EXPORT __declspec(dllimport)
#      else
#			define ICEXTRACT_EXPORT
#      endif
#    endif
#  endif

#  ifndef ICEXTRACT_NO_EXPORT
#    define ICEXTRACT_NO_EXPORT 
#  endif
#endif

#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/FORMAT/MzXMLFile.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/ANALYSIS/OPENSWATH/OPENSWATHALGO/DATAACCESS/SwathMap.h>
#include <OpenMS/FORMAT/SwathFile.h>
#include <OpenMS/CONCEPT/ClassTest.h>

struct MassScan
{ 
public:
	Eigen::VectorXd rts;
	Eigen::VectorXd precursors;
	std::vector<Eigen::MatrixXd> scans;
	std::vector<std::shared_ptr<MassScan> > childs;
};

class ICEXTRACT_EXPORT DDAMS
{

public:
	DDAMS() {}
	~DDAMS() {}

	//std::vector<OpenSwath::SwathMap> maps;
	std::vector<MassScan> maps;

	void loadfile(const std::string& filename);

	Eigen::MatrixXd getSpectra(int mapi, int scan);
	Eigen::MatrixXd getSpectra(int mapi, double rt);

	std::vector<Eigen::MatrixXd> getMS2Spectras(int scan);
	Eigen::VectorXd getPrecursors(int scan);

	Eigen::MatrixXd getTIC(int mapi);
	Eigen::VectorXd getRTs(int mapi);

	std::vector<Eigen::Vector4d> getRegion(int mapi, double rt_begin, double rt_end, double mz_begin, double mz_end);

};

class ICEXTRACT_EXPORT DIAMS  
{

public:
	DIAMS() {}
	~DIAMS() {}

	//std::vector<OpenSwath::SwathMap> maps;
	std::vector<MassScan> maps;

	Eigen::MatrixXd Swaths;

	void loadfile(const std::string& filename);

	Eigen::MatrixXd getSpectra(int mapi, int scan);
	Eigen::MatrixXd getSpectra(int mapi, double rt);

	Eigen::MatrixXd getTIC(int mapi);
	
	Eigen::VectorXd getRTs(int mapi);

	std::vector<Eigen::Vector4d> getRegion(int mapi, double rt_begin, double rt_end, double mz_begin, double mz_end);

	void saveSwath(std::string& path); //std::string& path


};

void ICEXTRACT_EXPORT CollectSWATH(MassScan &map, Eigen::MatrixXd &Scans, std::vector<Eigen::VectorXi> &Pos, std::vector<int> &Nums,  std::vector<double> &rts, std::vector<double> &ints);
void ICEXTRACT_EXPORT CollectSWATHByRT(MassScan &map, float s_rt, float e_rt, Eigen::MatrixXd &Scans, std::vector<Eigen::VectorXi> &Pos, std::vector<int> &Nums, std::vector<double> &ints, std::vector<double> &rts);
std::vector<Eigen::MatrixXd> ICEXTRACT_EXPORT PIC(MassScan &map, double EER, int MaxMissedScan, int min_lens, float min_int);
std::vector<Eigen::MatrixXd> ICEXTRACT_EXPORT PICByRT(MassScan &map, float s_rt, float e_rt, double EER, int MaxMissedScan, int min_lens, float min_int);
void ICEXTRACT_EXPORT ionPIC(Eigen::MatrixXd &Scans, std::set<int> &includes, std::vector<Eigen::VectorXi> &Pos, std::vector<int> &Nums, std::vector<double> &rts, int index, double EER, int MaxMissedScan, std::list<Eigen::VectorXd> &PeakCurve, std::set<int> &Sincludes);
void ICEXTRACT_EXPORT FindCloseIon(Eigen::MatrixXd &Scans, int &st, int &en,  Eigen::VectorXd &targetion, Eigen::VectorXd &closeion, int &mzpos, double &EER);
void ICEXTRACT_EXPORT FillMissValue(Eigen::MatrixXd &pic);
std::vector<Eigen::MatrixXd>  ICEXTRACT_EXPORT DDA_PIC(DDAMS &ddams, double EER1, float min_int1, int MaxMissedScan, int min_lens);
std::vector<Eigen::MatrixXd>  ICEXTRACT_EXPORT DIA_PIC(DIAMS &diams, double EER1, float min_int1, int MaxMissedScan, int min_lens);
std::vector<Eigen::MatrixXd>  ICEXTRACT_EXPORT MS2_PIC(DIAMS &diams, int map_ind, float s_rt, float e_rt, double EER, float min_int, int MaxMissedScan, int min_lens);
