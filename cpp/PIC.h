#ifdef DECDIA_STATIC_DEFINE
#  define DECDIA_EXPORT
#  define DECDIA_NO_EXPORT
#else
#  ifndef DECDIA_EXPORT
#    ifdef decDIA_EXPORTS
/* We are building this library */
#      ifdef _MSC_VER
#      		define DECDIA_EXPORT __declspec(dllexport)
#      else
#			define DECDIA_EXPORT
#      endif
#    else
/* We are using this library */
#      ifdef _MSC_VER
#      define DECDIA_EXPORT __declspec(dllimport)
#      else
#			define DECDIA_EXPORT
#      endif
#    endif
#  endif

#  ifndef DECDIA_NO_EXPORT
#    define DECDIA_NO_EXPORT 
#  endif
#endif

//#pragma once
//#include <iostream>
//#include <map>
//#include <functional>
//#include <vector>
//#include <memory>
//#include <string>
//#include <numeric>
//#include <Eigen/Core>
//#include <decdia_export.h>

#include <OpenMS\KERNEL\MSSpectrum.h>
#include <OpenMS\FORMAT\MzXMLFile.h>
#include <OpenMS\FORMAT\MzMLFile.h>
#include <OpenMS/ANALYSIS/OPENSWATH/OPENSWATHALGO/DATAACCESS/SwathMap.h>
#include <OpenMS/FORMAT/SwathFile.h>
#include <OpenMS/CONCEPT/ClassTest.h>

#include "decdia_export.h"

struct DECDIA_EXPORT MassScan
{
public:
	Eigen::VectorXd rts;
	Eigen::VectorXd precursors;
	std::vector<Eigen::MatrixXd> scans;
	//Eigen::VectorXf id;
	//float precursor_mz;
	//float RT;
	//float BIC;
	//float TIC;
	std::vector<std::shared_ptr<MassScan> > childs;
};

class DECDIA_EXPORT DIAMS
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

class DECDIA_EXPORT DDAMS
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

//void DECDIA_EXPORT CollectData(OpenMS::MSExperiment &exp, std::vector<Eigen::MatrixXd> &Scans, std::vector<Eigen::VectorXi> &Pos, std::vector<int> &Nums, std::vector<double> &ints, std::vector<double> &rts);
void DECDIA_EXPORT CollectSWATH(MassScan &map, Eigen::MatrixXd &Scans, std::vector<Eigen::VectorXi> &Pos, std::vector<int> &Nums,  std::vector<double> &rts, std::vector<double> &ints);
void DECDIA_EXPORT CollectSWATHByRT(MassScan &map, float s_rt, float e_rt, Eigen::MatrixXd &Scans, std::vector<Eigen::VectorXi> &Pos, std::vector<int> &Nums, std::vector<double> &ints, std::vector<double> &rts);
std::vector<Eigen::MatrixXd> DECDIA_EXPORT PIC(MassScan &map, double EER, int MaxMissedScan, int min_lens, float min_int);
std::vector<Eigen::MatrixXd> DECDIA_EXPORT PICByRT(MassScan &map, float s_rt, float e_rt, double EER, int MaxMissedScan, int min_lens, float min_int);
void DECDIA_EXPORT ionPIC(Eigen::MatrixXd &Scans, std::set<int> &includes, std::vector<Eigen::VectorXi> &Pos, std::vector<int> &Nums, std::vector<double> &rts, int index, double EER, int MaxMissedScan, std::list<Eigen::VectorXd> &PeakCurve, std::set<int> &Sincludes);
//void DECDIA_EXPORT FindCloseIon(Eigen::MatrixXd &closeScan, Eigen::VectorXd &closeion, Eigen::VectorXd &targetion, int &mzpos, double EER);
void DECDIA_EXPORT FindCloseIon(Eigen::MatrixXd &Scans, int &st, int &en,  Eigen::VectorXd &targetion, Eigen::VectorXd &closeion, int &mzpos, double &EER);
void DECDIA_EXPORT FillMissValue(Eigen::MatrixXd &pic);
std::vector<Eigen::MatrixXd>  DECDIA_EXPORT DIA_PIC(DIAMS &diams, double EER1, float min_int1, int MaxMissedScan, int min_lens);
std::vector<Eigen::MatrixXd>  DECDIA_EXPORT MS2_PIC(DIAMS &diams, int map_ind, float s_rt, float e_rt, double EER, float min_int, int MaxMissedScan, int min_lens);

std::vector<Eigen::MatrixXd>  DECDIA_EXPORT DDA_PIC(DDAMS &ddams, double EER1, float min_int1, int MaxMissedScan, int min_lens);
