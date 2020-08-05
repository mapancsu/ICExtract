#include <iostream>
#include <fstream>
#include <array>
#include <numeric>
#include "PIC.h"
#include <OpenMS\FORMAT\MzXMLFile.h>
#include <OpenMS\FORMAT\MzMLFile.h>
#include <OpenMS\KERNEL\MSExperiment.h>
#include <OpenMS\DATASTRUCTURES\String.h>

using namespace std;
using namespace Eigen;

int main() {
	string filename = "G:/DIA-Umpire-data/napedro_L120224_001.mzXML"; 
	string path = "G:/metabo/MTBLS417";
	DIAMS diams;
	diams.loadfile(filename);
	std::vector<Eigen::Vector4d> reg = diams.getRegion(0, 201.0, 220.0, 786.0, 792.0);
	std::vector<Eigen::MatrixXd> PIC1s = DIA_PIC(diams, 0.02, 300, 2, 5);
	std::vector<Eigen::MatrixXd> PIC2s = MS2_PIC(diams, 16, 200.0, 300.0, 0.02, 300, 3, 5);

	diams.saveSwath(path);
	Eigen::MatrixXd cc = diams.getSpectra(1, 5);
	Eigen::MatrixXd cc1 = diams.getSpectra(0, 101.5);
	Eigen::VectorXd rt1 = diams.getRTs(0);
	Eigen::MatrixXd tic = diams.getTIC(1);


	//string filename = "G:/MTBLS417/MS-DIAL/DDA/NEG/PH697338_neg_IDA.mzML";  
	//string filename = "E:/decPIC/Standards/1-8ng.mzXML";
	//DDAMS ddams;
	//ddams.loadfile(filename);
	//std::vector<Eigen::MatrixXd> PIC1s = DDA_PIC(ddams, 0.02, 300, 2, 5);

	//Eigen::MatrixXd cc = ddams.getSpectra(0, 5);
	//Eigen::MatrixXd cc1 = ddams.getSpectra(0, 101.5);
	//std::vector<Eigen::MatrixXd> cc2 = ddams.getMS2Spectras(20);
	//Eigen::VectorXd precurs = ddams.getPrecursors(20);
	//Eigen::VectorXd rt1 = ddams.getRTs(0);
	//Eigen::MatrixXd tic = ddams.getTIC(0);
	//std::vector<Eigen::Vector4d> reg = ddams.getRegion(0, 201.0, 220.0, 786.0, 792.0);

	return 0;
};