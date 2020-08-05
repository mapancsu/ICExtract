#pragma once
#include <iostream>
#include <direct.h>
#include <list>
#include <set>
#include "PIC.h"
#include <numeric>
using namespace std;

void DDAMS::loadfile(const std::string& filename)
{
	string suffixStr = filename.substr(filename.find_last_of('.') + 1);
	if (suffixStr == "mzML")
	{
		OpenMS::MzMLFile file;
		OpenMS::PeakMap exp;
		file.load(filename, exp);

		std::vector<int> Nr, Nu;
		for (int i = 0; i < exp.size(); i++)
		{
			OpenMS::MSSpectrum MSSpec = exp[i];
			if (MSSpec.getMSLevel() == 1)
			{
				Nr.push_back(i);
			}
			else if (MSSpec.getMSLevel() == 2)
			{
				Nu.push_back(i);
			}

		}
		std::vector < std::vector<int> > Ind;
		Ind.push_back(Nr);
		Ind.push_back(Nu);

		for (int i = 0; i < 2; i++)
		{
			std::vector<int> NN = Ind[i];

			int N = NN.size();
			Eigen::VectorXd rtss;
			std::vector<Eigen::MatrixXd> mss;
			mss.resize(N);
			rtss.resize(N);

			Eigen::VectorXd precurs;
			precurs.resize(N);

			MassScan ms;
			for (int j = 0; j < N; j++)
			{
				const OpenMS::MSSpectrum& MSSpec = exp[NN[j]];
				rtss[j] = MSSpec.getRT();
				if (i == 1)
				{
					precurs[j] = MSSpec.getPrecursors()[0].getMZ();
				}
				
				int NrPeak = MSSpec.size();
				Eigen::MatrixXd spec;
				spec.resize(NrPeak, 2);
				for (int k = 0; k < NrPeak; k++)
				{
					spec(k, 0) = MSSpec[k].getMZ();
					spec(k, 1) = MSSpec[k].getIntensity();

				}
				mss[j] = spec;
				std::vector<double> mzz(spec.col(0).data(), spec.col(0).data() + spec.rows());
			}

			ms.rts = rtss;
			ms.precursors = precurs;
			ms.scans = mss;
			maps.push_back(ms);
		}
	}
	else if (suffixStr == "mzXML")
	{
		OpenMS::MzXMLFile file;
		OpenMS::PeakMap exp;
		file.load(filename, exp);

		std::vector<int> Nr, Nu;
		for (int i = 0; i < exp.size(); i++)
		{
			OpenMS::MSSpectrum MSSpec = exp[i];
			if (MSSpec.getMSLevel() == 1)
			{
				Nr.push_back(i);
			}
			else if (MSSpec.getMSLevel() == 2)
			{
				Nu.push_back(i);
			}

		}
		std::vector < std::vector<int> > Ind;
		Ind.push_back(Nr);
		Ind.push_back(Nu);

		for (int i = 0; i < 2; i++)
		{
			std::vector<int> NN = Ind[i];

			int N = NN.size();
			Eigen::VectorXd rtss;
			std::vector<Eigen::MatrixXd> mss;
			mss.resize(N);
			rtss.resize(N);

			Eigen::VectorXd precurs;
			precurs.resize(N);

			MassScan ms;
			for (int j = 0; j < NN.size(); j++)
			{
				OpenMS::MSSpectrum MSSpec = exp[NN[j]];
				rtss[j] = MSSpec.getRT();
				if (i == 1)
				{
					precurs[j] = MSSpec.getPrecursors()[0].getMZ();
				}
				int NrPeak = MSSpec.size();
				Eigen::MatrixXd spec;
				spec.resize(NrPeak, 2);
				for (int k = 0; k < NrPeak; k++)
				{
					spec(k, 0) = MSSpec[k].getMZ();
					spec(k, 1) = MSSpec[k].getIntensity();

				}
				mss[j] = spec;
			}

			ms.rts = rtss;
			ms.precursors = precurs;
			ms.scans = mss;
			maps.push_back(ms);
		}
	}
}

Eigen::MatrixXd DDAMS::getSpectra(int mapi, int scan)
{
	Eigen::MatrixXd spec;

	if (0 <= mapi < maps.size())
	{
		if (0 <= scan < maps[mapi].rts.size())
		{
			spec = maps[mapi].scans[scan];
		}
	}
	return spec;
}

Eigen::MatrixXd DDAMS::getSpectra(int mapi, double rt)
{
	Eigen::MatrixXd spec;
	if (0 <= mapi < maps.size())
	{
		MassScan map = maps[mapi];
		Eigen::VectorXd rts = getRTs(mapi);

		int scan = std::upper_bound(rts.data(), rts.data() + rts.size(), rt, [](const double & a, const double & b)->bool
		{
			return a <= b;
		}) - rts.data();
		if (0 <= scan - 1 < rts.size())
		{
			spec = map.scans[scan - 1];
		}
	}
	return spec;
}

std::vector<Eigen::MatrixXd> DDAMS::getMS2Spectras(int scan)
{
	std::vector<Eigen::MatrixXd> MS2Specs;

	Eigen::VectorXd MS1rts = maps[0].rts;
	Eigen::VectorXd rts = maps[1].rts;
	double rt_begin = MS1rts[scan];
	double rt_end;
	if (scan == MS1rts.size() - 1)
	{
		rt_end = rts[rts.size() - 1];
	}
	else
	{
		rt_end = MS1rts[scan + 1];
	}

	int lower = std::upper_bound(rts.data(), rts.data() + rts.size(), rt_begin, [](const double & a, const double & b)->bool
	{
		return a <= b;
	}) - rts.data();

	int upper = std::upper_bound(rts.data(), rts.data() + rts.size(), rt_end, [](const double & a, const double & b)->bool
	{
		return a <= b;
	}) - rts.data();

	for (int i = lower; i < upper; i++)
	{
		Eigen::MatrixXd spec = getSpectra(1, i);
		MS2Specs.push_back(spec);
	}

	return MS2Specs;
}

Eigen::VectorXd DDAMS::getPrecursors(int scan)
{
	Eigen::VectorXd Precursors;

	Eigen::VectorXd MS1rts = maps[0].rts;
	Eigen::VectorXd rts = maps[1].rts;
	Eigen::VectorXd precurs = maps[1].precursors;
	double rt_begin = MS1rts[scan];
	double rt_end;
	if (scan == MS1rts.size() - 1)
	{
		rt_end = rts[rts.size() - 1];
	}
	else
	{
		rt_end = MS1rts[scan + 1];
	}

	int lower = std::upper_bound(rts.data(), rts.data() + rts.size(), rt_begin, [](const double & a, const double & b)->bool
	{
		return a <= b;
	}) - rts.data();

	int upper = std::upper_bound(rts.data(), rts.data() + rts.size(), rt_end, [](const double & a, const double & b)->bool
	{
		return a <= b;
	}) - rts.data();

	Precursors.resize(upper-lower);

	for (int i = lower; i < upper; i++)
	{
		Precursors[i-lower] = precurs[i];
	}

	return Precursors;
}

Eigen::VectorXd DDAMS::getRTs(int mapi)
{
	Eigen::VectorXd rts;
	if (0 <= mapi < maps.size())
	{
		rts = maps[mapi].rts;
	}
	return rts;
}

Eigen::MatrixXd DDAMS::getTIC(int mapi)
{
	Eigen::MatrixXd tic;
	if (0 <= mapi < maps.size())
	{
		MassScan map = maps[mapi];
		int scans = map.rts.size();

		tic.resize(scans, 2);
		for (int i = 0; i < scans; i++)
		{
			tic(i, 0) = map.rts(i);
			tic(i, 1) = accumulate(map.scans[i].col(1).data(), map.scans[i].col(1).data() + map.scans[i].rows(), 0);
		}
	}
	return tic;
}

std::vector<Eigen::Vector4d> DDAMS::getRegion(int mapi, double rt_begin, double rt_end, double mz_begin, double mz_end)
{
	std::vector<Eigen::Vector4d> ret;
	if (0 <= mapi < maps.size())
	{
		MassScan map = maps[mapi];
		Eigen::VectorXd rts = map.rts;

		//std::vector<double> rtt(rts.data(), rts.data() + rts.size());

		int lower = std::upper_bound(rts.data(), rts.data() + rts.size(), rt_begin, [](const double & a, const double & b)->bool
		{
			return a <= b;
		}) - rts.data();

		int upper = std::upper_bound(rts.data(), rts.data() + rts.size(), rt_end, [](const double & a, const double & b)->bool
		{
			return a <= b;
		}) - rts.data();

		for (int i = lower; i < upper; i++)
		{
			if (0 <= i < rts.size())
			{
				Eigen::Vector4d point;

				//std::vector<double> mzz(map.scans[i].col(0).data(), map.scans[i].col(0).data() + map.scans[i].rows());

				int lower1 = std::upper_bound(map.scans[i].col(0).data(), map.scans[i].col(0).data() + map.scans[i].col(0).rows(), mz_begin, [](const double & a, const double & b)->bool
				{
					return a <= b;
				}) - map.scans[i].col(0).data();

				int upper1 = std::upper_bound(map.scans[i].col(0).data(), map.scans[i].col(0).data() + map.scans[i].col(0).rows(), mz_end, [](const double & a, const double & b)->bool
				{
					return a <= b;
				}) - map.scans[i].col(0).data();

				for (int j = lower1; j < upper1; j++)
				{
					if (0 <= j < map.scans[i].rows())
					{
						point(0) = map.rts(i);
						point(1) = map.scans[i](j, 0);
						point(2) = map.scans[i](j, 1);
						point(3) = i;
						ret.push_back(point);
					}

				}
			}
		}
	}
	return ret;
}

std::vector<Eigen::MatrixXd> DDA_PIC(DDAMS &ddams, double EER1, float min_int1, int MaxMissedScan, int min_lens)
{
	MassScan  map = ddams.maps[0];
	std::vector<Eigen::MatrixXd> PICs = PIC(map, EER1, MaxMissedScan, min_lens, min_int1);
	return PICs;
}




void DIAMS::loadfile(const std::string& filename)
{
	string suffixStr = filename.substr(filename.find_last_of('.') + 1);
	if (suffixStr == "mzML")
	{
		boost::shared_ptr<OpenMS::ExperimentalSettings> meta = boost::shared_ptr<OpenMS::ExperimentalSettings>(new OpenMS::ExperimentalSettings());
		std::vector<OpenSwath::SwathMap> swathmaps = OpenMS::SwathFile().loadMzML(filename, "./", meta);
		Swaths.resize(swathmaps.size(), 3);
		MassScan ms;
		for (int i = 0; i < swathmaps.size(); i++)
		{
			Swaths(i, 0) = swathmaps[i].lower;
			Swaths(i, 1) = swathmaps[i].center;
			Swaths(i, 2) = swathmaps[i].upper;
			
			if (swathmaps[i].ms1 == true)
			{
				Swaths(i, 0) = -1;
				Swaths(i, 1) = -1;
				Swaths(i, 2) = -1;
			}

			int N = swathmaps[i].sptr->getNrSpectra();
			Eigen::VectorXd rtss;
			std::vector<Eigen::MatrixXd> mss;
			mss.resize(N);
			rtss.resize(N);

			Eigen::VectorXd precurs;
			precurs.resize(N);

			for (int j = 0; j < N; j++)
			{
				OpenSwath::SpectrumPtr spt = swathmaps[i].sptr->getSpectrumById(j);
				OpenSwath::SpectrumMeta meta = swathmaps[i].sptr->getSpectrumMetaById(j);
				rtss[j] = meta.RT;

				Eigen::MatrixXd spec;
				int NrPeak = spt->getMZArray()->data.size();
				spec.resize(NrPeak, 2);
				for (int k = 0; k < NrPeak; k++)
				{
					spec(k, 0) = spt->getMZArray()->data[k];
					spec(k, 1) = spt->getIntensityArray()->data[k];
				}
				mss[j] = spec;
			}
			ms.rts = rtss;
			ms.precursors = precurs;
			ms.scans = mss;
			maps.push_back(ms);
		};
	}
	else if (suffixStr == "mzXML")
	{
		boost::shared_ptr<OpenMS::ExperimentalSettings> meta = boost::shared_ptr<OpenMS::ExperimentalSettings>(new OpenMS::ExperimentalSettings());
		std::vector<OpenSwath::SwathMap> swathmaps = OpenMS::SwathFile().loadMzXML(filename, "./", meta);
		Swaths.resize(swathmaps.size(), 3);
		//maps.resize(swathmaps.size());
		MassScan ms;
		for (int i = 0; i < swathmaps.size(); i++)
		{
			Swaths(i, 0) = swathmaps[i].lower;
			Swaths(i, 1) = swathmaps[i].center;
			Swaths(i, 2) = swathmaps[i].upper;

			if (swathmaps[i].ms1 == true)
			{
				Swaths(i, 0) = -1;
				Swaths(i, 1) = -1;
				Swaths(i, 2) = -1;
			}

			int N = swathmaps[i].sptr->getNrSpectra();
			Eigen::VectorXd rtss;
			std::vector<Eigen::MatrixXd> mss;
			mss.resize(N);
			rtss.resize(N);

			Eigen::VectorXd precurs;
			precurs.resize(N);

			for (int j = 0; j < N; j++)
			{
				OpenSwath::SpectrumPtr spt = swathmaps[i].sptr->getSpectrumById(j);
				OpenSwath::SpectrumMeta meta = swathmaps[i].sptr->getSpectrumMetaById(j);
				rtss[j] = meta.RT;

				Eigen::MatrixXd spec;
				int NrPeak = spt->getMZArray()->data.size();
				spec.resize(NrPeak, 2);
				for (int k = 0; k < NrPeak; k++)
				{
					spec(k, 0) = spt->getMZArray()->data[k];
					spec(k, 1) = spt->getIntensityArray()->data[k];
				}
				mss[j] = spec;
			}
			ms.rts = rtss;
			ms.precursors = precurs;
			ms.scans = mss;
			maps.push_back(ms);
		}
	}
}

void DIAMS::saveSwath(std::string &path)
{
	//char path[1000];
	//_getcwd(path, 1000);
	path += "/swath.txt";
	const char* cpc = path.c_str();
	char* pc = new char[100];
	strcpy(pc, cpc);
	ofstream fout(pc);
	if (fout)
	{
		fout << "SWATH_no " << "Lower " << "Center " << "Upper " << endl;
		for (int i = 0; i < maps.size(); i++)
		{
			fout << i  << " " << "SWATH " << Swaths(i, 0) << " " << Swaths(i, 1) << " " << Swaths(i, 2) << endl;
		};
		fout.close();
	}
}

Eigen::MatrixXd DIAMS::getSpectra(int mapi, int scan)
{
	Eigen::MatrixXd spec;
	
	if (0 <= mapi < maps.size())
	{
		if (0 <= scan < maps[mapi].rts.size())
		{
			spec = maps[mapi].scans[scan];
		}
	}
	return spec;
}

Eigen::MatrixXd DIAMS::getSpectra(int mapi, double rt)
{
	Eigen::MatrixXd spec;
	if (0 <= mapi < maps.size())
	{
		MassScan map = maps[mapi];
		Eigen::VectorXd rts = getRTs(mapi);

		int scan = std::upper_bound(rts.data(), rts.data() + rts.size(), rt, [](const double & a, const double & b)->bool
		{
			return a <= b;
		}) - rts.data();
		if (0 <= scan - 1 < rts.size())
		{
			spec = map.scans[scan - 1];
		}
	}
	return spec;
}

Eigen::VectorXd DIAMS::getRTs(int mapi)
{
	Eigen::VectorXd rts;
	if (0 <= mapi < maps.size())
	{
		rts = maps[mapi].rts;
	}
	return rts;
}

Eigen::MatrixXd DIAMS::getTIC(int mapi)
{
	Eigen::MatrixXd tic;
	if (0 <= mapi < maps.size())
	{
		MassScan map = maps[mapi];
		int scans = map.rts.size();

		tic.resize(scans, 2);
		for (int i = 0; i < scans; i++)
		{
			tic(i, 0) = map.rts(i);
			tic(i, 1) = accumulate(map.scans[i].col(1).data(), map.scans[i].col(1).data() + map.scans[i].rows(), 0);
		}
	}
	return tic;
}

std::vector<Eigen::Vector4d> DIAMS::getRegion(int mapi, double rt_begin, double rt_end, double mz_begin, double mz_end)
{
	std::vector<Eigen::Vector4d> ret;
	if (0 <= mapi < maps.size())
	{
		MassScan map = maps[mapi];
		Eigen::VectorXd rts = map.rts;

		//std::vector<double> rtt(rts.data(), rts.data() + rts.size());
		
		int lower = std::upper_bound(rts.data(), rts.data() + rts.size(), rt_begin, [](const double & a, const double & b)->bool
		{
			return a <= b;
		}) - rts.data();

		int upper = std::upper_bound(rts.data(), rts.data() + rts.size(), rt_end, [](const double & a, const double & b)->bool
		{
			return a <= b;
		}) - rts.data();

		for (int i = lower; i < upper; i++)
		{
			if (0 <= i < rts.size())
			{
				Eigen::Vector4d point;

				//std::vector<double> mzz(map.scans[i].col(0).data(), map.scans[i].col(0).data() + map.scans[i].rows());

				int lower1 = std::upper_bound(map.scans[i].col(0).data(), map.scans[i].col(0).data() + map.scans[i].col(0).rows(), mz_begin, [](const double & a, const double & b)->bool
				{
					return a <= b;
				}) - map.scans[i].col(0).data();

				int upper1 = std::upper_bound(map.scans[i].col(0).data(), map.scans[i].col(0).data() + map.scans[i].col(0).rows(), mz_end, [](const double & a, const double & b)->bool
				{
					return a <= b;
				}) - map.scans[i].col(0).data();

				for (int j = lower1; j < upper1; j++)
				{
					if (0 <= j < map.scans[i].rows())
					{
						point(0) = map.rts(i);
						point(1) = map.scans[i](j, 0);
						point(2) = map.scans[i](j, 1);
						point(3) = i;
						ret.push_back(point);
					}

				}
			}
		}
	}
	return ret;
}

std::vector<Eigen::MatrixXd> DIA_PIC(DIAMS &diams, double EER1, float min_int1, int MaxMissedScan, int min_lens)
{
	MassScan  map = diams.maps[0];
	std::vector<Eigen::MatrixXd> PICs = PIC(map, EER1, MaxMissedScan, min_lens, min_int1);
	return PICs;
}

std::vector<Eigen::MatrixXd> MS2_PIC(DIAMS &diams, int map_ind, float s_rt, float e_rt, double EER, float min_int, int MaxMissedScan, int min_lens)
{
	MassScan  map = diams.maps[map_ind];
	std::vector<Eigen::MatrixXd> PICs = PICByRT(map,s_rt, e_rt, EER, MaxMissedScan, min_lens, min_int);
	return PICs;
}



void CollectSWATH(MassScan &map, Eigen::MatrixXd &Scans, std::vector<Eigen::VectorXi> &Pos, std::vector<int> &Nums, std::vector<double> &rts, std::vector<double> &ints)
{
	int NrScan = map.rts.size();
	int ions_num = 0;
	for (int i = 0; i < NrScan; i++)
	{
		Nums.push_back(map.scans[i].rows());
		ions_num = ions_num + map.scans[i].rows();
		rts.push_back(map.rts(i));
	};

	Scans.resize(ions_num, 3);
	for (int i = 0; i < NrScan; i++)
	{
		for (int j = 0; j < Nums[i]; j++)
		{
			int ind = accumulate(Nums.begin(), Nums.begin() + i, j);
			Scans(ind, 0) = map.rts(i);
			Scans(ind, 1) = map.scans[i](j, 0);
			Scans(ind, 2) = map.scans[i](j, 1);
			Eigen::Vector2i add(i, j);
			Pos.push_back(add);
			ints.push_back(map.scans[i](j, 1));
		};
	};
};

void CollectSWATHByRT(MassScan &map, float s_rt, float e_rt, Eigen::MatrixXd &Scans, std::vector<Eigen::VectorXi> &Pos, std::vector<int> &Nums, std::vector<double> &ints, std::vector<double> &rts)
{

	int NrScan = map.rts.size();
	std::vector<int> indes;
	int ions_num = 0;
	for (int i = 0; i < NrScan; i++)
	{
		double rt = map.rts(i);
		if (s_rt <= rt && e_rt >= rt)
		{
			Nums.push_back(map.scans[i].rows());
			rts.push_back(rt);
			indes.push_back(i);
			ions_num = ions_num + map.scans[i].rows();
		}
	};
	//Scans.resize(Nums.size());
	Scans.resize(ions_num, 3);
	for (int i = 0; i < indes.size(); i++)
	{
		double rt = rts[i];
		int spec_ind = indes[i];
		for (int j = 0; j < Nums[i]; j++)
		{
			int ind = accumulate(Nums.begin(), Nums.begin() + i, j);
			Eigen::Vector2i add(i, j);
			Scans(ind, 0) = rt;
			Scans(ind, 1) = map.scans[spec_ind](j, 0);
			Scans(ind, 2) = map.scans[spec_ind](j, 1);
			Pos.push_back(add);
			ints.push_back(map.scans[spec_ind](j, 1));
		}
	}
}

std::vector<Eigen::MatrixXd> PIC(MassScan &map, double EER, int MaxMissedScan, int min_lens, float min_int)
{
	std::vector<Eigen::MatrixXd> PICs;
	std::set<int> includes;

	Eigen::MatrixXd Scans;
	std::vector<Eigen::VectorXi> Pos;
	std::vector<int> Nums;
	std::vector<double> rts, ints;
	CollectSWATH(map, Scans, Pos, Nums, rts, ints);

	std::vector<int> s_index(ints.size());
	iota(s_index.begin(), s_index.end(), 0);
	std::sort(s_index.begin(), s_index.end(), [&ints](int i1, int i2) {return ints[i1] > ints[i2]; });

//#pragma omp parallel for 
	for (int idx = 0; idx < Scans.rows(); idx++)
	{
		std::list<Eigen::VectorXd> PeakCurve;
		std::set<int> Sincludes;

		if (Scans(s_index[idx], 2) <= min_int) { break; };

		ionPIC(Scans, includes, Pos, Nums, rts, s_index[idx], EER, MaxMissedScan, PeakCurve, Sincludes);
		if (PeakCurve.size() < min_lens) { continue; };

		//if (PICs.size() >= 900)
		//{
		//	std::cout << idx << '\n';
		//	min_lens = 5;
		//}


		Eigen::MatrixXd pic(PeakCurve.size(), 3);
		std::vector<double> ret(PeakCurve.size());
		int jj = 0;
		std::for_each(PeakCurve.begin(), PeakCurve.end(), [&pic, &ret, &jj](const Eigen::VectorXd &v) {
			pic.row(jj) = v;
			ret[jj] = v(2);
			jj++; });

		//Eigen::VectorXd idata = pic.col(2);
		//std::vector<double> tars(idata.data(), idata.data() + idata.size());

		FillMissValue(pic);
		//idata = pic.col(2);
		//std::vector<double> tars1(idata.data(), idata.data() + idata.size());

//#pragma omp critical (result)
		{
			includes.insert(Sincludes.begin(), Sincludes.end());
			PICs.push_back(pic);
		};

		//std::cout << idx << '\n';
	};
	return PICs;
};

std::vector<Eigen::MatrixXd> PICByRT(MassScan &map, float s_rt, float e_rt, double EER, int MaxMissedScan, int min_lens, float min_int)
{
	std::vector<Eigen::MatrixXd> PICs;
	std::set<int> includes;

	Eigen::MatrixXd Scans;
	std::vector<Eigen::VectorXi> Pos;
	std::vector<int> Nums;
	std::vector<double> ints;
	std::vector<double> rts;
	CollectSWATHByRT(map, s_rt, e_rt, Scans, Pos, Nums, ints, rts);

	std::vector<int> s_index(ints.size());
	iota(s_index.begin(), s_index.end(), 0);
	std::sort(s_index.begin(), s_index.end(), [&ints](int i1, int i2) {return ints[i1] > ints[i2]; });

	//#pragma omp parallel for 
	for (int idx = 0; idx < ints.size(); idx++)
	{
		std::list<Eigen::VectorXd> PeakCurve;
		std::set<int> Sincludes;

		if (ints[s_index[idx]] <= min_int) { break; };
		ionPIC(Scans, includes, Pos, Nums, rts, s_index[idx], EER, MaxMissedScan, PeakCurve, Sincludes);

		if (PeakCurve.size() < min_lens) { continue; };

		Eigen::MatrixXd pic(PeakCurve.size(), 3);
		std::vector<double> ret(PeakCurve.size());
		int jj = 0;
		std::for_each(PeakCurve.begin(), PeakCurve.end(), [&pic, &ret, &jj](const Eigen::VectorXd &v) {
			pic.row(jj) = v;
			ret[jj] = v(2);
			jj++; });
		
		FillMissValue(pic);

		//#pragma omp critical (result)
		{
			includes.insert(Sincludes.begin(), Sincludes.end());
			PICs.push_back(pic);
		};
		//std::cout << idx << '\n';
	};
	return PICs;
}

void FillMissValue(Eigen::MatrixXd &pic)
{
	int left = -1, right;
	std::set<int> Fish;
	for (int i = 0; i < pic.rows(); i++)
	{
		std::vector<int> MissP;
		right = -1;
		if (pic(i, 2) > 0.0)
		{
			left = i;
			continue;
		}
		else
		{
			if (Fish.find(i) != Fish.end())
			{
				continue;
			}
		};

		MissP.push_back(i);
		for (int j = i+1; j < pic.rows(); j++)
		{
			if (pic(j, 2) > 0.0)
			{
				right = j;
				break;
			}
			MissP.push_back(j);
		};
		if (left < 0 && right < 0)
		{
			return;
		}
		else if (left < 0 && right >0)
		{
			for (int it = 0; it < MissP.size(); it++)
			{
				pic(MissP[it], 2) = pic(right,2);
				Fish.insert(MissP[it]);
			}
		}
		else if (left > 0 && right < 0)
		{
			for (int it = 0; it < MissP.size(); it++)
			{
				pic(MissP[it], 2) = pic(left,2);
				Fish.insert(MissP[it]);
			}
		}
		else
		{
			double k, b;
			k = (pic(right, 2) - pic(left, 2)) / (MissP.size() + 1);
			b = pic(left, 2) - k*left;
			for (int it = 0; it < MissP.size(); it++)
			{
				pic(MissP[it], 2) = k*MissP[it] + b;
				Fish.insert(MissP[it]);
			}
		}
	}
};

void ionPIC(Eigen::MatrixXd &Scans, std::set<int> &includes, std::vector<Eigen::VectorXi> &Pos, std::vector<int> &Nums, std::vector<double> &rts, int index, double EER, int MaxMissedScan, std::list<Eigen::VectorXd> &PeakCurve, std::set<int> &Sincludes)
{
	//int ids = accumulate(Nums.begin(), Nums.begin() + Pos[index](0),0) + Pos[index](1);
	//std::cout << ids << '\n';

	if (includes.find(index) == includes.end())
	{
		Eigen::VectorXd ion = Scans.row(index);
		Sincludes.insert(index);

		Eigen::VectorXd add;
		add.resize(3);
		PeakCurve.push_back(ion);
		Eigen::VectorXd targetion, closeion, targetion1, closeion1;

		targetion = ion;
		int lmissedScan = 0, rmissedScan = 0, tmissedScan = 0, mzpos = -1;
		double bk = 0.0;
		for (int i = Pos[index](0) - 1; i >= 0; i--)
		{
			mzpos = -1;
			int st = accumulate(Nums.begin(), Nums.begin() + i, 0);
			int en = accumulate(Nums.begin(), Nums.begin() + i + 1, 0);
			FindCloseIon(Scans, st, en, targetion, closeion, mzpos, EER);

			if (mzpos == -1 || st == en)
			{
				add[0] = rts[i]; add[1] = targetion(1); add[2] = bk;
				PeakCurve.push_front(add);
				targetion = add;
				lmissedScan++;
			}
			else
			{
				int mzpos1 = -1;
				st = accumulate(Nums.begin(), Nums.begin() + i + 1, 0);
				en = accumulate(Nums.begin(), Nums.begin() + i + 2, 0);
				FindCloseIon(Scans, st, en, closeion, targetion1, mzpos1, EER);

				if ((mzpos1 == -1 && targetion(2) == 0) || targetion == targetion1)
				{
					if (includes.find(mzpos) == includes.end())
					{
						PeakCurve.push_front(closeion);
						targetion = closeion;
						Sincludes.insert(mzpos);
						lmissedScan = 0;
					}
					else
					{
						add[0] = closeion(0); add[1] = targetion(1); add[2] = bk;
						PeakCurve.push_front(add);
						targetion = add;
						lmissedScan++;
					};
				}
				else
				{
					add[0] = closeion(0); add[1] = targetion(1); add[2] = bk;
					PeakCurve.push_front(add);
					targetion = add;
					lmissedScan++;
				};
			};
			if (lmissedScan == MaxMissedScan)
			{
				break;
			};
		};
		for (int k = 0; k < lmissedScan; k++)
		{
			PeakCurve.pop_front();
		};

		targetion = ion;
		for (int j = Pos[index](0) + 1; j < Nums.size(); j++)
		{
			mzpos = -1;
			int st = accumulate(Nums.begin(), Nums.begin() + j, 0);
			int en = accumulate(Nums.begin(), Nums.begin() + j + 1, 0);
			FindCloseIon(Scans, st, en, targetion, closeion, mzpos, EER);

			if (mzpos == -1 || st==en)
			{
				add[0] = rts[j]; add[1] = targetion(1); add[2] = bk;
				PeakCurve.push_back(add);
				targetion = add;
				rmissedScan++;
			}
			else
			{
				int mzpos1 = -1;
				st = accumulate(Nums.begin(), Nums.begin() + j - 1, 0);
				en = accumulate(Nums.begin(), Nums.begin() + j, 0);
				FindCloseIon(Scans, st, en, closeion, targetion1, mzpos, EER);
				if ((mzpos1 == -1 && targetion(2) == 0) || targetion == targetion1)
				{
					if (includes.find(mzpos) == includes.end())
					{
						PeakCurve.push_back(closeion);
						Sincludes.insert(mzpos);
						targetion = closeion;
						rmissedScan = 0;
					}
					else
					{
						add[0] = closeion(0); add[1] = targetion(1); add[2] = bk;
						PeakCurve.push_back(add);
						targetion = add;
						rmissedScan++;
					};
				}
				else
				{
					add[0] = closeion(0); add[1] = targetion(1); add[2] = bk;
					PeakCurve.push_back(add);
					targetion = add;
					rmissedScan++;
				};
			};
			if (rmissedScan == MaxMissedScan)
			{
				break;
			};
		};
		for (int k = 0; k <rmissedScan; k++)
		{
			PeakCurve.pop_back();
		};
	};
};

void FindCloseIon(Eigen::MatrixXd &Scans, int &st, int &en, Eigen::VectorXd &targetion, Eigen::VectorXd &closeion, int &mzpos, double &EER)
{

	double mz_ds;
	mz_ds = EER;
	for (int ii = st; ii < en; ii++)
	{
		if (targetion[1] - Scans(ii, 1) > EER)
		{
			continue;
		};

		if (Scans(ii, 1) - targetion[1] > EER)
		{
			break;
		};

		if (abs(Scans(ii, 1) - targetion[1]) <= mz_ds)
		{
			mzpos = ii;
			closeion = Scans.row(ii);
			mz_ds = abs(Scans(ii, 1) - targetion[1]);
		};
	};
};
