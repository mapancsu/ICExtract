# ICExtract
Package for analyzing MS with Python. It provides an  ion-chromatogram-extraction method based on  Greed Nearest Neighbor algorithm  for raw LC-MS dataset (including DDA-MS data and DIA-MS data) effectively and quickly.

# Install

## Required Dependencies

* [Visual Studio Community 2015 with Update 3](http://download.microsoft.com/download/b/e/d/bedddfc4-55f4-4748-90a8-ffe38a40e89f/vs2015.3.com_enu.iso)
* [Anaconda Python 3.6.0 64bit](https://repo.continuum.io/archive/Anaconda3-4.3.1-Windows-x86_64.exe)
* [SWIG 3.0.10](https://sourceforge.net/projects/swig/files/swigwin/swigwin-3.0.10/)
* [CMake 3.7.1](https://cmake.org/files/v3.7/cmake-3.7.1-win64-x64.msi)
* [Eigen 3.3.3](http://bitbucket.org/eigen/eigen/get/3.3.3.zip) 
* [OpenMS 2.3.0](https://github.com/OpenMS)
	
## Download

* Download [ICExtract](https://github.com/mapancsu/ICExtract/archive/master.zip).
* Unzip it into ICExtract directory.

## Compile
* step 1 Compile OpenMS by [OpenMS  Guide](https://github.com/OpenMS/OpenMS/wiki/Building-OpenMS).
* step 2 Compile ICExtract. Open "VS2015 x64 Native Tools Command Prompt" ; Run following commands in the prompt.

	```shell
	cd ICExtract
	mkdir build
	cd build
	cmake .. -G "NMake Makefiles" -DCMAKE_BUILD_TYPE=Release
	nmake
	nmake install
	```

# Usage

* Download ICExtract_python.zip file from [url](https://github.com/mapancsu/ICExtract/releases/tag/ICExtract).
* Upzip ICExtract_python.zip file and go to /python directory.
* Run following Python code fragment to extract ion chromatogram from mzXML or mzML file.

	```python

	from _ICExtract import DDAMS, DDA_PIC, DIAMS, DIA_PIC, MS2_PIC
	import sys
	### for DDA-MS data
	mzfile="DDAfile.mzxml" ## dda-ms file
	mzfile=mzfile.encode(sys.getfilesystemencoding())
	ddams=DDAMS()
	ddams.load(mzfile)
	ics_ms1 = DDA_PIC(ddams, 0.025, 300.0, 10.0, 10.0, 5.0)
	
	### for DIA-MS data
	mzfile="DIAfile.mzxml" ## dia-ms file
	mzfile=mzfile.encode(sys.getfilesystemencoding())
	diams=DIAMS()
	diams.load(mzfile)
	ics_ms1 = DIA_PIC(diams, 0, 0.025, 300.0, 10.0, 10.0, 5.0) ## MS1
	ics_ms2 = MS2_PIC(diams, 1, 0.025, 300.0, 10.0, 10.0, 5.0) ## MS2 
	```

# Contact

For any questions, please contact:  [mapan_spy@163.com](mailto:mapan_spy@163.com)
