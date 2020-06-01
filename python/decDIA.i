// Tell swig the name of the module we're creating
%module decDIA

// Pull in the headers from Python itself and from our library
%{
#define SWIG_FILE_WITH_INIT
#define SWIG_PYTHON_STRICT_BYTE_CHAR
#include <Python.h>
#include "PIC.h"
//#include "decdia_export.h"
%}

%include <typemaps.i>
%include <std_string.i>
%include <std_vector.i>
%include <std_list.i>


namespace std {
   %template(BoolVector) vector<bool>;
   %template(IntVector) vector<int>;
   %template(DoubleVector) vector<double>;
   %template(StringVector) vector<string>;
   %template(ConstCharVector) vector<const char*>;
   %template(MassScanVector) vector<MassScan>;
   %template(VV3f) vector<Eigen::Vector3d>;
   %template(VV4f) vector<Eigen::Vector4d>;
   %template(PICVec) vector<Eigen::MatrixXd>;
   %template(PCList) list<Eigen::VectorXd>;
   %template(SWATHPIC) vector<vector<Eigen::MatrixXd>>;
}


// Eigen matrices into Numpy arrays.
%include <eigen.i>
//%include <eigen_add.i>
%eigen_typemaps(vector<vector<Eigen::MatrixXd>>)
%eigen_typemaps(Eigen::VectorXd)
%eigen_typemaps(Eigen::MatrixXd)
%eigen_typemaps(Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>)
%eigen_typemaps(Eigen::Vector3d)
%eigen_typemaps(Eigen::Vector4d)
%eigen_typemaps(Eigen::VectorXd)
%eigen_typemaps(Eigen::MatrixXd)
%eigen_typemaps(Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>)
%eigen_typemaps(Eigen::VectorXi)
%eigen_typemaps(Eigen::Matrix<int, Eigen::Dynamic, 1>)



// Tell swig to build bindings for everything in our library
%include <windows.i>
%include "PIC.h"
//%include "decdia_export.h"

