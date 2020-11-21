/*
* arrayProduct.cpp - example in MATLAB External Interfaces
*
* Multiplies an input scalar (multiplier)
* times a MxN matrix (inMatrix)
* and outputs a MxN matrix (outMatrix)
*
* Usage : from MATLAB
*         >> outMatrix = arrayProduct(multiplier, inMatrix)
*
* This is a C++ MEX-file for MATLAB.
* Copyright 2017 The MathWorks, Inc.
*
*/

#include "mex.hpp"
#include "mexAdapter.hpp"
#include <iostream>
using namespace std;

class MexFunction : public matlab::mex::Function {
public:
	void operator()(matlab::mex::ArgumentList outputs, matlab::mex::ArgumentList inputs) {
		matlab::data::TypedArray<double> Mex = std::move(inputs[0]);
		matlab::data::TypedArray<double> Path = std::move(inputs[1]);

		matlab::data::Array Result = calculateStoppingIndex(Mex, Path);
		outputs[0] = std::move(Result);
	}

	matlab::data::Array calculateStoppingIndex(matlab::data::TypedArray<double> Mex, matlab::data::TypedArray<double> Path) {
		double tmpResult=1;
		vector<double> MexXPoint, MexYPoint, PathXPoint, PathYPoint;
		matlab::data::ArrayFactory factory;

		size_t n = Mex.getNumberOfElements();
		n = n / 2;
		for (auto i : Mex)
		{
			if (n > 0) {
				MexXPoint.push_back(i);
				n--;
			}
			else
			{
				MexYPoint.push_back(i);
			}
		}
		n = Path.getNumberOfElements();
		n = n / 2;
		for (auto j : Path)
		{

			if (n > 0) {
				PathXPoint.push_back(j);
				n--;
			}
			else
			{
				PathYPoint.push_back(j);
			}
		}

		for (size_t j = 0; j < Path.getNumberOfElements() / 2 - 1; j++)
		{
			for (size_t i = 0; i < Mex.getNumberOfElements() / 2 - 1; i++)
			{
				double x1 = PathXPoint[j];
				double y1 = PathYPoint[j];
				double x2 = PathXPoint[j + 1];
				double y2 = PathYPoint[j + 1];

				double x3 = MexXPoint[i];
				double y3 = MexYPoint[i];
				double x4 = MexXPoint[i + 1];
				double y4 = MexYPoint[i + 1];

				double P = (x4 - x2) * (y4 - y3) - (x4 - x3) * (y4 - y2);
				double D = (x1 - x2) * (y4 - y3) - (x4 - x3) * (y1 - y2);
				double Q = (x1 - x2) * (y4 - y2) - (x4 - x2) * (y1 - y2);
				if (D != 0) {
					double lam = P / D;
					double eta = Q / D;
					if (!(lam >= 0 && lam <= 1 && eta >= 0 && eta <= 1))
					{
						continue;
					}
					tmpResult = lam * x1 + (1 - lam) * x2;
					goto Output;
				}
				if (P != 0 || Q != 0)
				{
					continue;
				}

				vector<double> t1, t2;
				t1.push_back(x1);
				t1.push_back(x2);
				sort(t1.begin(), t1.end());
				t2.push_back(x3);
				t2.push_back(x4);
				sort(t2.begin(), t2.end());

				if (t1[1] < t2[0] || t2[1] < t1[0])
				{
					continue;
				}
				if (t1[0] > t2[0])
				{
					tmpResult = t1[0];
					goto Output;
				}
				else
				{
					tmpResult = t2[0];
					goto Output;
				}
			}

		}

	Output:matlab::data::TypedArray<double>  Result = factory.createArray({ 1 },
		{ tmpResult });
	return Result;
	}
};
