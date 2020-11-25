#include "mex.hpp"
#include "mexAdapter.hpp"
#include <iostream>
using namespace std;

class MexFunction : public matlab::mex::Function {
public:
	void operator()(matlab::mex::ArgumentList outputs, matlab::mex::ArgumentList inputs) {

		matlab::data::TypedArray<double> MexX = move(inputs[0]);
		matlab::data::TypedArray<double> MexY = move(inputs[1]);

		matlab::data::TypedArray<double> PathX = move(inputs[2]);
		matlab::data::TypedArray<double> PathY = move(inputs[3]);

		matlab::data::ArrayFactory factory;
		matlab::data::TypedArray<double>  Result = factory.createArray({ 1 }, { calculateStoppingIndex(MexX, MexY, PathX, PathY) });
		outputs[0] = move(Result);
	}

	double calculateStoppingIndex(matlab::data::TypedArray<double> MexX, matlab::data::TypedArray<double> MexY,
		matlab::data::TypedArray<double> PathX, matlab::data::TypedArray<double> PathY) {
		vector<double> vectorPathX, vectorPathY, vectorMexX, vectorMexY;

		for (auto value : MexX)
		{
			vectorMexX.push_back(value);
		}
		for (auto value : MexY)
		{
			vectorMexY.push_back(value);
		}
		for (auto value : PathX)
		{
			vectorPathX.push_back(value);
		}
		for (auto value : PathY)
		{
			vectorPathY.push_back(value);
		}

		for (int i = 0; i < vectorPathX.size() - 1; i++)
		{
			double x1 = vectorPathX[i], y1 = vectorPathY[i];
			double x2 = vectorPathX[i + 1], y2 = vectorPathY[i + 1];

			for (int j = 0; j < vectorMexX.size() - 1; j++)
			{
				double x3 = vectorMexX[j], y3 = vectorMexY[j];
				double x4 = vectorMexX[j + 1], y4 = vectorMexY[j + 1];

				double P = (x4 - x2) * (y4 - y3) - (x4 - x3) * (y4 - y2);
				double D = (x1 - x2) * (y4 - y3) - (x4 - x3) * (y1 - y2);
				double Q = (x1 - x2) * (y4 - y2) - (x4 - x2) * (y1 - y2);
				//intersect not parallel
				if (D != 0) {
					double lam = P / D, eta = Q / D;
					if (!(lam >= 0 && lam <= 1 && eta >= 0 && eta <= 1)) continue;
					return lam * x1 + (1 - lam) * x2;
				}
				//not intersect
				if (P != 0 || Q != 0) continue;
				//parallel
				double t1 = min(std::max(y1, y2), max(y3, y4)),
					t2 = max(min(y1, y2), min(y3, y4)),
					tx2 = max(min(x1, x2), min(x3, x4));
				if (t1 >= t2) return tx2;
			}
		}
		return 1;
	}

};