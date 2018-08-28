#include "ScreenGenerator.h"
#include <json/json.h>
#include <fstream>
#include <iostream>
#include <vector>
#define _USE_MATH_DEFINES
#include <math.h>
#include <omp.h>

using namespace std;

// nasanen psf unit: cyc/deg
// TODO: 确保这个函数的正确性
double nasanenPSF(double x, double y) {
	// 参数
	static const double c = 0.525;
	static const double L = 250;
	static const double d = 3.91;

	double k = 1 / (c*log(L) + d);
	return 2 * M_PI*k / pow(k*k + 4 * M_PI*M_PI*(x*x + y*y), 1.5);
}

int main() {
	// setup screen info
	double labValue[3 * 16] = { 100, 58, 49, 19, 94, 52, 49, 19, 10,  0,  0,  0, 10,  1,  0,  0,
		                          0,-41, 82, 29, -8,-82, 75, -7,  5,-12, 33, 10,  0,-20, 26,  3,
					              0,-54, -4,-59,105, 29, 59, -4,  4,-21, -1,-14, 28,  9, 11,  1 };

	vector<double> labW = { 1, 0., 0. };

	vector<int> colorCombW = { 0,0,0,1,
	                           0,1,1,1,
	                           0,1,1,0,
	                           1,0,0,1};

	ScreenInfo info;
	info.dpi_ = 300;
	info.viewDistance_ = 9.5;
	info.colorNbr_ = 4; // CMYK
	info.filterSize_ = 21;
	info.screenSize_ = 64;
	info.screenLevel_ = 256;
	info.psf_ = nasanenPSF;
	info.colorCombW_ = colorCombW;
	info.labW_ = labW;
	info.labValue_ = labValue;

	ScreenGenerator screenGenerator(info);
	const vector<vector<int>>& screen = screenGenerator.buildScreen();

	// 保存到JSON文件
	ofstream fout;
	fout.open("bnScreenNew.json");
	vector<string> s = { "W","C","M","CM","Y","CY","MY","CMY","K","CK","MK","CMK","YK","CYK","MYK","CMYK" };
	Json::Value root;
	root["DPI"] = info.dpi_;
	root["screen size"]["height"] = info.screenSize_;
	root["screen size"]["width"] = info.screenSize_;
	root["screen level"] = info.screenLevel_;
	root["color mode"] = "CMYK";
	root["HVS"] = "nasanen";
	for (int i = 0; i < s.size(); ++i) {
		root["color combination weight"][s[i]] = colorCombW[i];
		root["color combination colorimetric values"][s[i]]["L"] = labValue[i];
		root["color combination colorimetric values"][s[i]]["a"] = labValue[s.size() + i];
		root["color combination colorimetric values"][s[i]]["b"] = labValue[2*s.size() + i];
	}
	root["Lab weight"]["L"] = info.labW_[0];
	root["Lab weight"]["a"] = info.labW_[1];
	root["Lab weight"]["b"] = info.labW_[2];

	vector<string> color = { "c","m","y","k" };
	for (int i = 0; i < color.size(); ++i) {
		for (int j = 0; j < screen[i].size(); ++j) {
			root["screen data"][color[i]][j] = screen[i][j];
		}
	}

	fout << root.toStyledString() << endl;

	cout << "\nScreen saved, done!" << endl;
	getchar();
}