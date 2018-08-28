#pragma once

#include <vector>
#include <functional>
// generate dbs screen
using namespace std;

class ScreenInfo {
public:
	int dpi_;
	double viewDistance_;
	int colorNbr_;
	int filterSize_;
	int screenSize_;
	int screenLevel_;
	vector<int> colorCombW_;
	vector<double> labW_;
	double *labValue_;
	function<double(double, double)> psf_;   // point spread function
};

class ScreenGenerator
{
public:
	ScreenGenerator(const ScreenInfo& info, int maxIterTimes=50);
	~ScreenGenerator();

	vector<vector<int>> buildScreen();

private:
	int colorNbr_;
	vector<double> colorCombW_;   // ����ɫ��ϵ�Ȩ��
	vector<double> labW_;         // l,a,b��Ȩ��
	vector<bool> finished_;
	vector<int*> dotProfiles_;
	const ScreenInfo& screenInfo_;
	int maxIterTimes_;            // ����������
	double **cpp_;

	void generateLevel(int l);
};

