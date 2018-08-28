#include "ScreenGenerator.h"
#include <memory>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <stdio.h>
#include <assert.h>
#include <string.h>

#define _USE_MATH_DEFINES
#include <math.h>

using namespace std;

#define WRAP_AROUND(v,s) if(v<0) v=v+s; if(v>s-1) v=v-s
#define USE_CORRECTION

ScreenGenerator::ScreenGenerator(const ScreenInfo& info, int maxIterTimes):
	maxIterTimes_(maxIterTimes),
	screenInfo_(info)
{
	// 产生随机数种子
	srand((unsigned)time(NULL));
}

ScreenGenerator::~ScreenGenerator()
{
	// 释放内存
	delete[] cpp_;
	for (int* dp : dotProfiles_) {
		delete[] dp;
	}
}

vector<vector<int>> ScreenGenerator::buildScreen()
{
	int fs = screenInfo_.filterSize_;
	int colorCombNbr = int(pow(2, screenInfo_.colorNbr_));
	// sanity check
	assert(fs % 2 != 0);                 // filter size 必须是奇数
	assert(screenInfo_.colorCombW_.size() == colorCombNbr);

	int sz = screenInfo_.screenSize_;

	// 准备工作
	// 首先生成0和screenLevel_对应的等级的dot profile
	finished_.resize(screenInfo_.screenLevel_ + 1);  // set up finished flag for each level
	std::fill(finished_.begin(), finished_.end(), false);
	dotProfiles_.resize(screenInfo_.screenLevel_ + 1);

	int *firstProfile = new int[sz*sz];
	memset(firstProfile, 0, sizeof(int)*sz*sz);
	dotProfiles_[0] = firstProfile;
	finished_[0] = true;

	int *lastProfile = new int[sz*sz];
	int maxLevelValue = int(pow(2, screenInfo_.colorNbr_)) - 1;
	for (int i = 0; i < sz*sz; ++i) {
		lastProfile[i] = maxLevelValue;
	}
	dotProfiles_[screenInfo_.screenLevel_] = lastProfile;
	finished_[screenInfo_.screenLevel_] = true;

	// 生成filter
	printf("\nConstructing filter...\n");
	vector<double> filter(fs*fs);
	int centerIdx = fs / 2;
	double scale = 180 / (M_PI*screenInfo_.viewDistance_);
	double scale2 = scale * scale;
	for (int i = 0; i < fs; ++i) {
		for (int j = 0; j < fs; ++j) {
			// PSF的单位是cyc/deg, 需要转化为 cyc/inch
			double x = (i - centerIdx) / double(screenInfo_.dpi_) * scale;  // 为了得到正确的结果，强制类型转换
			double y = (j - centerIdx) / double(screenInfo_.dpi_) * scale;
			filter[i*fs + j] = scale2 * screenInfo_.psf_(x, y) / screenInfo_.dpi_ / screenInfo_.dpi_; // 一定要除以dpi^2
		}
	}
	printf("done!\n"
		   "filter[center]=%f,\n"
		   "filter[0][0]=%f,\n"
		   "filter[filterSize-1][filterSize-1]=%f\n", filter[fs*fs / 2], filter[0], filter[fs*fs - 1]);

	// 计算cpp
	printf("\nCalculation Cpp...\n");
	int cppSize = fs * 2 - 1;
	double *cppFlat = new double[cppSize*cppSize];
	memset(cppFlat, 0, cppSize*cppSize*sizeof(double)); 	// 全部初始化为0
	cpp_ = new double* [cppSize];
	for (int i = 0; i < cppSize; ++i) {
		cpp_[i] = cppFlat + i * cppSize;
	}
	for (int i = 0; i < cppSize; ++i) {
		for (int j = 0; j < cppSize; ++j) {
			for (int k = 0; k < fs; k++) {
				for (int l = 0; l < fs; l++) {
					int nk = fs - 1 - i + k;
					int nl = fs - 1 - j + l;
					if (nk < fs && nk >= 0 && nl < fs && nl >= 0) {
						cpp_[i][j] += filter[k*fs+l]*filter[nk*fs+nl];
					}
				}
			}
		}
	}
	// 验证cpp的正确性
	printf("done!\n"
		"Cpp[center]=%f,\n"
		"Cpp[0][0]=%f,\n"
		"Cpp[CppSize-1][CppSize-1]=%f\n", cpp_[fs-1][fs-1], cpp_[0][0], cpp_[cppSize-1][cppSize-1]);

	printf("\nGenerating random process sequence...\n");
	// 灰度等级随机选择
	vector<int> randomSeq(screenInfo_.screenLevel_ - 1);
	for (int i = 0; i < randomSeq.size(); ++i)
		randomSeq[i] = i+1;
	for (int i = randomSeq.size()-1; i >= 1; --i) {
		swap(randomSeq[i], randomSeq[rand() % i]);
	}
	printf("done!\n");

	printf("\nStart dot profiles constructing...\n");
	for (int i = 0; i < screenInfo_.screenLevel_-1; ++i) {
		generateLevel(randomSeq[i]);
	}
	printf("done!\n");

	//  build screen from dot profiles
	printf("\nBuilding screen from dot profiles...\n");
	vector<vector<int>> screen;
	for (int i = 0; i < screenInfo_.colorNbr_; ++i) {
		int mask = 1 << i;
		vector<int> monoScreen(sz*sz, screenInfo_.screenLevel_);
		for (int j = 0; j < sz*sz; ++j) {
			for (int k = 0; k < dotProfiles_.size(); ++k) {
				if (dotProfiles_[k][j] & mask) {
					monoScreen[j] = k;
					break;
				}
			}
		}
		screen.push_back(monoScreen);
	}
	printf("done!\n");

	delete[] cppFlat;
	return screen;
}

void ScreenGenerator::generateLevel(int currentLevel)
{
	// 查找是否有前一等级和后一等级已经完成
	int pre, next;
	for (pre = currentLevel - 1; pre >= 0; --pre) {
		if (finished_[pre])
			break;
	}
	for (next = currentLevel + 1; next < screenInfo_.screenLevel_+1; ++next) {
		if (finished_[next])
			break;
	}

	printf("current level is:%d\n", currentLevel);
	printf("pre is:%d\n", pre);
	printf("next is:%d\n", next);

	int sz = screenInfo_.screenSize_;
	int cppSize = screenInfo_.filterSize_ * 2 - 1;
	double cppCenterValue = cpp_[cppSize / 2][cppSize / 2];

	// 计算本次需要增加的墨点数（单种颜色）
	int sc = sz*sz;
	int dotsPerLevel = sc / screenInfo_.screenLevel_;
	int newDots = dotsPerLevel * (currentLevel - pre);

	// 随机生成初始墨点分布
	int *dotProfile = new int[sc];
	memcpy(dotProfile, dotProfiles_[pre], sizeof(int)*sc);
	dotProfiles_[currentLevel] = dotProfile;
	for (int c = 0; c < screenInfo_.colorNbr_; ++c) {
		int mask = 1 << c;
		int i = 0;
		while (i < newDots) {
			int pos = rand() % sc;
			//  stack constraint
			if (!(dotProfile[pos] & mask) &&
				(dotProfiles_[next][pos] & mask)) {
				dotProfile[pos] |= mask;
				i++;
			}
		}
	}

	printf("generate Init dot profile done!\n");

	// 提高访问速度
	int **dotProfile2D = new int*[sz];
	for (int c = 0; c < sz; ++c) {
		dotProfile2D[c] = dotProfile + c * sz;
	}
	// 计算cpe
	int colorCombNbr = int(pow(2, screenInfo_.colorNbr_));
	double *cpe = new double[colorCombNbr * 3 * sc];
	vector<double> error(sc, 0);

	for (int i = 0; i < colorCombNbr; ++i){
		for (int j = 0; j < 3; ++j) {
			double *lrValue = screenInfo_.labValue_+j*colorCombNbr;
			double *p = cpe + (i * 3 + j)*sc;
			double sum = 0;
			for (int x = 0; x < sz; ++x) {
				for (int y = 0; y < sz; ++y) {
					sum += error[x*sz+y] = lrValue[dotProfile2D[x][y] & i];
				}
			}
			double mean = sum / sc;
			for (int x = 0; x < sz; ++x) {
				for (int y = 0; y < sz; ++y) {
					error[x*sz+y] -= mean;
				}
			}

			// TODO: 确保这里是对的
			for (int k = 0; k < sz; ++k)
				for (int l = 0; l < sz; ++l) {
					double acc = 0;
					for (int m = 0; m < cppSize; ++m)
						for (int n = 0; n < cppSize; ++n) {
							int ix = k - cppSize / 2 + m;
							int iy = l - cppSize / 2 + n;
							WRAP_AROUND(ix, sz);
							WRAP_AROUND(iy, sz);
							acc += error[ix*sz+iy] * cpp_[m][n];						
						}
					p[k*sz+l] = acc;
				}
			// cpe 计算完成
		}
	}
	// 设置修正项
	vector<double> correction(3 * colorCombNbr, 0);

	printf("begin swapping!\n");

	// 开始迭代优化
	for (int i = 0; i < maxIterTimes_; ++i) {
		int benefitPoints = 0;
		// 交换操作
		for (int c = 0; c < screenInfo_.colorNbr_; ++c) {
			int mask = 1 << c;
			for (int x = 0; x < sz; ++x)
				for (int y = 0; y < sz; ++y) {
					int p0 = dotProfile2D[x][y];
					int flat0 = x * sz + y;
					int maxK = -1;
					int maxL = -1;
					double maxD = 0;
					for (int k = 0; k < sz; ++k)
						for (int l = 0; l < sz; ++l) {
							int p1 = dotProfile2D[k][l];
							int flat1 = k * sz + l;
							//  stack constraint
							if ((p0&mask) != (p1&mask) &&
								!(dotProfiles_[pre][flat0] & mask) && !(dotProfiles_[pre][flat1] & mask) &&
								dotProfiles_[next][flat0] & mask && dotProfiles_[next][flat1] & mask)
							{
								int _p0 = p0 ^ mask;
								int _p1 = p1 ^ mask;
								// wrap around
								int n0 = abs(k - x);
								int n1 = abs(l - y);
								//if (n0 > cppSize / 2) {
								//	n0 = sz - n0;
								//}
								//if (n1 > cppSize / 2) {
								//	n1 = sz - n1;
								//}
								double cppn0n1 = 0;
								if (n0 <= cppSize / 2 && n1<= cppSize/2) {
									cppn0n1 = cpp_[n0+cppSize/2][n1+cppSize/2];
								}
								double totalD = 0;
								for (int ccb = 0; ccb < colorCombNbr; ++ccb) {
									if (screenInfo_.colorCombW_[ccb] > 0 && ccb&mask) {
										double colorD = 0;
										int p0m = p0 & ccb;
										int p1m = p1 & ccb;
										int _p0m = _p0 & ccb;
										int _p1m = _p1 & ccb;
										for (int lr = 0; lr < 3; ++lr) {
											double *cpeTemp = cpe + (ccb * 3 + lr)*sc;
											double *lrValue = screenInfo_.labValue_ + lr * colorCombNbr;
											double a0 = lrValue[_p0m] - lrValue[p0m];
											double a1 = lrValue[_p1m] - lrValue[p1m];
											double d = (a0*a0 + a1*a1)*cppCenterValue + 2 * a0*a1*cppn0n1 + 
												           2 * (a0*cpeTemp[x*sz + y] + a1 * cpeTemp[k*sz +l]);
#ifdef USE_CORRECTION
											//计算修正项
											double change = a0 + a1;
											double corr = 2 * change*correction[ccb * 3 + lr] + change * change;
											corr = corr / sc;
											d -= corr;
#endif // USE_CORRECTION
											colorD += d * screenInfo_.labW_[lr];
										}
										totalD += colorD * screenInfo_.colorCombW_[ccb];
									}
								}
								if (totalD < maxD) {
									maxD = totalD;
									maxK = k;
									maxL = l;
								}
							}
						}
					// 存在交换操作可以减少色彩波动，更新cpe
					if (maxD < 0) {
						++benefitPoints;
						int nx = maxK;
						int ny = maxL;
						WRAP_AROUND(nx, sz);
						WRAP_AROUND(ny, sz);
						int p1 = dotProfile2D[nx][ny];

						// 执行交换操作
						dotProfile2D[x][y] ^= mask;
						dotProfile2D[nx][ny] ^= mask;

						int _p0 = p0 ^ mask;
						int _p1 = p1 ^ mask;
						for (int ccb = 0; ccb<colorCombNbr; ++ccb) {
							if (screenInfo_.colorCombW_[ccb] > 0 && ccb&mask) {
								int p0m = p0 & ccb;
								int p1m = p1 & ccb;
								int _p0m = _p0 & ccb;
								int _p1m = _p1 & ccb;
								for (int lr = 0; lr < 3; ++lr) {
									double *cpeTemp = cpe + (ccb * 3 + lr)*	sc;
									double *lrValue = screenInfo_.labValue_+lr*colorCombNbr;
									double a0 = lrValue[_p0m] - lrValue[p0m];
									double a1 = lrValue[_p1m] - lrValue[p1m];
									// 更新修正项
									correction[ccb * 3 + lr] += (a0 + a1);
									for (int ix = 0; ix < cppSize; ++ix) {
										for (int iy = 0; iy < cppSize; ++iy) {
											int p0x = x - cppSize / 2 + ix;
											int p0y = y - cppSize / 2 + iy;
											int p1x = nx - cppSize / 2 + ix;
											int p1y = ny - cppSize / 2 + iy;
											WRAP_AROUND(p0x, sz);
											WRAP_AROUND(p0y, sz);
											WRAP_AROUND(p1x, sz);
											WRAP_AROUND(p1y, sz);
											cpeTemp[p0x*sz + p0y] += a0*cpp_[ix][iy];
											cpeTemp[p1x*sz + p1y] += a1*cpp_[ix][iy];							
										}
									}
								}
							}
						}
					}
				}
		}

		printf("level[%d] %d th iteration finished.\n"
			"Total benefit points: %d\n", currentLevel, i, benefitPoints);
		if (benefitPoints == 0) {
			break;
		}
	}
	// 设置完成标志
	finished_[currentLevel] = true;
	delete[] dotProfile2D;
	delete[] cpe;
}
