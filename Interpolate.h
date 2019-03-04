#pragma once

#ifdef AFX_CLASS  
#define AFX_EX_CLASS _declspec(dllexport)  
#else  
#define AFX_EX_CLASS _declspec(dllimport)  
#endif 

class AFX_EX_CLASS CInterpolate
{
public:
	//
	// 构造与析构
	//

	CInterpolate();
	virtual ~CInterpolate();

	//
	// 将字符串转换成结点值
	//

	static int GetNodesFromString(string s, int n, double dblNodes[], const string& sDelim = " ");

	//
	// 插值算法函数
	//

	// 一元全区间不等距插值
	static double GetValueLagrange(int n, double x[], double y[], double t);
	// 一元全区间等距插值
	static double GetValueLagrange(int n, double x0, double xStep, double y[], double t);
	// 一元三点不等距插值
	static double GetValueLagrange3(int n, double x[], double y[], double t);
	// 一元三点等距插值
	static double GetValueLagrange3(int n, double x0, double xStep, double y[], double t);
	// 连分式不等距插值
	static double GetValuePqs(int n, double x[], double y[], double t);
	// 连分式等距插值
	static double GetValuePqs(int n, double x0, double xStep, double y[], double t);
	// 埃尔米特不等距插值
	static double GetValueHermite(int n, double x[], double y[], double dy[], double t);
	// 埃尔米特等距插值
	static double GetValueHermite(int n, double x0, double xStep, double y[], double dy[], double t);
	// 埃特金不等距逐步插值
	static double GetValueAitken(int n, double x[], double y[], double t, double eps = 0.000001);
	// 埃特金等距逐步插值
	static double GetValueAitken(int n, double x0, double xStep, double y[], double t, double eps = 0.000001);
	// 光滑不等距插值
	static double GetValueAkima(int n, double x[], double y[], double t, double s[], int k = -1);
	// 光滑等距插值
	static double GetValueAkima(int n, double x0, double xStep, double y[], double t, double s[], int k = -1);
	// 第一种边界条件的三次样条函数插值、微商与积分
	static double GetValueSpline1(int n, double x[], double y[], double dy[], double ddy[],
		int m, double t[], double z[], double dz[], double ddz[]);
	// 第二种边界条件的三次样条函数插值、微商与积分
	static double GetValueSpline2(int n, double x[], double y[], double dy[], double ddy[],
		int m, double t[], double z[], double dz[], double ddz[]);
	// 第三种边界条件的三次样条函数插值、微商与积分
	static double GetValueSpline3(int n, double x[], double y[], double dy[], double ddy[],
		int m, double t[], double z[], double dz[], double ddz[]);
	// 二元三点插值
	static double GetValueTqip(int n, double x[], int m, double y[], double z[], double u, double v);
	// 二元全区间插值
	static double GetValueLagrange2(int n, double x[], int m, double y[], double z[], double u, double v);

public:	
	//三次样条插值
	static double Spline(double x[], double y[], int n, double ddy1, double ddyn, double t[], int m, double z[]);
	//拉格朗日插值
	static double Lagrange(double *x, double *y, double xx, int n);

};

