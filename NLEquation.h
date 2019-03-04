#pragma once

#ifdef AFX_CLASS  
#define AFX_EX_CLASS _declspec(dllexport)  
#else  
#define AFX_EX_CLASS _declspec(dllimport)  
#endif 

#include "LEquation.h"
#include "Matrix.h"
#include "Complex.h"

class AFX_EX_CLASS CNLEquation
{
	//
	// 公有接口函数
	//
public:

	//
	// 构造与析构
	//

	CNLEquation();
	virtual ~CNLEquation();

	//
	// 虚函数：计算方程左端函数值，必须在引申类中覆盖该类函数
	//

	virtual double Func(double x)
	{
		return 0.0;
	}
	virtual double Func(int n, double x[])
	{
		return 0.0;
	}
	virtual void Func(double x, double y[])
	{
	}
	virtual double Func(double x, double y)
	{
		return 0.0;
	}
	virtual double Func(double x[], double y[])
	{
		return 0.0;
	}
	virtual void FuncMJ(int n, double x[], double p[])
	{
	}

	//
	// 非线性方程求解算法
	//

	// 求非线性方程实根的对分法
	int GetRootBisect(int nNumRoots, double x[], double xStart, double xEnd, double dblStep, double eps = 0.000001);
	// 求非线性方程一个实根的牛顿法
	bool GetRootNewton(double* x, int nMaxIt = 60, double eps = 0.000001);
	// 求非线性方程一个实根的埃特金迭代法
	bool GetRootAitken(double* x, int nMaxIt = 60, double eps = 0.000001);
	// 求非线性方程一个实根的连分式解法
	bool GetRootPq(double* x, double eps = 0.000001);
	// 求实系数代数方程全部根的QR方法
	bool GetRootQr(int n, double dblCoef[], double xr[], double xi[], int nMaxIt = 60, double eps = 0.000001);
	// 求实系数代数方程全部根的牛顿下山法
	bool GetRootNewtonDownHill(int n, double dblCoef[], double xr[], double xi[]);
	// 求复系数代数方程全部根的牛顿下山法
	bool GetRootNewtonDownHill(int n, double ar[], double ai[], double xr[], double xi[]);
	// 求非线性方程一个实根的蒙特卡洛法
	void GetRootMonteCarlo(double* x, double xStart, int nControlB, double eps = 0.000001);
	// 求实函数或复函数方程一个复根的蒙特卡洛法
	void GetRootMonteCarlo(double* x, double* y, double xStart, int nControlB, double eps = 0.000001);

	//
	// 非线性方程组求解算法
	//

	// 求非线性方程组一组实根的梯度法
	bool GetRootsetGrad(int n, double x[], int nMaxIt = 500, double eps = 0.000001);
	// 求非线性方程组一组实根的拟牛顿法
	bool GetRootsetNewton(int n, double x[], double t, double h, int nMaxIt = 500, double eps = 0.000001);
	// 求非线性方程组最小二乘解的广义逆法
	bool GetRootsetGinv(int m, int n, double x[], double eps1 = 0.000001, double eps2 = 0.000001);
	// 求非线性方程组一组实根的蒙特卡洛法
	void GetRootsetMonteCarlo(int n, double x[], double xStart, int nControlB, double eps = 0.000001);

	//
	// 内部函数
	//
private:
	void g60(double* t, double* x, double* y, double* x1, double* y1, double* dx, double* dy, double* p, double* q, int* k, int* it);
	void g90(double xr[], double xi[], double dblCoef[], double* x, double* y, double* p, double* q, double* w, int* k);
	void g65(double* x, double* y, double* x1, double* y1, double* dx, double* dy, double* dd, double* dc, double* c, int* k, int* is, int* it);
	void g60c(double* t, double* x, double* y, double* x1, double* y1, double* dx, double* dy, double* p, double* q, int* k, int* it);
	void g90c(double xr[], double xi[], double ar[], double ai[], double* x, double* y, double* p, double* w, int* k);
	void g65c(double* x, double* y, double* x1, double* y1, double* dx, double* dy, double* dd, double* dc, double* c, int* k, int* is, int* it);
	double rnd(double* r);
};

