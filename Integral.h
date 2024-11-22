#pragma once

#ifdef AFX_CLASS  
#define AFX_EX_CLASS _declspec(dllexport)  
#else  
#define AFX_EX_CLASS _declspec(dllimport)  
#endif 

namespace PersonalMethod {

	class AFX_EX_CLASS CIntegral
	{
	public:
		//
		// 构造与析构
		//

		CIntegral();
		virtual ~CIntegral();

		//
		// 虚函数：计算积分函数值，必须在派生类中覆盖该函数
		//

		virtual double Func(double x)
		{
			return 0.0;
		}

		//
		// 数值积分求解算法
		//

		// 变步长梯形求积法
		double GetValueTrapezia(double a, double b, double eps = 0.000001);
		// 变步长辛卜生求积法
		double GetValueSimpson(double a, double b, double eps = 0.000001);
		// 自适应梯形求积法
		double GetValueATrapezia(double a, double b, double d, double eps = 0.000001);
		// 龙贝格求积法
		double GetValueRomberg(double a, double b, double eps = 0.000001);
		// 计算一维积分的连分式法
		double GetValuePq(double a, double b, double eps = 0.000001);
		// 高振荡函数求积法
		double GetValuePart(double a, double b, int m, int n, double fa[], double fb[], double s[]);
		// 勒让德－高斯求积法
		double GetValueLegdGauss(double a, double b, double eps = 0.000001);
		// 拉盖尔－高斯求积法
		double GetValueLgreGauss();
		// 埃尔米特－高斯求积法
		double GetValueHermiteGauss();

		//
		// 内部函数
		//
	private:
		void ppp(double x0, double x1, double h, double f0, double f1, double t0, double eps, double d, double t[]);

	};

}