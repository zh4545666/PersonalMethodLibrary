#pragma once

#ifdef AFX_CLASS  
#define AFX_EX_CLASS _declspec(dllexport)  
#else  
#define AFX_EX_CLASS _declspec(dllimport)  
#endif 

#include "Matrix.h"

class AFX_EX_CLASS CLEquation
{
	//
	// 公有接口函数
	//
public:

	//
	// 构造与析构
	//

	CLEquation();				// 默认构造函数
	// 指定系数和常数构造函数
	CLEquation(const CMatrix& mtxCoef, const CMatrix& mtxConst);
	virtual ~CLEquation();		// 析构函数
	// 初始化
	bool Init(const CMatrix& mtxCoef, const CMatrix& mtxConst);

	//
	// 属性
	//

	CMatrix GetCoefMatrix() const;	// 获取系数矩阵
	CMatrix GetConstMatrix() const;	// 获取常数矩阵
	int	GetNumEquations() const;	// 获取方程个数
	int	GetNumUnknowns() const;		// 获取未知数个数

	//
	// 线性方程组求解算法
	//

	// 全选主元高斯消去法
	bool GetRootsetGauss(CMatrix& mtxResult);
	// 全选主元高斯－约当消去法
	bool GetRootsetGaussJordan(CMatrix& mtxResult);
	// 复系数方程组的全选主元高斯消去法
	bool GetRootsetGauss(const CMatrix& mtxCoefImag, const CMatrix& mtxConstImag, CMatrix& mtxResult, CMatrix& mtxResultImag);
	// 复系数方程组的全选主元高斯－约当消去法
	bool GetRootsetGaussJordan(const CMatrix& mtxCoefImag, const CMatrix& mtxConstImag, CMatrix& mtxResult, CMatrix& mtxResultImag);
	// 求解三对角线方程组的追赶法
	bool GetRootsetTriDiagonal(CMatrix& mtxResult);
	// 一般带型方程组的求解
	bool GetRootsetBand(int nBandWidth, CMatrix& mtxResult);
	// 求解对称方程组的分解法
	bool GetRootsetDjn(CMatrix& mtxResult);
	// 求解对称正定方程组的平方根法
	bool GetRootsetCholesky(CMatrix& mtxResult);
	// 求解大型稀疏方程组的全选主元高斯－约去消去法
	bool GetRootsetGgje(CMatrix& mtxResult);
	// 求解托伯利兹方程组的列文逊方法
	bool GetRootsetTlvs(CMatrix& mtxResult);
	// 高斯－赛德尔迭代法
	bool GetRootsetGaussSeidel(CMatrix& mtxResult, double eps = 0.000001);
	// 求解对称正定方程组的共轭梯度法
	void GetRootsetGrad(CMatrix& mtxResult, double eps = 0.000001);
	// 求解线性最小二乘问题的豪斯荷尔德变换法
	bool GetRootsetMqr(CMatrix& mtxResult, CMatrix& mtxQ, CMatrix& mtxR);
	// 求解线性最小二乘问题的广义逆法
	bool GetRootsetGinv(CMatrix& mtxResult, CMatrix& mtxAP, CMatrix& mtxU, CMatrix& mtxV, double eps = 0.000001);
	// 病态方程组的求解
	bool GetRootsetMorbid(CMatrix& mtxResult, int nMaxIt = 60, double eps = 0.000001);

	//
	// 保护性数据成员
	//
protected:
	CMatrix	m_mtxCoef;		// 系数矩阵
	CMatrix m_mtxConst;		// 常数矩阵

};

