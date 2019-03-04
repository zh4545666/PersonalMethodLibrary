#pragma once

#ifdef AFX_CLASS  
#define AFX_EX_CLASS _declspec(dllexport)  
#else  
#define AFX_EX_CLASS _declspec(dllimport)  
#endif 

#include "Matrix.h"

class AFX_EX_CLASS CIntegralTransform
{
public:
	//
	// 构造与析构
	//

	CIntegralTransform();
	~CIntegralTransform();

public:
	//
	// 积分变换方法
	//

	//快速傅里叶变换 r=log2(N) 
	static void FFT(double* TD, CComplex* FD, int r);
	static void FFT(CComplex* TD, CComplex* FD, int r);
	static void IFFT(CComplex* FD, double* TD, int r);
	static void IFFT(CComplex* FD, CComplex* TD, int r);

	//傅里叶变换
	static void DFT(double* TD, CComplex* FD, size_t N);
	static void DFT(CComplex* TD, CComplex* FD, size_t N);
	static void IDFT(CComplex* FD, double* TD, size_t N);
	static void IDFT(CComplex* FD, CComplex* TD, size_t N);

	//希尔伯特
	static void HilBert(double* TD,CComplex* FD,size_t N);
	static void HilBert(CComplex* TD, CComplex* FD, size_t N);

	//S变换 (FD must be preallocated, with N/2+1 rows and N columns)
	static void STransform(double* TD,CComplex* FD,size_t N);
	//S逆变换(与原曲线形态相同，幅值不一致)
	static void ISTransform(CComplex* FD,double* TD,size_t N);

	//离散余弦变换（一维）
	static void DCT(double* TD, double* FD, size_t N);
	static void IDCT(double* FD, double* TD, size_t N);
	//离散余弦变换（二维，矩阵须为方阵）
	static void DCT2(CMatrix* TD, CMatrix* FD);
	static void IDCT2(CMatrix* FD, CMatrix* TD);

private:
	static inline double gauss(int n, int m);
};

