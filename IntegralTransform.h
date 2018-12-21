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
	// ����������
	//

	CIntegralTransform();
	~CIntegralTransform();

public:
	//
	// ���ֱ任����
	//

	//���ٸ���Ҷ�任 r=log2(N) 
	static void FFT(double* TD, CComplex* FD, int r);
	static void FFT(CComplex* TD, CComplex* FD, int r);
	static void IFFT(CComplex* FD, double* TD, int r);
	static void IFFT(CComplex* FD, CComplex* TD, int r);

	//����Ҷ�任
	static void DFT(double* TD, CComplex* FD, size_t N);
	static void DFT(CComplex* TD, CComplex* FD, size_t N);
	static void IDFT(CComplex* FD, double* TD, size_t N);
	static void IDFT(CComplex* FD, CComplex* TD, size_t N);

	//ϣ������
	static void HilBert(double* TD,CComplex* FD,size_t N);
	static void HilBert(CComplex* TD, CComplex* FD, size_t N);

	//S�任 (FD must be preallocated, with N/2+1 rows and N columns)
	static void STransform(double* TD,CComplex* FD,size_t N);
	//S��任(��ԭ������̬��ͬ����ֵ��һ��)
	static void ISTransform(CComplex* FD,double* TD,size_t N);

	//��ɢ���ұ任��һά��
	static void DCT(double* TD, double* FD, size_t N);
	static void IDCT(double* FD, double* TD, size_t N);
	//��ɢ���ұ任����ά��������Ϊ����
	static void DCT2(CMatrix* TD, CMatrix* FD);
	static void IDCT2(CMatrix* FD, CMatrix* TD);

private:
	static inline double gauss(int n, int m);
};

