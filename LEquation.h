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
	// ���нӿں���
	//
public:

	//
	// ����������
	//

	CLEquation();				// Ĭ�Ϲ��캯��
	// ָ��ϵ���ͳ������캯��
	CLEquation(const CMatrix& mtxCoef, const CMatrix& mtxConst);
	virtual ~CLEquation();		// ��������
	// ��ʼ��
	bool Init(const CMatrix& mtxCoef, const CMatrix& mtxConst);

	//
	// ����
	//

	CMatrix GetCoefMatrix() const;	// ��ȡϵ������
	CMatrix GetConstMatrix() const;	// ��ȡ��������
	int	GetNumEquations() const;	// ��ȡ���̸���
	int	GetNumUnknowns() const;		// ��ȡδ֪������

	//
	// ���Է���������㷨
	//

	// ȫѡ��Ԫ��˹��ȥ��
	bool GetRootsetGauss(CMatrix& mtxResult);
	// ȫѡ��Ԫ��˹��Լ����ȥ��
	bool GetRootsetGaussJordan(CMatrix& mtxResult);
	// ��ϵ���������ȫѡ��Ԫ��˹��ȥ��
	bool GetRootsetGauss(const CMatrix& mtxCoefImag, const CMatrix& mtxConstImag, CMatrix& mtxResult, CMatrix& mtxResultImag);
	// ��ϵ���������ȫѡ��Ԫ��˹��Լ����ȥ��
	bool GetRootsetGaussJordan(const CMatrix& mtxCoefImag, const CMatrix& mtxConstImag, CMatrix& mtxResult, CMatrix& mtxResultImag);
	// ������Խ��߷������׷�Ϸ�
	bool GetRootsetTriDiagonal(CMatrix& mtxResult);
	// һ����ͷ���������
	bool GetRootsetBand(int nBandWidth, CMatrix& mtxResult);
	// ���ԳƷ�����ķֽⷨ
	bool GetRootsetDjn(CMatrix& mtxResult);
	// ���Գ������������ƽ������
	bool GetRootsetCholesky(CMatrix& mtxResult);
	// ������ϡ�跽�����ȫѡ��Ԫ��˹��Լȥ��ȥ��
	bool GetRootsetGgje(CMatrix& mtxResult);
	// ����в����ȷ����������ѷ����
	bool GetRootsetTlvs(CMatrix& mtxResult);
	// ��˹�����¶�������
	bool GetRootsetGaussSeidel(CMatrix& mtxResult, double eps = 0.000001);
	// ���Գ�����������Ĺ����ݶȷ�
	void GetRootsetGrad(CMatrix& mtxResult, double eps = 0.000001);
	// ���������С��������ĺ�˹�ɶ��±任��
	bool GetRootsetMqr(CMatrix& mtxResult, CMatrix& mtxQ, CMatrix& mtxR);
	// ���������С��������Ĺ����淨
	bool GetRootsetGinv(CMatrix& mtxResult, CMatrix& mtxAP, CMatrix& mtxU, CMatrix& mtxV, double eps = 0.000001);
	// ��̬����������
	bool GetRootsetMorbid(CMatrix& mtxResult, int nMaxIt = 60, double eps = 0.000001);

	//
	// ���������ݳ�Ա
	//
protected:
	CMatrix	m_mtxCoef;		// ϵ������
	CMatrix m_mtxConst;		// ��������

};

