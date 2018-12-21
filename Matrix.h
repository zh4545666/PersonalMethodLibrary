#pragma once

#ifdef AFX_CLASS  
#define AFX_EX_CLASS _declspec(dllexport)  
#else  
#define AFX_EX_CLASS _declspec(dllimport)  
#endif 

#include "StringBuffer.h"

class AFX_EX_CLASS CMatrix
{
	//
	// ���нӿں���
	//
public:

	//
	// ����������
	//

	CMatrix();										// �������캯��
	CMatrix(int nRows, int nCols);					// ָ�����й��캯��
	CMatrix(int nRows, int nCols, double value[]);	// ָ�����ݹ��캯��
	CMatrix(int nSize);								// �����캯��
	CMatrix(int nSize, double value[]);				// ָ�����ݷ����캯��
	CMatrix(const CMatrix& other);					// �������캯��
	bool	Init(int nRows, int nCols);				// ��ʼ������	
	bool	MakeUnitMatrix(int nSize);				// �������ʼ��Ϊ��λ����
	virtual ~CMatrix();								// ��������

	//
	// ��������ʾ
	//

	// ���ַ���ת��Ϊ��������
	bool FromString(string s, const string& sDelim = " ", bool bLineBreak = true);
	// ������ת��Ϊ�ַ���
	string ToString(const string& sDelim = " ", bool bLineBreak = true) const;
	// �������ָ����ת��Ϊ�ַ���
	string RowToString(int nRow, const string& sDelim = " ") const;
	// �������ָ����ת��Ϊ�ַ���
	string ColToString(int nCol, const string& sDelim = " ") const;

	//
	// Ԫ����ֵ����
	//

	bool	SetElement(int nRow, int nCol, double value);	// ����ָ��Ԫ�ص�ֵ
	double	GetElement(int nRow, int nCol) const;			// ��ȡָ��Ԫ�ص�ֵ
	void    SetData(double value[]);						// ���þ����ֵ(value���ȱ������m_nNumColumns*m_nNumRows)
	int		GetNumColumns() const;							// ��ȡ���������
	int		GetNumRows() const;								// ��ȡ���������
	int     GetRowVector(int nRow, double* pVector) const;	// ��ȡ�����ָ���о���
	int     GetColVector(int nCol, double* pVector) const;	// ��ȡ�����ָ���о���
	double* GetData() const;								// ��ȡ�����ֵ

	//
	// ��ѧ����
	//

	CMatrix& operator=(const CMatrix& other);
	bool operator==(const CMatrix& other) const;
	bool operator!=(const CMatrix& other) const;
	CMatrix	operator+(const CMatrix& other) const;
	CMatrix	operator-(const CMatrix& other) const;
	CMatrix	operator*(double value) const;
	CMatrix	operator*(const CMatrix& other) const;
	// ������˷�
	bool CMul(const CMatrix& AR, const CMatrix& AI, const CMatrix& BR, const CMatrix& BI, CMatrix& CR, CMatrix& CI) const;
	// �����ת��
	CMatrix Transpose() const;

	//
	// �㷨
	//

	// ʵ���������ȫѡ��Ԫ��˹��Լ����
	bool InvertGaussJordan();
	// �����������ȫѡ��Ԫ��˹��Լ����
	bool InvertGaussJordan(CMatrix& mtxImag);
	// �Գ��������������
	bool InvertSsgj();
	// �в����Ⱦ�������İ����ط���
	bool InvertTrench();
	// ������ʽֵ��ȫѡ��Ԫ��˹��ȥ��
	double DetGauss();
	// ������ȵ�ȫѡ��Ԫ��˹��ȥ��
	int RankGauss();
	// �Գ��������������˹���ֽ�������ʽ����ֵ
	bool DetCholesky(double* dblDet);
	// ��������Ƿֽ�
	bool SplitLU(CMatrix& mtxL, CMatrix& mtxU);
	// һ��ʵ�����QR�ֽ�
	bool SplitQR(CMatrix& mtxQ);
	// һ��ʵ���������ֵ�ֽ�
	bool SplitUV(CMatrix& mtxU, CMatrix& mtxV, double eps = 0.000001);
	// ������������ֵ�ֽⷨ
	bool GInvertUV(CMatrix& mtxAP, CMatrix& mtxU, CMatrix& mtxV, double eps = 0.000001);
	// Լ���Գƾ���Ϊ�Գ����Խ���ĺ�˹�ɶ��±任��
	bool MakeSymTri(CMatrix& mtxQ, CMatrix& mtxT, double dblB[], double dblC[]);
	// ʵ�Գ����Խ����ȫ������ֵ�����������ļ���
	bool SymTriEigenv(double dblB[], double dblC[], CMatrix& mtxQ, int nMaxIt = 60, double eps = 0.000001);
	// Լ��һ��ʵ����Ϊ���겮�����ĳ������Ʊ任��
	void MakeHberg();
	// ����겮�����ȫ������ֵ��QR����
	bool HBergEigenv(double dblU[], double dblV[], int nMaxIt = 60, double eps = 0.000001);
	// ��ʵ�Գƾ�������ֵ�������������ſɱȷ�
	bool JacobiEigenv(double dblEigenValue[], CMatrix& mtxEigenVector, int nMaxIt = 60, double eps = 0.000001);
	// ��ʵ�Գƾ�������ֵ�������������ſɱȹ��ط�
	bool JacobiEigenv2(double dblEigenValue[], CMatrix& mtxEigenVector, double eps = 0.000001);

	//
	// ���������ݳ�Ա
	//
protected:
	int	m_nNumColumns;			// ��������
	int	m_nNumRows;				// ��������
	double*	m_pData;			// �������ݻ�����

	//
	// �ڲ�����
	//
private:
	void ppp(double a[], double e[], double s[], double v[], int m, int n);
	void sss(double fg[2], double cs[2]);

};


#include "Complex.h"
class AFX_EX_CLASS CComplexMatrix
{
public:
	//
	// ����������
	//

	CComplexMatrix();										// �������캯��
	CComplexMatrix(int nRows, int nCols);					// ָ�����й��캯��
	CComplexMatrix(int nRows, int nCols, CComplex value[]);	// ָ�����ݹ��캯��
	CComplexMatrix(int nSize);								// �����캯��
	CComplexMatrix(int nSize, CComplex value[]);			// ָ�����ݷ����캯��
	CComplexMatrix(const CComplexMatrix& other);			// �������캯��
	bool	Init(int nRows, int nCols);				// ��ʼ������	
	bool	MakeUnitMatrix(int nSize);				// �������ʼ��Ϊ��λ����
	virtual ~CComplexMatrix();

	//
	// ��������ʾ
	//

	// ���ַ���ת��Ϊ��������
	bool FromString(string s, const string& sDelim = " ", bool bLineBreak = true);
	// ������ת��Ϊ�ַ���
	string ToString(const string& sDelim = " ", bool bLineBreak = true) const;
	// �������ָ����ת��Ϊ�ַ���
	string RowToString(int nRow, const string& sDelim = " ") const;
	// �������ָ����ת��Ϊ�ַ���
	string ColToString(int nCol, const string& sDelim = " ") const;

	//
	// Ԫ����ֵ����
	//

	bool		SetElement(int nRow, int nCol, CComplex& pValue);	// ����ָ��Ԫ�ص�ֵ
	CComplex	GetElement(int nRow, int nCol) const;			// ��ȡָ��Ԫ�ص�ֵ
	void		SetData(CComplex value[]);						// ���þ����ֵ(value���ȱ������m_nNumColumns*m_nNumRows)
	int			GetNumColumns() const;							// ��ȡ���������
	int			GetNumRows() const;								// ��ȡ���������
	int			GetRowVector(int nRow, CComplex* pVector) const;	// ��ȡ�����ָ���о���
	int			GetColVector(int nCol, CComplex* pVector) const;	// ��ȡ�����ָ���о���
	CComplex*	GetData() const;								// ��ȡ�����ֵ

	//
	// ��ѧ����
	//

	CComplexMatrix& operator=(const CComplexMatrix& other);
	bool operator==(const CComplexMatrix& other) const;
	bool operator!=(const CComplexMatrix& other) const;
	CComplexMatrix	operator+(const CComplexMatrix& other) const;
	CComplexMatrix	operator-(const CComplexMatrix& other) const;
	CComplexMatrix	operator*(double value) const;
	CComplexMatrix	operator*(const CComplexMatrix& other) const;
	CComplexMatrix Transpose() const;

private:


	//
	// ���������ݳ�Ա
	//
protected:
	int	m_nNumColumns;			// ��������
	int	m_nNumRows;				// ��������
	CComplex*	m_pData;			// �������ݻ�����

};


