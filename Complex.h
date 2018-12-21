#pragma once

#ifdef AFX_CLASS  
#define AFX_EX_CLASS _declspec(dllexport)  
#else  
#define AFX_EX_CLASS _declspec(dllimport)  
#endif 

#include "StringBuffer.h"

class AFX_EX_CLASS CComplex
{
	//
	// ���нӿں���
	//
public:

	//
	// ����������
	//

	CComplex();							// �������캯��
	CComplex(double dblX, double dblY);	// ָ��ֵ���캯��
	CComplex(const CComplex& other);	// �������캯��
	virtual ~CComplex() {};				// ��������

	//
	// ��������ʾ
	//

	void SetData(double dblX, double dblY);
	void SetReal(double dblX);	// ָ��������ʵ��
	void SetImag(double dblY);	// ָ���������鲿
	double GetReal();			// ȡ������ʵ��
	double GetImag();			// ȡ�������鲿
	string ToString(bool bflag = true) const;	// ������ת��Ϊ"a+bj"��ʽ���ַ���
	// ��"a,b"��ʽ���ַ�����ת��Ϊ��������aΪ������ʵ����bΪ�������鲿(sDelim����Ϊ'i'��'j')
	void FromString(string s, const string& sDelim = "");

	//
	// ��ѧ����
	//

	bool operator==(const CComplex& cpxX) const;
	bool operator!=(const CComplex& cpxX) const;
	CComplex& operator=(const CComplex& cpxX);
	CComplex operator+(const CComplex& cpxX) const;
	CComplex operator-(const CComplex& cpxX) const;
	CComplex operator*(const CComplex& cpxX) const;
	CComplex operator*(const double& cpxX) const;
	CComplex operator/(const CComplex& cpxX) const;
	double Abs() const;	// ������ģ

	//
	// ��������
	//

	void Root(int n, CComplex cpxR[]) const;		// �����ĸ�
	CComplex Pow(double dblW) const;				// ������ʵ��ָ��
	CComplex Pow(CComplex cpxW, int n = 0) const;	// �����ĸ���ָ��
	CComplex Log() const;							// �����Ķ���
	CComplex Sin() const;							// ����������
	CComplex Cos() const;							// ����������
	CComplex Tan() const;							// ����������

	//
	// ����������
	//
protected:
	double	m_dblX;		// ������ʵ��
	double	m_dblY;		// �������鲿
};

