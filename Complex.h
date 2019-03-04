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
	// 公有接口函数
	//
public:

	//
	// 构造与析构
	//

	CComplex();							// 基本构造函数
	CComplex(double dblX, double dblY);	// 指定值构造函数
	CComplex(const CComplex& other);	// 拷贝构造函数
	virtual ~CComplex() {};				// 析构函数

	//
	// 输入与显示
	//

	void SetData(double dblX, double dblY);
	void SetReal(double dblX);	// 指定复数的实部
	void SetImag(double dblY);	// 指定复数的虚部
	double GetReal();			// 取复数的实部
	double GetImag();			// 取复数的虚部
	string ToString(bool bflag = true) const;	// 将复数转化为"a+bj"形式的字符串
	// 将"a,b"形式的字符串对转化为复数，以a为复数的实部，b为复数的虚部(sDelim不能为'i'或'j')
	void FromString(string s, const string& sDelim = "");

	//
	// 数学运算
	//

	bool operator==(const CComplex& cpxX) const;
	bool operator!=(const CComplex& cpxX) const;
	CComplex& operator=(const CComplex& cpxX);
	CComplex operator+(const CComplex& cpxX) const;
	CComplex operator-(const CComplex& cpxX) const;
	CComplex operator*(const CComplex& cpxX) const;
	CComplex operator*(const double& cpxX) const;
	CComplex operator/(const CComplex& cpxX) const;
	double Abs() const;	// 复数的模

	//
	// 函数运算
	//

	void Root(int n, CComplex cpxR[]) const;		// 复数的根
	CComplex Pow(double dblW) const;				// 复数的实幂指数
	CComplex Pow(CComplex cpxW, int n = 0) const;	// 复数的复幂指数
	CComplex Log() const;							// 复数的对数
	CComplex Sin() const;							// 复数的正弦
	CComplex Cos() const;							// 复数的余弦
	CComplex Tan() const;							// 复数的正切

	//
	// 保护性数据
	//
protected:
	double	m_dblX;		// 复数的实部
	double	m_dblY;		// 复数的虚部
};

