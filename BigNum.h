#pragma once

#include "stdafx.h"

#ifdef AFX_CLASS  
#define AFX_EX_CLASS _declspec(dllexport)  
#else  
#define AFX_EX_CLASS _declspec(dllimport)  
#endif 

class AFX_EX_CLASS CBigNum
{
private:
	int a[500];    //可以控制大数的位数 
	int len;       //大数长度
public:
	CBigNum(){ len = 1; memset(a, 0, sizeof(a)); }   //构造函数
	CBigNum(const int);       //将一个int类型的变量转化为大数
	CBigNum(const char*, bool bIgnoreNonnumeric = true);     //将一个字符串类型的变量转化为大数
	CBigNum(const CBigNum &);  //拷贝构造函数
	CBigNum &operator=(const CBigNum &);   //重载赋值运算符，大数之间进行赋值运算
	CBigNum operator+(const CBigNum &) const;   //重载加法运算符，两个大数之间的相加运算 
	CBigNum operator-(const CBigNum &) const;   //重载减法运算符，两个大数之间的相减运算 
	CBigNum operator*(const CBigNum &) const;   //重载乘法运算符，两个大数之间的相乘运算 
	CBigNum operator/(const int   &) const;    //重载除法运算符，大数对一个整数进行相除运算
	CBigNum operator^(const int  &) const;    //大数的n次方运算
	int    operator%(const int  &) const;    //大数对一个int类型的变量进行取模运算    
	bool   operator>(const CBigNum & T)const;   //大数和另一个大数的大小比较
	bool   operator>(const int & t)const;      //大数和一个int类型的变量的大小比较
	string ToString();

	friend istream& operator>>(istream&, CBigNum&);   //重载输入运算符
	friend ostream& operator<<(ostream&, CBigNum&);   //重载输出运算符
};

