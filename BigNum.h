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
	int a[500];    //���Կ��ƴ�����λ�� 
	int len;       //��������
public:
	CBigNum(){ len = 1; memset(a, 0, sizeof(a)); }   //���캯��
	CBigNum(const int);       //��һ��int���͵ı���ת��Ϊ����
	CBigNum(const char*, bool bIgnoreNonnumeric = true);     //��һ���ַ������͵ı���ת��Ϊ����
	CBigNum(const CBigNum &);  //�������캯��
	CBigNum &operator=(const CBigNum &);   //���ظ�ֵ�����������֮����и�ֵ����
	CBigNum operator+(const CBigNum &) const;   //���ؼӷ����������������֮���������� 
	CBigNum operator-(const CBigNum &) const;   //���ؼ������������������֮���������� 
	CBigNum operator*(const CBigNum &) const;   //���س˷����������������֮���������� 
	CBigNum operator/(const int   &) const;    //���س����������������һ�����������������
	CBigNum operator^(const int  &) const;    //������n�η�����
	int    operator%(const int  &) const;    //������һ��int���͵ı�������ȡģ����    
	bool   operator>(const CBigNum & T)const;   //��������һ�������Ĵ�С�Ƚ�
	bool   operator>(const int & t)const;      //������һ��int���͵ı����Ĵ�С�Ƚ�
	string ToString();

	friend istream& operator>>(istream&, CBigNum&);   //�������������
	friend ostream& operator<<(ostream&, CBigNum&);   //������������
};

