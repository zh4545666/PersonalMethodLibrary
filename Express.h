#pragma once

#include "stdafx.h"

#define EXP_PTR_NUM 11  //���������
#define EXP_FUN_NUM 12  //��������

#ifdef AFX_CLASS  
#define AFX_EX_CLASS _declspec(dllexport)  
#else  
#define AFX_EX_CLASS _declspec(dllimport)  
#endif 

/* ���ʷ������� */
typedef enum{
	TKT_NUMBER,    //���� 
	TKT_OPERATOR,  //����� 
	TKT_FUNCTION,  //���� 

	TKT_ENDSIGN,   //������ 
	TKT_UNKNOW     //δ֪���� 
}TokenType_Express;


class AFX_EX_CLASS CExpress
{
public:
	CExpress(const string &expression);
	CExpress();
	~CExpress();

public:
	int getVal(double &res, const string &expression);
	int getVal(double &res);

private:
	string expression;        //���ʽ 
	string token;             //ÿ�ζ�ȡ�ĵ��� 
	TokenType_Express tkType; //��ȡ�ĵ������� 
	int pos, length;          //��ȡ��λ�úͳ��� 

	/* ��̬���� */
	static string ptrList[];   //֧�ֵ������
	static int ptrArgCnt[];    //�������������

	static string funList[];   //֧�ֵĺ����б� 
	static int funArgCnt[];    //������������� 
	static int preceMap[][EXP_PTR_NUM];      //�������ȼ��� 

	/* ������� */
	void debug();

	/* ��ȡ��һ������ */
	void readToken();

	/* �Ƚ���������������� */
	int comparePrece(const string &ptr1, const std::string &ptr2);

	/* ������������� */
	double calculate(const string &ptr, double arg[]);

	/* ������������ */
	double callFun(const string &fun, double arg[]);

	/* ��ȡ�������� */
	int getPtrIndex(const string &ptr);

	/* ��ȡ������� */
	int getFunIndex(const string &fun);

	/* �Ӳ�����ջ��ȡn������ */
	bool getArg(stack<double> &opnd, double arg[], int n);
};

