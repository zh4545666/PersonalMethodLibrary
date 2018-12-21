#pragma once

#ifdef AFX_CLASS  
#define AFX_EX_CLASS _declspec(dllexport)  
#else  
#define AFX_EX_CLASS _declspec(dllimport)  
#endif 

class AFX_EX_CLASS CIntegral
{
public:
	//
	// ����������
	//

	CIntegral();
	virtual ~CIntegral();

	//
	// �麯����������ֺ���ֵ���������������и��Ǹú���
	//

	virtual double Func(double x)
	{
		return 0.0;
	}

	//
	// ��ֵ��������㷨
	//

	// �䲽�����������
	double GetValueTrapezia(double a, double b, double eps = 0.000001);
	// �䲽�������������
	double GetValueSimpson(double a, double b, double eps = 0.000001);
	// ����Ӧ���������
	double GetValueATrapezia(double a, double b, double d, double eps = 0.000001);
	// �����������
	double GetValueRomberg(double a, double b, double eps = 0.000001);
	// ����һά���ֵ�����ʽ��
	double GetValuePq(double a, double b, double eps = 0.000001);
	// ���񵴺��������
	double GetValuePart(double a, double b, int m, int n, double fa[], double fb[], double s[]);
	// ���õ£���˹�����
	double GetValueLegdGauss(double a, double b, double eps = 0.000001);
	// ���Ƕ�����˹�����
	double GetValueLgreGauss();
	// �������أ���˹�����
	double GetValueHermiteGauss();

	//
	// �ڲ�����
	//
private:
	void ppp(double x0, double x1, double h, double f0, double f1, double t0, double eps, double d, double t[]);

};

