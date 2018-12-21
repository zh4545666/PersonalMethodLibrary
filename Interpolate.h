#pragma once

#ifdef AFX_CLASS  
#define AFX_EX_CLASS _declspec(dllexport)  
#else  
#define AFX_EX_CLASS _declspec(dllimport)  
#endif 

class AFX_EX_CLASS CInterpolate
{
public:
	//
	// ����������
	//

	CInterpolate();
	virtual ~CInterpolate();

	//
	// ���ַ���ת���ɽ��ֵ
	//

	static int GetNodesFromString(string s, int n, double dblNodes[], const string& sDelim = " ");

	//
	// ��ֵ�㷨����
	//

	// һԪȫ���䲻�Ⱦ��ֵ
	static double GetValueLagrange(int n, double x[], double y[], double t);
	// һԪȫ����Ⱦ��ֵ
	static double GetValueLagrange(int n, double x0, double xStep, double y[], double t);
	// һԪ���㲻�Ⱦ��ֵ
	static double GetValueLagrange3(int n, double x[], double y[], double t);
	// һԪ����Ⱦ��ֵ
	static double GetValueLagrange3(int n, double x0, double xStep, double y[], double t);
	// ����ʽ���Ⱦ��ֵ
	static double GetValuePqs(int n, double x[], double y[], double t);
	// ����ʽ�Ⱦ��ֵ
	static double GetValuePqs(int n, double x0, double xStep, double y[], double t);
	// �������ز��Ⱦ��ֵ
	static double GetValueHermite(int n, double x[], double y[], double dy[], double t);
	// �������صȾ��ֵ
	static double GetValueHermite(int n, double x0, double xStep, double y[], double dy[], double t);
	// ���ؽ𲻵Ⱦ��𲽲�ֵ
	static double GetValueAitken(int n, double x[], double y[], double t, double eps = 0.000001);
	// ���ؽ�Ⱦ��𲽲�ֵ
	static double GetValueAitken(int n, double x0, double xStep, double y[], double t, double eps = 0.000001);
	// �⻬���Ⱦ��ֵ
	static double GetValueAkima(int n, double x[], double y[], double t, double s[], int k = -1);
	// �⻬�Ⱦ��ֵ
	static double GetValueAkima(int n, double x0, double xStep, double y[], double t, double s[], int k = -1);
	// ��һ�ֱ߽���������������������ֵ��΢�������
	static double GetValueSpline1(int n, double x[], double y[], double dy[], double ddy[],
		int m, double t[], double z[], double dz[], double ddz[]);
	// �ڶ��ֱ߽���������������������ֵ��΢�������
	static double GetValueSpline2(int n, double x[], double y[], double dy[], double ddy[],
		int m, double t[], double z[], double dz[], double ddz[]);
	// �����ֱ߽���������������������ֵ��΢�������
	static double GetValueSpline3(int n, double x[], double y[], double dy[], double ddy[],
		int m, double t[], double z[], double dz[], double ddz[]);
	// ��Ԫ�����ֵ
	static double GetValueTqip(int n, double x[], int m, double y[], double z[], double u, double v);
	// ��Ԫȫ�����ֵ
	static double GetValueLagrange2(int n, double x[], int m, double y[], double z[], double u, double v);

public:	
	//����������ֵ
	static double Spline(double x[], double y[], int n, double ddy1, double ddyn, double t[], int m, double z[]);
	//�������ղ�ֵ
	static double Lagrange(double *x, double *y, double xx, int n);

};

