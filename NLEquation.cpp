#include "stdafx.h"
#include "NLEquation.h"


//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////
// �������캯��
//////////////////////////////////////////////////////////////////////
CNLEquation::CNLEquation()
{
}

//////////////////////////////////////////////////////////////////////
// ��������
//////////////////////////////////////////////////////////////////////
CNLEquation::~CNLEquation()
{
}

//////////////////////////////////////////////////////////////////////
// ������Է���ʵ���ĶԷַ�
//
// ����ʱ���븲�Ǽ��㷽����˺���f(x)ֵ���麯��double Func(double x)
//
// ������
// 1. int nNumRoots - ��[xStart, xEnd]��ʵ��������Ԥ��ֵ
// 2. double x[] - һά���飬����Ϊm������������[xStart, xEnd]��������
//                 ��ʵ����ʵ�������ɺ���ֵ����
// 3. double xStart - ����������˵�
// 4. double xEnd - ���������Ҷ˵�
// 5. double dblStep - �������ʱ���õĲ���
// 6. double eps - ���ȿ��Ʋ�����Ĭ��ֵΪ0.000001
//
// ����ֵ��int �ͣ���õ�ʵ������Ŀ
//////////////////////////////////////////////////////////////////////
int CNLEquation::GetRootBisect(int nNumRoots, double x[], double xStart, double xEnd, double dblStep, double eps /*= 0.000001*/)
{
	int n, js;
	double z, y, z1, y1, z0, y0;

	// ���ĸ�����0
	n = 0;

	// ����˵㿪ʼ����
	z = xStart;
	y = Func(z);

	// ѭ�����
	while ((z <= xEnd + dblStep / 2.0) && (n != nNumRoots))
	{
		if (fabs(y)<eps)
		{
			n = n + 1;
			x[n - 1] = z;
			z = z + dblStep / 2.0;
			y = Func(z);
		}
		else
		{
			z1 = z + dblStep;
			y1 = Func(z1);

			if (fabs(y1)<eps)
			{
				n = n + 1;
				x[n - 1] = z1;
				z = z1 + dblStep / 2.0;
				y = Func(z);
			}
			else if (y*y1>0.0)
			{
				y = y1;
				z = z1;
			}
			else
			{
				js = 0;
				while (js == 0)
				{
					if (fabs(z1 - z)<eps)
					{
						n = n + 1;
						x[n - 1] = (z1 + z) / 2.0;
						z = z1 + dblStep / 2.0; y = Func(z);
						js = 1;
					}
					else
					{
						z0 = (z1 + z) / 2.0;
						y0 = Func(z0);
						if (fabs(y0)<eps)
						{
							x[n] = z0;
							n = n + 1;
							js = 1;
							z = z0 + dblStep / 2.0;
							y = Func(z);
						}
						else if ((y*y0)<0.0)
						{
							z1 = z0;
							y1 = y0;
						}
						else
						{
							z = z0;
							y = y0;
						}
					}
				}
			}
		}
	}

	// ����ʵ������Ŀ
	return(n);
}

//////////////////////////////////////////////////////////////////////
// ������Է���һ��ʵ����ţ�ٷ�
//
// ����ʱ���븲�Ǽ��㷽����˺���f(x)����һ�׵���f'(x)ֵ���麯��
//                  void Func(double x, double y[])
//         y(0) ����f(x)��ֵ
//         y(1) ����f'(x)��ֵ
//
// ������
// 1. double *x - ���������ֵ���²�⣩��������������õ�һ��ʵ��
// 2. int nMaxIt - �ݹ������Ĭ��ֵΪ60
// 3. double eps - ���ȿ��Ʋ�����Ĭ��ֵΪ0.000001
//
// ����ֵ��bool �ͣ�����Ƿ�ɹ�
//////////////////////////////////////////////////////////////////////
bool CNLEquation::GetRootNewton(double* x, int nMaxIt /*= 60*/, double eps /*= 0.000001*/)
{
	int l;
	double y[2], d, p, x0, x1;

	// ����ֵ
	l = nMaxIt;
	x0 = *x;
	Func(x0, y);

	// ��⣬���ƾ���
	d = eps + 1.0;
	while ((d >= eps) && (l != 0))
	{
		if (y[1] == 0.0)
			return false;

		x1 = x0 - y[0] / y[1];
		Func(x1, y);

		d = fabs(x1 - x0);
		p = fabs(y[0]);
		if (p>d)
			d = p;
		x0 = x1;
		l = l - 1;
	}

	*x = x1;

	return true;
}

//////////////////////////////////////////////////////////////////////
// ������Է���һ��ʵ���İ��ؽ������
//
// ����ʱ���븲�Ǽ��㷽����˺���f(x)ֵ���麯��double Func(double x)
//
// ������
// 1. double *x - ���������ֵ���²�⣩��������������õ�һ��ʵ��
// 2. int nMaxIt - �ݹ������Ĭ��ֵΪ60
// 3. double eps - ���ȿ��Ʋ�����Ĭ��ֵΪ0.000001
//
// ����ֵ��bool �ͣ�����Ƿ�ɹ�
//////////////////////////////////////////////////////////////////////
bool CNLEquation::GetRootAitken(double* x, int nMaxIt /*= 60*/, double eps /*= 0.000001*/)
{
	int flag, l;
	double u, v, x0;

	// �������
	l = 0;
	x0 = *x;
	flag = 0;

	// �������
	while ((flag == 0) && (l != nMaxIt))
	{
		l = l + 1;
		u = Func(x0);
		v = Func(u);
		if (fabs(u - v)<eps)
		{
			x0 = v;
			flag = 1;
		}
		else
			x0 = v - (v - u)*(v - u) / (v - 2.0*u + x0);
	}

	*x = x0;

	// �Ƿ���ָ���ĵ��������ڴﵽ��⾫��
	return (nMaxIt > l);
}

//////////////////////////////////////////////////////////////////////
// ������Է���һ��ʵ��������ʽ�ⷨ
//
// ����ʱ���븲�Ǽ��㷽����˺���f(x)ֵ���麯��double Func(double x)
//
// ������
// 1. double *x - ���������ֵ���²�⣩��������������õ�һ��ʵ��
// 2. double eps - ���ȿ��Ʋ�����Ĭ��ֵΪ0.000001
//
// ����ֵ��bool �ͣ�����Ƿ�ɹ�
//////////////////////////////////////////////////////////////////////
bool CNLEquation::GetRootPq(double* x, double eps /*= 0.000001*/)
{
	int i, j, m, it, l;
	double a[10], y[10], z, h, x0, q;

	// �������
	l = 10;
	q = 1.0e+35;
	x0 = *x;
	h = 0.0;

	// ����ʽ���
	while (l != 0)
	{
		l = l - 1;
		j = 0;
		it = l;
		while (j <= 7)
		{
			if (j <= 2)
				z = x0 + 0.1*j;
			else
				z = h;

			y[j] = Func(z);
			h = z;
			if (j == 0)
				a[0] = z;
			else
			{
				m = 0;
				i = 0;
				while ((m == 0) && (i <= j - 1))
				{
					if (fabs(h - a[i]) + 1.0 == 1.0)
						m = 1;
					else
						h = (y[j] - y[i]) / (h - a[i]);

					i = i + 1;
				}
				a[j] = h;
				if (m != 0)
					a[j] = q;
				h = 0.0;
				for (i = j - 1; i >= 0; i--)
				{
					if (fabs(a[i + 1] + h) + 1.0 == 1.0)
						h = q;
					else
						h = -y[i] / (a[i + 1] + h);
				}

				h = h + a[0];
			}

			if (fabs(y[j]) >= eps)
				j = j + 1;
			else
			{
				j = 10;
				l = 0;
			}
		}

		x0 = h;
	}

	*x = h;

	// �Ƿ���10������ʽ�����ʵ����
	return (10>it);
}

//////////////////////////////////////////////////////////////////////
// ��ʵϵ����������ȫ������QR����
//
// ������
// 1. int n - ����ʽ���̵Ĵ���
// 2. double dblCoef[] - һά���飬����Ϊn+1�������ݴ������δ��n�ζ���ʽ���̵�n+1��ϵ��
// 3. double xr[] - һά���飬����Ϊn������n������ʵ��
// 4. double xi[] - һά���飬����Ϊn������n�������鲿
// 5. int nMaxIt - ����������Ĭ��ֵΪ 60
// 6. double eps - ���ȿ��Ʋ�����Ĭ��ֵΪ 0.000001
//
// ����ֵ��bool �ͣ�����Ƿ�ɹ�
//////////////////////////////////////////////////////////////////////
bool CNLEquation::GetRootQr(int n, double dblCoef[], double xr[], double xi[], int nMaxIt /*= 60*/, double eps /*= 0.000001*/)
{
	// ��ʼ������
	CMatrix mtxQ;
	mtxQ.Init(n, n);
	double *q = mtxQ.GetData();

	// ������겮�����
	int j;
	for (j = 0; j <= n - 1; j++)
		q[j] = -dblCoef[n - j - 1] / dblCoef[n];

	for (j = n; j <= n*n - 1; j++)
		q[j] = 0.0;

	for (int i = 0; i <= n - 2; i++)
		q[(i + 1)*n + i] = 1.0;

	// ����겮����������ֵ��������������Ϊ���̵Ľ�
	if (mtxQ.HBergEigenv(xr, xi, nMaxIt, eps))
		return true;

	return false;
}

//////////////////////////////////////////////////////////////////////
// ��ʵϵ����������ȫ������ţ����ɽ��
//
// ������
// 1. int n - ����ʽ���̵Ĵ���
// 2. double dblCoef[] - һά���飬����Ϊn+1�������ݴ������δ��n�ζ���ʽ���̵�n+1��ϵ��
// 3. double xr[] - һά���飬����Ϊn������n������ʵ��
// 4. double xi[] - һά���飬����Ϊn������n�������鲿
//
// ����ֵ��bool �ͣ�����Ƿ�ɹ�
//////////////////////////////////////////////////////////////////////
bool CNLEquation::GetRootNewtonDownHill(int n, double dblCoef[], double xr[], double xi[])
{
	int m, i, jt, k, is, it;
	double t, x, y, x1, y1, dx, dy, p, q, w, dd, dc, c;
	double g, u, v, pq, g1, u1, v1;

	// ��ʼ�ж�
	m = n;
	while ((m>0) && (fabs(dblCoef[m]) + 1.0 == 1.0))
		m = m - 1;

	// ���ʧ��
	if (m <= 0)
		return false;

	for (i = 0; i <= m; i++)
		dblCoef[i] = dblCoef[i] / dblCoef[m];

	for (i = 0; i <= m / 2; i++)
	{
		w = dblCoef[i];
		dblCoef[i] = dblCoef[m - i];
		dblCoef[m - i] = w;
	}

	// �������
	k = m;
	is = 0;
	w = 1.0;
	jt = 1;
	while (jt == 1)
	{
		pq = fabs(dblCoef[k]);
		while (pq<1.0e-12)
		{
			xr[k - 1] = 0.0;
			xi[k - 1] = 0.0;
			k = k - 1;
			if (k == 1)
			{
				xr[0] = -dblCoef[1] * w / dblCoef[0];
				xi[0] = 0.0;

				return true;
			}

			pq = fabs(dblCoef[k]);
		}

		q = log(pq);
		q = q / (1.0*k);
		q = exp(q);
		p = q;
		w = w*p;
		for (i = 1; i <= k; i++)
		{
			dblCoef[i] = dblCoef[i] / q;
			q = q*p;
		}

		x = 0.0001;
		x1 = x;
		y = 0.2;
		y1 = y;
		dx = 1.0;
		g = 1.0e+37;
	l40:
		u = dblCoef[0]; v = 0.0;
		for (i = 1; i <= k; i++)
		{
			p = u*x1;
			q = v*y1;
			pq = (u + v)*(x1 + y1);
			u = p - q + dblCoef[i];
			v = pq - p - q;
		}

		g1 = u*u + v*v;
		if (g1 >= g)
		{
			if (is != 0)
			{
				it = 1;
				g65(&x, &y, &x1, &y1, &dx, &dy, &dd, &dc, &c, &k, &is, &it);
				if (it == 0)
					goto l40;
			}
			else
			{
				g60(&t, &x, &y, &x1, &y1, &dx, &dy, &p, &q, &k, &it);
				if (t >= 1.0e-03)
					goto l40;

				if (g>1.0e-18)
				{
					it = 0;
					g65(&x, &y, &x1, &y1, &dx, &dy, &dd, &dc, &c, &k, &is, &it);
					if (it == 0)
						goto l40;
				}
			}

			g90(xr, xi, dblCoef, &x, &y, &p, &q, &w, &k);
		}
		else
		{
			g = g1;
			x = x1;
			y = y1;
			is = 0;
			if (g <= 1.0e-22)
				g90(xr, xi, dblCoef, &x, &y, &p, &q, &w, &k);
			else
			{
				u1 = k*dblCoef[0];
				v1 = 0.0;
				for (i = 2; i <= k; i++)
				{
					p = u1*x;
					q = v1*y;
					pq = (u1 + v1)*(x + y);
					u1 = p - q + (k - i + 1)*dblCoef[i - 1];
					v1 = pq - p - q;
				}

				p = u1*u1 + v1*v1;
				if (p <= 1.0e-20)
				{
					it = 0;
					g65(&x, &y, &x1, &y1, &dx, &dy, &dd, &dc, &c, &k, &is, &it);
					if (it == 0)
						goto l40;

					g90(xr, xi, dblCoef, &x, &y, &p, &q, &w, &k);
				}
				else
				{
					dx = (u*u1 + v*v1) / p;
					dy = (u1*v - v1*u) / p;
					t = 1.0 + 4.0 / k;
					g60(&t, &x, &y, &x1, &y1, &dx, &dy, &p, &q, &k, &it);
					if (t >= 1.0e-03)
						goto l40;

					if (g>1.0e-18)
					{
						it = 0;
						g65(&x, &y, &x1, &y1, &dx, &dy, &dd, &dc, &c, &k, &is, &it);
						if (it == 0)
							goto l40;
					}

					g90(xr, xi, dblCoef, &x, &y, &p, &q, &w, &k);
				}
			}
		}

		if (k == 1)
			jt = 0;
		else
			jt = 1;
	}

	return true;
}


//////////////////////////////////////////////////////////////////////
// �ڲ�����
//////////////////////////////////////////////////////////////////////
void CNLEquation::g60(double* t, double* x, double* y, double* x1, double* y1, double* dx, double* dy, double* p, double* q, int* k, int* it)
{
	*it = 1;
	while (*it == 1)
	{
		*t = *t / 1.67;
		*it = 0;
		*x1 = *x - (*t)*(*dx);
		*y1 = *y - (*t)*(*dy);
		if (*k >= 50)
		{
			*p = sqrt((*x1)*(*x1) + (*y1)*(*y1));
			*q = exp(85.0 / (*k));
			if (*p >= *q)
				*it = 1;
		}
	}
}

//////////////////////////////////////////////////////////////////////
// �ڲ�����
//////////////////////////////////////////////////////////////////////
void CNLEquation::g90(double xr[], double xi[], double dblCoef[], double* x, double* y, double* p, double* q, double* w, int* k)
{
	int i;

	if (fabs(*y) <= 1.0e-06)
	{
		*p = -(*x);
		*y = 0.0;
		*q = 0.0;
	}
	else
	{
		*p = -2.0*(*x);
		*q = (*x)*(*x) + (*y)*(*y);
		xr[*k - 1] = (*x)*(*w);
		xi[*k - 1] = -(*y)*(*w);
		*k = *k - 1;
	}

	for (i = 1; i <= *k; i++)
	{
		dblCoef[i] = dblCoef[i] - dblCoef[i - 1] * (*p);
		dblCoef[i + 1] = dblCoef[i + 1] - dblCoef[i - 1] * (*q);
	}

	xr[*k - 1] = (*x)*(*w);
	xi[*k - 1] = (*y)*(*w);
	*k = *k - 1;
	if (*k == 1)
	{
		xr[0] = -dblCoef[1] * (*w) / dblCoef[0];
		xi[0] = 0.0;
	}
}

//////////////////////////////////////////////////////////////////////
// �ڲ�����
//////////////////////////////////////////////////////////////////////
void CNLEquation::g65(double* x, double* y, double* x1, double* y1, double* dx, double* dy, double* dd, double* dc, double* c, int* k, int* is, int* it)
{
	if (*it == 0)
	{
		*is = 1;
		*dd = sqrt((*dx)*(*dx) + (*dy)*(*dy));
		if (*dd>1.0)
			*dd = 1.0;
		*dc = 6.28 / (4.5*(*k));
		*c = 0.0;
	}

	while (true)
	{
		*c = *c + (*dc);
		*dx = (*dd)*cos(*c);
		*dy = (*dd)*sin(*c);
		*x1 = *x + *dx;
		*y1 = *y + *dy;
		if (*c <= 6.29)
		{
			*it = 0;
			return;
		}

		*dd = *dd / 1.67;
		if (*dd <= 1.0e-07)
		{
			*it = 1;
			return;
		}

		*c = 0.0;
	}
}

//////////////////////////////////////////////////////////////////////
// ��ϵ����������ȫ������ţ����ɽ��
//
// ������
// 1. int n - ����ʽ���̵Ĵ���
// 2. double ar[] - һά���飬����Ϊn+1�������ݴ������δ��n�ζ���ʽ���̵�n+1��ϵ����ʵ��
// 3. double ai[] - һά���飬����Ϊn+1�������ݴ������δ��n�ζ���ʽ���̵�n+1��ϵ�����鲿
// 4. double xr[] - һά���飬����Ϊn������n������ʵ��
// 5. double xi[] - һά���飬����Ϊn������n�������鲿
//
// ����ֵ��bool �ͣ�����Ƿ�ɹ�
//////////////////////////////////////////////////////////////////////
bool CNLEquation::GetRootNewtonDownHill(int n, double ar[], double ai[], double xr[], double xi[])
{
	int m, i, jt, k, is, it;
	double t, x, y, x1, y1, dx, dy, p, q, w, dd, dc, c;
	double g, u, v, pq, g1, u1, v1;

	// ��ʼ�ж�
	m = n;
	p = sqrt(ar[m] * ar[m] + ai[m] * ai[m]);
	while ((m>0) && (p + 1.0 == 1.0))
	{
		m = m - 1;
		p = sqrt(ar[m] * ar[m] + ai[m] * ai[m]);
	}

	// ���ʧ��
	if (m <= 0)
		return false;

	for (i = 0; i <= m; i++)
	{
		ar[i] = ar[i] / p;
		ai[i] = ai[i] / p;
	}

	for (i = 0; i <= m / 2; i++)
	{
		w = ar[i];
		ar[i] = ar[m - i];
		ar[m - i] = w;
		w = ai[i];
		ai[i] = ai[m - i];
		ai[m - i] = w;
	}

	// �������
	k = m;
	is = 0;
	w = 1.0;
	jt = 1;
	while (jt == 1)
	{
		pq = sqrt(ar[k] * ar[k] + ai[k] * ai[k]);
		while (pq<1.0e-12)
		{
			xr[k - 1] = 0.0;
			xi[k - 1] = 0.0;
			k = k - 1;
			if (k == 1)
			{
				p = ar[0] * ar[0] + ai[0] * ai[0];
				xr[0] = -w*(ar[0] * ar[1] + ai[0] * ai[1]) / p;
				xi[0] = w*(ar[1] * ai[0] - ar[0] * ai[1]) / p;

				return true;
			}

			pq = sqrt(ar[k] * ar[k] + ai[k] * ai[k]);
		}

		q = log(pq);
		q = q / (1.0*k);
		q = exp(q);
		p = q;
		w = w*p;
		for (i = 1; i <= k; i++)
		{
			ar[i] = ar[i] / q;
			ai[i] = ai[i] / q;
			q = q*p;
		}

		x = 0.0001;
		x1 = x;
		y = 0.2;
		y1 = y;
		dx = 1.0;
		g = 1.0e+37;
	l40:
		u = ar[0];
		v = ai[0];
		for (i = 1; i <= k; i++)
		{
			p = u*x1;
			q = v*y1;
			pq = (u + v)*(x1 + y1);
			u = p - q + ar[i];
			v = pq - p - q + ai[i];
		}

		g1 = u*u + v*v;
		if (g1 >= g)
		{
			if (is != 0)
			{
				it = 1;
				g65c(&x, &y, &x1, &y1, &dx, &dy, &dd, &dc, &c, &k, &is, &it);
				if (it == 0)
					goto l40;
			}
			else
			{
				g60c(&t, &x, &y, &x1, &y1, &dx, &dy, &p, &q, &k, &it);
				if (t >= 1.0e-03)
					goto l40;

				if (g>1.0e-18)
				{
					it = 0;
					g65c(&x, &y, &x1, &y1, &dx, &dy, &dd, &dc, &c, &k, &is, &it);
					if (it == 0)
						goto l40;
				}
			}

			g90c(xr, xi, ar, ai, &x, &y, &p, &w, &k);
		}
		else
		{
			g = g1;
			x = x1;
			y = y1;
			is = 0;
			if (g <= 1.0e-22)
				g90c(xr, xi, ar, ai, &x, &y, &p, &w, &k);
			else
			{
				u1 = k*ar[0];
				v1 = ai[0];
				for (i = 2; i <= k; i++)
				{
					p = u1*x;
					q = v1*y;
					pq = (u1 + v1)*(x + y);
					u1 = p - q + (k - i + 1)*ar[i - 1];
					v1 = pq - p - q + (k - i + 1)*ai[i - 1];
				}

				p = u1*u1 + v1*v1;
				if (p <= 1.0e-20)
				{
					it = 0;
					g65c(&x, &y, &x1, &y1, &dx, &dy, &dd, &dc, &c, &k, &is, &it);
					if (it == 0)
						goto l40;

					g90c(xr, xi, ar, ai, &x, &y, &p, &w, &k);
				}
				else
				{
					dx = (u*u1 + v*v1) / p;
					dy = (u1*v - v1*u) / p;
					t = 1.0 + 4.0 / k;
					g60c(&t, &x, &y, &x1, &y1, &dx, &dy, &p, &q, &k, &it);
					if (t >= 1.0e-03)
						goto l40;

					if (g>1.0e-18)
					{
						it = 0;
						g65c(&x, &y, &x1, &y1, &dx, &dy, &dd, &dc, &c, &k, &is, &it);
						if (it == 0)
							goto l40;
					}

					g90c(xr, xi, ar, ai, &x, &y, &p, &w, &k);
				}
			}
		}

		if (k == 1)
			jt = 0;
		else
			jt = 1;
	}

	return true;
}

//////////////////////////////////////////////////////////////////////
// �ڲ�����
//////////////////////////////////////////////////////////////////////
void CNLEquation::g60c(double* t, double* x, double* y, double* x1, double* y1, double* dx, double* dy, double* p, double* q, int* k, int* it)
{
	*it = 1;
	while (*it == 1)
	{
		*t = *t / 1.67;
		*it = 0;
		*x1 = *x - *t*(*dx);
		*y1 = *y - *t*(*dy);
		if (*k >= 30)
		{
			*p = sqrt(*x1*(*x1) + *y1*(*y1));
			*q = exp(75.0 / (*k));
			if (*p >= *q)
				*it = 1;
		}
	}
}

//////////////////////////////////////////////////////////////////////
// �ڲ�����
//////////////////////////////////////////////////////////////////////
void CNLEquation::g90c(double xr[], double xi[], double ar[], double ai[], double* x, double* y, double* p, double* w, int* k)
{
	int i;
	for (i = 1; i <= *k; i++)
	{
		ar[i] = ar[i] + ar[i - 1] * (*x) - ai[i - 1] * (*y);
		ai[i] = ai[i] + ar[i - 1] * (*y) + ai[i - 1] * (*x);
	}

	xr[*k - 1] = *x*(*w);
	xi[*k - 1] = *y*(*w);
	*k = *k - 1;
	if (*k == 1)
	{
		*p = ar[0] * ar[0] + ai[0] * ai[0];
		xr[0] = -*w*(ar[0] * ar[1] + ai[0] * ai[1]) / (*p);
		xi[0] = *w*(ar[1] * ai[0] - ar[0] * ai[1]) / (*p);
	}
}

//////////////////////////////////////////////////////////////////////
// �ڲ�����
//////////////////////////////////////////////////////////////////////
void CNLEquation::g65c(double* x, double* y, double* x1, double* y1, double* dx, double* dy, double* dd, double* dc, double* c, int* k, int* is, int* it)
{
	if (*it == 0)
	{
		*is = 1;
		*dd = sqrt(*dx*(*dx) + *dy*(*dy));
		if (*dd>1.0)
			*dd = 1.0;
		*dc = 6.28 / (4.5*(*k));
		*c = 0.0;
	}

	while (true)
	{
		*c = *c + *dc;
		*dx = *dd*cos(*c);
		*dy = *dd*sin(*c);
		*x1 = *x + *dx;
		*y1 = *y + *dy;
		if (*c <= 6.29)
		{
			*it = 0;
			return;
		}

		*dd = *dd / 1.67;
		if (*dd <= 1.0e-07)
		{
			*it = 1;
			return;
		}

		*c = 0.0;
	}
}

//////////////////////////////////////////////////////////////////////
// ������Է���һ��ʵ�������ؿ��巨
//
// ����ʱ���븲�Ǽ��㷽����˺���f(x)ֵ���麯��double Func(double x)
//
// ������
// 1. double* x - �����ֵ���²�⣩��������õ�ʵ��
// 2. double xStart - ���ȷֲ��Ķ˵��ֵ
// 3. int nControlB - ���Ʋ���
// 4. double eps - ���ƾ��ȣ�Ĭ��ֵΪ0.000001
//
// ����ֵ����
//////////////////////////////////////////////////////////////////////
void CNLEquation::GetRootMonteCarlo(double* x, double xStart, int nControlB, double eps /*= 0.000001*/)
{
	int k;
	double xx, a, y, x1, y1, r;

	// �������
	a = xStart;
	k = 1;
	r = 1.0;

	// ��ֵ
	xx = *x;
	y = Func(xx);

	// ���ȿ������
	while (a >= eps)
	{
		x1 = rnd(&r);

		x1 = -a + 2.0*a*x1;
		x1 = xx + x1;
		y1 = Func(x1);

		k = k + 1;
		if (fabs(y1) >= fabs(y))
		{
			if (k>nControlB)
			{
				k = 1;
				a = a / 2.0;
			}
		}
		else
		{
			k = 1;
			xx = x1;
			y = y1;
			if (fabs(y)<eps)
			{
				*x = xx;
				return;
			}
		}
	}

	*x = xx;
}

//////////////////////////////////////////////////////////////////////
// �ڲ�����
//////////////////////////////////////////////////////////////////////
double CNLEquation::rnd(double* r)
{
	int m;
	double s, u, v, p;

	s = 65536.0;
	u = 2053.0;
	v = 13849.0;
	m = (int)(*r / s);
	*r = *r - m*s;
	*r = u*(*r) + v;
	m = (int)(*r / s);
	*r = *r - m*s; p = *r / s;

	return(p);
}

//////////////////////////////////////////////////////////////////////
// ��ʵ�����򸴺�������һ�����������ؿ��巨
//
// ����ʱ���븲�Ǽ��㷽����˺�����ģֵ||f(x, y)||���麯��
//           double Func(double x, double y)
//
// ������
// 1. double* x - �����ֵ���²�⣩��ʵ����������õĸ���ʵ��
// 2. double* y - �����ֵ���²�⣩���鲿��������õĸ����鲿
// 3. double xStart - ���ȷֲ��Ķ˵��ֵ
// 4. int nControlB - ���Ʋ���
// 5. double eps - ���ƾ��ȣ�Ĭ��ֵΪ0.000001
// 6. double xi[] - һά���飬����Ϊn������n�������鲿
//
// ����ֵ����
//////////////////////////////////////////////////////////////////////
void CNLEquation::GetRootMonteCarlo(double* x, double* y, double xStart, int nControlB, double eps /*= 0.000001*/)
{
	int k;
	double xx, yy, a, r, z, x1, y1, z1;

	// ����������ֵ
	a = xStart;
	k = 1;
	r = 1.0;
	xx = *x;
	yy = *y;
	z = Func(xx, yy);

	// ���ȿ������
	while (a >= eps)
	{
		x1 = -a + 2.0*a*rnd(&r);
		x1 = xx + x1;
		y1 = -a + 2.0*a*rnd(&r);
		y1 = yy + y1;

		z1 = Func(x1, y1);

		k = k + 1;
		if (z1 >= z)
		{
			if (k>nControlB)
			{
				k = 1;
				a = a / 2.0;
			}
		}
		else
		{
			k = 1;
			xx = x1;
			yy = y1;
			z = z1;
			if (z<eps)
			{
				*x = xx;
				*y = yy;
				return;
			}
		}
	}

	*x = xx;
	*y = yy;
}

//////////////////////////////////////////////////////////////////////
// ������Է�����һ��ʵ�����ݶȷ�
//
// ����ʱ���븲�Ǽ��㷽����˺���f(x)ֵ����ƫ����ֵ���麯��
//         double Func(double x[], double y[])
//
// ������
// 1. int n - ���̵ĸ�����Ҳ��δ֪���ĸ���
// 2. double x[] - һά���飬����Ϊn�����һ���ֵx0, x1, ��, xn-1��
//                 ����ʱ��ŷ������һ��ʵ��
// 3. int nMaxIt - ����������Ĭ��ֵΪ60
// 4. double eps - ���ƾ��ȣ�Ĭ��ֵΪ0.000001
//
// ����ֵ��bool �ͣ�����Ƿ�ɹ�
//////////////////////////////////////////////////////////////////////
bool CNLEquation::GetRootsetGrad(int n, double x[], int nMaxIt /*= 60*/, double eps /*= 0.000001*/)
{
	int l, j;
	double f, d, s, *y;

	y = new double[n];

	l = nMaxIt;
	f = Func(x, y);

	// ���ƾ��ȣ��������
	while (f >= eps)
	{
		l = l - 1;
		if (l == 0)
		{
			delete[] y;
			return true;
		}

		d = 0.0;
		for (j = 0; j <= n - 1; j++)
			d = d + y[j] * y[j];
		if (d + 1.0 == 1.0)
		{
			delete[] y;
			return false;
		}

		s = f / d;
		for (j = 0; j <= n - 1; j++)
			x[j] = x[j] - s*y[j];

		f = Func(x, y);
	}

	delete[] y;

	// �Ƿ�����Ч���������ڴﵽ����
	return (nMaxIt>l);
}

//////////////////////////////////////////////////////////////////////
// ������Է�����һ��ʵ������ţ�ٷ�
//
// ����ʱ���븲�Ǽ��㷽����˺���f(x)ֵ����ƫ����ֵ���麯��
//         double Func(double x[], double y[])
//
// ������
// 1. int n - ���̵ĸ�����Ҳ��δ֪���ĸ���
// 2. double x[] - һά���飬����Ϊn�����һ���ֵx0, x1, ��, xn-1��
//                 ����ʱ��ŷ������һ��ʵ��
// 3. double t - ����h��С�ı�����0<t<1
// 4. double h - ������ֵ
// 3. int nMaxIt - ����������Ĭ��ֵΪ60
// 4. double eps - ���ƾ��ȣ�Ĭ��ֵΪ0.000001
//
// ����ֵ��bool �ͣ�����Ƿ�ɹ�
//////////////////////////////////////////////////////////////////////
bool CNLEquation::GetRootsetNewton(int n, double x[], double t, double h, int nMaxIt /*= 60*/, double eps /*= 0.000001*/)
{
	int i, j, l;
	double am, z, beta, d, *y, *a, *b;

	y = new double[n];

	// �������
	CMatrix mtxCoef(n, n);
	CMatrix mtxConst(n, 1);
	a = mtxCoef.GetData();
	b = mtxConst.GetData();

	// �������
	l = nMaxIt;
	am = 1.0 + eps;
	while (am >= eps)
	{
		Func(x, b);

		am = 0.0;
		for (i = 0; i <= n - 1; i++)
		{
			z = fabs(b[i]);
			if (z>am)
				am = z;
		}

		if (am >= eps)
		{
			l = l - 1;
			if (l == 0)
			{
				delete[] y;
				return false;
			}

			for (j = 0; j <= n - 1; j++)
			{
				z = x[j];
				x[j] = x[j] + h;

				Func(x, y);

				for (i = 0; i <= n - 1; i++)
					a[i*n + j] = y[i];

				x[j] = z;
			}

			// ����ȫѡ��Ԫ��˹��Ԫ��
			CLEquation leqs(mtxCoef, mtxConst);
			CMatrix mtxResult;
			if (!leqs.GetRootsetGauss(mtxResult))
			{
				delete[] y;
				return false;
			}

			mtxConst = mtxResult;
			b = mtxConst.GetData();

			beta = 1.0;
			for (i = 0; i <= n - 1; i++)
				beta = beta - b[i];

			if (beta == 0.0)
			{
				delete[] y;
				return false;
			}

			d = h / beta;
			for (i = 0; i <= n - 1; i++)
				x[i] = x[i] - d*b[i];

			h = t*h;
		}
	}

	delete[] y;

	// �Ƿ�����Ч���������ڴﵽ����
	return (nMaxIt>l);
}

//////////////////////////////////////////////////////////////////////
// ������Է�������С���˽�Ĺ����淨
//
// ����ʱ��1. �븲�Ǽ��㷽����˺���f(x)ֵ����ƫ����ֵ���麯��
//             double Func(double x[], double y[])
//         2. �븲�Ǽ����ſɱȾ��������麯��
//             double FuncMJ(int n, double x[], double y[])
//
// ������
// 1. int m - ���̵ĸ���
// 2. int n - δ֪���ĸ���
// 3. double x[] - һά���飬����Ϊn�����һ���ֵx0, x1, ��, xn-1��Ҫ��
//                 ��ȫΪ0������ʱ��ŷ��������С���˽⣬��m=nʱ������
//                 �����Է�����Ľ�
// 4. double eps1 - ��С���˽�ľ��ȿ��ƾ��ȣ�Ĭ��ֵΪ0.000001
// 5. double eps2 - ����ֵ�ֽ�ľ��ȿ��ƾ��ȣ�Ĭ��ֵΪ0.000001
//
// ����ֵ��bool �ͣ�����Ƿ�ɹ�
//////////////////////////////////////////////////////////////////////
bool CNLEquation::GetRootsetGinv(int m, int n, double x[], double eps1 /*= 0.000001*/, double eps2 /*= 0.000001*/)
{
	int i, j, k, l, kk, jt;
	double y[10], b[10], alpha, z, h2, y1, y2, y3, y0, h1;
	double *p, *d, *dx, *w;

	// ���Ʋ���
	int ka = max(m, n) + 1;
	w = new double[ka];

	// �趨��������Ϊ60���������
	l = 60;
	alpha = 1.0;
	while (l>0)
	{
		CMatrix mtxP(m, n), mtxD(m, 1);
		p = mtxP.GetData();
		d = mtxD.GetData();

		Func(x, d);
		FuncMJ(n, x, p);

		// �������Է�����
		CLEquation leqs(mtxP, mtxD);
		// ��ʱ����
		CMatrix mtxAP, mtxU, mtxV;
		// �����
		CMatrix mtxDX;
		// ���ڹ��������С���˽�
		if (!leqs.GetRootsetGinv(mtxDX, mtxAP, mtxU, mtxV, eps2))
		{
			delete[] w;
			return false;
		}

		dx = mtxDX.GetData();

		j = 0;
		jt = 1;
		h2 = 0.0;
		while (jt == 1)
		{
			jt = 0;
			if (j <= 2)
				z = alpha + 0.01*j;
			else
				z = h2;

			for (i = 0; i <= n - 1; i++)
				w[i] = x[i] - z*dx[i];

			Func(w, d);

			y1 = 0.0;
			for (i = 0; i <= m - 1; i++)
				y1 = y1 + d[i] * d[i];
			for (i = 0; i <= n - 1; i++)
				w[i] = x[i] - (z + 0.00001)*dx[i];

			Func(w, d);

			y2 = 0.0;
			for (i = 0; i <= m - 1; i++)
				y2 = y2 + d[i] * d[i];

			y0 = (y2 - y1) / 0.00001;

			if (fabs(y0)>1.0e-10)
			{
				h1 = y0; h2 = z;
				if (j == 0)
				{
					y[0] = h1;
					b[0] = h2;
				}
				else
				{
					y[j] = h1;
					kk = 0;
					k = 0;
					while ((kk == 0) && (k <= j - 1))
					{
						y3 = h2 - b[k];
						if (fabs(y3) + 1.0 == 1.0)
							kk = 1;
						else
							h2 = (h1 - y[k]) / y3;

						k = k + 1;
					}

					b[j] = h2;
					if (kk != 0)
						b[j] = 1.0e+35;

					h2 = 0.0;
					for (k = j - 1; k >= 0; k--)
						h2 = -y[k] / (b[k + 1] + h2);

					h2 = h2 + b[0];
				}

				j = j + 1;
				if (j <= 7)
					jt = 1;
				else
					z = h2;
			}
		}

		alpha = z;
		y1 = 0.0;
		y2 = 0.0;
		for (i = 0; i <= n - 1; i++)
		{
			dx[i] = -alpha*dx[i];
			x[i] = x[i] + dx[i];
			y1 = y1 + fabs(dx[i]);
			y2 = y2 + fabs(x[i]);
		}

		// ���ɹ�
		if (y1<eps1*y2)
		{
			delete[] w;
			return true;
		}

		l = l - 1;
	}

	delete[] w;

	// ���ʧ��
	return false;
}

//////////////////////////////////////////////////////////////////////
// ������Է�����һ��ʵ�������ؿ��巨
//
// ����ʱ���븲�Ǽ��㷽�����ģ����ֵ||F||���麯��
//         double Func(int n, double x[])
//         �䷵��ֵΪSqr(f1*f1 + f2*f2 + �� + fn*fn)
//
// ������
// 1. double x[] - һά���飬����Ϊn�����һ���ֵ������ʱ��ŷ��̵�һ��ʵ��
// 2. double xStart - ���ȷֲ��Ķ˵��ֵ
// 3. int nControlB - ���Ʋ���
// 4. double eps - ���ƾ��ȣ�Ĭ��ֵΪ0.000001
//
// ����ֵ����
//////////////////////////////////////////////////////////////////////
void CNLEquation::GetRootsetMonteCarlo(int n, double x[], double xStart, int nControlB, double eps /*= 0.000001*/)
{
	int k, i;
	double a, r, *y, z, z1;

	y = new double[n];

	// ��ֵ
	a = xStart;
	k = 1;
	r = 1.0;

	z = Func(n, x);

	// �þ��ȿ��Ƶ������
	while (a >= eps)
	{
		for (i = 0; i <= n - 1; i++)
			y[i] = -a + 2.0*a*rnd(&r) + x[i];

		z1 = Func(n, y);

		k = k + 1;
		if (z1 >= z)
		{
			if (k>nControlB)
			{
				k = 1;
				a = a / 2.0;
			}
		}
		else
		{
			k = 1;
			for (i = 0; i <= n - 1; i++)
				x[i] = y[i];

			// ���ɹ�
			z = z1;
			if (z<eps)
			{
				delete[] y;
				return;
			}
		}
	}

	delete[] y;
}
