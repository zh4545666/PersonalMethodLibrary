#include "stdafx.h"
#include "NLEquation.h"

using namespace PersonalMethod;

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////
// 基本构造函数
//////////////////////////////////////////////////////////////////////
CNLEquation::CNLEquation()
{
}

//////////////////////////////////////////////////////////////////////
// 析构函数
//////////////////////////////////////////////////////////////////////
CNLEquation::~CNLEquation()
{
}

//////////////////////////////////////////////////////////////////////
// 求非线性方程实根的对分法
//
// 调用时，须覆盖计算方程左端函数f(x)值的虚函数double Func(double x)
//
// 参数：
// 1. int nNumRoots - 在[xStart, xEnd]内实根个数的预估值
// 2. double x[] - 一维数组，长度为m。返回在区间[xStart, xEnd]内搜索到
//                 的实根，实根个数由函数值返回
// 3. double xStart - 求根区间的左端点
// 4. double xEnd - 求根区间的右端点
// 5. double dblStep - 搜索求根时采用的步长
// 6. double eps - 精度控制参数，默认值为0.000001
//
// 返回值：int 型，求得的实根的数目
//////////////////////////////////////////////////////////////////////
int CNLEquation::GetRootBisect(int nNumRoots, double x[], double xStart, double xEnd, double dblStep, double eps /*= 0.000001*/)
{
	int n, js;
	double z, y, z1, y1, z0, y0;

	// 根的个数清0
	n = 0;

	// 从左端点开始搜索
	z = xStart;
	y = Func(z);

	// 循环求解
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

	// 返回实根的数目
	return(n);
}

//////////////////////////////////////////////////////////////////////
// 求非线性方程一个实根的牛顿法
//
// 调用时，须覆盖计算方程左端函数f(x)及其一阶导数f'(x)值的虚函数
//                  void Func(double x, double y[])
//         y(0) 返回f(x)的值
//         y(1) 返回f'(x)的值
//
// 参数：
// 1. double *x - 传入迭代初值（猜测解），返回在区间求得的一个实根
// 2. int nMaxIt - 递归次数，默认值为60
// 3. double eps - 精度控制参数，默认值为0.000001
//
// 返回值：bool 型，求解是否成功
//////////////////////////////////////////////////////////////////////
bool CNLEquation::GetRootNewton(double* x, int nMaxIt /*= 60*/, double eps /*= 0.000001*/)
{
	int l;
	double y[2], d, p, x0, x1;

	// 条件值
	l = nMaxIt;
	x0 = *x;
	Func(x0, y);

	// 求解，控制精度
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
// 求非线性方程一个实根的埃特金迭代法
//
// 调用时，须覆盖计算方程左端函数f(x)值的虚函数double Func(double x)
//
// 参数：
// 1. double *x - 传入迭代初值（猜测解），返回在区间求得的一个实根
// 2. int nMaxIt - 递归次数，默认值为60
// 3. double eps - 精度控制参数，默认值为0.000001
//
// 返回值：bool 型，求解是否成功
//////////////////////////////////////////////////////////////////////
bool CNLEquation::GetRootAitken(double* x, int nMaxIt /*= 60*/, double eps /*= 0.000001*/)
{
	int flag, l;
	double u, v, x0;

	// 求解条件
	l = 0;
	x0 = *x;
	flag = 0;

	// 迭代求解
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

	// 是否在指定的迭代次数内达到求解精度
	return (nMaxIt > l);
}

//////////////////////////////////////////////////////////////////////
// 求非线性方程一个实根的连分式解法
//
// 调用时，须覆盖计算方程左端函数f(x)值的虚函数double Func(double x)
//
// 参数：
// 1. double *x - 传入迭代初值（猜测解），返回在区间求得的一个实根
// 2. double eps - 精度控制参数，默认值为0.000001
//
// 返回值：bool 型，求解是否成功
//////////////////////////////////////////////////////////////////////
bool CNLEquation::GetRootPq(double* x, double eps /*= 0.000001*/)
{
	int i, j, m, it, l;
	double a[10], y[10], z, h, x0, q;

	// 求解条件
	l = 10;
	q = 1.0e+35;
	x0 = *x;
	h = 0.0;

	// 连分式求解
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

	// 是否在10阶连分式内求的实根？
	return (10>it);
}

//////////////////////////////////////////////////////////////////////
// 求实系数代数方程全部根的QR方法
//
// 参数：
// 1. int n - 多项式方程的次数
// 2. double dblCoef[] - 一维数组，长度为n+1，按降幂次序依次存放n次多项式方程的n+1个系数
// 3. double xr[] - 一维数组，长度为n，返回n个根的实部
// 4. double xi[] - 一维数组，长度为n，返回n个根的虚部
// 5. int nMaxIt - 迭代次数，默认值为 60
// 6. double eps - 精度控制参数，默认值为 0.000001
//
// 返回值：bool 型，求解是否成功
//////////////////////////////////////////////////////////////////////
bool CNLEquation::GetRootQr(int n, double dblCoef[], double xr[], double xi[], int nMaxIt /*= 60*/, double eps /*= 0.000001*/)
{
	// 初始化矩阵
	CMatrix mtxQ;
	mtxQ.Init(n, n);
	double *q = mtxQ.GetData();

	// 构造赫申伯格矩阵
	int j;
	for (j = 0; j <= n - 1; j++)
		q[j] = -dblCoef[n - j - 1] / dblCoef[n];

	for (j = n; j <= n*n - 1; j++)
		q[j] = 0.0;

	for (int i = 0; i <= n - 2; i++)
		q[(i + 1)*n + i] = 1.0;

	// 求赫申伯格矩阵的特征值和特征向量，即为方程的解
	if (mtxQ.HBergEigenv(xr, xi, nMaxIt, eps))
		return true;

	return false;
}

//////////////////////////////////////////////////////////////////////
// 求实系数代数方程全部根的牛顿下山法
//
// 参数：
// 1. int n - 多项式方程的次数
// 2. double dblCoef[] - 一维数组，长度为n+1，按降幂次序依次存放n次多项式方程的n+1个系数
// 3. double xr[] - 一维数组，长度为n，返回n个根的实部
// 4. double xi[] - 一维数组，长度为n，返回n个根的虚部
//
// 返回值：bool 型，求解是否成功
//////////////////////////////////////////////////////////////////////
bool CNLEquation::GetRootNewtonDownHill(int n, double dblCoef[], double xr[], double xi[])
{
	int m, i, jt, k, is, it;
	double t, x, y, x1, y1, dx, dy, p, q, w, dd, dc, c;
	double g, u, v, pq, g1, u1, v1;

	// 初始判断
	m = n;
	while ((m>0) && (fabs(dblCoef[m]) + 1.0 == 1.0))
		m = m - 1;

	// 求解失败
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

	// 迭代求解
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
// 内部函数
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
// 内部函数
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
// 内部函数
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
// 求复系数代数方程全部根的牛顿下山法
//
// 参数：
// 1. int n - 多项式方程的次数
// 2. double ar[] - 一维数组，长度为n+1，按降幂次序依次存放n次多项式方程的n+1个系数的实部
// 3. double ai[] - 一维数组，长度为n+1，按降幂次序依次存放n次多项式方程的n+1个系数的虚部
// 4. double xr[] - 一维数组，长度为n，返回n个根的实部
// 5. double xi[] - 一维数组，长度为n，返回n个根的虚部
//
// 返回值：bool 型，求解是否成功
//////////////////////////////////////////////////////////////////////
bool CNLEquation::GetRootNewtonDownHill(int n, double ar[], double ai[], double xr[], double xi[])
{
	int m, i, jt, k, is, it;
	double t, x, y, x1, y1, dx, dy, p, q, w, dd, dc, c;
	double g, u, v, pq, g1, u1, v1;

	// 初始判断
	m = n;
	p = sqrt(ar[m] * ar[m] + ai[m] * ai[m]);
	while ((m>0) && (p + 1.0 == 1.0))
	{
		m = m - 1;
		p = sqrt(ar[m] * ar[m] + ai[m] * ai[m]);
	}

	// 求解失败
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

	// 迭代求解
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
// 内部函数
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
// 内部函数
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
// 内部函数
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
// 求非线性方程一个实根的蒙特卡洛法
//
// 调用时，须覆盖计算方程左端函数f(x)值的虚函数double Func(double x)
//
// 参数：
// 1. double* x - 传入初值（猜测解），返回求得的实根
// 2. double xStart - 均匀分布的端点初值
// 3. int nControlB - 控制参数
// 4. double eps - 控制精度，默认值为0.000001
//
// 返回值：无
//////////////////////////////////////////////////////////////////////
void CNLEquation::GetRootMonteCarlo(double* x, double xStart, int nControlB, double eps /*= 0.000001*/)
{
	int k;
	double xx, a, y, x1, y1, r;

	// 求解条件
	a = xStart;
	k = 1;
	r = 1.0;

	// 初值
	xx = *x;
	y = Func(xx);

	// 精度控制求解
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
// 内部函数
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
// 求实函数或复函数方程一个复根的蒙特卡洛法
//
// 调用时，须覆盖计算方程左端函数的模值||f(x, y)||的虚函数
//           double Func(double x, double y)
//
// 参数：
// 1. double* x - 传入初值（猜测解）的实部，返回求得的根的实部
// 2. double* y - 传入初值（猜测解）的虚部，返回求得的根的虚部
// 3. double xStart - 均匀分布的端点初值
// 4. int nControlB - 控制参数
// 5. double eps - 控制精度，默认值为0.000001
// 6. double xi[] - 一维数组，长度为n，返回n个根的虚部
//
// 返回值：无
//////////////////////////////////////////////////////////////////////
void CNLEquation::GetRootMonteCarlo(double* x, double* y, double xStart, int nControlB, double eps /*= 0.000001*/)
{
	int k;
	double xx, yy, a, r, z, x1, y1, z1;

	// 求解条件与初值
	a = xStart;
	k = 1;
	r = 1.0;
	xx = *x;
	yy = *y;
	z = Func(xx, yy);

	// 精度控制求解
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
// 求非线性方程组一组实根的梯度法
//
// 调用时，须覆盖计算方程左端函数f(x)值及其偏导数值的虚函数
//         double Func(double x[], double y[])
//
// 参数：
// 1. int n - 方程的个数，也是未知数的个数
// 2. double x[] - 一维数组，长度为n，存放一组初值x0, x1, …, xn-1，
//                 返回时存放方程组的一组实根
// 3. int nMaxIt - 迭代次数，默认值为60
// 4. double eps - 控制精度，默认值为0.000001
//
// 返回值：bool 型，求解是否成功
//////////////////////////////////////////////////////////////////////
bool CNLEquation::GetRootsetGrad(int n, double x[], int nMaxIt /*= 60*/, double eps /*= 0.000001*/)
{
	int l, j;
	double f, d, s, *y;

	y = new double[n];

	l = nMaxIt;
	f = Func(x, y);

	// 控制精度，迭代求解
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

	// 是否在有效迭代次数内达到精度
	return (nMaxIt>l);
}

//////////////////////////////////////////////////////////////////////
// 求非线性方程组一组实根的拟牛顿法
//
// 调用时，须覆盖计算方程左端函数f(x)值及其偏导数值的虚函数
//         double Func(double x[], double y[])
//
// 参数：
// 1. int n - 方程的个数，也是未知数的个数
// 2. double x[] - 一维数组，长度为n，存放一组初值x0, x1, …, xn-1，
//                 返回时存放方程组的一组实根
// 3. double t - 控制h大小的变量，0<t<1
// 4. double h - 增量初值
// 3. int nMaxIt - 迭代次数，默认值为60
// 4. double eps - 控制精度，默认值为0.000001
//
// 返回值：bool 型，求解是否成功
//////////////////////////////////////////////////////////////////////
bool CNLEquation::GetRootsetNewton(int n, double x[], double t, double h, int nMaxIt /*= 60*/, double eps /*= 0.000001*/)
{
	int i, j, l;
	double am, z, beta, d, *y, *a, *b;

	y = new double[n];

	// 构造矩阵
	CMatrix mtxCoef(n, n);
	CMatrix mtxConst(n, 1);
	a = mtxCoef.GetData();
	b = mtxConst.GetData();

	// 迭代求解
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

			// 调用全选主元高斯消元法
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

	// 是否在有效迭代次数内达到精度
	return (nMaxIt>l);
}

//////////////////////////////////////////////////////////////////////
// 求非线性方程组最小二乘解的广义逆法
//
// 调用时，1. 须覆盖计算方程左端函数f(x)值及其偏导数值的虚函数
//             double Func(double x[], double y[])
//         2. 须覆盖计算雅可比矩阵函数的虚函数
//             double FuncMJ(int n, double x[], double y[])
//
// 参数：
// 1. int m - 方程的个数
// 2. int n - 未知数的个数
// 3. double x[] - 一维数组，长度为n，存放一组初值x0, x1, …, xn-1，要求
//                 不全为0，返回时存放方程组的最小二乘解，当m=n时，即是
//                 非线性方程组的解
// 4. double eps1 - 最小二乘解的精度控制精度，默认值为0.000001
// 5. double eps2 - 奇异值分解的精度控制精度，默认值为0.000001
//
// 返回值：bool 型，求解是否成功
//////////////////////////////////////////////////////////////////////
bool CNLEquation::GetRootsetGinv(int m, int n, double x[], double eps1 /*= 0.000001*/, double eps2 /*= 0.000001*/)
{
	int i, j, k, l, kk, jt;
	double y[10], b[10], alpha, z, h2, y1, y2, y3, y0, h1;
	double *p, *d, *dx, *w;

	// 控制参数
	int ka = max(m, n) + 1;
	w = new double[ka];

	// 设定迭代次数为60，迭代求解
	l = 60;
	alpha = 1.0;
	while (l>0)
	{
		CMatrix mtxP(m, n), mtxD(m, 1);
		p = mtxP.GetData();
		d = mtxD.GetData();

		Func(x, d);
		FuncMJ(n, x, p);

		// 构造线性方程组
		CLEquation leqs(mtxP, mtxD);
		// 临时矩阵
		CMatrix mtxAP, mtxU, mtxV;
		// 解矩阵
		CMatrix mtxDX;
		// 基于广义逆的最小二乘解
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

		// 求解成功
		if (y1<eps1*y2)
		{
			delete[] w;
			return true;
		}

		l = l - 1;
	}

	delete[] w;

	// 求解失败
	return false;
}

//////////////////////////////////////////////////////////////////////
// 求非线性方程组一组实根的蒙特卡洛法
//
// 调用时，须覆盖计算方程左端模函数值||F||的虚函数
//         double Func(int n, double x[])
//         其返回值为Sqr(f1*f1 + f2*f2 + … + fn*fn)
//
// 参数：
// 1. double x[] - 一维数组，长度为n，存放一组初值，返回时存放方程的一组实根
// 2. double xStart - 均匀分布的端点初值
// 3. int nControlB - 控制参数
// 4. double eps - 控制精度，默认值为0.000001
//
// 返回值：无
//////////////////////////////////////////////////////////////////////
void CNLEquation::GetRootsetMonteCarlo(int n, double x[], double xStart, int nControlB, double eps /*= 0.000001*/)
{
	int k, i;
	double a, r, *y, z, z1;

	y = new double[n];

	// 初值
	a = xStart;
	k = 1;
	r = 1.0;

	z = Func(n, x);

	// 用精度控制迭代求解
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

			// 求解成功
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
