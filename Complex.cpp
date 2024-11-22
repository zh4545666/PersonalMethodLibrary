#include "stdafx.h"
#include "Complex.h"

using namespace PersonalMethod;

//////////////////////////////////////////////////////////////////////
// 基本构造函数
//////////////////////////////////////////////////////////////////////
CComplex::CComplex()
{
	m_dblX = 0.0;
	m_dblY = 0.0;
}

//////////////////////////////////////////////////////////////////////
// 指定值构造函数
//
// 参数：
// 1. double dblX - 指定的实部
// 2. double dblY - 指定的虚部
//////////////////////////////////////////////////////////////////////
CComplex::CComplex(double dblX, double dblY)
{
	m_dblX = dblX;
	m_dblY = dblY;
}

//////////////////////////////////////////////////////////////////////
// 拷贝构造函数
//
// 参数：
// 1. const CComplex& other - 源复数
//////////////////////////////////////////////////////////////////////
CComplex::CComplex(const CComplex& other)
{
	m_dblX = other.m_dblX;
	m_dblY = other.m_dblY;
}

//////////////////////////////////////////////////////////////////////
// 指定复数的实部和虚部
//
// 参数：
// 1. double dblX - 复数的实部
// 2. double dblY - 复数的虚部
//////////////////////////////////////////////////////////////////////
void CComplex::SetData(double dblX, double dblY)
{
	m_dblX = dblX;
	m_dblY = dblY;
}

//////////////////////////////////////////////////////////////////////
// 指定复数的实部
//
// 参数：
// 1. double dblX - 复数的实部
//////////////////////////////////////////////////////////////////////
void CComplex::SetReal(double dblX)
{
	m_dblX = dblX;
}

//////////////////////////////////////////////////////////////////////
// 指定复数的虚部
//
// 参数：
// 1. double dblX - 复数的虚部
//////////////////////////////////////////////////////////////////////
void CComplex::SetImag(double dblY)
{
	m_dblY = dblY;
}

//////////////////////////////////////////////////////////////////////
// 取复数的实部
//
// 参数：  无
//
// 返回值：double 型，复数的实部
//////////////////////////////////////////////////////////////////////
double CComplex::GetReal()
{
	return m_dblX;
}

//////////////////////////////////////////////////////////////////////
// 取复数的虚部
//
// 参数：  无
//
// 返回值：double 型，复数的虚部
//////////////////////////////////////////////////////////////////////
double CComplex::GetImag()
{
	return m_dblY;
}

//////////////////////////////////////////////////////////////////////
// 将复数转化为"a+bj"形式的字符串
//
// 参数：  无
//
// 返回值：string 对象，"a+bj"形式的字符串
//////////////////////////////////////////////////////////////////////
string CComplex::ToString(bool bflag) const
{
	string s;
	ostringstream ostr;
	if (m_dblX != 0.0)
	{
		if (bflag)
		{
			if (m_dblY > 0)
				ostr << m_dblX << " + " << m_dblY << "j";
			else if (m_dblY < 0)
				ostr << m_dblX << " - " << abs(m_dblY) << "j";
			else
				ostr << m_dblX;
		}
		else
		{
			if (m_dblY > 0)
				ostr << m_dblX << "+" << m_dblY << "j";
			else if (m_dblY < 0)
				ostr << m_dblX << "-" << abs(m_dblY) << "j";
			else
				ostr << m_dblX;
		}
	}
	else
	{
		if (m_dblY > 0)
			ostr << m_dblY << "j";
		else if (m_dblY < 0)
			ostr << "-" << abs(m_dblY) << "j";
		else
			ostr << m_dblX;
	}
	s = ostr.str();
	return s;
}

//////////////////////////////////////////////////////////////////////
// 将"a,b"形式的字符串转化为复数，以a为复数的实部，b为复数的虚部
//
// 参数：
// 1. string s - "a,b"形式的字符串，a为复数的实部，b为复数的虚部
// 2. const string& sDelim - a, b之间的分隔符，默认为空格
//
// 返回值：无
//////////////////////////////////////////////////////////////////////
void CComplex::FromString(string s, const string& sDelim /*= " "*/)
{
	m_dblX = 0;
	m_dblY = 0;
	if (sDelim == "i" || sDelim == "j")
		return;

	if (sDelim.empty())
	{
		int pos = s.find('+');
		if (pos < 0)
		{
			pos = s.rfind('-');
			if (pos == 0)
				pos = -1;
		}

		if (pos < 0)
		{
			string s1 = s;
			double dblX = atof(s1.c_str());

			int pos1 = s.find('i');
			int pos2 = s.find('j');
			if (pos1 >= 0 || pos2 >= 0)
				m_dblY += dblX;
			else
				m_dblX += dblX;
		}
		else
		{
			string s1 = s.substr(0, pos);
			string s2 = s.substr(pos + 1, s.length());
			CTokenizer::TrimRightSpace(s1);
			CTokenizer::TrimLeftSpace(s1);
			double dblX = atof(s1.c_str());
			CTokenizer::TrimRightSpace(s2);
			CTokenizer::TrimLeftSpace(s2);
			double dblY = atof(s2.c_str());

			if (s[pos] == '-')
				dblY *= -1;

			int pos1 = s1.find('i');
			int pos2 = s1.find('j');
			if (pos1 >= 0 || pos2 >= 0)
				m_dblY += dblX;
			else
				m_dblX += dblX;
			pos1 = s2.find('i');
			pos2 = s2.find('j');
			if (pos1 >= 0 || pos2 >= 0)
				m_dblY += dblY;
			else
				m_dblX += dblY;
		}
	}
	else
	{
		int pos = s.find(sDelim);
		if (pos >= 0)
		{
			string s1 = s.substr(0, pos);
			string s2 = s.substr(pos + sDelim.size(), s.length());
			CTokenizer::TrimRightSpace(s1);
			CTokenizer::TrimLeftSpace(s1);
			m_dblX = atof(s1.c_str());
			CTokenizer::TrimRightSpace(s2);
			CTokenizer::TrimLeftSpace(s2);
			m_dblY = atof(s2.c_str());
		}
		else
		{
			CTokenizer::TrimRightSpace(s);
			CTokenizer::TrimLeftSpace(s);
			m_dblX = atof(s.c_str());
		}
	}

}

//////////////////////////////////////////////////////////////////////
// 重载运算符==，比较两个复数是否相等
//
// 参数：
// 1. const CComplex& cpxX - 用于比较的复数
//
// 返回值：bool型，相等则为true，否则为false
//////////////////////////////////////////////////////////////////////
bool CComplex::operator==(const CComplex& cpxX) const
{
	return (m_dblX == cpxX.m_dblX && m_dblY == cpxX.m_dblY);
}

//////////////////////////////////////////////////////////////////////
// 重载运算符!=，比较两个复数是否不等
//
// 参数：
// 1. const CComplex& cpxX - 用于比较的复数
//
// 返回值：bool型，不相等则为true，相等为false
//////////////////////////////////////////////////////////////////////
bool CComplex::operator!=(const CComplex& cpxX) const
{
	return (m_dblX != cpxX.m_dblX || m_dblY != cpxX.m_dblY);
}

//////////////////////////////////////////////////////////////////////
// 重载运算符=，给复数赋值
//
// 参数：
// 1. const CComplex& cpxX - 用于给复数赋值的源复数
//
// 返回值：CComplex型的引用，所引用的复数与cpxX相等
//////////////////////////////////////////////////////////////////////
CComplex& CComplex::operator=(const CComplex& cpxX)
{
	m_dblX = cpxX.m_dblX;
	m_dblY = cpxX.m_dblY;

	return *this;
}

//////////////////////////////////////////////////////////////////////
// 重载运算符+，实现复数的加法
//
// 参数：
// 1. const CComplex& cpxX - 与指定复数相加的复数
//
// 返回值：CComplex型，指定复数与cpxX相加之和
//////////////////////////////////////////////////////////////////////
CComplex CComplex::operator+(const CComplex& cpxX) const
{
	double x = m_dblX + cpxX.m_dblX;
	double y = m_dblY + cpxX.m_dblY;

	return CComplex(x, y);
}

//////////////////////////////////////////////////////////////////////
// 重载运算符-，实现复数的减法
//
// 参数：
// 1. const CComplex& cpxX - 与指定复数相减的复数
//
// 返回值：CComplex型，指定复数减去cpxX之差
//////////////////////////////////////////////////////////////////////
CComplex CComplex::operator-(const CComplex& cpxX) const
{
	double x = m_dblX - cpxX.m_dblX;
	double y = m_dblY - cpxX.m_dblY;

	return CComplex(x, y);
}

//////////////////////////////////////////////////////////////////////
// 重载运算符*，实现复数的乘法
//
// 参数：
// 1. const CComplex& cpxX - 与指定复数相乘的复数
//
// 返回值：CComplex型，指定复数与cpxX相乘之积
//////////////////////////////////////////////////////////////////////
CComplex CComplex::operator*(const CComplex& cpxX) const
{
	double x = m_dblX * cpxX.m_dblX - m_dblY * cpxX.m_dblY;
	double y = m_dblX * cpxX.m_dblY + m_dblY * cpxX.m_dblX;

	return CComplex(x, y);
}

//////////////////////////////////////////////////////////////////////
// 重载运算符*，实现复数的乘法
//
// 参数：
// 1. const double& cpxX - 与指定复数相乘的复数
//
// 返回值：CComplex型，指定复数与cpxX相乘之积
//////////////////////////////////////////////////////////////////////
CComplex CComplex::operator*(const double& cpxX) const
{
	double x = m_dblX * cpxX;
	double y = m_dblY * cpxX;

	return CComplex(x, y);
}

//////////////////////////////////////////////////////////////////////
// 重载运算符/，实现复数的除法
//
// 参数：
// 1. const CComplex& cpxX - 与指定复数相除的复数
//
// 返回值：CComplex型，指定复数除与cpxX之商
//////////////////////////////////////////////////////////////////////
CComplex CComplex::operator/(const CComplex& cpxX) const
{
	double e, f, x, y;

	if (fabs(cpxX.m_dblX) >= fabs(cpxX.m_dblY))
	{
		e = cpxX.m_dblY / cpxX.m_dblX;
		f = cpxX.m_dblX + e * cpxX.m_dblY;

		x = (m_dblX + m_dblY * e) / f;
		y = (m_dblY - m_dblX * e) / f;
	}
	else
	{
		e = cpxX.m_dblX / cpxX.m_dblY;
		f = cpxX.m_dblY + e * cpxX.m_dblX;

		x = (m_dblX * e + m_dblY) / f;
		y = (m_dblY * e - m_dblX) / f;
	}

	return CComplex(x, y);
}

//////////////////////////////////////////////////////////////////////
// 计算复数的模
//
// 参数：无
//
// 返回值：double型，指定复数的模
//////////////////////////////////////////////////////////////////////
double CComplex::Abs() const
{
	// 求取实部和虚部的绝对值
	double x = fabs(m_dblX);
	double y = fabs(m_dblY);

	if (m_dblX == 0)
		return y;
	if (m_dblY == 0)
		return x;


	// 计算模
	if (x > y)
		return (x * sqrt(1 + (y / x) * (y / x)));

	return (y * sqrt(1 + (x / y) * (x / y)));
}

//////////////////////////////////////////////////////////////////////
// 计算复数的根
//
// 参数：
// 1. int n - 待求根的根次
// 2. CComplex cpxR[] - CComplex型数组，长度为n，返回复数的所有根
//
// 返回值：无
//////////////////////////////////////////////////////////////////////
void CComplex::Root(int n, CComplex cpxR[]) const
{
	if (n<1)
		return;

	double q = atan2(m_dblY, m_dblX);
	double r = sqrt(m_dblX*m_dblX + m_dblY*m_dblY);
	if (r != 0)
	{
		r = (1.0 / n)*log(r);
		r = exp(r);
	}

	for (int k = 0; k <= n - 1; k++)
	{
		double t = (2.0*k*3.1415926 + q) / n;
		cpxR[k].m_dblX = r*cos(t);
		cpxR[k].m_dblY = r*sin(t);
	}
}

//////////////////////////////////////////////////////////////////////
// 计算复数的实幂指数
//
// 参数：
// 1. double dblW - 待求实幂指数的幂次
//
// 返回值：CComplex型，复数的实幂指数值
//////////////////////////////////////////////////////////////////////
CComplex CComplex::Pow(double dblW) const
{
	// 常量
	const double PI = 3.14159265358979;

	// 局部变量
	double r, t;

	// 特殊值处理
	if ((m_dblX == 0) && (m_dblY == 0))
		return CComplex(0, 0);

	// 幂运算公式中的三角函数运算
	if (m_dblX == 0)
	{
		if (m_dblY > 0)
			t = 1.5707963268;
		else
			t = -1.5707963268;
	}
	else
	{
		if (m_dblX > 0)
			t = atan2(m_dblY, m_dblX);
		else
		{
			if (m_dblY >= 0)
				t = atan2(m_dblY, m_dblX) + PI;
			else
				t = atan2(m_dblY, m_dblX) - PI;
		}
	}

	// 模的幂
	r = exp(dblW * log(sqrt(m_dblX * m_dblX + m_dblY * m_dblY)));

	// 复数的实幂指数
	return CComplex(r * cos(dblW * t), r * sin(dblW * t));
}

//////////////////////////////////////////////////////////////////////
// 计算复数的复幂指数
//
// 参数：
// 1. CComplex cpxW - 待求复幂指数的幂次
// 2. int n - 控制参数，默认值为0。当n=0时，求得的结果为复幂指数的主值。
//
// 返回值：CComplex型，复数的复幂指数值
//////////////////////////////////////////////////////////////////////
CComplex CComplex::Pow(CComplex cpxW, int n /*= 0*/) const
{
	// 常量
	const double PI = 3.14159265358979;
	// 局部变量
	double r, s, u, v;

	// 特殊值处理
	if (m_dblX == 0)
	{
		if (m_dblY == 0)
			return CComplex(0, 0);

		s = 1.5707963268 * (fabs(m_dblY) / m_dblY + 4 * n);
	}
	else
	{
		s = 2 * PI * n + atan2(m_dblY, m_dblX);

		if (m_dblX < 0)
		{
			if (m_dblY > 0)
				s = s + PI;
			else
				s = s - PI;
		}
	}

	// 求幂运算公式
	r = 0.5 * log(m_dblX * m_dblX + m_dblY * m_dblY);
	v = cpxW.m_dblX * r + cpxW.m_dblY * s;
	u = exp(cpxW.m_dblX * r - cpxW.m_dblY * s);

	return CComplex(u * cos(v), u * sin(v));
}

//////////////////////////////////////////////////////////////////////
// 计算复数的自然对数
//
// 参数：无
//
// 返回值：CComplex型，复数的自然对数值
//////////////////////////////////////////////////////////////////////
CComplex CComplex::Log() const
{
	double p = log(sqrt(m_dblX*m_dblX + m_dblY*m_dblY));
	return CComplex(p, atan2(m_dblY, m_dblX));
}

//////////////////////////////////////////////////////////////////////
// 计算复数的正弦
//
// 参数：无
//
// 返回值：CComplex型，复数的正弦值
//////////////////////////////////////////////////////////////////////
CComplex CComplex::Sin() const
{
	int i;
	double x, y, y1, br, b1, b2, c[6];

	// 切比雪夫公式的常数系数
	c[0] = 1.13031820798497;
	c[1] = 0.04433684984866;
	c[2] = 0.00054292631191;
	c[3] = 0.00000319843646;
	c[4] = 0.00000001103607;
	c[5] = 0.00000000002498;

	y1 = exp(m_dblY);
	x = 0.5 * (y1 + 1 / y1);
	if (fabs(m_dblY) >= 1)
		y = 0.5 * (y1 - 1 / y1);
	else
	{
		b1 = 0;
		b2 = 0;
		y1 = 2 * (2 * m_dblY * m_dblY - 1);
		for (i = 5; i >= 0; --i)
		{
			br = y1 * b1 - b2 - c[i];
			if (i != 0)
			{
				b2 = b1;
				b1 = br;
			}
		}

		y = m_dblY * (br - b1);
	}

	// 组合计算结果
	x = x * sin(m_dblX);
	y = y * cos(m_dblX);

	return CComplex(x, y);
}

//////////////////////////////////////////////////////////////////////
// 计算复数的余弦
//
// 参数：无
//
// 返回值：CComplex型，复数的余弦值
//////////////////////////////////////////////////////////////////////
CComplex CComplex::Cos() const
{
	int i;
	double x, y, y1, br, b1, b2, c[6];

	// 切比雪夫公式的常数系数
	c[0] = 1.13031820798497;
	c[1] = 0.04433684984866;
	c[2] = 0.00054292631191;
	c[3] = 0.00000319843646;
	c[4] = 0.00000001103607;
	c[5] = 0.00000000002498;

	y1 = exp(m_dblY);
	x = 0.5 * (y1 + 1 / y1);
	if (fabs(m_dblY) >= 1)
		y = 0.5 * (y1 - 1 / y1);
	else
	{
		b1 = 0;
		b2 = 0;
		y1 = 2 * (2 * m_dblY * m_dblY - 1);
		for (i = 5; i >= 0; --i)
		{
			br = y1 * b1 - b2 - c[i];
			if (i != 0)
			{
				b2 = b1;
				b1 = br;
			}
		}

		y = m_dblY * (br - b1);
	}

	// 组合计算结果
	x = x * cos(m_dblX);
	y = -y * sin(m_dblX);

	return CComplex(x, y);
}

//////////////////////////////////////////////////////////////////////
// 计算复数的正切
//
// 参数：无
//
// 返回值：CComplex型，复数的正切值
//////////////////////////////////////////////////////////////////////
CComplex CComplex::Tan() const
{
	return Sin() / Cos();
}
