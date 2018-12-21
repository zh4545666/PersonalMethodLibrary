#include "stdafx.h"
#include "Complex.h"


//////////////////////////////////////////////////////////////////////
// �������캯��
//////////////////////////////////////////////////////////////////////
CComplex::CComplex()
{
	m_dblX = 0.0;
	m_dblY = 0.0;
}

//////////////////////////////////////////////////////////////////////
// ָ��ֵ���캯��
//
// ������
// 1. double dblX - ָ����ʵ��
// 2. double dblY - ָ�����鲿
//////////////////////////////////////////////////////////////////////
CComplex::CComplex(double dblX, double dblY)
{
	m_dblX = dblX;
	m_dblY = dblY;
}

//////////////////////////////////////////////////////////////////////
// �������캯��
//
// ������
// 1. const CComplex& other - Դ����
//////////////////////////////////////////////////////////////////////
CComplex::CComplex(const CComplex& other)
{
	m_dblX = other.m_dblX;
	m_dblY = other.m_dblY;
}

//////////////////////////////////////////////////////////////////////
// ָ��������ʵ�����鲿
//
// ������
// 1. double dblX - ������ʵ��
// 2. double dblY - �������鲿
//////////////////////////////////////////////////////////////////////
void CComplex::SetData(double dblX, double dblY)
{
	m_dblX = dblX;
	m_dblY = dblY;
}

//////////////////////////////////////////////////////////////////////
// ָ��������ʵ��
//
// ������
// 1. double dblX - ������ʵ��
//////////////////////////////////////////////////////////////////////
void CComplex::SetReal(double dblX)
{
	m_dblX = dblX;
}

//////////////////////////////////////////////////////////////////////
// ָ���������鲿
//
// ������
// 1. double dblX - �������鲿
//////////////////////////////////////////////////////////////////////
void CComplex::SetImag(double dblY)
{
	m_dblY = dblY;
}

//////////////////////////////////////////////////////////////////////
// ȡ������ʵ��
//
// ������  ��
//
// ����ֵ��double �ͣ�������ʵ��
//////////////////////////////////////////////////////////////////////
double CComplex::GetReal()
{
	return m_dblX;
}

//////////////////////////////////////////////////////////////////////
// ȡ�������鲿
//
// ������  ��
//
// ����ֵ��double �ͣ��������鲿
//////////////////////////////////////////////////////////////////////
double CComplex::GetImag()
{
	return m_dblY;
}

//////////////////////////////////////////////////////////////////////
// ������ת��Ϊ"a+bj"��ʽ���ַ���
//
// ������  ��
//
// ����ֵ��string ����"a+bj"��ʽ���ַ���
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
// ��"a,b"��ʽ���ַ���ת��Ϊ��������aΪ������ʵ����bΪ�������鲿
//
// ������
// 1. string s - "a,b"��ʽ���ַ�����aΪ������ʵ����bΪ�������鲿
// 2. const string& sDelim - a, b֮��ķָ�����Ĭ��Ϊ�ո�
//
// ����ֵ����
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
// ���������==���Ƚ����������Ƿ����
//
// ������
// 1. const CComplex& cpxX - ���ڱȽϵĸ���
//
// ����ֵ��bool�ͣ������Ϊtrue������Ϊfalse
//////////////////////////////////////////////////////////////////////
bool CComplex::operator==(const CComplex& cpxX) const
{
	return (m_dblX == cpxX.m_dblX && m_dblY == cpxX.m_dblY);
}

//////////////////////////////////////////////////////////////////////
// ���������!=���Ƚ����������Ƿ񲻵�
//
// ������
// 1. const CComplex& cpxX - ���ڱȽϵĸ���
//
// ����ֵ��bool�ͣ��������Ϊtrue�����Ϊfalse
//////////////////////////////////////////////////////////////////////
bool CComplex::operator!=(const CComplex& cpxX) const
{
	return (m_dblX != cpxX.m_dblX || m_dblY != cpxX.m_dblY);
}

//////////////////////////////////////////////////////////////////////
// ���������=����������ֵ
//
// ������
// 1. const CComplex& cpxX - ���ڸ�������ֵ��Դ����
//
// ����ֵ��CComplex�͵����ã������õĸ�����cpxX���
//////////////////////////////////////////////////////////////////////
CComplex& CComplex::operator=(const CComplex& cpxX)
{
	m_dblX = cpxX.m_dblX;
	m_dblY = cpxX.m_dblY;

	return *this;
}

//////////////////////////////////////////////////////////////////////
// ���������+��ʵ�ָ����ļӷ�
//
// ������
// 1. const CComplex& cpxX - ��ָ��������ӵĸ���
//
// ����ֵ��CComplex�ͣ�ָ��������cpxX���֮��
//////////////////////////////////////////////////////////////////////
CComplex CComplex::operator+(const CComplex& cpxX) const
{
	double x = m_dblX + cpxX.m_dblX;
	double y = m_dblY + cpxX.m_dblY;

	return CComplex(x, y);
}

//////////////////////////////////////////////////////////////////////
// ���������-��ʵ�ָ����ļ���
//
// ������
// 1. const CComplex& cpxX - ��ָ����������ĸ���
//
// ����ֵ��CComplex�ͣ�ָ��������ȥcpxX֮��
//////////////////////////////////////////////////////////////////////
CComplex CComplex::operator-(const CComplex& cpxX) const
{
	double x = m_dblX - cpxX.m_dblX;
	double y = m_dblY - cpxX.m_dblY;

	return CComplex(x, y);
}

//////////////////////////////////////////////////////////////////////
// ���������*��ʵ�ָ����ĳ˷�
//
// ������
// 1. const CComplex& cpxX - ��ָ��������˵ĸ���
//
// ����ֵ��CComplex�ͣ�ָ��������cpxX���֮��
//////////////////////////////////////////////////////////////////////
CComplex CComplex::operator*(const CComplex& cpxX) const
{
	double x = m_dblX * cpxX.m_dblX - m_dblY * cpxX.m_dblY;
	double y = m_dblX * cpxX.m_dblY + m_dblY * cpxX.m_dblX;

	return CComplex(x, y);
}

//////////////////////////////////////////////////////////////////////
// ���������*��ʵ�ָ����ĳ˷�
//
// ������
// 1. const double& cpxX - ��ָ��������˵ĸ���
//
// ����ֵ��CComplex�ͣ�ָ��������cpxX���֮��
//////////////////////////////////////////////////////////////////////
CComplex CComplex::operator*(const double& cpxX) const
{
	double x = m_dblX * cpxX;
	double y = m_dblY * cpxX;

	return CComplex(x, y);
}

//////////////////////////////////////////////////////////////////////
// ���������/��ʵ�ָ����ĳ���
//
// ������
// 1. const CComplex& cpxX - ��ָ����������ĸ���
//
// ����ֵ��CComplex�ͣ�ָ����������cpxX֮��
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
// ���㸴����ģ
//
// ��������
//
// ����ֵ��double�ͣ�ָ��������ģ
//////////////////////////////////////////////////////////////////////
double CComplex::Abs() const
{
	// ��ȡʵ�����鲿�ľ���ֵ
	double x = fabs(m_dblX);
	double y = fabs(m_dblY);

	if (m_dblX == 0)
		return y;
	if (m_dblY == 0)
		return x;


	// ����ģ
	if (x > y)
		return (x * sqrt(1 + (y / x) * (y / x)));

	return (y * sqrt(1 + (x / y) * (x / y)));
}

//////////////////////////////////////////////////////////////////////
// ���㸴���ĸ�
//
// ������
// 1. int n - ������ĸ���
// 2. CComplex cpxR[] - CComplex�����飬����Ϊn�����ظ��������и�
//
// ����ֵ����
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
// ���㸴����ʵ��ָ��
//
// ������
// 1. double dblW - ����ʵ��ָ�����ݴ�
//
// ����ֵ��CComplex�ͣ�������ʵ��ָ��ֵ
//////////////////////////////////////////////////////////////////////
CComplex CComplex::Pow(double dblW) const
{
	// ����
	const double PI = 3.14159265358979;

	// �ֲ�����
	double r, t;

	// ����ֵ����
	if ((m_dblX == 0) && (m_dblY == 0))
		return CComplex(0, 0);

	// �����㹫ʽ�е����Ǻ�������
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

	// ģ����
	r = exp(dblW * log(sqrt(m_dblX * m_dblX + m_dblY * m_dblY)));

	// ������ʵ��ָ��
	return CComplex(r * cos(dblW * t), r * sin(dblW * t));
}

//////////////////////////////////////////////////////////////////////
// ���㸴���ĸ���ָ��
//
// ������
// 1. CComplex cpxW - ������ָ�����ݴ�
// 2. int n - ���Ʋ�����Ĭ��ֵΪ0����n=0ʱ����õĽ��Ϊ����ָ������ֵ��
//
// ����ֵ��CComplex�ͣ������ĸ���ָ��ֵ
//////////////////////////////////////////////////////////////////////
CComplex CComplex::Pow(CComplex cpxW, int n /*= 0*/) const
{
	// ����
	const double PI = 3.14159265358979;
	// �ֲ�����
	double r, s, u, v;

	// ����ֵ����
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

	// �������㹫ʽ
	r = 0.5 * log(m_dblX * m_dblX + m_dblY * m_dblY);
	v = cpxW.m_dblX * r + cpxW.m_dblY * s;
	u = exp(cpxW.m_dblX * r - cpxW.m_dblY * s);

	return CComplex(u * cos(v), u * sin(v));
}

//////////////////////////////////////////////////////////////////////
// ���㸴������Ȼ����
//
// ��������
//
// ����ֵ��CComplex�ͣ���������Ȼ����ֵ
//////////////////////////////////////////////////////////////////////
CComplex CComplex::Log() const
{
	double p = log(sqrt(m_dblX*m_dblX + m_dblY*m_dblY));
	return CComplex(p, atan2(m_dblY, m_dblX));
}

//////////////////////////////////////////////////////////////////////
// ���㸴��������
//
// ��������
//
// ����ֵ��CComplex�ͣ�����������ֵ
//////////////////////////////////////////////////////////////////////
CComplex CComplex::Sin() const
{
	int i;
	double x, y, y1, br, b1, b2, c[6];

	// �б�ѩ��ʽ�ĳ���ϵ��
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

	// ��ϼ�����
	x = x * sin(m_dblX);
	y = y * cos(m_dblX);

	return CComplex(x, y);
}

//////////////////////////////////////////////////////////////////////
// ���㸴��������
//
// ��������
//
// ����ֵ��CComplex�ͣ�����������ֵ
//////////////////////////////////////////////////////////////////////
CComplex CComplex::Cos() const
{
	int i;
	double x, y, y1, br, b1, b2, c[6];

	// �б�ѩ��ʽ�ĳ���ϵ��
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

	// ��ϼ�����
	x = x * cos(m_dblX);
	y = -y * sin(m_dblX);

	return CComplex(x, y);
}

//////////////////////////////////////////////////////////////////////
// ���㸴��������
//
// ��������
//
// ����ֵ��CComplex�ͣ�����������ֵ
//////////////////////////////////////////////////////////////////////
CComplex CComplex::Tan() const
{
	return Sin() / Cos();
}
