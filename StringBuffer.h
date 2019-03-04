#pragma once

#include <bitset>

class CTokenizer
{
public:
	CTokenizer(const string& cs, const string& csDelim) : m_cs(cs), m_nCurPos(0)
	{
		SetDelimiters(csDelim);
	}
	void SetDelimiters(const string& csDelim)
	{
		for (size_t i = 0; i < csDelim.size(); ++i)
			m_sDelimiter.set(static_cast<BYTE>(csDelim[i]));
	}

	bool Next(string& cs)
	{
		cs.clear();

		while (m_nCurPos < m_cs.size() && m_sDelimiter[static_cast<BYTE>(m_cs[m_nCurPos])])
			++m_nCurPos;

		if (m_nCurPos >= m_cs.size())
			return false;

		int nStartPos = m_nCurPos;
		while (m_nCurPos < m_cs.size() && !m_sDelimiter[static_cast<BYTE>(m_cs[m_nCurPos])])
			++m_nCurPos;

		cs = m_cs.substr(nStartPos, m_nCurPos - nStartPos);

		return true;
	}

	string	Tail() const
	{
		size_t nCurPos = m_nCurPos;

		while (nCurPos < m_cs.size() && m_sDelimiter[static_cast<BYTE>(m_cs[nCurPos])])
			++nCurPos;

		string csResult;

		if (nCurPos < m_cs.size())
			csResult = m_cs.substr(nCurPos);

		return csResult;
	}

public:
	static void TrimRightSpace(string &pString)
	{
		int     strLength = 0;

		if (pString.empty())
		{
			return;
		}

		strLength = pString.size();

		for (int i = strLength - 1; i >= 0; --i)
		{
			if (!isspace(pString[i]))
			{
				pString = pString.substr(0, i + 1);
				break;
			}
		}
		return;
	}
	static void  TrimLeftSpace(string &pString)
	{
		if (pString.empty())
		{
			return;
		}

		unsigned long i = 0;
		while (i<pString.size())
		{
			if (!isspace(pString.at(i)))
			{
				if (i>1)
					i -= 1;
				pString = pString.substr(i, pString.size());
				break;
			}
			++i;
		}
		return;
	}

private:
	string m_cs;
	std::bitset<256> m_sDelimiter;
	size_t m_nCurPos;
};