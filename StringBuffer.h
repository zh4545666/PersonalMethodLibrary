#pragma once

#include <bitset>
#include <string>

namespace PersonalMethod {

	typedef unsigned char BYTE;
	typedef unsigned __int64 size_t;

	class CTokenizer
	{
	public:
		CTokenizer(const std::string& cs, const std::string& csDelim) : m_cs(cs), m_nCurPos(0)
		{
			SetDelimiters(csDelim);
		}
		void SetDelimiters(const std::string& csDelim)
		{
			for (size_t i = 0; i < csDelim.size(); ++i)
				m_sDelimiter.set(static_cast<BYTE>(csDelim[i]));
		}

		bool Next(std::string& cs)
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

		std::string	Tail() const
		{
			size_t nCurPos = m_nCurPos;

			while (nCurPos < m_cs.size() && m_sDelimiter[static_cast<BYTE>(m_cs[nCurPos])])
				++nCurPos;

			std::string csResult;

			if (nCurPos < m_cs.size())
				csResult = m_cs.substr(nCurPos);

			return csResult;
		}

	public:
		static void TrimRightSpace(std::string& pString)
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
		static void  TrimLeftSpace(std::string& pString)
		{
			if (pString.empty())
			{
				return;
			}

			unsigned long i = 0;
			while (i < pString.size())
			{
				if (!isspace(pString.at(i)))
				{
					if (i > 1)
						i -= 1;
					pString = pString.substr(i, pString.size());
					break;
				}
				++i;
			}
			return;
		}

	private:
		std::string m_cs;
		std::bitset<256> m_sDelimiter;
		size_t m_nCurPos;
	};

}