#include "stdafx.h"
#include "DataAnalysis.h"


CDataAnalysis::CDataAnalysis()
{
	make_crc_table();
}


CDataAnalysis::~CDataAnalysis()
{
}

void CDataAnalysis::GetSTD(double* pData, size_t N, double pResult[3])
{
	size_t i;
	pResult[0] = 0.0;
	for (i = 0; i < N; i++) pResult[0] += pData[i] / N;
	pResult[1] = 0.0;
	for (i = 0; i < N; i++)
		pResult[1] += (pData[i] - pResult[0])*(pData[i] - pResult[0]);
	pResult[1] = pResult[1] / N;
	pResult[2] = sqrt(pResult[1]);
}

/* Make the table for a fast CRC. */
void CDataAnalysis::make_crc_table(void)
{
	unsigned long c;
	int n, k;

	for (n = 0; n < 256; n++) {
		c = (unsigned long)n;
		for (k = 0; k < 8; k++) {
			if (c & 1)
				c = 0xedb88320L ^ (c >> 1);
			else
				c = c >> 1;
		}
		crc_table[n] = c;
	}
	crc_table_computed = 1;
}

/* Update a running CRC with the bytes buf[0..len-1]--the CRC
should be initialized to all 1's, and the transmitted value
is the 1's complement of the final running CRC (see the
crc() routine below). */

unsigned long CDataAnalysis::update_crc(unsigned long crc, unsigned char *buf,
	int len)
{
	unsigned long c = crc;
	int n;

	if (!crc_table_computed)
		make_crc_table();
	for (n = 0; n < len; n++) {
		c = crc_table[(c ^ buf[n]) & 0xff] ^ (c >> 8);
	}
	return c;
}

/* Return the CRC of the bytes buf[0..len-1]. */
unsigned long CDataAnalysis::crc(unsigned char *buf, int len)
{
	return update_crc(0xffffffffL, buf, len) ^ 0xffffffffL;
}

bool CDataAnalysis::huffman_node_compare(const huffman_node* other, const huffman_node* other2)
{
	if (other == NULL || other2 == NULL)
		return false;
	return other->weight < other2->weight;
}

CDataAnalysis::huffman_node* CDataAnalysis::huffman_tree_create(map<unsigned char, int> &word, map<unsigned char, string> &res)
{
	huffman_node* root = NULL;
	vector<huffman_node*> huffman_tree_node;

	map<unsigned char, int>::iterator it;
	for (it = word.begin(); it != word.end(); it++)
	{
		huffman_node* pnode = new huffman_node();
		pnode->c = it->first;
		pnode->weight = it->second;
		pnode->left = NULL;
		pnode->right = NULL;
		pnode->bValue = true;
		huffman_tree_node.push_back(pnode);
	}

	while (huffman_tree_node.size() > 0)
	{
		sort(huffman_tree_node.begin(), huffman_tree_node.end(), huffman_node_compare);
		if (huffman_tree_node.size() == 1)
		{
			root = huffman_tree_node[0];
			huffman_tree_node.erase(huffman_tree_node.begin());
		}
		else
		{
			huffman_node* node_1 = huffman_tree_node[0];
			huffman_node* node_2 = huffman_tree_node[1];
			huffman_tree_node.erase(huffman_tree_node.begin());
			huffman_tree_node.erase(huffman_tree_node.begin());
			huffman_node* node = new huffman_node();
			node->weight = node_1->weight + node_2->weight;
			(node_1->weight < node_2->weight) ? (node->left = node_1, node->right = node_2) : (node->left = node_2, node->right = node_1);
			huffman_tree_node.push_back(node);
		}
	}
	get_huffman_code(root);

	res.clear();
	print_huffman_pre(root, res);

	return root;
}

void CDataAnalysis::huffman_tree_del(huffman_node* pnode)
{
	if (pnode == NULL)
		return;

	huffman_tree_del(pnode->left);
	huffman_tree_del(pnode->right);
	delete pnode;
	pnode = NULL;
}

int CDataAnalysis::get_huffman_code(huffman_node* pnode)
{
	if (pnode == NULL)
		return 1;
	huffman_node* p = pnode;
	queue<huffman_node*> q;
	q.push(p);
	while (q.size() > 0)
	{
		p = q.front();
		q.pop();
		if (p->left != NULL)
		{
			q.push(p->left);
			memcpy((p->left)->hufman_code, p->hufman_code, sizeof(char)*LEN);
			char* ptr = p->left->hufman_code;
			while (*ptr != '\0')
				ptr++;
			*ptr = '0';
		}
		if (p->right != NULL)
		{
			q.push(p->right);
			memcpy((p->right)->hufman_code, p->hufman_code, sizeof(char)*LEN);
			char* ptr = p->right->hufman_code;
			while (*ptr != '\0')
				ptr++;
			*ptr = '1';
		}
	}

	return 0;

}

void CDataAnalysis::print_huffman_pre(huffman_node* pnode, map<unsigned char, string>& res)
{
	if (pnode != NULL)
	{
		if (pnode->bValue)
			res[pnode->c] = pnode->hufman_code;
		print_huffman_pre(pnode->left, res);
		print_huffman_pre(pnode->right, res);
	}
}

void CDataAnalysis::printHex(std::ifstream& ifs, std::ostream& ostream){
	using namespace std;
	ostream << setfill('0') << hex << uppercase;

	unsigned char byte;
	unsigned long count = 0;
	while (true){
		ostream << setw(8) << count << "    ";
		for (int i = 0; i<8; ++i){
			if (ifs.read((char*)&byte, 1))
				ostream << setw(2) << (int)byte << " ";
			else
			{
				ostream << setfill(' ') << dec;
				return;
			}
		}
		ostream << " ";
		for (int i = 0; i<8; ++i){
			if (ifs.read((char*)&byte, 1))
				ostream << setw(2) << (int)byte << " ";
			else
			{
				ostream << setfill(' ') << dec;
				return;
			}
		}
		ostream << endl;
		count += 16;
	}
}

int CDataAnalysis::UnsignedCharToInt(unsigned char c1, unsigned char c2)
{
	int res;
	res = int(c1) * 256 + int(c2);
	return res;
}

void CDataAnalysis::IntToUnsignedChar(int val, unsigned char res[2])
{
	memset(res, 0, sizeof(unsigned char)* 2);
	res[0] = val / 0xff;
	val %= 0xff;
	res[1] = val;
}

string CDataAnalysis::UnsignedCharToString(unsigned char c)
{
	//输出8位二进制值
	int value = (int)c;
	string s;
	while (value > 0)
	{
		if (value % 2 == 0)
			s = '0' + s;
		else
			s = '1' + s;
		value /= 2;
	}

	while (s.length() < 8)
		s = '0' + s;

	return s;
}

int CDataAnalysis::KmpMatch(char* s, char* p)
{
	int sLen = strlen(s);
	int pLen = strlen(p);

	int* next = new int[pLen];
	GetNextval(p, next);

	int i = 0, j = 0;
	while (i < sLen && j < pLen)
	{
		if (s[i] == p[j])
		{
			i++; j++;
		}
		else
		{ 
			j = next[j];
		}
	}
	if (j == pLen)
		return i - j;
	else
		return -1;
}

void CDataAnalysis::GetNextval(char* p, int next[])
{
	int pLen = strlen(p);
	next[0] = -1;
	int k = -1;
	int j = 0;
	while (j < pLen - 1)
	{
		//p[k]表示前缀，p[j]表示后缀    
		if (k == -1 || p[j] == p[k])
		{
			++j;
			++k;
			if (p[j] != p[k])
				next[j] = k;
			else
				//因为不能出现p[j] = p[ next[j ]]，所以当出现时需要继续递归，k = next[k] = next[next[k]]  
				next[j] = next[k];
		}
		else
		{
			k = next[k];
		}
	}
}
