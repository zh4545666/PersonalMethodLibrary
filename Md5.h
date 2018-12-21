#ifndef CMD5_H  
#define CMD5_H  

/* Type define */
typedef unsigned char byte;
typedef unsigned int uint32;

#ifdef AFX_CLASS  
#define AFX_EX_CLASS _declspec(dllexport)  
#else  
#define AFX_EX_CLASS _declspec(dllimport)  
#endif

/* MD5 declaration. */
class AFX_EX_CLASS CMD5 {
public:
	//
	// 4�����캯��
	//

	CMD5();                                           //Ĭ�Ϲ��캯��
	CMD5(const void *input, size_t length);           //�����ڴ��ַ�볤����Ϣ�Ĺ��캯��
	CMD5(const string &str);                          //�����ַ����Ĺ��캯��
	CMD5(ifstream &in);                               //�������Ĺ��캯��

	//
	// 3��Update����
	//

	void update(const void *input, size_t length);    //��CMD5����������ڴ��
	void update(const string &str);                   //����ַ���
	void update(ifstream &in);                        //�����

	//
	// �������
	//

	const byte* digest();                             //����MD5��,������ָ������ָ��
	string toString();                                //����MD5��,���������Ӧ���ַ���
	void reset();                                     //����
private:
	void update(const byte *input, size_t length);
	void final();
	void transform(const byte block[64]);
	void encode(const uint32 *input, byte *output, size_t length);
	void decode(const byte *input, uint32 *output, size_t length);
	string bytesToHexString(const byte *input, size_t length);

	/* class uncopyable */
	CMD5(const CMD5&);
	CMD5& operator=(const CMD5&);
private:
	uint32 _state[4];   /* state (ABCD) */
	uint32 _count[2];   /* number of bits, modulo 2^64 (low-order word first) */
	byte _buffer[64];   /* input buffer */
	byte _digest[16];   /* message digest */
	bool _finished;     /* calculate finished ? */

	static const byte PADDING[64];  /* padding for calculate */
	static const char HEX[16];
	static const size_t BUFFER_SIZE = 1024;
};

#endif/*CMD5_H*/ 
