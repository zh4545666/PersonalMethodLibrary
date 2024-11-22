#pragma once  
  
#ifdef AFX_CLASS  
#define AFX_EX_CLASS _declspec(dllexport)  
#else  
#define AFX_EX_CLASS _declspec(dllimport)  
#endif 

namespace PersonalMethod {

#define AES_KEY_ROW_NUMBER 4  
#define AES_KEY_COLUMN_NUMBER 4  
#define AES_ROUND_COUNT 10  

    typedef unsigned char BYTE;
    typedef unsigned __int64 size_t;

    class AFX_EX_CLASS CAES
    {
    public:
        CAES() = delete;

        //初始化AES密钥，默认为4*4
        CAES(unsigned char* key, size_t N = 16);
        virtual ~CAES(void);

        //原始数据，加密/解密后数据，长度
        void Encrypt(BYTE*, BYTE*, size_t);
        void Decrypt(BYTE*, BYTE*, size_t);
        //原始数据（加密/解密后数据），长度
        void Encrypt(BYTE*, size_t);
        void Decrypt(BYTE*, size_t);
        //输入文件名，输出文件名
        bool Encrypt(char*, char*);
        bool Decrypt(char*, char*);

    private:

        BYTE swapbox[11][4][4];

        BYTE* Cipher(BYTE* input);
        BYTE* InvCipher(BYTE* input);

        BYTE* Cipher(void* input, size_t length);
        BYTE* InvCipher(void* input, size_t length);

        void KeyExpansion(BYTE* key, BYTE w[][4][AES_KEY_COLUMN_NUMBER]);
        BYTE FFmul(BYTE a, BYTE b);

        void SubBytes(BYTE state[][AES_KEY_COLUMN_NUMBER]);
        void ShiftRows(BYTE state[][AES_KEY_COLUMN_NUMBER]);
        void MixColumns(BYTE state[][AES_KEY_COLUMN_NUMBER]);
        void AddRoundKey(BYTE state[][AES_KEY_COLUMN_NUMBER], BYTE k[][AES_KEY_COLUMN_NUMBER]);

        void InvSubBytes(BYTE state[][AES_KEY_COLUMN_NUMBER]);
        void InvShiftRows(BYTE state[][AES_KEY_COLUMN_NUMBER]);
        void InvMixColumns(BYTE state[][AES_KEY_COLUMN_NUMBER]);

    };

}