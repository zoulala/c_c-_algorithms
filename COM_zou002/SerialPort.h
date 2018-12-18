//////////////////////////////////////////////////////////////////////
// COPYRIGHT NOTICE  
// Copyright (c) 2016, �������Ϣ�������޹�˾Algorithm Group  ����Ȩ������  
// All rights reserved.  
//
// @file    SerialPort.h    
// @brief   ����ͨ��ͷ�ļ�  
//
// ����������
//
//
// @version 1.0     
// @author  LingWei zou    
// @E-mail��zlw2008ok@126.com  
// @date    2016/11/17  
//  
//
// �޶�˵���� 
////////////////////////////////////////////////////////////////////

#ifndef SerialPort_h
#define SerialPort_h


#include <Windows.h>  

class SerialPort
{
public:
    SerialPort(void);  
    ~SerialPort(void); 

public:
	bool OpenPort(UINT PortNo );
	bool InitPort(UINT BaudRate);
	bool OpenListenThread();	//���������߳�
    static UINT WINAPI ListenThread(void* pParam);//�߳�ִ�в���
	bool CloseListenThread();	//�رռ����߳�
	bool ReadPort(char *str);
	bool WritePort(unsigned char* WData, unsigned int length);
	bool ClearPort();
	bool ClosePort();

private:
	HANDLE hcom;

	/** �߳��˳���־���� */   
   static bool s_bExit;  
 
    /** �߳̾�� */   
    volatile HANDLE    m_hListenThread;  

};



#endif