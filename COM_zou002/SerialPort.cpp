//////////////////////////////////////////////////////////////////////
// COPYRIGHT NOTICE  
// Copyright (c) 2016, �������Ϣ�������޹�˾Algorithm Group  ����Ȩ������  
// All rights reserved.  
//
// @file    SerialPort.cpp    
// @brief   ����ͨ�ź��������ļ�  
//
// �򿪡����á���д�͹رյȲ���
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


#include "StdAfx.h"  
#include "SerialPort.h"  
#include <process.h>  
#include <iostream>  

#define buff_n 1   //��ȡ����buff��С
/** �߳��˳���־ */   
bool SerialPort::s_bExit = false;  

SerialPort::SerialPort(void)
{
	//m_hComm = INVALID_HANDLE_VALUE;  
    m_hListenThread = INVALID_HANDLE_VALUE;  
}
SerialPort::~SerialPort(void)
{
	CloseListenThread();  
    ClosePort(); 
}


bool SerialPort::OpenPort(UINT PortNo )
{
    /** �Ѵ��ڵı��ת��Ϊ�豸�� */   
    char szPort[50];  
    sprintf_s(szPort, "COM%d", PortNo);  
 
   //  Open a handle to the specified com port.
   hcom = CreateFileA( szPort,
                      GENERIC_READ | GENERIC_WRITE,
                      0,      //  must be opened with exclusive-access
                      NULL,   //  default security attributes
                      OPEN_EXISTING, //  must use OPEN_EXISTING
                      0,//FILE_ATTRIBUTE_NORMAL|FILE_FLAG_OVERLAPPED      //  not overlapped I/O
                      NULL ); //  hTemplate must be NULL for comm devices

   if (hcom == INVALID_HANDLE_VALUE) 
   {
       //  Handle the error.
       //printf ("CreateFile failed with error %d.\n", GetLastError());
       return false;
   } 
   return true;
}

bool SerialPort::InitPort(UINT BaudRate)
{
	DCB Dcb;
	bool fSuccess;
	//DCB Dcbb;
	fSuccess = GetCommState(hcom, &Dcb);	

	Dcb.BaudRate = BaudRate;     //  baud rate
	Dcb.ByteSize = 8;             //  data size, xmit and rcv
	Dcb.Parity   = NOPARITY;      //  parity bit
	Dcb.StopBits = ONESTOPBIT;    //  stop bit

	fSuccess = SetCommState(hcom, &Dcb);
	if (!fSuccess) 
	{
		//  Handle the error.
		printf ("SetCommState failed with error %d.\n", GetLastError());
		return false;
	}
	return true;
}

bool SerialPort::OpenListenThread()	//���������߳�
{

    /** ����߳��Ƿ��Ѿ������� */   
    if (m_hListenThread != INVALID_HANDLE_VALUE)  
    {  
        /** �߳��Ѿ����� */   
        return false;  
    }  
 
    s_bExit = false;  
    /** �߳�ID */   
    UINT threadId;  
    /** �����������ݼ����߳� */   
    m_hListenThread = (HANDLE)_beginthreadex(NULL, 0, ListenThread, this, 0, &threadId);  
    if (!m_hListenThread)  
    {  
        return false;  
    }  

    /** �����̵߳����ȼ�,������ͨ�߳� */   
    if (!SetThreadPriority(m_hListenThread, THREAD_PRIORITY_ABOVE_NORMAL))  
    {  
        return false;  
    }  

	return true;
}
UINT SerialPort::ListenThread(void* pParam)
{
	 /** �õ������ָ�� */   
	SerialPort *pSerialPort = reinterpret_cast<SerialPort*>(pParam);
    char str;
	while(!pSerialPort->s_bExit)
	{
		pSerialPort->ReadPort(&str);
 
		printf("%c",str);
	}

	return 0;
}
bool SerialPort::CloseListenThread()	//�رռ����߳�
{
	if (m_hListenThread != INVALID_HANDLE_VALUE)  
    {  
        /** ֪ͨ�߳��˳� */   
        s_bExit = true;  
 
        /** �ȴ��߳��˳� */   
        Sleep(10);  

        /** ���߳̾����Ч */   
        CloseHandle( m_hListenThread );  
        m_hListenThread = INVALID_HANDLE_VALUE;  
    }
	return true;
}
bool SerialPort::ReadPort(char *str)//
{
	//char str[100];
	DWORD wCount;
	bool bReadStat;
	bReadStat=ReadFile(hcom,str,buff_n,&wCount,NULL);//���������� Ϊÿ�ζ�ȡ���ֽ�������ȡ���˾ͷ��غ�����û�ж�ȡ��һֱ�ȴ�����������ȡ�������ݴ���str��
	
	if(!bReadStat)
	{
		printf("��ȡ����ʧ�ܣ�");
		return false;
	}

	return true;
}

bool SerialPort::WritePort(unsigned char* WData, unsigned int length)
{
	return true;
}
bool SerialPort::ClearPort()
{
	return PurgeComm(hcom, PURGE_RXCLEAR | PURGE_RXABORT);  //������ڻ�����
}
bool SerialPort::ClosePort()
{
	return true;
}
