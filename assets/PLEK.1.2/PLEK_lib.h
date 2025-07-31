///////////////////////////////////////////////////////////////////
//                                                               
//  PLEK - predictor of lncRNAs and mRNAs based on k-mer scheme  
//  Authors: Aimin Li, Junying Zhang                             
//  Contacts: LiAiminMail@gmail.com, jyzhang@mail.xidian.edu.cn  
//  Webcite: https://sourceforge.net/projects/plek/                      
//  Updated on: Mar 22 2014                                                 
//                                                               
///////////////////////////////////////////////////////////////////

#include <string.h>
#include <math.h>
#include <ctype.h>

char * str_lefttrim(char *dst,char *src, int lentrim)
{
    char *p_src = src;
    char *q_dst = dst;
    int len = strlen(src);
    if(lentrim>len) lentrim = len;
    while(lentrim--) *(q_dst++) = *(p_src++);
    *(q_dst++)='\0';
    return dst;
}

char * str_substring(char *dst,char *src, int sublength,int start_position) 
{
    char *p = src;
    char *q = dst;
    int len = strlen(src);
    if(sublength>len) sublength = len-start_position;   
    if(start_position<0) start_position=0;    
    if(start_position>len) return NULL;
    p += start_position;
    while(sublength--) *(q++) = *(p++);
    *(q++)='\0'; 
    return dst;
}

char * str_righttrim(char *dst,char *src, int lentrim)
{
    char *p_src = src;
    char *q_dst = dst;
    int len = strlen(src);
    if(lentrim>len) lentrim = len;
    p_src += (len-lentrim);   
    while(lentrim)
	{
		*(q_dst++) = *(p_src++); 
		lentrim--; 
	}
    return dst;
}


char *str_reverse(char *str, size_t len)
{
	
    char  *start = str;
    char  *end = str + len - 1;
    char  ch;
	
    if (str != NULL)
    {
        while (start < end)
        {
            ch = *start;
            *start++ = *end;
            *end-- = ch;
        }
    }
    return str;
}


/* get the position of the ordinal-th ch in the str */
long str_search_char(char str[], char ch, long ordinal)
{
	/* get the position of the ordinal-th ch in the str */
	long identic=0,j=0;
	char *p;
	p=str;
	while (*p!='\0')
	{
		
		if (*p==ch)
		{
			identic++;
			if (identic==ordinal)
			{
				return j;
			}
		}
		j++;
		p++;
	}
	return -1;
}

long min_long(long a,long b,long c)
{
	long t;
	t=a>b?b:a;
	return t>c?c:t;
}
