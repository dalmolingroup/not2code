///////////////////////////////////////////////////////////////////
//                                                               
//  PLEK - predictor of lncRNAs and mRNAs based on k-mer scheme  
//  Authors: Aimin Li, Junying Zhang                             
//  Contacts: LiAiminMail@gmail.com, jyzhang@mail.xidian.edu.cn  
//  Webcite: https://sourceforge.net/projects/plek/                      
//  Updated on: June 13, 2014                                                 
//                                                               
///////////////////////////////////////////////////////////////////

#include <ctype.h>
#include <time.h>
#include <math.h>
#include "PLEK_lib.h"

#ifndef MAX_LINE_LENGTH  
#define MAX_LINE_LENGTH 1000000  
#endif 

#ifndef MAX_SAMPLE_COUNT  
#define MAX_SAMPLE_COUNT 2000000  
#endif 

#ifndef REMOVE_TEMP_FILES  
#define REMOVE_TEMP_FILES 0  
#endif 

void check_two_line_fasta_file(char filename[],FILE *logs, long outmsg)
{
	FILE *fid_origin;
	long rowid=0,seq_count=0;

	char *buff_str_desc;
	char *buff_str_seq;
	
	
	buff_str_desc=(char*)malloc(MAX_LINE_LENGTH*sizeof(char));
	buff_str_seq=(char*)malloc(MAX_LINE_LENGTH*sizeof(char));
	
	
	fid_origin=fopen(filename,"r"); 
	if(fid_origin==NULL)
	{
		printf("ERROR: Cannot open %s\n",filename);
		fprintf(logs, "ERROR: Cannot open %s\n",filename);
		fflush(logs);
		return  ;
	}
	
	if(outmsg) printf("Check the format of input fasta file ... \n");
	fprintf(logs,"Check the format of input fasta file ... \n");
	fflush(logs);
	
	while(!feof(fid_origin))
	{
		/* read two lines  */
		if(fgets(buff_str_desc,MAX_LINE_LENGTH,fid_origin)==NULL)
		{
			//printf("read failed!\n");	
			break;		 
		}
		if(fgets(buff_str_seq,MAX_LINE_LENGTH,fid_origin)==NULL)
		{
			//printf("read failed!\n");	
			break;		 
		}
		if (buff_str_desc[0]==10 || buff_str_desc[0]==13)
		{
			break;
		}
		
		if (buff_str_desc[0]!='>')  /* description rows*/
		{
			printf("ERROR: INVALID SYMBOL %c, file: %s, row id: %ld\n", buff_str_desc[0],filename,rowid);
			fprintf(logs,"ERROR: INVALID SYMBOL %c, file: %s, row id: %ld\n", buff_str_desc[0],filename,rowid);
			return;
		}

		if( buff_str_seq[0]!='A' && buff_str_seq[0]!='C' && buff_str_seq[0]!='G'&& buff_str_seq[0]!='T' 
			 &&	buff_str_seq[0]!='a' && buff_str_seq[0]!='c' && buff_str_seq[0]!='g'&& buff_str_seq[0]!='t'
				&& buff_str_seq[0]!='N'&& buff_str_seq[0]!='n')
		{
			printf("ERROR: INVALID SYMBOL %c, file: %s, row id: %ld\n", buff_str_seq[0],filename,rowid);
			fprintf(logs,"ERROR: INVALID SYMBOL %c, file: %s, row id: %ld\n", buff_str_seq[0],filename,rowid);
			return;
		}

	  
		seq_count++;
	 
		
		rowid++; 
		rowid++;
		if(rowid%10000==0)
		{
			if(outmsg) printf("   %s, row id: %ld\n", filename,rowid);
			fprintf(logs,"   %s, row id: %ld\n", filename,rowid);
			fflush(logs);
		}
		
		
	}
	if(outmsg) printf("   Total lines: %ld \n", rowid );
	fprintf(logs,"   Total lines: %ld \n", rowid );
	if(outmsg) printf("   Sequence count: %ld\n\n", seq_count );
	fprintf(logs,"   Sequence count: %ld\n\n", seq_count );
	fflush(logs);
	fclose(fid_origin);
 
	free(buff_str_desc); /* free */
	free(buff_str_seq); /* free */
 
}




void extract_str_from_fasta(char file_name_origin[],long is_desc_line,long start_pos,long end_pos,char **p,long *validcount)
{
	/* for comparing two files with same ID lines */

	long rowid=0;
	//long max_line_count=100000 ;
	FILE *fid_origin ;	
	char *buff_str_desc,temp[20000];
	*validcount=0;
	


	/* line by line */	
	buff_str_desc=(char*)malloc(MAX_LINE_LENGTH*sizeof(char));	
	fid_origin=fopen(file_name_origin,"r"); 
	if(fid_origin==NULL)
	{
		printf("ERROR: Cannot open %s\n",file_name_origin);
		return ;
	}
	
	*validcount=0;
	while(!feof(fid_origin))
	{
		if(fgets(buff_str_desc,MAX_LINE_LENGTH,fid_origin)==NULL)
		{
			break;		 
		}
		
		rowid++;
		if(rowid%10000==0)
			printf("extract_str_from_fasta, row id: %ld, %s\n", rowid,file_name_origin);
		if (is_desc_line && buff_str_desc[0]=='>')
		{	
			strcpy(p[*validcount], str_substring(temp,buff_str_desc,end_pos-start_pos+1,start_pos));
			p[*validcount][end_pos]='\0';
			(*validcount)++;
		}

		if (!is_desc_line && buff_str_desc[0]!='>')
		{			
			strcpy(p[*validcount], str_substring(temp,buff_str_desc,end_pos-start_pos+1,start_pos));
			p[*validcount][end_pos]='\0';
			(*validcount)++;
		}
	 
		
	}
	printf("   Total lines: %ld\n\n", rowid );
	fclose(fid_origin);
	 
	free(buff_str_desc); /* free */
	//if (REMOVE_TEMP_FILES) remove(file_name_origin);
}


void extract_seq_with_identical_id(char file_name_seq_desc[],char file_name_valid_id[],char file_name_new[],long start_pos,long end_pos,FILE *logs)
{

	FILE *fid_origin, *file_new;
	long rowid=0,seq_count=0,j,max_line_count=200000;
 
	char *buff_str_desc;
	char *buff_str_seq;

	char **p;
	long *validcount,val_count;
	validcount=&val_count;

	p = (char **)malloc(sizeof(char *)*max_line_count); 
	for (j=0; j<max_line_count; j++)
	{
		p[j] = (char *)malloc(sizeof(char)*(end_pos-start_pos+2));  // '\0'
	}

	extract_str_from_fasta(file_name_valid_id,1,start_pos,end_pos,p,validcount);


	buff_str_desc=(char*)malloc(MAX_LINE_LENGTH*sizeof(char));
	buff_str_seq=(char*)malloc(MAX_LINE_LENGTH*sizeof(char));
	
	file_new=fopen(file_name_new,"w");
	if(file_new==NULL )
	{
		printf("ERROR: Cannot open %s\n",file_name_new);
		fprintf(logs,"ERROR: Cannot open %s\n",file_name_new);
		fflush(logs);
		return;
	}	
 
	fid_origin=fopen(file_name_seq_desc,"r"); 
	if(fid_origin==NULL)
	{
		printf("ERROR: Cannot open %s\n",file_name_seq_desc);
		fprintf(logs, "ERROR: Cannot open %s\n",file_name_seq_desc);
		fflush(logs);
		return;
	}

 
	while(!feof(fid_origin))
	{
		/* read two lines  */
		if(fgets(buff_str_desc,MAX_LINE_LENGTH,fid_origin)==NULL)
		{
			//printf("read failed!\n");	
			break;		 
		}
		if(fgets(buff_str_seq,MAX_LINE_LENGTH,fid_origin)==NULL)
		{
			//printf("read failed!\n");	
			break;		 
		}
	
		if (strlen(buff_str_desc)>12)
		{
			for (j=0;j<val_count;j++)
			{
				if (strstr(buff_str_desc,p[j]))
				{
					fprintf(file_new,"%s",buff_str_desc); 
					fprintf(file_new,"%s",buff_str_seq); 
					
					seq_count++;
				}
			}


		
		}
		
		rowid++;
		if(rowid%10000==0)
		{
			printf("filtering seq, %s, row id: %ld\n", file_name_seq_desc,rowid);
			fprintf(logs, "filtering seq, %s, row id: %ld\n", file_name_seq_desc,rowid);
			fflush(logs);
		}
		 
	}
	printf("   Total lines: %ld \n", rowid );
	fprintf(logs, "   Total lines: %ld \n", rowid );
	printf("   Sequence count: %ld\n\n", seq_count );
	fprintf(logs, "   Sequence count: %ld\n\n", seq_count );
	fflush(logs);
	fclose(fid_origin); 
	fclose(file_new); 
	free(buff_str_desc); /* free */
	free(buff_str_seq); /* free */
	if (REMOVE_TEMP_FILES) remove(file_name_seq_desc);
}


void snap_str_from_file(char file_name_origin[],char file_name_new[],char ch,long start_ordinal,long end_ordinal)
{

	/* line by line */
	FILE *fid_origin, *file_new;
	long rowid=0;
	char temp[500];
	
	char *buff_str_desc;	
	
	buff_str_desc=(char*)malloc(MAX_LINE_LENGTH*sizeof(char));	
	
	file_new=fopen(file_name_new,"w");
	if(file_new==NULL )
	{
		printf("ERROR: Cannot open %s\n",file_name_new);
		return;
	}	
	
	fid_origin=fopen(file_name_origin,"r"); 
	if(fid_origin==NULL)
	{
		printf("ERROR: Cannot open %s\n",file_name_origin);
		return;
	}
	
	while(!feof(fid_origin))
	{
		if(fgets(buff_str_desc,MAX_LINE_LENGTH,fid_origin)==NULL)
		{
			break;		 
		}
		
		rowid++;
		if(rowid%10000==0)
			printf("extracting seq, row id: %ld, %s\n", rowid,file_name_origin);
		if (str_search_char(buff_str_desc,ch,start_ordinal)>0 && str_search_char(buff_str_desc,ch,end_ordinal)>0) 
		{			
			fprintf(file_new,"%s\n", str_substring(temp,buff_str_desc,
				str_search_char(buff_str_desc,ch,end_ordinal)-str_search_char(buff_str_desc,ch,start_ordinal)-1,
				str_search_char(buff_str_desc,ch,start_ordinal)+1)); 
		}
		 
		
	}
	printf("   Total lines: %ld\n\n", rowid );
	fclose(fid_origin);
	fclose(file_new); 
	free(buff_str_desc); /* free */
	//if (REMOVE_TEMP_FILES) remove(file_name_origin);

}

void snap_str_from_file2(char file_name_origin[],char file_name_new[], long type)
{
	/*   */
	/* line by line */
	FILE *fid_origin, *file_new;
	long rowid=0;
	char temp1[50000],temp2[50000];
	
	char *buff_str_desc;	
	
	buff_str_desc=(char*)malloc(MAX_LINE_LENGTH*sizeof(char));	
	
	file_new=fopen(file_name_new,"w");
	if(file_new==NULL )
	{
		printf("ERROR: Cannot open %s\n",file_name_new);
		return;
	}	
	
	fid_origin=fopen(file_name_origin,"r"); 
	if(fid_origin==NULL)
	{
		printf("ERROR: Cannot open %s\n",file_name_origin);
		return;
	}
	
	while(!feof(fid_origin))
	{
		if(fgets(buff_str_desc,MAX_LINE_LENGTH,fid_origin)==NULL)
		{
			break;		 
		}
		
		rowid++;
		if(rowid%10000==0)
			printf("extracting seq, row id: %ld, %s\n", rowid,file_name_origin);
		
	
		if (strlen(buff_str_desc)>10 && type==1) 
		{	
			long start_char=1;
			long end_char=2;
			char ch='|';
			long position=4;
			long start_pos=str_search_char(buff_str_desc,ch,position);
			str_substring(temp1,buff_str_desc,end_char-start_char+1,start_char-1);
			str_substring(temp2,buff_str_desc,strlen(buff_str_desc)-start_pos,start_pos);
			fprintf(file_new,"%s\n",strcat(temp1,temp2)); 
		}

		
		if (strlen(buff_str_desc)>10 && type==2) 
		{		
			char ch='|';
			long start_char=1;
			long end_char=2;
			long start_ch=3;
			long end_ch=5;
			long start_pos=str_search_char(buff_str_desc,ch,start_ch);
			long end_pos=str_search_char(buff_str_desc,ch,end_ch);
			str_substring(temp1,buff_str_desc,end_char-start_char+1,start_char-1);
			str_substring(temp2,buff_str_desc,end_pos-start_pos-1,start_pos+1);
			fprintf(file_new,"%s\n",strcat(temp1,temp2)); 
		}
	
		if (strlen(buff_str_desc)>10 && type==3) 
		{	
			char ch='|';
			long start_ch=1;
			long end_ch=2;
			long start_ch2=2;
			long end_ch2=3;
			long start_pos=str_search_char(buff_str_desc,ch,start_ch);
			long end_pos=str_search_char(buff_str_desc,ch,end_ch);
			long start_pos2=str_search_char(buff_str_desc,ch,start_ch2);
			long end_pos2=str_search_char(buff_str_desc,ch,end_ch2);
			str_substring(temp1,buff_str_desc,end_pos-start_pos-1,start_pos+1);
			str_substring(temp2,buff_str_desc,end_pos2-start_pos2-1,start_pos2+1);
			
			fprintf(file_new,"%s\n",strcat(temp1,temp2)); 
		}
		
		
	}
	printf("   Total lines: %ld\n\n", rowid );
	fclose(fid_origin);
	fclose(file_new); 
	free(buff_str_desc); /* free */
	//if (REMOVE_TEMP_FILES) remove(file_name_origin);
	
}


void split_seq_and_desc(char file_name_fasta[],char file_name_seq[],
					char file_name_desc[] ,FILE * logs, long outmsg, long rmtempfile)
{


	long rowid=0,seq_count=0;
	char *buff_data;
	char *seq;
	FILE *fid_seq,*fid_desc, *fid_reads; 
 
	long is_desc_line=1;
	long is_seq_line=0;

	buff_data=(char*)malloc(MAX_LINE_LENGTH*sizeof(char));
	seq=(char*)malloc(MAX_LINE_LENGTH*sizeof(char));
	fid_reads=fopen(file_name_fasta,"r"); 
	if(fid_reads==NULL)
	{
		printf("ERROR: Cannot open %s\n",file_name_fasta);
		fprintf(logs, "ERROR: Cannot open %s\n",file_name_fasta);
		fflush(logs);
		return;
	}
	fid_seq=fopen(file_name_seq,"w");
	fid_desc=fopen(file_name_desc,"w");
	if(fid_seq==NULL || fid_desc==NULL)
	{
		printf("ERROR: Cannot open %s\n",file_name_seq);
		fprintf(logs, "ERROR: Cannot open %s\n",file_name_seq);
		fflush(logs);
		return;
	}	 

	if(outmsg) printf("Split define lines and nucleotide lines ... \n");
	fprintf(logs,"Split define lines and nucleotide lines ... \n");
	fflush(logs);


	rowid=0;
	while(!feof(fid_reads))
	{
		if(fgets(buff_data,MAX_LINE_LENGTH,fid_reads)==NULL)
		{
			//printf("read file failed ( in %s )!\n",file_name_fasta);	
			break;
		}
		else
		{
			if(buff_data[0]==10 || buff_data[0]==13 )
				break;
		
			if(buff_data[0]=='>') /* description */
			{
			   if(is_seq_line==1) /*   */
			   {
		      		fprintf(fid_seq,"%s",seq); 
					seq_count++;
					is_seq_line=0;
			   }
			   else
			   {
			   
			   }
			   fprintf(fid_desc,"%s",buff_data); 
			   is_desc_line=1;
			}
			else
			{
				
			
				if(is_desc_line==1)
				{  
					strcpy(seq,buff_data);
					is_seq_line=1;
				}
				else
				
				{
					long slenth=strlen(seq);
					if(seq[slenth-1]==10) 
					{ 
						seq[slenth-1]='\0'; 
					}
					strcat(seq,buff_data);
					is_seq_line=1;
				}
				is_desc_line=0;

			}

			if(rowid%10000==0)
			{
				if(outmsg) printf("   %s, row index: %ld\n", file_name_fasta,rowid);
				fprintf(logs,"   %s, row index: %ld\n", file_name_fasta,rowid);
				fflush(logs);
			}
			
			rowid++; /* row++ */
		
		}
	}
	if(is_seq_line==1)
	{
		fprintf(fid_seq,"%s",seq); 
		seq_count++;
		is_seq_line=0;
	}
	if(outmsg) printf("   Total lines: %ld\n",rowid);
	fprintf(logs,"   Total lines: %ld\n",rowid);
	if(outmsg) printf("   Sequence count: %ld\n\n",seq_count);
	fprintf(logs,"   Sequence count: %ld\n\n",seq_count);
	
	fflush(logs);
	fclose(fid_reads);
	fclose(fid_seq);
	fclose(fid_desc); 
	free(buff_data);   /* free */
	free(seq);   /* free */ 
	if (rmtempfile) remove(file_name_fasta);
}




void to_one_line(char file_name_origin[],char new_file_name[],FILE *logs, long outmsg)
{

	long rowid=0,seq_count=0;
	char *buff_data, tempchar;
	char *seq;
	FILE *fid_new,*fid_origin;	 
	long is_desc_line=1; /* is > lines/rows  */
	long is_seq_line=0; /* is sequence data */
	long slenth;

	buff_data=(char*)malloc(MAX_LINE_LENGTH*sizeof(char));
	seq=(char*)malloc(MAX_LINE_LENGTH*sizeof(char));
	fid_origin=fopen(file_name_origin,"r"); 
	if(fid_origin==NULL)
	{
		printf("WARNING: Cannot open file: %s.\n",file_name_origin);
		fprintf(logs,"WARNING: Cannot open file: %s.\n",file_name_origin);
		fflush(logs);
		return;
	}

	fid_new=fopen(new_file_name,"w");
	if(fid_new==NULL )
	{
		puts("WARNING: Cannot open file\n");
		fprintf(logs,"WARNING: Cannot open file\n");
		fflush(logs);
		return;
	} 

	if(outmsg) printf("Change the format of the input fasta file ... \n");
	fprintf(logs,"Change the format of the input fasta file ... \n");
	fflush(logs);


	while(!feof(fid_origin))
	{
		if(fgets(buff_data,MAX_LINE_LENGTH,fid_origin)==NULL)
		{
			//printf("read failed !\n");	
			break;
		}
		else
		{
			tempchar=buff_data[0];
			if(tempchar==10 || tempchar==13) /* is carriage return or line break */
				break;
			
			if(tempchar=='>' || (tempchar!='A' && tempchar!='C' && tempchar!='G'&& tempchar!='T' 
				&& tempchar!='a' && tempchar!='c' && tempchar!='g'&& tempchar!='t' 
				&& tempchar!='N' && tempchar!='n' )) /* description rows*/
			{

				if(tempchar!='>') buff_data[0]='>';

			   if(is_seq_line==1) /* if the last row is seq , then write */
			   {
		      		fprintf(fid_new,"%s",seq); /* write seq */
					seq_count++;
					is_seq_line=0;
			   }
			   else
			   {
			   
			   }
			   fprintf(fid_new,"%s",buff_data);  /* write desc */
			   is_desc_line=1;
			}
			else
			{
				
				
				if(is_desc_line==1)
				{  
					strcpy(seq,buff_data);
					is_seq_line=1;
				}
				else
					
				{
					slenth=strlen(seq);
					if(seq[slenth-1]==10 || seq[slenth-1]==13) 
					{ 
						seq[slenth-1]='\0'; 
					}
					strcat(seq,buff_data);
					is_seq_line=1;
				}
				is_desc_line=0;

			}

			// printf("to one line, %s, row id: %ld\n", file_name_origin,rowid);
			
			 
			if(rowid%100000==0)
			{
				if(outmsg) printf("   %s, row id: %ld\n", file_name_origin,rowid);
				fprintf(logs,"   %s, row id: %ld\n", file_name_origin,rowid);
				fflush(logs);
			}
			 
			rowid++; /* row++ */
		}
	}
	if(is_seq_line==1)
	{
		fprintf(fid_new,"%s",seq); 
		seq_count++;
		is_seq_line=0;
	}
	if(outmsg) printf("   Total lines: %ld\n", rowid );
	fprintf(logs, "   Total lines: %ld\n", rowid );
	if(outmsg) printf("   Sequence count: %ld\n\n", seq_count );
	fprintf(logs, "   Sequence count: %ld\n\n", seq_count );

	
	fflush(logs);
	fclose(fid_origin);
	fclose(fid_new);
	free(buff_data); /* free */
	free(seq); /* free */
	fflush(logs);
	/* remove(file_name_origin);	*/
}


void merge_two_file_to_one(char file_name1[],char file_name2[],char file_name_new[])
{

	FILE *fid_1,*fid_2,*file_new;
	char *buff_data;

	printf("merging %s and %s ... \n",file_name1,file_name2);

	buff_data=(char*)malloc(MAX_LINE_LENGTH*sizeof(char));
	file_new=fopen(file_name_new,"w");
	if(file_name_new==NULL )
	{
		puts("ERROR: Cannot open files\n");
		return;
	}

	/* fid_1 */
	fid_1=fopen(file_name1,"r"); 
	if(fid_1==NULL)
	{
		puts("ERROR: Cannot open  \n");
		return;
	}
	
	while(!feof(fid_1))
	{
		if(fgets(buff_data,MAX_LINE_LENGTH,fid_1)==NULL)
		{
			break;		 
		}
		else
		{
		 	fprintf(file_new,"%s",buff_data); 
		}
	}
	fclose(fid_1);

	/* fid_2 */
	fid_2=fopen(file_name2,"r"); 
	if(fid_2==NULL)
	{
		puts("ERROR: Cannot open.\n");
		return;
	}
	while(!feof(fid_1))
	{
		if(fgets(buff_data,MAX_LINE_LENGTH,fid_1)==NULL)
		{
			break;		 
		}
		else
		{
		 	fprintf(file_new,"%s",buff_data); 
		}
	}

	free(buff_data); /* free */
	fclose(fid_2);
	fclose(file_new); 
	if (REMOVE_TEMP_FILES) remove(file_name1);
	if (REMOVE_TEMP_FILES) remove(file_name2);
}

void fasta_file_get_some_rows(char file_name1[],char file_name_new[],long rowcout)
{

	FILE *fid_1,*file_new; 
	long i=0;
	char *buff_data;
	
	buff_data=(char*)malloc(MAX_LINE_LENGTH*sizeof(char));
	file_new=fopen(file_name_new,"w");
	if(file_name_new==NULL )
	{
		puts("ERROR: Cannot open files\n");
		return;
	}
	
	/* fid_1 */
	fid_1=fopen(file_name1,"r"); 
	if(fid_1==NULL)
	{
		puts("ERROR: Cannot open  \n");
		return;
	}
	
	while(!feof(fid_1))
	{
		if(fgets(buff_data,MAX_LINE_LENGTH,fid_1)==NULL)
		{
			//puts("read a line from file failed (in posit)!\n");	
			break;		 
		}
		else
		{
			i++;
			fprintf(file_new,"%s",buff_data); 
			if (i==rowcout)
			{
				break;
			}

		}
	}
	free(buff_data); /* free */
	fclose(fid_1);
	fclose(file_new); 
}



long remove_short_sequence(char file_name_origin[],char file_name_new[],
						   size_t length, FILE *logs, long outmsg, long rmtempfile, 
						   long is_known_classlabel)
{

	FILE *fid_origin, *file_new;
	long rowid=0,seq_count=0;
	long max_line_length=0;
	long len, removed_count=0;
 
	char *buff_str_desc;
	char *buff_str_seq;


	buff_str_desc=(char*)malloc(MAX_LINE_LENGTH*sizeof(char));
	buff_str_seq=(char*)malloc(MAX_LINE_LENGTH*sizeof(char));
	
	file_new=fopen(file_name_new,"w");
	if(file_new==NULL )
	{
		printf("ERROR: Cannot open %s\n",file_name_new);
		fprintf(logs, "ERROR: Cannot open %s\n",file_name_new);
		fflush(logs);
		return -1;
	}	
 
	fid_origin=fopen(file_name_origin,"r"); 
	if(fid_origin==NULL)
	{
		printf("ERROR: Cannot open %s\n",file_name_origin);
		fprintf(logs, "ERROR: Cannot open %s\n",file_name_origin);
		fflush(logs);
		return -1;
	}

	if(outmsg) printf("Remove sequences whose lengths are not more than minlength  ... \n");
	fprintf(logs,"Remove sequences whose lengths are not more than minlength   ... \n");
	fflush(logs);
	
	while(!feof(fid_origin))
	{
		/* read two lines  */
		if(fgets(buff_str_desc,MAX_LINE_LENGTH,fid_origin)==NULL)
		{
			//printf("read failed!\n");	
			break;		 
		}
		if(fgets(buff_str_seq,MAX_LINE_LENGTH,fid_origin)==NULL)
		{
			//printf("read failed!\n");	
			break;		 
		}
		if (buff_str_desc[0]==10 || buff_str_desc[0]==13)
		{
			break;
		}
	
		len=strlen(buff_str_seq);
		if (len>max_line_length)
		{
			max_line_length=len;
		}
		
		if (len>(long)length)
		{
			fprintf(file_new,"%s",buff_str_desc); 
			fprintf(file_new,"%s",buff_str_seq); 
			seq_count++;
		}
		else
		{
			removed_count++;
		}
		
		rowid++; rowid++;
		if(rowid%10000==0)
		{
			if(outmsg) printf("   %s, row id: %ld\n", file_name_origin,rowid);
			fprintf(logs,"   %s, row id: %ld\n", file_name_origin,rowid);
			fflush(logs);
		}

		 
	}
	
	if(seq_count>MAX_SAMPLE_COUNT)
	{
		printf("Please change \'#define MAX_SAMPLE_COUNT 2000000\' to \'#define MAX_SAMPLE_COUNT 10000000\' in the following files:\nPLEK_kmer.h\nPLEK_spsn.c\nPLEK_main.c\nPLEK_fastafile.h\n");
	}
	
	if(removed_count>0) /* some transcripts were removed. */
	{
		if(is_known_classlabel==1)
			if(strstr(file_name_origin, "pos_")>0) printf("	WARNING: mRNAs, %ld short sequences (<=%ld nt) were removed. %ld retained.\n", removed_count, length, seq_count );
			else printf("	WARNING: ncRNAs, %ld short sequences (<=%ld nt) were removed. %ld retained.\n", removed_count, length, seq_count );
		else
			printf("	WARNING: %ld short sequences (<=%ld nt) were removed. %ld retained.\n", removed_count, length, seq_count );
	}
	else
	{
		if(is_known_classlabel==1)
			if(strstr(file_name_origin, "pos_")>0) printf("	mRNAs, %ld sequences.\n", seq_count );
			else printf("	ncRNAs, %ld sequences.\n", seq_count);
		else
			printf("	%ld sequences.\n", seq_count );
	}
		
	
	if(outmsg) printf("   Total lines: %ld \n", rowid );
	fprintf(logs,"   Total lines: %ld \n", rowid );
	if(outmsg) printf("   Sequence count: %ld\n\n", seq_count );
	fprintf(logs,"   Sequence count: %ld\n\n", seq_count );
	fflush(logs);
	fclose(fid_origin);
	fclose(file_new); 
	free(buff_str_desc); /* free */
	free(buff_str_seq); /* free */
	if (rmtempfile)  remove(file_name_origin);

	return max_line_length;
}



void extract_seq_with_special_str(char file_name_origin[],char file_name_new[],char str[],FILE * logs)
{


	FILE *fid_origin, *file_new;
	long rowid=0,seq_count=0;
 
	char *buff_str_desc;
	char *buff_str_seq;

	buff_str_desc=(char*)malloc(MAX_LINE_LENGTH*sizeof(char));
	buff_str_seq=(char*)malloc(MAX_LINE_LENGTH*sizeof(char));
	
	file_new=fopen(file_name_new,"w");
	if(file_new==NULL )
	{
		printf("ERROR: Cannot open %s\n",file_name_new);
		fprintf(logs,"ERROR: Cannot open %s\n",file_name_new);
		fflush(logs);
		return;
	}	
 
	fid_origin=fopen(file_name_origin,"r"); 
	if(fid_origin==NULL)
	{
		printf("ERROR: Cannot open %s\n",file_name_origin);
		fprintf(logs, "ERROR: Cannot open %s\n",file_name_origin);
		fflush(logs);
		return;
	}

	 
	printf("filtering seq ... key: %s\n",str);
	fprintf(logs, "filtering seq ... key: %s\n",str);
	while(!feof(fid_origin))
	{
		/* read two lines  */
		if(fgets(buff_str_desc,MAX_LINE_LENGTH,fid_origin)==NULL)
		{
			//printf("read failed!\n");	
			break;		 
		}
		if(fgets(buff_str_seq,MAX_LINE_LENGTH,fid_origin)==NULL)
		{
			//printf("read failed!\n");	
			break;		 
		}
	
		if (strstr(buff_str_desc,str))
		{
			fprintf(file_new,"%s",buff_str_desc); 
			fprintf(file_new,"%s",buff_str_seq); 
			seq_count++;
		}
		
		rowid++;
		if(rowid%10000==0)
		{
			printf("filtering seq, %s, row id: %ld\n", file_name_origin,rowid);
			fprintf(logs, "filtering seq, %s, row id: %ld\n", file_name_origin,rowid);
			fflush(logs);
		}
		 
	}
	printf("   Total lines: %ld \n", rowid );
	fprintf(logs, "   Total lines: %ld \n", rowid );
	printf("   Sequence count: %ld\n\n", seq_count );
	fprintf(logs, "   Sequence count: %ld\n\n", seq_count );
	fflush(logs);
	fclose(fid_origin); 
	fclose(file_new); 
	free(buff_str_desc); /* free */
	free(buff_str_seq); /* free */
	if (REMOVE_TEMP_FILES) remove(file_name_origin);
}


void extract_seq_with_special_str2(char file_name_origin[],char file_name_origin2[],char file_name_new[],FILE * logs)
{
 


	FILE *fid_origin, *file_new,*fid_origin2;
	long rowid=0,seq_count=0;
 
	char *buff_str_desc;
	char *buff_str_seq;
	char *buff_str_desc2;
	
	char str_geneid[50];

	buff_str_desc=(char*)malloc(MAX_LINE_LENGTH*sizeof(char));
	buff_str_seq=(char*)malloc(MAX_LINE_LENGTH*sizeof(char));
		buff_str_desc2=(char*)malloc(MAX_LINE_LENGTH*sizeof(char));
	
	file_new=fopen(file_name_new,"w");
	if(file_new==NULL )
	{
		printf("ERROR: Cannot open %s\n",file_name_new);
		fprintf(logs,"ERROR: Cannot open %s\n",file_name_new);
		fflush(logs);
		return;
	}	
 
	fid_origin=fopen(file_name_origin,"r"); 
	if(fid_origin==NULL)
	{
		printf("ERROR: Cannot open %s\n",file_name_origin);
		fprintf(logs, "ERROR: Cannot open %s\n",file_name_origin);
		fflush(logs);
		return;
	}
	fid_origin2=fopen(file_name_origin2,"r"); 
	if(fid_origin2==NULL)
	{
		printf("ERROR: Cannot open %s\n",file_name_origin2);
		fprintf(logs, "ERROR: Cannot open %s\n",file_name_origin2);
		fflush(logs);
		return;
	}
	 
 
	while(!feof(fid_origin))
	{
		/* read two lines  */
		if(fgets(buff_str_desc,MAX_LINE_LENGTH,fid_origin)==NULL)
		{
			//printf("read failed!\n");	
			break;		 
		}
		if(fgets(buff_str_seq,MAX_LINE_LENGTH,fid_origin)==NULL)
		{
			//printf("read failed!\n");	
			break;		 
		}
	
		//>ENST00000473358.1|ENSG00000243485.1|OTTHUMG00000000959.2|OTTHUMT00000002840.1|MIR1302-2-001|MIR1302-2|712|
		// get ENST00000473358.1 , then determine this exist in file 2.

		str_substring(str_geneid,buff_str_desc,17,1);  

		rewind(fid_origin2) ;
		while(!feof(fid_origin2))
		{
			/* read two lines  */
			if(fgets(buff_str_desc2,MAX_LINE_LENGTH,fid_origin2)==NULL)
			{
				//printf("read failed!\n");	
				break;		 
			}
			if (strstr(buff_str_desc2,str_geneid))
			{
				fprintf(file_new,"%s",buff_str_desc); 
				fprintf(file_new,"%s",buff_str_seq); 
				seq_count++;
				break;

			}

		}

 
		
		rowid++;
		if(rowid%1000==0)
		{
			printf("filtering seq, %s, row id: %ld\n", file_name_origin,rowid);
			fprintf(logs, "filtering seq, %s, row id: %ld\n", file_name_origin,rowid);
			fflush(logs);
		}
		 
	}
	printf("   Total lines: %ld \n", rowid );
	fprintf(logs, "   Total lines: %ld \n", rowid );
	printf("   Sequence count: %ld\n\n", seq_count );
	fprintf(logs, "   Sequence count: %ld\n\n", seq_count );
	fflush(logs);
	fclose(fid_origin); fclose(fid_origin2); 
	fclose(file_new); 
	free(buff_str_desc); /* free */
	free(buff_str_seq); /* free */
	free(buff_str_desc2);
	if (REMOVE_TEMP_FILES) remove(file_name_origin);
}


void extract_seq_without_special_str(char file_name_origin[],char file_name_new[],char str[],FILE * logs)
{
	

	FILE *fid_origin, *file_new;
	long rowid=0,seq_count=0;
 
	char *buff_str_desc;
	char *buff_str_seq;

	buff_str_desc=(char*)malloc(MAX_LINE_LENGTH*sizeof(char));
	buff_str_seq=(char*)malloc(MAX_LINE_LENGTH*sizeof(char));
	
	file_new=fopen(file_name_new,"w");
	if(file_new==NULL )
	{
		printf("ERROR: Cannot open %s\n",file_name_new);
		fprintf(logs,"ERROR: Cannot open %s\n",file_name_new);
		fflush(logs);
		return;
	}	
 
	fid_origin=fopen(file_name_origin,"r"); 
	if(fid_origin==NULL)
	{
		printf("ERROR: Cannot open %s\n",file_name_origin);
		fprintf(logs, "ERROR: Cannot open %s\n",file_name_origin);
		fflush(logs);
		return;
	}

	 
	printf("filtering seq ... key: %s\n",str);
	fprintf(logs, "filtering seq ... key: %s\n",str);
	while(!feof(fid_origin))
	{
		/* read two lines  */
		if(fgets(buff_str_desc,MAX_LINE_LENGTH,fid_origin)==NULL)
		{
			//printf("read failed!\n");	
			break;		 
		}
		if(fgets(buff_str_seq,MAX_LINE_LENGTH,fid_origin)==NULL)
		{
			//printf("read failed!\n");	
			break;		 
		}
	
		if (strstr(buff_str_desc,str))
		{
			
		}
		else
		{
			fprintf(file_new,"%s",buff_str_desc); 
			fprintf(file_new,"%s",buff_str_seq); 
			seq_count++;
		}
		
		rowid++;
		if(rowid%10000==0)
		{
			printf("filtering seq, %s, row id: %ld\n", file_name_origin,rowid);
			fprintf(logs, "filtering seq, %s, row id: %ld\n", file_name_origin,rowid);
			fflush(logs);
		}
		 
	}
	printf("   Total lines: %ld \n", rowid );
	fprintf(logs, "   Total lines: %ld \n", rowid );
	printf("   Sequence count: %ld\n\n", seq_count );
	fprintf(logs, "   Sequence count: %ld\n\n", seq_count );
	fflush(logs);
	fclose(fid_origin); 
	fclose(file_new); 
	free(buff_str_desc); /* free */
	free(buff_str_seq); /* free */
	if (REMOVE_TEMP_FILES) remove(file_name_origin);
}

void extract_seq_with_special_seqstr(char file_name_origin[],char file_name_new[],char seq_startedwithstr[],FILE * logs)
{


	FILE *fid_origin, *file_new;
	long rowid=0,seq_count=0;
 
	char *buff_str_desc;
	char *buff_str_seq;

	buff_str_desc=(char*)malloc(MAX_LINE_LENGTH*sizeof(char));
	buff_str_seq=(char*)malloc(MAX_LINE_LENGTH*sizeof(char));
	
	file_new=fopen(file_name_new,"w");
	if(file_new==NULL )
	{
		printf("ERROR: Cannot open %s\n",file_name_new);
		fprintf(logs,"ERROR: Cannot open %s\n",file_name_new);
		fflush(logs);
		return;
	}	
 
	fid_origin=fopen(file_name_origin,"r"); 
	if(fid_origin==NULL)
	{
		printf("ERROR: Cannot open %s\n",file_name_origin);
		fprintf(logs, "ERROR: Cannot open %s\n",file_name_origin);
		fflush(logs);
		return;
	}

	printf("filtering seq ... key: %s\n",seq_startedwithstr);
	fprintf(logs, "filtering seq ... key: %s\n",seq_startedwithstr);
	while(!feof(fid_origin))
	{
		/* read two lines  */
		if(fgets(buff_str_desc,MAX_LINE_LENGTH,fid_origin)==NULL)
		{
			//printf("read failed!\n");	
			break;		 
		}
		if(fgets(buff_str_seq,MAX_LINE_LENGTH,fid_origin)==NULL)
		{
			//printf("read failed!\n");	
			break;		 
		}
		rowid++;
		if(rowid%10000==0)
		{
			printf("filtering seq, %s, row id: %ld\n",file_name_origin, rowid);
			fprintf(logs, "filtering seq, %s, row id: %ld\n" ,file_name_origin, rowid);
			fflush(logs);
		}
		if ( buff_str_seq[0]==seq_startedwithstr[0] && buff_str_seq[1]==seq_startedwithstr[1] && buff_str_seq[2]==seq_startedwithstr[2] )
		{
			fprintf(file_new,"%s",buff_str_desc); 
			fprintf(file_new,"%s",buff_str_seq); 
			seq_count++;
		}			
		 
	}
	printf("   Total lines: %ld \n", rowid );
	fprintf(logs, "   Total lines: %ld \n", rowid );
	printf("   Sequence count: %ld\n\n", seq_count );
	fprintf(logs, "   Sequence count: %ld\n\n", seq_count );
	fflush(logs);
	fclose(fid_origin); 
	fclose(file_new); 
	free(buff_str_desc); /* free */
	free(buff_str_seq); /* free */
	if (REMOVE_TEMP_FILES) remove(file_name_origin);
}

//
void extract_seq_remove_two_ends(char file_name_origin[],char file_name_new[],FILE * logs)
{

	
	FILE *fid_origin, *file_new;
	long rowid=0,seq_count=0;
	
	char *buff_str_desc;
	char *buff_str_seq,*temp;
	
	//start_coden=['ATG'], stop_coden=['TAG','TAA','TGA']
	long ATG_postion=0;
	long TAG_postion=0;
	long TAA_postion=0;
	long TGA_postion=0;
	long start=0,end=0;
	long i=0,seq_len=0,min_len=1000000;

	buff_str_desc=(char*)malloc(MAX_LINE_LENGTH*sizeof(char));
	buff_str_seq=(char*)malloc(MAX_LINE_LENGTH*sizeof(char));
	temp=(char*)malloc(MAX_LINE_LENGTH*sizeof(char));
	
	file_new=fopen(file_name_new,"w");
	if(file_new==NULL )
	{
		printf("ERROR: Cannot open %s\n",file_name_new);
		fprintf(logs,"ERROR: Cannot open %s\n",file_name_new);
		fflush(logs);
		return;
	}	
	
	fid_origin=fopen(file_name_origin,"r"); 
	if(fid_origin==NULL)
	{
		printf("ERROR: Cannot open %s\n",file_name_origin);
		fprintf(logs, "ERROR: Cannot open %s\n",file_name_origin);
		fflush(logs);
		return;
	}
	
	printf("filtering seq ... two ends\n");
	fprintf(logs, "filtering seq ... two ends\n");
	while(!feof(fid_origin))
	{
		i=0;
		seq_len=0;
		min_len=1000000;

		/* read two lines  */
		if(fgets(buff_str_desc,MAX_LINE_LENGTH,fid_origin)==NULL)
		{
			//printf("read failed!\n");	
			break;		 
		}
		if(fgets(buff_str_seq,MAX_LINE_LENGTH,fid_origin)==NULL)
		{
			//printf("read failed!\n");	
			break;		 
		}
		rowid++;
		if(rowid%10000==0)
		{
			printf("filtering seq, %s, row id: %ld\n",file_name_origin, rowid);
			fprintf(logs, "filtering seq, %s, row id: %ld\n" ,file_name_origin, rowid);
			fflush(logs);
		}

		//start_coden=['ATG'], stop_coden=['TAG','TAA','TGA']
		i=0;
		ATG_postion=0;
		while (buff_str_seq[i+2]!='\0')
		{
			if (buff_str_seq[i]=='A' && buff_str_seq[i+1]=='T' && buff_str_seq[i+2]=='G')
			{
				ATG_postion=i;
				break;
			}
			i++;
		}
		i=2;
		seq_len=strlen(buff_str_seq);
		TAG_postion=seq_len-2;
		TAA_postion=seq_len-2;
		TGA_postion=seq_len-2;
		while (seq_len-i>2)
		{
			if (buff_str_seq[seq_len-i-2]=='T' && buff_str_seq[seq_len-i-1]=='A' && buff_str_seq[seq_len-i]=='G')
			{
				TAG_postion=seq_len-i;
				break;
			}
			if (buff_str_seq[seq_len-i-2]=='T' && buff_str_seq[seq_len-i-1]=='A' && buff_str_seq[seq_len-i]=='A')
			{
				TAA_postion=seq_len-i;
				break;
			}
			if (buff_str_seq[seq_len-i-2]=='T' && buff_str_seq[seq_len-i-1]=='G' && buff_str_seq[seq_len-i]=='A')
			{
				TGA_postion=seq_len-i;
				break;
			}
			i++;
		}

		start=ATG_postion;
		end=min_long(TAG_postion,TAA_postion,TGA_postion);
		if (end-start<min_len)
		{
			min_len=end-start;		
		}
		//150
		if (end-start<150 && start>0)
		{
			start=0;
		}
		if (end-start<150 && start==0)
		{
			end=seq_len-2;
		}

		fprintf(logs, "seq_len=%ld, new_len=%ld, diff_len=%ld\n",seq_len, end-start,seq_len-(end-start));

		fprintf(file_new,"%s",buff_str_desc); 
		fprintf(file_new,"%s\n",str_substring(temp, buff_str_seq,end-start+1,start)); 
		fflush(file_new);
		seq_count++;
		 		
		
	}
	fprintf(logs, "min_len=%ld\n", min_len);
	printf("   Total lines: %ld \n", rowid );
	fprintf(logs, "   Total lines: %ld \n", rowid );
	printf("   Sequence count: %ld\n\n", seq_count );
	fprintf(logs, "   Sequence count: %ld\n\n", seq_count );
	fflush(logs);
	fclose(fid_origin); 
	fclose(file_new); 
	free(temp);
	free(buff_str_desc); /* free */
	free(buff_str_seq); /* free */
	if (REMOVE_TEMP_FILES) remove(file_name_origin);
}


void extract_lines_without_special_str(char file_name_origin[],char file_name_new[],char str[])
{
	FILE *fid_origin, *file_new;
	long rowid=0; 
	char *buff_str_desc; 

	buff_str_desc=(char*)malloc(MAX_LINE_LENGTH*sizeof(char)); 
	
	file_new=fopen(file_name_new,"w");
	if(file_new==NULL )
	{
		printf("ERROR: Cannot open %s\n",file_name_new);
		return;
	}	
 
	fid_origin=fopen(file_name_origin,"r"); 
	if(fid_origin==NULL)
	{
		printf("ERROR: Cannot open %s\n",file_name_origin);
		return;
	}
	
	while(!feof(fid_origin))
	{
		if(fgets(buff_str_desc,MAX_LINE_LENGTH,fid_origin)==NULL)
		{
			break;		 
		}
	 
		rowid++;
		if(rowid%10000==0)
				printf("filtering seq, row id: %ld, %s\n", rowid,file_name_origin);
		if (strstr(buff_str_desc,str))
		{
		
		 
		}
		else	
		{
			fprintf(file_new,"%s",buff_str_desc); 

		}
		 
	}
	printf("   Total lines: %ld\n\n", rowid );
	fclose(fid_origin);
	fclose(file_new); 
	free(buff_str_desc);  /* free */
	//if (REMOVE_TEMP_FILES) remove(file_name_origin);
}

void extract_lines_with_special_str(char file_name_origin[],char file_name_new[],char str[])
{
	/* line by line */
	FILE *fid_origin, *file_new;
	long rowid=0;
	
	char *buff_str_desc;	
	
	buff_str_desc=(char*)malloc(MAX_LINE_LENGTH*sizeof(char));	
	
	file_new=fopen(file_name_new,"w");
	if(file_new==NULL )
	{
		printf("ERROR: Cannot open %s\n",file_name_new);
		return;
	}	
	
	fid_origin=fopen(file_name_origin,"r"); 
	if(fid_origin==NULL)
	{
		printf("ERROR: Cannot open %s\n",file_name_origin);
		return;
	}
	
	while(!feof(fid_origin))
	{
		if(fgets(buff_str_desc,MAX_LINE_LENGTH,fid_origin)==NULL)
		{
			break;		 
		}
		
		rowid++;
		if(rowid%10000==0)
			printf("filtering seq, row id: %ld, %s\n", rowid,file_name_origin);
		if (strstr(buff_str_desc,str))
		{			
				fprintf(file_new,"%s",buff_str_desc); 
		}
		else	
		{	
			
		}
		
	}
	printf("   Total lines: %ld\n\n", rowid );
	fclose(fid_origin);
	fclose(file_new); 
	free(buff_str_desc); /* free */
	//if (REMOVE_TEMP_FILES) remove(file_name_origin);
}



void extract_first_lines(char file_name_origin[],char file_name_new[],long linecount)
{
	
	FILE *fid_origin, *file_new;
	long rowid=0;
	
	char *buff_str;
	buff_str=(char*)malloc(MAX_LINE_LENGTH*sizeof(char));
	
	
	file_new=fopen(file_name_new,"w");
	if(file_new==NULL )
	{
		printf("ERROR: Cannot open %s\n",file_name_new);
		return;
	}	
	
	fid_origin=fopen(file_name_origin,"r"); 
	if(fid_origin==NULL)
	{
		printf("ERROR: Cannot open %s\n",file_name_origin);
		return;
	}
	
	while(!feof(fid_origin))
	{
		if(fgets(buff_str,MAX_LINE_LENGTH,fid_origin)==NULL)
		{
			break;		 
		}	
		rowid++;	
		fprintf(file_new,"%s",buff_str); 
		if (rowid>=linecount)
		{
			break;	
		}
	}
	
	fclose(fid_origin);
	fclose(file_new); 
	free(buff_str); /* free */
}

void get_longest_transcript(char file_name_origin[],char file_name_new[],char splitter[],long position,FILE * logs)
{


	FILE *fid_origin, *file_new;
	long rowid=0,i,j,maxlength,maxlenpos;
	
	char *buff_str_desc;
	char *buff_str_seq;
	char *p,**gene_name; /*==============*/
	long *length,*ismax,long_transcript_count=0;
	
	gene_name = (char**)malloc(sizeof(char*)*MAX_SAMPLE_COUNT);
	buff_str_desc=(char*)malloc(MAX_LINE_LENGTH*sizeof(char));
	buff_str_seq=(char*)malloc(MAX_LINE_LENGTH*sizeof(char));	
	length=(long*)malloc(MAX_SAMPLE_COUNT*sizeof(long));
	ismax=(long*)malloc(MAX_SAMPLE_COUNT*sizeof(long));

	for (i=0;i<MAX_SAMPLE_COUNT;i++)
	{
		ismax[i]=-1; /* 0,1 */
	}
	
	file_new=fopen(file_name_new,"w");
	if(file_new==NULL )
	{
		printf("ERROR: Cannot open %s\n",file_name_new);
		fprintf(logs,"ERROR: Cannot open %s\n",file_name_new);
		fflush(logs);
		return;
	}	
	
	fid_origin=fopen(file_name_origin,"r"); 
	if(fid_origin==NULL)
	{
		printf("ERROR: Cannot open %s\n",file_name_origin);
		fprintf(logs, "ERROR: Cannot open %s\n",file_name_origin);
		fflush(logs);
		return;
	}

	printf("longest transcript ... \n");
	fprintf(logs,"longest transcript ... \n");
	fflush(logs);


	/* get gene_name and seq length */
	while(!feof(fid_origin))
	{
		/* read two lines  */
		if(fgets(buff_str_desc,MAX_LINE_LENGTH,fid_origin)==NULL)
		{
			//printf("read failed!\n");	
			break;		 
		}
		if(fgets(buff_str_seq,MAX_LINE_LENGTH,fid_origin)==NULL)
		{
			//printf("read failed!\n");	
			break;		 
		}
	
		if(rowid%10000==0)
		{
			printf("longest transcript, %s, row id: %ld\n", file_name_origin,rowid);
			fprintf(logs,"longest transcript, %s, row id: %ld\n", file_name_origin,rowid);
			fflush(logs);
	}
		gene_name[rowid]=(char*)malloc(30*sizeof(char));
		p=strtok(buff_str_desc,splitter);
		for (i=0;i<position;i++)
		{
			p=strtok(NULL,splitter);
		}
	 
		strcpy(gene_name[rowid],p);
		//puts(gene_name[rowid]);

		length[rowid]=strlen(buff_str_seq);

		rowid++;		
	}
	/* compare and mark max */
	for (i=0;i<rowid;i++)
	{
		
		if (ismax[i]==-1)/* has not been compared */
		{
			maxlength=length[i];
			maxlenpos=i;
			ismax[i]=1;
			for (j=i+1;j<rowid;j++)
			{
				if (strcmp(gene_name[i],gene_name[j])==0)
				{
					if (maxlength>length[j])
					{					
						ismax[j]=0;
					}
					else
					{
						ismax[maxlenpos]=0;
						maxlength=length[j];
						maxlenpos=j;						
					}
					
				}
			}
		}
	}


	/* extract  */
	rewind(fid_origin);
	rowid=0;
	while(!feof(fid_origin))
	{
		/* read two lines  */
		if(fgets(buff_str_desc,MAX_LINE_LENGTH,fid_origin)==NULL)
		{
			//printf("read failed!\n");	
			break;		 
		}
		if(fgets(buff_str_seq,MAX_LINE_LENGTH,fid_origin)==NULL)
		{
			//printf("read failed!\n");	
			break;		 
		}
		
		if (ismax[rowid]==1)
		{
			fprintf(file_new,"%s",buff_str_desc); 
			fprintf(file_new,"%s",buff_str_seq); 
			long_transcript_count++;
		}
 
		
		rowid++;		
	}


	printf("   Total lines: %ld\n", rowid );
	fprintf(logs, "   Total lines: %ld\n", rowid );
	printf("long_transcript_count: %ld\n\n", long_transcript_count );
	fprintf(logs, "long_transcript_count: %ld\n\n", long_transcript_count );
	
	fflush(logs);
	fclose(fid_origin);
	fclose(file_new); 
	free(buff_str_desc);/* free */
	free(buff_str_seq);/* free */
	free(length);/* free */
	free(ismax);/* free */
	for (i=0;i<rowid;i++)
	{
		free(gene_name[i]);/* free */
	}
	free(gene_name);/* free */
	if (REMOVE_TEMP_FILES) remove(file_name_origin);
}

void copy_file(char file_name_origin[],char file_name_new[])
{

	FILE *fid_origin, *file_new;
	long rowid=0;

	char *buff_data;
	buff_data=(char*)malloc(MAX_LINE_LENGTH*sizeof(char));
 
	file_new=fopen(file_name_new,"w");
	if(file_new==NULL )
	{
		printf("ERROR: Cannot open %s\n",file_name_new);
		return;
	} 
	fid_origin=fopen(file_name_origin,"r"); 
	if(fid_origin==NULL)
	{
		printf("ERROR: Cannot open %s\n",file_name_origin);
		return;
	}

	while(!feof(fid_origin))
	{
		if(fgets(buff_data,MAX_LINE_LENGTH,fid_origin)==NULL)
		{
			break;		 
		}	
		fprintf(file_new,"%s",buff_data); 
		rowid++;
		if(rowid%10000==0)
				printf("copying file, row id: %ld, %s\n", rowid,file_name_origin);	 		
		 
	}
	printf("   Total lines: %ld\n\n", rowid );
	fclose(fid_origin);
	fclose(file_new); 
	free(buff_data);/* free */
	if (REMOVE_TEMP_FILES) remove(file_name_origin);
}

void samples_add_label(char file_name_origin[],char file_name_new[],long label)
{
	
	FILE *fid_origin, *file_new;
	long rowid=0;
	
	char *buff_data;
	buff_data=(char*)malloc(MAX_LINE_LENGTH*sizeof(char));
	
	file_new=fopen(file_name_new,"w");
	if(file_new==NULL )
	{
		printf("ERROR: Cannot open %s\n",file_name_new);
		return;
	} 
	fid_origin=fopen(file_name_origin,"r"); 
	if(fid_origin==NULL)
	{
		printf("ERROR: Cannot open %s\n",file_name_origin);
		return;
	}
	
	while(!feof(fid_origin))
	{
		if(fgets(buff_data,MAX_LINE_LENGTH,fid_origin)==NULL)
		{
			break;		 
		}
		fprintf(file_new,"%ld ",label); 
		fprintf(file_new,"%s",buff_data); 
		rowid++;
		if(rowid%10000==0)
			printf("adding label, row id: %ld, %s\n", rowid,file_name_origin);	 		
		
	}
	printf("   Total lines: %ld\n\n", rowid );
	fclose(fid_origin);
	fclose(file_new); 
	free(buff_data);/* free */
	if (REMOVE_TEMP_FILES) remove(file_name_origin);
}


void chars_to_upper_case(char file_name_origin[],char file_name_new[],
						 FILE * logs, long outmsg, long rmtempfile)
{

	FILE *fid_origin, *file_new;
	long rowid=0,i;
 	char *p;
	char *buff_data;

	buff_data=(char*)malloc(MAX_LINE_LENGTH*sizeof(char));
 
	file_new=fopen(file_name_new,"w");
	if(file_new==NULL )
	{
		printf("ERROR: Cannot open %s\n",file_name_new);
		fprintf(logs,"ERROR: Cannot open %s\n",file_name_new);
		fflush(logs);
		return;
	} 
	fid_origin=fopen(file_name_origin,"r"); 
	if(fid_origin==NULL)
	{
		printf("ERROR: Cannot open %s\n",file_name_origin);
		fprintf(logs, "ERROR: Cannot open %s\n",file_name_origin);
		fflush(logs);
		return;
	}	

	if(outmsg) printf("Nucleotides to upper case ... \n");
	fprintf(logs,"Nucleotides to upper case ... \n");
	fflush(logs);

	rowid=0;
	while(!feof(fid_origin))
	{
		/* read two lines  */
		if(fgets(buff_data,MAX_LINE_LENGTH,fid_origin)==NULL)
		{
			//printf("read failed!\n");	
			break;		 
		}
		if(buff_data[0]==10 || buff_data[0]==13 )
		{
			break;	
		}
		rowid++;
		if(rowid%10000==0)
		{
			if(outmsg) printf("   %s, row id: %ld\n", file_name_origin,rowid);
			fprintf(logs,"   %s, row id: %ld\n", file_name_origin,rowid);
			fflush(logs);
		}
		
		p=buff_data;
		i=0;
		while(*p)
		{
			buff_data[i]=toupper(*p);
			i++;
			p++;
		}
		
		fprintf(file_new,"%s",buff_data); 
	 		
		 
	}
	if(outmsg) printf("   Total lines: %ld\n\n", rowid );
	fprintf(logs,"   Total lines: %ld\n\n", rowid );
	fflush(logs);
	fclose(fid_origin);
	fclose(file_new); 
	free(buff_data); /* free */
	if (rmtempfile) remove(file_name_origin);
}


long fasta_file_line_count(char file_name[],FILE * logs)
{
 
 	char *buff_data;
	FILE *fid_1;
	long row=0;

	buff_data=(char*)malloc(MAX_LINE_LENGTH*sizeof(char));
 
	fid_1=fopen(file_name,"r"); 
	if(fid_1==NULL)
	{
		printf("ERROR: Cannot open %s\n",file_name);
		fprintf(logs,"ERROR: Cannot open %s\n",file_name);
		fflush(logs);
		return -1;
	}

	while(!feof(fid_1))
	{
		if(fgets(buff_data,MAX_LINE_LENGTH,fid_1)==NULL)
		{
			//puts("read failed!\n");	
				break;		 
		}
		row++;
	}
	fclose(fid_1);
	free(buff_data);/* free */
	return row;
}


long fasta_file_string_count(char file_name[],char strs[])
{
 
	FILE *fid_1;
	long row=0; 
	char *buff_data;
	buff_data=(char*)malloc(MAX_LINE_LENGTH*sizeof(char));
	fid_1=fopen(file_name,"r"); 
	if(fid_1==NULL)
	{
		printf("ERROR: Cannot open %s\n",file_name);
		return -1;
	}

	while(!feof(fid_1))
	{
		if(fgets(buff_data,MAX_LINE_LENGTH,fid_1)==NULL)
		{
			//puts("read a line from file failed !\n");	
			break;		 
		}
		if(strstr(buff_data,strs))
		{
			row++;
		}

	}
	fclose(fid_1);
	free(buff_data); /* free */
	return row;
}

long fasta_file_string_count2(char file_name[],char strs[],char strs2[])
{
	/* fasta file, strs count in the file_name */
 
		char *buff_data;
	long row=0;
	FILE *fid_1;
 	buff_data=(char*)malloc(MAX_LINE_LENGTH*sizeof(char));
	fid_1=fopen(file_name,"r"); 

	if(fid_1==NULL)
	{
		printf("ERROR: Cannot open %s\n",file_name);
		return -1;
	}

	while(!feof(fid_1))
	{
		if(fgets(buff_data,MAX_LINE_LENGTH,fid_1)==NULL)
		{
			//puts("read a line from file failed !\n");	
			break;		 
		}
		if(strstr(buff_data,strs) && strstr(buff_data,strs2))
		{
			row++;
		}

	}
	fclose(fid_1);
	free(buff_data);/* free */
	return row;
}





void find_short_sequence(char  file_name[])
{

	long i,min=100000000,x=0;
	char *buff_data;
	FILE *fid2,*fid;
	fid=fopen(file_name,"r");
	if (fid==NULL)
	{
		printf("ERROR: Cannot open file");
	}
	fid2=fopen("short_cDNA","w");
	if (fid2==NULL)
	{
		printf("ERROR: Cannot open file");
	}

	buff_data=(char*)malloc(MAX_LINE_LENGTH*sizeof(char));

	while(!feof(fid))
	{
		if(fgets(buff_data,MAX_LINE_LENGTH,fid))
		{
			x++;
			i=strlen(buff_data);
			if(i<min)
				min=i;
			if(i<50)
			{	
				//printf("========================================================================");
				printf("id=%ld ;  len=%ld  ;%s\n",x,i,buff_data);
			 	fprintf(fid2,"%s",buff_data);	
			}
			
		}
	}
	printf("min=%ld  \n",min);

	fclose(fid);
	fclose(fid2);
	free(buff_data);/* free */

}

void file_row_count(char  file_name[])
{
	/*   */
	long i=0,rowcount=0;
	FILE *fid;
	char *buff_data;
	buff_data=(char*)malloc(MAX_LINE_LENGTH*sizeof(char));
	fid=fopen(file_name,"r");
	if (fid==NULL)
	{
		printf("ERROR: Cannot open file");
	}
 
	
	while(!feof(fid))
	{
		if(fgets(buff_data,MAX_LINE_LENGTH,fid))
		{
		 
			i=strlen(buff_data);
			if(i>=2) rowcount++; 
			
		}
	}
	printf("rowcount=%ld  \n",rowcount);
	
	fclose(fid);
 	free(buff_data);/* free */
	
}
  

void randomlize(long *a, long n)
{
 /* permute the elements of a */	
	
	long i = 0,j = 0, k = 0;	
	srand(time(0));  

	for(i = 0; i < n; i++)		
	{		
		j = rand()%n;		
		k = a[i];		
		a[i] = a[j];		
		a[j] = k;		
	}	
}

void randomlize_evennumber(long *a, long n)
{

	
	long i = 0,j = 0, k = 0,t;	
	long *newa;
	srand(time(0)); 

	for(i = 0; i < n; i++)		
	{		
		j = rand()%n;		
		k = a[i];		
		a[i] = a[j];		
		a[j] = k;		
	}

	newa=(long *)malloc(n*sizeof(long));
	t=0;
	for(i = 0; i < n; i++)		
	{
		if(a[i]%2==0)
		{
			newa[t]=a[i];
			t++;
		}
	}
 	for(i = 0; i < n; i++)		
	{
		a[i]=newa[i];
	}
	free(newa);
}


void SelectionSort(long *a, long n)
{

    long i, j, index, value;
	
    for (i = 0; i < n - 1; i ++) 
	{
        index = i;
        value = a[i];
        for (j = i + 1; j < n; j ++)
            if (value > a[j]) 
			{
                index = j;
                value = a[j];
            }
		a[index] = a[i];
		a[i] = value;
    }
}


void random_select_rows(long total_row_count,long new_count,char fname[],char new_fname[],FILE * logs)
{
	
	FILE *fid_fname,*fid_new_file;
	long *all_ids;
	long i,j,k;
	char *buff_data; /* read line from sequence file */
	all_ids=(long *)malloc(total_row_count*sizeof(long));
	buff_data=(char*)malloc(MAX_LINE_LENGTH*sizeof(char));
	for (i=0;i<total_row_count;i++)
	{
		all_ids[i]=i; /* 0 -- max_sample_size-1 */
	}
	randomlize(all_ids,total_row_count);
	SelectionSort(all_ids,new_count); /* sorting */
	
	
	fid_fname=fopen(fname,"r"); 
	if(fid_fname==NULL )
	{
		printf("ERROR: Cannot open file\n");
		return;
	}
	fid_new_file=fopen(new_fname,"w"); 
	if(fid_new_file==NULL )
	{
		printf("ERROR: Cannot open file\n");
		return;
	}
	
	k=0;
	for(j=0;j<total_row_count;j++)
	{
		if(j%10000==0)
		{
			printf("selecting rows %s, row index: %ld\n", fname,j);
			fprintf(logs,"selecting rows %s, row index: %ld\n", fname,j);
			fflush(logs);
		 
		}
		if (j==all_ids[k])
		{ 	
			fgets(buff_data,MAX_LINE_LENGTH,fid_fname);	
			fprintf(fid_new_file,"%s",buff_data); 
			k++;
			if (k>=new_count)
			{
				break;
			}
		}
		else
		{
			fgets(buff_data,MAX_LINE_LENGTH,fid_fname);	
		}
	} 
	
	fclose(fid_fname);
	fclose(fid_new_file);
	free(all_ids);
	free(buff_data);
}


void random_select_rows_with_desc(long total_row_count,long new_count,char fname[],
						char fname_desc[],char new_fname[],char new_fname_desc[],
						FILE *logs, long outmsg, long rmtempfile)
{

	FILE *fid_fname,*fid_new_file,*fid_fname_desc,*fid_new_desc;
	long *all_ids;
	long i,j,k;
	char *buff_data; /* read line from sequence file */
	all_ids=(long *)malloc(total_row_count*sizeof(long));
	buff_data=(char*)malloc(MAX_LINE_LENGTH*sizeof(char));
	for (i=0;i<total_row_count;i++)
	{
		all_ids[i]=i; /* 0 -- max_sample_size-1 */
	}
	randomlize(all_ids,total_row_count);
	SelectionSort(all_ids,new_count); /* sorting */
	
	
	fid_fname=fopen(fname,"r"); 
	if(fid_fname==NULL )
	{
		printf("ERROR: Cannot open file.\n");
		return;
	}
	fid_fname_desc=fopen(fname_desc,"r"); 
	if(fid_fname_desc==NULL )
	{
		printf("ERROR: Cannot open file.\n");
		return;
	}


	fid_new_file=fopen(new_fname,"w"); 
	if(fid_new_file==NULL )
	{
		printf("ERROR: Cannot open file.\n");
		return;
	}
	fid_new_desc=fopen(new_fname_desc,"w"); 
	if(fid_new_desc==NULL )
	{
		printf("ERROR: Cannot open file.\n");
		return;
	}

	if(outmsg) printf("\nRandomly select sequences from file ... \n");
	fprintf(logs,"\nRandomly select sequences from file ... \n");
	fflush(logs);
	
	k=0;
	for(j=0;j<total_row_count;j++)
	{
		if(j%10000==0)
		{
			 if(outmsg) printf("   %s, row index: %ld\n", fname,j);
			 fprintf(logs,"   %s, row index: %ld\n", fname,j);
		}
		if (j==all_ids[k])
		{ 	
			fgets(buff_data,MAX_LINE_LENGTH,fid_fname);	
		    fprintf(fid_new_file,"%s",buff_data); 

			fgets(buff_data,MAX_LINE_LENGTH,fid_fname_desc);	
		    fprintf(fid_new_desc,"%s",buff_data); 

			k++;
			if (k>=new_count)
			{
				break;
			}
		}
		else
		{
			fgets(buff_data,MAX_LINE_LENGTH,fid_fname);	
			fgets(buff_data,MAX_LINE_LENGTH,fid_fname_desc);
		}
	} 
	
	fclose(fid_fname);
	fclose(fid_new_file);
	fclose(fid_fname_desc);	
	fclose(fid_new_desc);

	free(all_ids);
	free(buff_data);
	if(rmtempfile)
	{
		remove(fname);
		remove(fname_desc);
	}
}

void random_select_rows_train_test(long total_row_count,long train_count,char fname[],char train_file[],char test_file[])
{

	FILE *fid_fname,*fid_train_file,*fid_test_file;
	long *all_ids;
	long i,j,k;
	char *buff_data; /* read line from sequence file */
	all_ids=(long *)malloc(total_row_count*sizeof(long));
	buff_data=(char*)malloc(MAX_LINE_LENGTH*sizeof(char));
	for (i=0;i<total_row_count;i++)
	{
		all_ids[i]=i; /* 0 -- max_sample_size-1 */
	}
	randomlize(all_ids,total_row_count);
	SelectionSort(all_ids,train_count); /* sorting */
	
	
	fid_fname=fopen(fname,"r"); 
	if(fid_fname==NULL )
	{
		printf("ERROR: Cannot open file\n");
		return;
	}
	fid_train_file=fopen(train_file,"w"); 
	if(fid_train_file==NULL )
	{
		printf("ERROR: Cannot open file\n");
		return;
	}
	fid_test_file=fopen(test_file,"w"); 
	if(fid_test_file==NULL )
	{
		printf("ERROR: Cannot open file\n");
		return;
	}
	
	
	k=0;
	for(j=0;j<total_row_count;j++)
	{
		if(j%10000==0)
			 printf("selecting rows %s, row index: %ld\n", fname,j);
		if (j==all_ids[k])
		{ 	
			fgets(buff_data,MAX_LINE_LENGTH,fid_fname);	
		    fprintf(fid_train_file,"%s",buff_data); 
			k++;
			if (k>=train_count)
			{
				break;
			}
		}
		else
		{
			fgets(buff_data,MAX_LINE_LENGTH,fid_fname);
			 fprintf(fid_test_file,"%s",buff_data); 
		}
	} 
	
	fclose(fid_fname);
	fclose(fid_train_file);
	fclose(fid_test_file);	
	free(all_ids);
	free(buff_data);
}

void random_select_rows_twoadjacent(long total_row_count,long new_count,char fname[],char new_fname[])
{

	FILE *fid_fname,*fid_new_file;

	long *all_ids;
	long i,j,k;
	char *buff_data; /* read line from sequence file */
	all_ids=(long *)malloc(total_row_count*sizeof(long));
	buff_data=(char*)malloc(MAX_LINE_LENGTH*sizeof(char));
	for (i=0;i<total_row_count;i++)
	{
		all_ids[i]=i; /* 0 -- max_sample_size-1 */
	}
	randomlize_evennumber(all_ids,total_row_count); 
	SelectionSort(all_ids,new_count/2); /* sorting */  
	fid_fname=fopen(fname,"r"); 
	if(fid_fname==NULL )
	{
		printf("ERROR: Cannot open file\n");
		return;
	}
	fid_new_file=fopen(new_fname,"w"); 
	if(fid_new_file==NULL )
	{
		printf("ERROR: Cannot open file\n");
		return;
	}
	
	k=0;
	for(j=0;j<total_row_count;j++)
	{
		if(j%10000==0)
			 printf("selecting rows %s, row index: %ld\n", fname,j);
		if (j==all_ids[k])
		{			
			if(fgets(buff_data,MAX_LINE_LENGTH,fid_fname))
			fprintf(fid_new_file,"%s",buff_data); 

			if(fgets(buff_data,MAX_LINE_LENGTH,fid_fname))
			fprintf(fid_new_file,"%s",buff_data); 
			 				 
					
			k++;
			if (k>=new_count/2)
			{
				break;
			}
		}
		else
		{
			fgets(buff_data,MAX_LINE_LENGTH,fid_fname);
			fgets(buff_data,MAX_LINE_LENGTH,fid_fname);
		}
		j++; /* j=j+2 */
	} 
	
	fclose(fid_fname);
	fclose(fid_new_file);
	free(all_ids);
	free(buff_data);
}

void random_select_rows_twoadjacent2(long total_row_count,long new_count,char fname[],char new_fname[],char new_fname2[])
{

	FILE *fid_fname,*fid_new_file,*fid_new_file2;
	
	long *all_ids;
	long i,j,k;
	char *buff_data; /* read line from sequence file */
	all_ids=(long *)malloc(total_row_count*sizeof(long));
	buff_data=(char*)malloc(MAX_LINE_LENGTH*sizeof(char));
	for (i=0;i<total_row_count;i++)
	{
		all_ids[i]=i; /* 0 -- max_sample_size-1 */
	}
	randomlize_evennumber(all_ids,total_row_count); 
	SelectionSort(all_ids,new_count/2); /* sorting */
	
	fid_fname=fopen(fname,"r"); 
	if(fid_fname==NULL )
	{
		printf("ERROR: Cannot open file\n");
		return;
	}
	fid_new_file=fopen(new_fname,"w"); 
	if(fid_new_file==NULL )
	{
		printf("ERROR: Cannot open file\n");
		return;
	}
	fid_new_file2=fopen(new_fname2,"w"); 
	if(fid_new_file2==NULL )
	{
		printf("ERROR: Cannot open file\n");
		return;
	}
	
	k=0;
	for(j=0;j<total_row_count;j++)
	{
		if(j%10000==0)
			printf("selecting rows %s, row index: %ld\n", fname,j);
		if (j==all_ids[k])
		{			
			if(fgets(buff_data,MAX_LINE_LENGTH,fid_fname))
				fprintf(fid_new_file,"%s",buff_data); 
			
			if(fgets(buff_data,MAX_LINE_LENGTH,fid_fname))
				fprintf(fid_new_file,"%s",buff_data); 
			
			
			k++;
			if (k>=new_count/2)
			{
				break;
			}
		}
		else
		{
			fgets(buff_data,MAX_LINE_LENGTH,fid_fname);
			fprintf(fid_new_file2,"%s",buff_data); 
			fgets(buff_data,MAX_LINE_LENGTH,fid_fname);
			fprintf(fid_new_file2,"%s",buff_data); 
		}
		j++; /* j=j+2 */
	} 
	
	fclose(fid_fname);
	fclose(fid_new_file);
	fclose(fid_new_file2);
	free(all_ids);
	free(buff_data);
}


void save_matrix_to_file(double **m,char new_file[],long rowcount,long columncount )
  {
	 
	  FILE   *file_new;  
	  long i,j;
	   
	  
	  file_new=fopen(new_file,"w");
	  if(file_new==NULL )
	  {
		  printf("ERROR: Cannot open %s\n",new_file);
		  return;
	  }
	  
	  for(i=0;i<rowcount;i++)
	  {
		for(j=0;j<columncount;j++)
		{
			fprintf(file_new,"%.12lf ",m[i][j]);	
		}
		fputc('\n',file_new);
	  }
	   
	  
	  
	  fclose(file_new); 
}

void save_1d_array_to_file(long *m,char new_file[], long columncount )
{

	FILE   *file_new;  
	long  j;
	   
	
	file_new=fopen(new_file,"w");
	if(file_new==NULL )
	{
		printf("ERROR: Cannot open %s\n",new_file);
		return;
	}
	
	for(j=0;j<columncount;j++)
	{
		fprintf(file_new,"%ld  ",*(m+j));	
	} 	
	
	fclose(file_new); 
}

void training_label_vector_forlibsvm(long plusfrom,long plusto,long minusfrom,long minusto,char filename[])
{
	  FILE   *file_new;  
	  long i;
	   
	  
	  file_new=fopen(filename,"w");
	  if(file_new==NULL )
	  {
		  printf("ERROR: Cannot open %s\n",filename);
		  return;
	  }
	  
	  for(i=0;i<plusto;i++)
	  {		 
		  fprintf(file_new,"%d\n",1);	
	  }
	  for(i=minusfrom-1;i<minusto;i++)
	  {		 
		  fprintf(file_new,"%d\n",2);	
	  }	  
	  
	  fclose(file_new); 
}

void print_1d_array_longint(long a[],long length)
{
   long i;
   for (i=0;i<length;i++)
   {
	   printf("%ld ",a[i]);
   }
   printf("\n");
}

void print_1d_array_double(double a[],long length)
{
	long i;
	for (i=0;i<length;i++)
	{
		printf("%lf ",a[i]);
	}
	printf("\n");
}

void split_train_test(char fname[],char train_file[],char test_file[],long groups[],long times)
{
	/*    */
	FILE *fid_fname,*fid_train_file,*fid_test_file;
	long rowid=0;
	char *buff_data;  
 
	buff_data=(char*)malloc(MAX_LINE_LENGTH*sizeof(char));
 	
	fid_fname=fopen(fname,"r"); 
	if(fid_fname==NULL )
	{
		printf("ERROR: Cannot open file\n");
		return;
	}
	fid_train_file=fopen(train_file,"w"); 
	if(fid_train_file==NULL )
	{
		printf("ERROR: Cannot open file\n");
		return;
	}
	fid_test_file=fopen(test_file,"w"); 
	if(fid_test_file==NULL )
	{
		printf("ERROR: Cannot open file\n");
		return;
	}
	
	
	while(!feof(fid_fname))
	{
		
		fgets(buff_data,MAX_LINE_LENGTH,fid_fname);	
		if (times==groups[rowid])
		{ 			
			 fprintf(fid_test_file,"%s",buff_data); 
		}
		else
		{	
			fprintf(fid_train_file,"%s",buff_data); 				
		}
		rowid++;
		if(rowid%10000==0)
			printf("spliting train and test, %s, row index: %ld, times: %ld\n", fname,rowid,times);
	} 
	
	fclose(fid_fname);
	fclose(fid_train_file);
	fclose(fid_test_file);	
	free(buff_data);
	//if (REMOVE_TEMP_FILES) remove(train_file);
}


void compare_numerical_file(char file1[],char file2[],char file_diff[])
{

	FILE  *fid_file_diff,*fid_file1,*fid_file2;
	double x,y;
	long i=0,diffcount=0; 
 	char *buff_data;

	fid_file1=fopen(file1,"r");
	if(fid_file1==NULL)
	{
		printf("ERROR: Cannot open %s\n",file1);
		return  ;
	}
	
	fid_file2=fopen(file2,"r");
	if(fid_file2==NULL)
	{
		printf("ERROR: Cannot open %s\n",file2);
		return  ;
	}
	fid_file_diff=fopen(file_diff,"w");
	if(fid_file_diff==NULL)
	{
		printf("ERROR: Cannot open %s\n",file_diff);
		return  ;
	}
	 

	i=0;
	buff_data=(char*)malloc(MAX_LINE_LENGTH*sizeof(char));
	while(!feof(fid_file1)  && !feof(fid_file2))
	{
		double t;
		fscanf(fid_file1,"%lf",&x);
		fscanf(fid_file2,"%lf",&y);
		t=fabs(x-y);
		if (t>0.000001)
		{
			diffcount++;
			fprintf(fid_file_diff,"%ld:%lf,%lf ",diffcount,x,y);

		}
		i++;
	 
		
	}

	fclose(fid_file1);
	fclose(fid_file2);
	fclose(fid_file_diff);
 
	free(buff_data);
 

}

char * mystrcat(char left[],char right[],char newstring[])
{ 
	strcpy(newstring,left);
	strcat(newstring,right);
	return newstring;
}

void get_longest_transcript_lncrna(char origin_file[],char new_file[],long type, FILE * logs)
{
 
	switch(type)
	{
		case 1: // 
			 
			break;
		case 2: //gencode
			get_longest_transcript(origin_file,new_file,"|",1,logs); // for genecode
			break;
		case 3: //copy,  noncode_lncrna
			printf("longest transcript ... \ncopying ... \n");
			fprintf(logs,"longest transcript ... \ncopying ... \n");
			fflush(logs);		
			copy_file(origin_file,new_file);
			break;
	}


}

void get_longest_transcript_mrna(char origin_file[],char new_file[],long type, FILE * logs)
{
	switch(type)
	{
		case 1: //ensembl,TAIR,  mouse, human, arabidopsis
			get_longest_transcript(origin_file,new_file," :",10,logs);
			break;
		case 2: // refseq
			get_longest_transcript(origin_file,new_file," |",3,logs); //>gi|52856423|ref|NM_001005338.1| Homo sapiens olfactory receptor, family 5, subfamily H, member 1 (OR5H1), mRNA
			break;
		case 3: //copy
			printf("longest transcript ... \ncopying ... \n");
			fprintf(logs,"longest transcript ... \ncopying ... \n");
			fflush(logs);
			copy_file(origin_file,new_file);
				break;
	 
	}
	
}
 

void filter_samples_lncrna(char origin_file[],char new_file[],long type, FILE * logs)
{
	
	switch(type)
	{
		case 1: //noncode_lncrna
			 
			break;
		case 2: //gencode_lncrna
			 
			break;
		case 3: //copy
			printf("filter_samples_lncrna ... \n");
			fprintf(logs,"filter_samples_lncrna ... \n");
			fflush(logs);
			copy_file(origin_file,new_file);
				break;
	}
	
	
}

void filter_samples_mrna(char origin_file[],char new_file[],long type, FILE * logs,char out_file[])
{
	char temp1[200],temp2[200];



	switch(type)
	{
		case 1: //ensembl_mrna, mouse, human
			extract_seq_with_special_str(origin_file
				,mystrcat(out_file,"neg_filt_seq1",temp2),"transcript_biotype:protein_coding",logs);
			
			extract_seq_with_special_str(mystrcat(out_file,"neg_filt_seq1",temp1)
				,mystrcat(out_file,"neg_filt_seq2",temp2),"gene_biotype:protein_coding",logs);
			
			extract_seq_with_special_str(mystrcat(out_file,"neg_filt_seq2",temp1)
				,new_file,"cdna:known",logs);
			break;
		case 2: // ensembl_TAIR_Arabidopsis
			extract_seq_with_special_str(origin_file
				,mystrcat(out_file,"neg_filt_seq1",temp2),"protein",logs);
			
			extract_seq_with_special_seqstr(mystrcat(out_file,"neg_filt_seq1",temp1)
				,new_file,"ATG",logs);

			break;
		case 3: //copy		 
			printf("filter_samples_mrna\ncopying ... \n");
			fprintf(logs,"filter_samples_mrna\ncopying ... \n");
			fflush(logs);
			copy_file(origin_file,new_file);
				break;

		case 4: // refseq
			extract_seq_with_special_str(origin_file
				,new_file,"mRNA",logs);
			
			//extract_seq_with_special_seqstr(mystrcat(out_file,"neg_filt_seq1",temp1)
			//	,new_file,"ATG",logs);	//2600 
			
			break;
	}
}
 
 
 
