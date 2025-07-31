///////////////////////////////////////////////////////////////////
//                                                               
//  PLEK - predictor of lncRNAs and mRNAs based on k-mer scheme  
//  Authors: Aimin Li, Junying Zhang                            
//  Contacts: LiAiminMail@gmail.com, jyzhang@mail.xidian.edu.cn  
//  Webcite: https://sourceforge.net/projects/plek/                      
//  Updated on: June 7, 2014                                                  
//                                                               
///////////////////////////////////////////////////////////////////

/*  complie on linux: # LANG=C gcc -g -Wall PLEK_main.c -o PLEK -lm  */ 
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include "PLEK_kmer.h" 
#include "PLEK_fastafile.h"

#ifndef _MAX_SAMPLE_COUNT
#define _MAX_SAMPLE_COUNT 2000000
#endif

int main(int argc,char *argv[])
{	
	long mer_count=5, steplength_is_1=1, denominator_type=5;  
	char positive_fasta[500]="", negative_fasta[500]="", prefix[500]="_";  

	long mustbebalanced=0, pos_class_label=1, neg_class_label=0, onlyk=0, is_known_classlabel=0; 
	char allsvm[]="allsvm", allsvmdesc[]="allsvmdesc"; 
	long pos_seq_count=0,  neg_seq_count=0,  balanced_count=0, i, outmsg=0, rmtempfile=1 ;
 
	char temp1[500], temp2[500], temp3[500], temp4[500], temp5[500], temp6[500];
	FILE *logs,  *fid_allsvm, *fid_allsvmdesc, *fid_positive_fasta, *fid_negative_fasta; 
	long max_seq_length=0, pos_max_seq_length=0, neg_max_seq_length=0, min_seq_len=200;
	
	
	clock_t start_time, end_time;  
	double  duration_time;  
	time_t rawtime;
	struct tm * timeinfo;
    start_time=clock();  
	
	if(argc>1)
	{
		for (i=0;i<argc;i++)
		{
			if (strcmp(argv[i],"-k")==0)
			{
				mer_count=atol(argv[i+1]);
			}
			if (strcmp(argv[i],"-s")==0)
			{
				steplength_is_1=atol(argv[i+1]);
			}
			if (strcmp(argv[i],"-d")==0)
			{
				denominator_type=atol(argv[i+1]);
			}
			if (strcmp(argv[i],"-p")==0)
			{
				strcpy(positive_fasta,argv[i+1]);
			}
			if (strcmp(argv[i],"-f")==0)
			{
				strcpy(positive_fasta,argv[i+1]);
			}
			if (strcmp(argv[i],"-n")==0)
			{
				strcpy(negative_fasta,argv[i+1]);
			}
			if (strcmp(argv[i],"-o")==0)
			{
				strcpy(prefix,argv[i+1]);
				strcat(prefix,"_");
			}
			if (strcmp(argv[i],"-b")==0)
			{
				mustbebalanced=1; //NOTE: -b , need to balance;  NO -b, not balance
			}
			if (strcmp(argv[i],"-l")==0)
			{
				min_seq_len=atol(argv[i+1]);
			}
			if (strcmp(argv[i],"-isoutmsg")==0)
			{
				outmsg=atol(argv[i+1]);
			}
			if (strcmp(argv[i],"-isrmtempfile")==0)
			{
				rmtempfile=atol(argv[i+1]);
			}
		 
		}
	}
	else
	{	
		printf("usage on linux/unix: ./PLEK -p positive_samples -n negative_samples -o output_file_prefix\n\n");
		printf("usage on linux/unix: ./PLEK -p positive_samples -n negative_samples -o output_file_prefix -b \n\n");
		printf("usage on linux/unix: ./PLEK -p positive_samples -o output_file_prefix\n\n");
		printf("usage on linux/unix: ./PLEK -n negative_samples -o output_file_prefix\n\n");
		printf("usage on linux/unix: ./PLEK -f samples -o output_file_prefix\n\n");
		return 1;
	}

	

	if (strlen(positive_fasta)>0)
	{ 
		fid_positive_fasta=fopen(positive_fasta,"r");
		if( fid_positive_fasta==NULL)
		{
			 printf("WARNING: CANNOT OPEN FILE %s\n",positive_fasta);  
			 return -1;
		}
		fclose(fid_positive_fasta);
	}
	if (strlen(negative_fasta)>0)
	{
		fid_negative_fasta=fopen(negative_fasta,"r");
		if( fid_negative_fasta==NULL)
		{
			printf("WARNING: CANNOT OPEN FILE %s\n",negative_fasta);  
			return -1;
		}
		fclose(fid_negative_fasta);
	} 
	
	if(strlen(positive_fasta)>0 && strlen(negative_fasta)>0)
		is_known_classlabel=1;
	 
 
	logs=fopen(mystrcat(prefix,"logs",temp1),"w");
	if(logs==NULL)
	{
		 logs=fopen(mystrcat(prefix,"_logs",temp1),"w");
	}
	time(&rawtime);timeinfo=localtime(&rawtime);
	 
	fprintf(logs, "%s\n", "==========================================================="); 
	fprintf(logs, "The current date/time is: %s\n", asctime (timeinfo) ); 
	fflush(logs);


	initialize_kmer(mer_count);

	fid_allsvm=fopen(mystrcat(prefix,allsvm,temp1),"w"); 
	if( fid_allsvm==NULL)
	{
		printf("cannot open file: %s\n",mystrcat(prefix,allsvm,temp1));
		return -1;
	}
	fclose(fid_allsvm);	

	fid_allsvmdesc=fopen(mystrcat(prefix,allsvmdesc,temp1),"w"); 
	if( fid_allsvmdesc==NULL)
	{
		printf("cannot open file: %s\n",mystrcat(prefix,allsvmdesc,temp1));
		return -1;
	}
	fclose(fid_allsvmdesc);	


	if (strlen(positive_fasta)>0)
	{ 
		/* positive  */       
		to_one_line(positive_fasta,mystrcat(prefix,"pos_1row",temp1),logs,outmsg);

		check_two_line_fasta_file(mystrcat(prefix,"pos_1row",temp1),logs, outmsg);
		
		pos_max_seq_length=	remove_short_sequence(mystrcat(prefix,"pos_1row",temp1)
			,mystrcat(prefix,"pos_1row_200",temp2),min_seq_len,logs, outmsg, rmtempfile,
			is_known_classlabel); //==========		
			
		split_seq_and_desc(mystrcat(prefix,"pos_1row_200",temp1),
			mystrcat(prefix,"pos_seq_data",temp2), 
			mystrcat(prefix,"pos_seq_desc",temp3) ,logs, outmsg, rmtempfile);

		chars_to_upper_case(mystrcat(prefix,"pos_seq_data",temp1)
			,mystrcat(prefix,"pos_seq_upcase",temp2),logs, outmsg, rmtempfile);
 
		pos_seq_count=fasta_file_line_count(mystrcat(prefix,"pos_seq_upcase",temp1),logs);

		if (pos_seq_count<=0)
		{
			printf("WARNING: No valid positive sequences.\n");
			fprintf(logs, "WARNING: No valid positive sequences.\n");
			fflush(logs);
			return -1;
		}	
	}
	

	if (strlen(negative_fasta)>0)
	{
		/* negative  */
		to_one_line(negative_fasta,mystrcat(prefix,"neg_1row",temp1),logs,outmsg);

		check_two_line_fasta_file(mystrcat(prefix,"neg_1row",temp1),logs, outmsg);

		neg_max_seq_length=remove_short_sequence(mystrcat(prefix,"neg_1row",temp1)
			,mystrcat(prefix,"neg_1row_200",temp2), min_seq_len,logs, outmsg, rmtempfile,
			is_known_classlabel); //==========

		split_seq_and_desc(mystrcat(prefix,"neg_1row_200",temp1),
			mystrcat(prefix,"neg_seq_data",temp2), 
			mystrcat(prefix,"neg_seq_desc",temp3) ,logs, outmsg, rmtempfile);	

		chars_to_upper_case(mystrcat(prefix,"neg_seq_data",temp1) 
			,mystrcat(prefix,"neg_seq_upcase",temp2),logs, outmsg, rmtempfile);
		 
		neg_seq_count=fasta_file_line_count(mystrcat(prefix,"neg_seq_upcase",temp1),logs);

		if (neg_seq_count<=0)
		{
			printf("WARNING: No valid negative sequences.\n");
			fprintf(logs, "WARNING: No valid negative sequences.\n");
			fflush(logs);
			return -1;
		}
	}

	if (mustbebalanced)
	{
		balanced_count=pos_seq_count>neg_seq_count?neg_seq_count:pos_seq_count;
		max_seq_length=pos_max_seq_length>neg_max_seq_length?pos_max_seq_length:neg_max_seq_length;
		if(outmsg) printf("Balanced sample number: %ld\n",balanced_count);
		fprintf(logs,"Balanced sample number: %ld\n",balanced_count);
		fflush(logs);
	}

 
	if (strlen(positive_fasta)>0)
	{ 
		/* positive , first*/
		if (mustbebalanced && balanced_count!=pos_seq_count)
		{ 
			random_select_rows_with_desc(pos_seq_count,balanced_count,
				mystrcat(prefix,"pos_seq_upcase",temp1),
				mystrcat(prefix,"pos_seq_desc",temp2),
				mystrcat(prefix,"pos_seq_balanced",temp3),
				mystrcat(prefix,"pos_seq_desc_new",temp4),logs, outmsg, rmtempfile);

			frequency_label_svm(mystrcat(prefix,"pos_seq_balanced",temp1),
				mystrcat(prefix,"pos_svm",temp2),
	 			mer_count,onlyk,pos_class_label ,balanced_count,steplength_is_1,
				mystrcat(prefix,"pos_seq_desc_new",temp4),
				mystrcat(prefix,allsvm,temp5),denominator_type,prefix,logs,
				mustbebalanced,mystrcat(prefix,allsvmdesc,temp6), outmsg, rmtempfile); 
		
		}
		else
		{
 			frequency_label_svm(mystrcat(prefix,"pos_seq_upcase",temp1),
				mystrcat(prefix,"pos_svm",temp2),
	 			mer_count,onlyk,pos_class_label ,balanced_count,steplength_is_1,
				mystrcat(prefix,"pos_seq_desc",temp4),
				mystrcat(prefix,allsvm,temp5),denominator_type,prefix,logs,
				mustbebalanced,mystrcat(prefix,allsvmdesc,temp6), outmsg, rmtempfile); 
		 
		}	
	}

	if (strlen(negative_fasta)>0)
	{

		/* negative ,second */
		if (mustbebalanced && balanced_count!=neg_seq_count)
		{
			random_select_rows_with_desc(neg_seq_count,balanced_count,
				mystrcat(prefix,"neg_seq_upcase",temp1),
				mystrcat(prefix,"neg_seq_desc",temp2),
				mystrcat(prefix,"neg_seq_balanced",temp3),
				mystrcat(prefix,"neg_seq_desc_new",temp4),logs, outmsg, rmtempfile);


			frequency_label_svm(mystrcat(prefix,"neg_seq_balanced",temp1),
				mystrcat(prefix,"neg_svm",temp2), 
	 			mer_count,onlyk,neg_class_label, balanced_count,steplength_is_1,
				mystrcat(prefix,"neg_seq_desc_new",temp4),
				mystrcat(prefix,allsvm,temp5),denominator_type,prefix,logs,
				mustbebalanced,mystrcat(prefix,allsvmdesc,temp6), outmsg, rmtempfile); 
			 
		}
		else
		{

 			frequency_label_svm(mystrcat(prefix,"neg_seq_upcase",temp1),
				mystrcat(prefix,"neg_svm",temp2), 
	 			mer_count,onlyk,neg_class_label, balanced_count,steplength_is_1,
				mystrcat(prefix,"neg_seq_desc",temp4),
				mystrcat(prefix,allsvm,temp5),denominator_type,prefix,logs,
				mustbebalanced,mystrcat(prefix,allsvmdesc,temp6), outmsg, rmtempfile); 
		 
		} 

	}

     
    
	free_mer(mer_count); 

	//if(outmsg) printf("The current date/time is: %s", asctime (timeinfo) ); 
 	
	end_time = clock();  
	duration_time = (double)(end_time-start_time) / CLOCKS_PER_SEC;  
	//if(outmsg) printf( "\007duration: %lf seconds\n", duration_time ); 
	

	fflush(logs);
	time(&rawtime);timeinfo=localtime(&rawtime);
	fprintf(logs," [date]: %s \n",asctime (timeinfo));
	fprintf(logs," steplength_is_1: %ld\n",steplength_is_1); 
	fprintf(logs," denominator_type: %ld\n",denominator_type); 	
	fprintf(logs," positive_fasta file: %s \n",positive_fasta);
	fprintf(logs," negative_fasta file: %s \n",negative_fasta); 
	fprintf(logs," max_seq_length: %ld\n",max_seq_length); 	
	fprintf(logs," pos_seq_count: %ld\n",pos_seq_count);
	fprintf(logs," neg_seq_count: %ld\n",neg_seq_count);
	fprintf(logs," balanced_count: %ld\n",balanced_count);
	fprintf(logs," mer_count: %ld\n",mer_count); 
	fprintf(logs," mustbebalanced: %ld\n",mustbebalanced); 
	fprintf(logs," min_seq_len: %ld\n",min_seq_len); 
	fprintf(logs," onlyk: %ld\n",onlyk);
	fprintf(logs, "duration: %lf seconds\n", duration_time ); 	
	fprintf(logs, "The current date/time is: %s", asctime (timeinfo) ); 
	fprintf(logs, "%s\n", "==========================================================="); 

	fclose(logs); 
	
	return 1;

	
}
