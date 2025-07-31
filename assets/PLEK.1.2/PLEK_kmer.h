///////////////////////////////////////////////////////////////////
//                                                               
//  PLEK - predictor of lncRNAs and mRNAs based on k-mer scheme  
//  Authors: Aimin Li, Junying Zhang                             
//  Contacts: LiAiminMail@gmail.com, jyzhang@mail.xidian.edu.cn  
//  Webcite: https://sourceforge.net/projects/plek/                      
//  Updated on: Mar 22 2014                                              
//                                                               
///////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <ctype.h>

#ifndef MAX_LINE_LENGTH  
	#define MAX_LINE_LENGTH 1000000 
#endif 

#ifndef MAX_SAMPLE_COUNT  
	#define MAX_SAMPLE_COUNT 2000000  
#endif 

#ifndef REMOVE_TEMP_FILES  
	#define REMOVE_TEMP_FILES 0  
#endif 


char mer_1[4];
char **mer_2,**mer_3,**mer_4,**mer_5,**mer_6;


long power_long(long x,long y)
{
	long i,r;
	if (y<0)
	{
		printf("power_long error.");
		return 1;
	}
	if (y==0)
	{
		return 1;
	}

	r=x;
	for (i=1;i<y;i++)
	{
		r=r*x;
	}
	return r;
}

double denominator_value(long kth,long len,long deno_type,long kmer_count)
{
	return (len-kth+1)*1.0*power_long(4,kmer_count-kth);
}
 
void initialize_kmer(long mer_count)
{

	long i,j,k,l,m,n  ;

	/*////////////////////////////////////////*/
	/* malloc memory */
	if (mer_count>=2)
	{
		mer_2 = (char **)malloc(sizeof(char *)*16); 
		for (j=0; j<16; j++)
		{
			mer_2[j] = (char *)malloc(sizeof(char)*2); 
		}
	}


	if (mer_count>=3)
	{
		mer_3 = (char **)malloc(sizeof(char *)*64); 
		for (j=0; j<64; j++)
		{
			mer_3[j] = (char *)malloc(sizeof(char)*3); 
		}
	}


	if (mer_count>=4)
	{
		mer_4 = (char **)malloc(sizeof(char *)*256); 
		for (j=0; j<256; j++)
		{
			mer_4[j] = (char *)malloc(sizeof(char)*4); 
		}
	}


	if (mer_count>=5)
	{
		mer_5 = (char **)malloc(sizeof(char *)*1024); 
		for (j=0; j<1024; j++)
		{
			mer_5[j] = (char *)malloc(sizeof(char)*5); 
		}
	}


	if (mer_count>=6)
	{
		mer_6 = (char **)malloc(sizeof(char *)*4096); 
		for (j=0; j<4096; j++)
		{
			mer_6[j] = (char *)malloc(sizeof(char)*6); 
		}
	}


	/*	////////////////////////////////////////////////// */
	/* initialize mer string */
	////////////////////////////////////////////////////////

	/* 1-mer */
	 mer_1[0]='A';
	 mer_1[1]='G';
	 mer_1[2]='C';
	 mer_1[3]='T';

	 if (mer_count>=2)
	 {
		 /* 2-mer */
		 for( i=0;i<4;++i)
		 {
			 for( j=0;j<4;++j)
			 {
				 long index=i*4+j;
				 mer_2[index][0]=mer_1[i];
				 mer_2[index][1]=mer_1[j];
			 }
		}
	}

	if (mer_count>=3)
	{
		/* 3-mer */
		for( i=0;i<4;++i)
		{
			for( j=0;j<4;++j)
			{
				for( k=0;k<4;++k)
				{
					long index=(i*4+j)*4+k;
					mer_3[index][0]=mer_1[i];
					mer_3[index][1]=mer_1[j];
					mer_3[index][2]=mer_1[k];
				}
			}
		}
	}

	
	if (mer_count>=4)
	{
		/* 4-mer */
		for(i=0;i<4;++i)
		{
			for( j=0;j<4;++j)
			{
				for( k=0;k<4;++k)
				{
					for( l=0;l<4;++l)
					{
						long index=((i*4+j)*4+k)*4+l;
						mer_4[index][0]=mer_1[i];
						mer_4[index][1]=mer_1[j];
						mer_4[index][2]=mer_1[k];
						mer_4[index][3]=mer_1[l];
					}
				}
			}
		}
	}


	if (mer_count>=5)
	{
		/* 5-mer */
		for(i=0;i<4;++i)
		{
			for( j=0;j<4;++j)
			{
				for( k=0;k<4;++k)
				{
					for( l=0;l<4;++l)
					{
						for( m=0;m<4;++m)
						{
							long index=(((i*4+j)*4+k)*4+l)*4+m;
							mer_5[index][0]=mer_1[i];
							mer_5[index][1]=mer_1[j];
							mer_5[index][2]=mer_1[k];
							mer_5[index][3]=mer_1[l];
							mer_5[index][4]=mer_1[m];
						}
					}
				}
			}
		}
	}



	if (mer_count>=6)
	{
		/* 6-mer */
		for(i=0;i<4;++i)
		{
			for( j=0;j<4;++j)
			{
				for( k=0;k<4;++k)
				{
					for( l=0;l<4;++l)
					{
						for( m=0;m<4;++m)
						{
							for( n=0;n<4;++n)
							{
								long index=((((i*4+j)*4+k)*4+l)*4+m)*4+n;
								mer_6[index][0]=mer_1[i];
								mer_6[index][1]=mer_1[j];
								mer_6[index][2]=mer_1[k];
								mer_6[index][3]=mer_1[l];
								mer_6[index][4]=mer_1[m];
								mer_6[index][5]=mer_1[n];
							}
						}
					}
				}
			}
		}
	}
	

	



} /* end of initialize_kmer */

void free_mer(long mer_count)
{
	long j;
	/*  free mer */
	if (mer_count>=2)
	{
		for (j=0; j<16; j++)
		{
			free(mer_2[j]);
		}
		free(mer_2);
	}	

	if (mer_count>=3)
	{
		for (j=0; j<64; j++)
		{
			free(mer_3[j]);
		}
		free(mer_3);
	}
	

	if (mer_count>=4)
	{
		for (j=0; j<256; j++)
		{
			free(mer_4[j]);
		}
		free(mer_4);
	}
	

	if (mer_count>=5)
	{
		for (j=0; j<1024; j++)
		{
			free(mer_5[j]);
		}
		free(mer_5);
	}
	

	if (mer_count>=6)
	{
		for (j=0; j<4096; j++)
		{
			free(mer_6[j]);
		}
		free(mer_6);
	}

} /* end of free_mer */


void write_kmer_to_file(long mer_count)
{

	//////////////////////////////////////////////////////////////	
	/*	write mer_string to mer_string 					*/
	//////////////////////////////////////////////////////////////

	long i ;

	FILE *file_id_mer_string=fopen("mer_string","w");
	if(file_id_mer_string==NULL)
	{
		puts("write failed\n");
		exit(1);
	}
	/* 1-mer */
	for(i=0;i<4;++i)
	{
		fprintf(file_id_mer_string,"%c ",mer_1[i]);
	}
	fputc('\n',file_id_mer_string);
	if (mer_count>=2)
	{
		/* 2-mer */
		for(i=0;i<16;++i)
		{
			fprintf(file_id_mer_string,"%c",mer_2[i][0]);
			fprintf(file_id_mer_string,"%c",mer_2[i][1]);
			fputc(' ',file_id_mer_string);
		}
	fputs("\n",file_id_mer_string);
	}


	if (mer_count>=3)
	{
		/* 3-mer */
		for(i=0;i<64;++i)
		{
			fprintf(file_id_mer_string,"%c",mer_3[i][0]);
			fprintf(file_id_mer_string,"%c",mer_3[i][1]);
			fprintf(file_id_mer_string,"%c",mer_3[i][2]);
			fputc(' ',file_id_mer_string);
		}
	fputs("\n",file_id_mer_string);	
	}


	if (mer_count>=4)
	{
		/* 4-mer */
		for(i=0;i<256;++i)
		{
			fprintf(file_id_mer_string,"%c",mer_4[i][0]);
			fprintf(file_id_mer_string,"%c",mer_4[i][1]);
			fprintf(file_id_mer_string,"%c",mer_4[i][2]);
			fprintf(file_id_mer_string,"%c",mer_4[i][3]);
			fputc(' ',file_id_mer_string);
		}
	fputs("\n",file_id_mer_string);	
	}


	if (mer_count>=5)
	{
		/* 5-mer */
		for(i=0;i<1024;++i)
		{
			fprintf(file_id_mer_string,"%c",mer_5[i][0]);
			fprintf(file_id_mer_string,"%c",mer_5[i][1]);
			fprintf(file_id_mer_string,"%c",mer_5[i][2]);
			fprintf(file_id_mer_string,"%c",mer_5[i][3]);
			fprintf(file_id_mer_string,"%c",mer_5[i][4]);
			fputc(' ',file_id_mer_string);
		}
	fputs("\n",file_id_mer_string);
	}


	if (mer_count>=6)
	{
		/* 6-mer */
		for(i=0;i<4096;++i)
		{
			fprintf(file_id_mer_string,"%c",mer_6[i][0]);
			fprintf(file_id_mer_string,"%c",mer_6[i][1]);
			fprintf(file_id_mer_string,"%c",mer_6[i][2]);
			fprintf(file_id_mer_string,"%c",mer_6[i][3]);
			fprintf(file_id_mer_string,"%c",mer_6[i][4]);
			fprintf(file_id_mer_string,"%c",mer_6[i][5]);
			fputc(' ',file_id_mer_string);
		}
	fputs("\n",file_id_mer_string);
	}




	fclose(file_id_mer_string);

}/* end of write_kmer_to_file() */



void frequency_vertor(char sequence_only_file[],char vector_file_name[], 
					  char mean_frequency[],long mer_count, long onlyk)

{

	/* Input: a seq file.  */
	/* Output: k-mer used frequency file, a row corresponding to a seq, */

	/* onlyk:0 -- all(1 to k); 1 -- only 1; 2 -- only 2.... */

	 
	double *vector_data; /*  frequency vector */
	double *total_average; /* positive average, negative average */
	char *buff_data;

	long i,len; /* length of sequence */
	long row_id=0; /* row counter */
	long vector_dimen;
 
	double count_mer1[4]={0.000};  /* the count of A,T,C,G; 1-mer */
	double *count_mer2; /* the count of AA,AT,AC,AG...; 2-mer */
	double *count_mer3;
	double *count_mer4;
	double *count_mer5;
	double *count_mer6;
 

	long i1,i2,i3,i4,i5,i6 ;
	long begin_mer1=0,begin_mer2=0,begin_mer3=0,begin_mer4=0,begin_mer5=0,begin_mer6=0	 ;			
	long a1,a2,a3,a4,a5,a6;
	long p1,p2,p3,p4,p5,p6;
	long q2,q3,q4,q5,q6;
	FILE *fid_sequence_file,*fid_vector,*fid_mean;


	if(mer_count==1) vector_dimen=4;
	if(mer_count==2) vector_dimen=20;
	if(mer_count==3) vector_dimen=84;
	if(mer_count==4) vector_dimen=340;
	if(mer_count==5) vector_dimen=1364;
	if(mer_count==6) vector_dimen=5460;


	buff_data=(char*)malloc(MAX_LINE_LENGTH*sizeof(char)); 

	vector_data=(double*)malloc(vector_dimen*sizeof(double));
	for(i=0;i<vector_dimen;i++)
		vector_data[i]=0.0;
	total_average=(double*)malloc(vector_dimen*sizeof(double));
	for(i=0;i<vector_dimen;i++)
		total_average[i]=0.0;

	if (mer_count>=2)
	{
		count_mer2=(double*)malloc(16*sizeof(double));
	}
	if (mer_count>=3)
	{
		count_mer3=(double*)malloc(64*sizeof(double));
	}	
	if (mer_count>=4)
	{
		count_mer4=(double*)malloc(256*sizeof(double));
	}	
	if (mer_count>=5)
	{
		count_mer5=(double*)malloc(1024*sizeof(double));
	}	
	if (mer_count>=6)
	{
		count_mer6=(double*)malloc(4096*sizeof(double));
	}	


	fid_sequence_file=fopen(sequence_only_file,"r");
	fid_vector=fopen(vector_file_name,"w");
	fid_mean=fopen(mean_frequency,"w"); 
	if(fid_vector==NULL  )
	{
		printf("Cannot open %s\n",vector_file_name);
		return;
	}
	if(fid_sequence_file==NULL )
	{
		printf("Cannot open %s\n",sequence_only_file);
		return;
	}
	if( fid_mean==NULL)
	{
		printf("Cannot open %s\n",mean_frequency);
		return;
	}

	while(!feof(fid_sequence_file))
	{
		if(fgets(buff_data,MAX_LINE_LENGTH,fid_sequence_file))
		{
			/* initial to zero */
			if (mer_count>=1)
			{
				for (i=0;i<4;i++)
				{
					count_mer1[i]=0.0;
				}
			}

			if (mer_count>=2)
			{
				for (i=0;i<16;i++)
				{
					count_mer2[i]=0.0;
				}
			}
			
			if (mer_count>=3)
			{
				for (i=0;i<64;i++)
				{
					count_mer3[i]=0.0;
				}
			}
			
			if (mer_count>=4)
			{
				for (i=0;i<256;i++)
				{
					count_mer4[i]=0.0;
				}
			}				
			
			if (mer_count>=5)
			{
				for (i=0;i<1024;i++)
				{
					count_mer5[i]=0.0;
				}
			}				
			
			if (mer_count>=6)
			{
				for (i=0;i<4096;i++)
				{
					count_mer6[i]=0.0;
				}
			}				
			
			

			if(row_id%1000==0)
				printf("frequency, %s, row index: %ld\n", vector_file_name,row_id);

			row_id++;

			/* computation */

			len=strlen(buff_data);
			if(buff_data[len-1]==10) /* carriage return */
			{
				len=len-1;
			}
			if(buff_data[len-1]==13) /* line break */
			{
				len=len-1;
			}


		 
			if(mer_count>=1 && (onlyk==0 || onlyk==1))
			{
				//1mer  
				begin_mer1=0;
				while(begin_mer1<len)
				{
					for(i1=0;i1<4;++i1)
					{
						if(buff_data[begin_mer1]==mer_1[i1])
						{
							count_mer1[i1]=count_mer1[i1]+1;
						}
					}
					begin_mer1=begin_mer1+1;
				}
				for(a1=0;a1<4;++a1)
				{		
					// fprintf(fid_vector_count,"%c:%lf ",mer_1[a1],count_mer1[a1]);
					count_mer1[a1]=count_mer1[a1]/denominator_value(1,len,5,mer_count);
				}
			}

			if(mer_count>=2 && (onlyk==0 || onlyk==2))
			{
				//2mer 
				begin_mer2=0;
				while(begin_mer2<len-1)
				{
					for(i2=0;i2<16;++i2)
					{			
						if(buff_data[begin_mer2]==mer_2[i2][0] &&
							buff_data[begin_mer2+1]==mer_2[i2][1])
						{
							count_mer2[i2]=count_mer2[i2]+1;					
							break; 
						}			
					}
					begin_mer2=begin_mer2+2;
				}
				for( a2=0;a2<16;++a2)
				{				
					// fprintf(fid_vector_count,"%c%c:%lf ",mer_2[a2][0],mer_2[a2][1],count_mer2[a2]);
					count_mer2[a2]=count_mer2[a2]/denominator_value(2,len,5,mer_count);
				}
			}

			if(mer_count>=3 && (onlyk==0 || onlyk==3))
			{
				//3mer
				begin_mer3=0;
				while(begin_mer3<len-2)
				{
					for(i3=0;i3<64;++i3)
					{				
						if(buff_data[begin_mer3]==mer_3[i3][0] &&
							buff_data[begin_mer3+1]==mer_3[i3][1] &&
							buff_data[begin_mer3+2]==mer_3[i3][2])
						{
							count_mer3[i3]=count_mer3[i3]+1;
							break;
						}
					}
					begin_mer3=begin_mer3+3;
				}
				for( a3=0;a3<64;++a3)
				{		
					// fprintf(fid_vector_count,"%c%c%c:%lf ",mer_3[a3][0],mer_3[a3][1],mer_3[a3][2],count_mer3[a3]);
					count_mer3[a3]=count_mer3[a3]/denominator_value(3,len,5,mer_count);
				}	
			}

			if(mer_count>=4 && (onlyk==0 || onlyk==4))
			{
				//4mer
				begin_mer4=0;
				while(begin_mer4<len-3)
				{
					for(i4=0;i4<256;++i4)
					{			
						if(buff_data[begin_mer4]==mer_4[i4][0] &&
							buff_data[begin_mer4+1]==mer_4[i4][1] &&
							buff_data[begin_mer4+2]==mer_4[i4][2] &&
							buff_data[begin_mer4+3]==mer_4[i4][3] )
						{
							count_mer4[i4]+=1;								
							break;
						}
								
					}
					begin_mer4=begin_mer4+4;
				}

				//write vector_data to array
				for( a4=0;a4<256;++a4)
				{
					// fprintf(fid_vector_count,"%c%c%c%c:%lf ",mer_4[a4][0],mer_4[a4][1],mer_4[a4][2],mer_4[a4][3],count_mer4[a4]);
					count_mer4[a4]=count_mer4[a4]/denominator_value(4,len,5,mer_count);				
				}
			}


			if(mer_count>=5 && (onlyk==0 || onlyk==5))
			{
				//5mer 
				begin_mer5=0;
				while(begin_mer5<len-4)
				{
					for(i5=0;i5<1024;++i5)
					{			
						if(buff_data[begin_mer5]==mer_5[i5][0] &&
							buff_data[begin_mer5+1]==mer_5[i5][1] &&
							buff_data[begin_mer5+2]==mer_5[i5][2] &&
							buff_data[begin_mer5+3]==mer_5[i5][3] &&
							buff_data[begin_mer5+4]==mer_5[i5][4] )
						{
							count_mer5[i5]+=1;								
							break;
						}
								
					}
					begin_mer5=begin_mer5+5;
				}
				//write vector_data to array
				for( a5=0;a5<1024;++a5)
				{
					// fprintf(fid_vector_count,"%c%c%c%c%c:%lf ",mer_5[a5][0],mer_5[a5][1],mer_5[a5][2],mer_5[a5][3],mer_5[a5][4],count_mer5[a5]);
					count_mer5[a5]=count_mer5[a5]/denominator_value(5,len,5,mer_count);				
				}
			}

			if(mer_count>=6 && (onlyk==0 || onlyk==6))
			{
				//6mer 
				begin_mer6=0;
				while(begin_mer6<len-5)
				{
					for(i6=0;i6<4096;++i6)
					{			
						if(buff_data[begin_mer6]==mer_6[i6][0] &&
							buff_data[begin_mer6+1]==mer_6[i6][1] &&
							buff_data[begin_mer6+2]==mer_6[i6][2] &&
							buff_data[begin_mer6+3]==mer_6[i6][3] &&
							buff_data[begin_mer6+4]==mer_6[i6][4] &&
							buff_data[begin_mer6+5]==mer_6[i6][5])
						{
							count_mer6[i6]+=1;								
							break;
						}
						
					}
					begin_mer6=begin_mer6+6;
				}
				//write vector_data to array
				for( a6=0;a6<4096;++a6)
				{
					count_mer6[a6]=count_mer6[a6]/denominator_value(6,len,5,mer_count);				
				}
			}

			

		
			/* combine all frequency together */

			if(mer_count>=1 && (onlyk==0 || onlyk==1))
			{
				for(p1=0;p1<4;p1++)
				{
					vector_data[p1]=count_mer1[p1];
					total_average[p1]+=vector_data[p1];  
					fprintf(fid_vector,"%.9lf ",vector_data[p1]);
				}
			}

			if(mer_count>=2 && (onlyk==0 || onlyk==2))
			{
				for( p2=4,q2=0;q2<16;p2++,q2++)
				{
					vector_data[p2]=count_mer2[q2];
					total_average[p2]+=vector_data[p2];  
					fprintf(fid_vector,"%.9lf ",vector_data[p2]);
				}
			}

			if(mer_count>=3 && (onlyk==0 || onlyk==3))
			{
				for( p3=20,q3=0;q3<64;p3++,q3++)
				{
					vector_data[p3]=count_mer3[q3];
					total_average[p3]+=vector_data[p3];  
					fprintf(fid_vector,"%.9lf ",vector_data[p3]);
				}
			}
		
			if(mer_count>=4 && (onlyk==0 || onlyk==4))
			{
				for( p4=84,q4=0;q4<256;p4++,q4++)
				{
					vector_data[p4]=count_mer4[q4];
					total_average[p4]+=vector_data[p4];  
					fprintf(fid_vector,"%.9lf ",vector_data[p4]);
				}			
			}

			if(mer_count>=5 && (onlyk==0 || onlyk==5))
			{ 
				for( p5=340,q5=0;q5<1024;p5++,q5++)
				{
					vector_data[p5]=count_mer5[q5];
					total_average[p5]+=vector_data[p5];  
					fprintf(fid_vector,"%.9lf ",vector_data[p5]);
				}
			}

			if(mer_count>=6 && (onlyk==0 || onlyk==6))
			{ 
				for( p6=1364,q6=0;q6<4096;p6++,q6++)
				{
					vector_data[p6]=count_mer6[q6];
					total_average[p6]+=vector_data[p6];  
					fprintf(fid_vector,"%.9lf ",vector_data[p6]);
				}
			}

			

		 
			fputc('\n',fid_vector);	
		 
		}
	}

	/* mean */
	if(onlyk==0)
	{
		for(p1=0;p1<vector_dimen;++p1)
		{
			total_average[p1]/=row_id;
			fprintf(fid_mean,"%.9lf ",total_average[p1]);
		}
	}
	else
	{
		switch(onlyk)
		{
		case 1:   
			for(p1=0;p1<4;p1++)
			{
				total_average[p1]/=row_id;
				fprintf(fid_mean,"%.9lf ",total_average[p1]);
			}
			break;
		case 2:   
			for( p2=4,q2=0;q2<16;p2++,q2++)
			{
				total_average[p2]/=row_id;
				fprintf(fid_mean,"%.9lf ",total_average[p2]);
			}
			break;
		case 3:   
			for( p3=20,q3=0;q3<64;p3++,q3++)
			{
				total_average[p3]/=row_id;
				fprintf(fid_mean,"%.9lf ",total_average[p3]);
			}
			break;
		case 4:   
			for( p4=84,q4=0;q4<256;p4++,q4++)
			{
				total_average[p4]/=row_id;
				fprintf(fid_mean,"%.9lf ",total_average[p4]);
			}
			break;
		case 5:   
			for( p5=340,q5=0;q5<1024;p5++,q5++)
			{
				total_average[p5]/=row_id;
				fprintf(fid_mean,"%.9lf ",total_average[p5]);
			}
			break;
		case 6:   
			for( p6=1364,q6=0;q6<4096;p6++,q6++)
			{
				total_average[p6]/=row_id;
				fprintf(fid_mean,"%.9lf ",total_average[p6]);
			}
			break;
		}
	}

	fputc('\n',fid_mean);
	fclose(fid_mean);
	fclose(fid_vector);
	fclose(fid_sequence_file); 
 	free(buff_data);/* free */
	free(vector_data);/* free */
	free(total_average);/* free */

	if (mer_count>=2) free(count_mer2);	/* free */
	if (mer_count>=3) free(count_mer3);	/* free */
	if (mer_count>=4) free(count_mer4);	/* free */	
	if (mer_count>=5) free(count_mer5);	/* free */
	if (mer_count>=6) free(count_mer6);	/* free */
 

	if (REMOVE_TEMP_FILES) remove(sequence_only_file);

    /* end of frequency */ 
}


  void frequency_vertor_with_class_label(char sequence_only_file[],char vector_file_name[], 
					  char mean_frequency[],long mer_count, long onlyk,long class_label,long steplength_is_1)

{
	 // frequency_vertor_with_class_label("positive_seq_balanced","pos_freq_vector","pos_freq_mean",
	//		mer_count,onlyk,pos_class_label);


	/* Input: a seq file.  */
	/* Output: k-mer used frequency file, a row corresponding to a seq, */

	/* onlyk:0 -- all(1 to k); 1 -- only 1; 2 -- only 2.... */

 
	double *vector_data; /*  frequency vector */
	double *total_average; /* positive average, negative average */
	char *buff_data;

	long i,len; /* length of sequence */
	long row_id=0; /* row counter */
	long vector_dimen;
 
	double count_mer1[4]={0.000};  /* the count of A,T,C,G; 1-mer */
	double *count_mer2; /* the count of AA,AT,AC,AG...; 2-mer */
	double *count_mer3;
	double *count_mer4;
	double *count_mer5;
	double *count_mer6;
 

	long i1,i2,i3,i4,i5,i6;
	long begin_mer1=0,begin_mer2=0,begin_mer3=0,begin_mer4=0,begin_mer5=0,begin_mer6=0;			
	long a1,a2,a3,a4,a5,a6;
	long p1,p2,p3,p4,p5,p6;
	long q2,q3,q4,q5,q6;
	FILE *fid_sequence_file,*fid_vector,*fid_mean;


	if(mer_count==1) vector_dimen=4;
	if(mer_count==2) vector_dimen=20;
	if(mer_count==3) vector_dimen=84;
	if(mer_count==4) vector_dimen=340;
	if(mer_count==5) vector_dimen=1364;
	if(mer_count==6) vector_dimen=5460;


	buff_data=(char*)malloc(MAX_LINE_LENGTH*sizeof(char)); 

	vector_data=(double*)malloc(vector_dimen*sizeof(double));
	for(i=0;i<vector_dimen;i++)
		vector_data[i]=0.0;
	total_average=(double*)malloc(vector_dimen*sizeof(double));
	for(i=0;i<vector_dimen;i++)
		total_average[i]=0.0;

	if (mer_count>=2)
	{
		count_mer2=(double*)malloc(16*sizeof(double));
	}
	if (mer_count>=3)
	{
		count_mer3=(double*)malloc(64*sizeof(double));
	}	
	if (mer_count>=4)
	{
		count_mer4=(double*)malloc(256*sizeof(double));
	}	
	if (mer_count>=5)
	{
		count_mer5=(double*)malloc(1024*sizeof(double));
	}	
	if (mer_count>=6)
	{
		count_mer6=(double*)malloc(4096*sizeof(double));
	}	
 

	fid_sequence_file=fopen(sequence_only_file,"r");
	fid_vector=fopen(vector_file_name,"w");
	fid_mean=fopen(mean_frequency,"w"); 
	if(fid_vector==NULL  )
	{
		printf("Cannot open %s\n",vector_file_name);
		return;
	}
	if(fid_sequence_file==NULL )
	{
		printf("Cannot open %s\n",sequence_only_file);
		return;
	}
	if( fid_mean==NULL)
	{
		printf("Cannot open %s\n",mean_frequency);
		return;
	}

	while(!feof(fid_sequence_file))
	{
		if(fgets(buff_data,MAX_LINE_LENGTH,fid_sequence_file))
		{
			/* initial to zero */
			if (mer_count>=1)
			{
				for (i=0;i<4;i++)
				{
					count_mer1[i]=0.0;
				}
			}

			if (mer_count>=2)
			{
				for (i=0;i<16;i++)
				{
					count_mer2[i]=0.0;
				}
			}
			
			if (mer_count>=3)
			{
				for (i=0;i<64;i++)
				{
					count_mer3[i]=0.0;
				}
			}
			
			if (mer_count>=4)
			{
				for (i=0;i<256;i++)
				{
					count_mer4[i]=0.0;
				}
			}				
			
			if (mer_count>=5)
			{
				for (i=0;i<1024;i++)
				{
					count_mer5[i]=0.0;
				}
			}				
			
			if (mer_count>=6)
			{
				for (i=0;i<4096;i++)
				{
					count_mer6[i]=0.0;
				}
			}				
			
			

			if(row_id%1000==0)
				printf("frequency, %s, row index: %ld\n", vector_file_name,row_id);

			row_id++;

			/* computation */

			len=strlen(buff_data);
			if(buff_data[len-1]==10) /* carriage return */
			{
				len=len-1;
			}
			if(buff_data[len-1]==13) /* line break */
			{
				len=len-1;
			}


		 
			if(mer_count>=1 && (onlyk==0 || onlyk==1))
			{
				//1mer  
				begin_mer1=0;
				while(begin_mer1<len)
				{
					for(i1=0;i1<4;++i1)
					{
						if(buff_data[begin_mer1]==mer_1[i1])
						{
							count_mer1[i1]=count_mer1[i1]+1;
						}
					}
					begin_mer1=begin_mer1+1;
				}
				for(a1=0;a1<4;++a1)
				{		
					// fprintf(fid_vector_count,"%c:%lf ",mer_1[a1],count_mer1[a1]);
					count_mer1[a1]=count_mer1[a1]/(double)len;
				}
			}

			if(mer_count>=2 && (onlyk==0 || onlyk==2))
			{
				//2mer 
				begin_mer2=0;
				while(begin_mer2<len-1)
				{
					for(i2=0;i2<16;++i2)
					{			
						if(buff_data[begin_mer2]==mer_2[i2][0] &&
							buff_data[begin_mer2+1]==mer_2[i2][1])
						{
							count_mer2[i2]=count_mer2[i2]+1;					
							break; 
						}			
					}
					begin_mer2=steplength_is_1?begin_mer2+1:begin_mer2+2;
				}
				for( a2=0;a2<16;++a2)
				{				
					// fprintf(fid_vector_count,"%c%c:%lf ",mer_2[a2][0],mer_2[a2][1],count_mer2[a2]);
					count_mer2[a2]=count_mer2[a2]/(len/2.0);
				}
			}

			if(mer_count>=3 && (onlyk==0 || onlyk==3))
			{
				//3mer
				begin_mer3=0;
				while(begin_mer3<len-2)
				{
					for(i3=0;i3<64;++i3)
					{				
						if(buff_data[begin_mer3]==mer_3[i3][0] &&
							buff_data[begin_mer3+1]==mer_3[i3][1] &&
							buff_data[begin_mer3+2]==mer_3[i3][2])
						{
							count_mer3[i3]=count_mer3[i3]+1;
							break;
						}
					}
					 begin_mer3=steplength_is_1?begin_mer3+1:begin_mer3+3;
				}
				for( a3=0;a3<64;++a3)
				{		
					// fprintf(fid_vector_count,"%c%c%c:%lf ",mer_3[a3][0],mer_3[a3][1],mer_3[a3][2],count_mer3[a3]);
					count_mer3[a3]=count_mer3[a3]/(len/3.0);
				}	
			}

			if(mer_count>=4 && (onlyk==0 || onlyk==4))
			{
				//4mer
				begin_mer4=0;
				while(begin_mer4<len-3)
				{
					for(i4=0;i4<256;++i4)
					{			
						if(buff_data[begin_mer4]==mer_4[i4][0] &&
							buff_data[begin_mer4+1]==mer_4[i4][1] &&
							buff_data[begin_mer4+2]==mer_4[i4][2] &&
							buff_data[begin_mer4+3]==mer_4[i4][3] )
						{
							count_mer4[i4]+=1;								
							break;
						}
								
					}
					 begin_mer4=steplength_is_1?begin_mer4+1:begin_mer4+4;
				}

				//write vector_data to array
				for( a4=0;a4<256;++a4)
				{
					// fprintf(fid_vector_count,"%c%c%c%c:%lf ",mer_4[a4][0],mer_4[a4][1],mer_4[a4][2],mer_4[a4][3],count_mer4[a4]);
					count_mer4[a4]=count_mer4[a4]/(len/4.0);				
				}
			}


			if(mer_count>=5 && (onlyk==0 || onlyk==5))
			{
				//5mer 
				begin_mer5=0;
				while(begin_mer5<len-4)
				{
					for(i5=0;i5<1024;++i5)
					{			
						if(buff_data[begin_mer5]==mer_5[i5][0] &&
							buff_data[begin_mer5+1]==mer_5[i5][1] &&
							buff_data[begin_mer5+2]==mer_5[i5][2] &&
							buff_data[begin_mer5+3]==mer_5[i5][3] &&
							buff_data[begin_mer5+4]==mer_5[i5][4] )
						{
							count_mer5[i5]+=1;								
							break;
						}
								
					}
					 begin_mer5=steplength_is_1?begin_mer5+1:begin_mer5+5;
				}
				//write vector_data to array
				for( a5=0;a5<1024;++a5)
				{
					// fprintf(fid_vector_count,"%c%c%c%c%c:%lf ",mer_5[a5][0],mer_5[a5][1],mer_5[a5][2],mer_5[a5][3],mer_5[a5][4],count_mer5[a5]);
					count_mer5[a5]=count_mer5[a5]/(len/5.0);				
				}
			}

			if(mer_count>=6 && (onlyk==0 || onlyk==6))
			{
				//6mer 
				begin_mer6=0;
				while(begin_mer6<len-5)
				{
					for(i6=0;i6<4096;++i6)
					{			
						if(buff_data[begin_mer6]==mer_6[i6][0] &&
							buff_data[begin_mer6+1]==mer_6[i6][1] &&
							buff_data[begin_mer6+2]==mer_6[i6][2] &&
							buff_data[begin_mer6+3]==mer_6[i6][3] &&
							buff_data[begin_mer6+4]==mer_6[i6][4] &&
							buff_data[begin_mer6+5]==mer_6[i6][5])
						{
							count_mer6[i6]+=1;								
							break;
						}
						
					}
					 begin_mer6=steplength_is_1?begin_mer6+1:begin_mer6+6;
				}
				//write vector_data to array
				for( a6=0;a6<4096;++a6)
				{
					count_mer6[a6]=count_mer6[a6]/(len/6.0);				
				}
			}

			

		
			/* combine all frequency together */
			
			fprintf(fid_vector,"%ld ",class_label);

			if(mer_count>=1 && (onlyk==0 || onlyk==1))
			{
				for(p1=0;p1<4;p1++)
				{
					vector_data[p1]=count_mer1[p1];
					total_average[p1]+=vector_data[p1];  
					fprintf(fid_vector,"%.9lf ",vector_data[p1]);
				}
			}

			if(mer_count>=2 && (onlyk==0 || onlyk==2))
			{
				for( p2=4,q2=0;q2<16;p2++,q2++)
				{
					vector_data[p2]=count_mer2[q2];
					total_average[p2]+=vector_data[p2];  
					fprintf(fid_vector,"%.9lf ",vector_data[p2]);
				}
			}

			if(mer_count>=3 && (onlyk==0 || onlyk==3))
			{
				for( p3=20,q3=0;q3<64;p3++,q3++)
				{
					vector_data[p3]=count_mer3[q3];
					total_average[p3]+=vector_data[p3];  
					fprintf(fid_vector,"%.9lf ",vector_data[p3]);
				}
			}
		
			if(mer_count>=4 && (onlyk==0 || onlyk==4))
			{
				for( p4=84,q4=0;q4<256;p4++,q4++)
				{
					vector_data[p4]=count_mer4[q4];
					total_average[p4]+=vector_data[p4];  
					fprintf(fid_vector,"%.9lf ",vector_data[p4]);
				}			
			}

			if(mer_count>=5 && (onlyk==0 || onlyk==5))
			{ 
				for( p5=340,q5=0;q5<1024;p5++,q5++)
				{
					vector_data[p5]=count_mer5[q5];
					total_average[p5]+=vector_data[p5];  
					fprintf(fid_vector,"%.9lf ",vector_data[p5]);
				}
			}

			if(mer_count>=6 && (onlyk==0 || onlyk==6))
			{ 
				for( p6=1364,q6=0;q6<4096;p6++,q6++)
				{
					vector_data[p6]=count_mer6[q6];
					total_average[p6]+=vector_data[p6];  
					fprintf(fid_vector,"%.9lf ",vector_data[p6]);
				}
			}

		

		 
			fputc('\n',fid_vector);	
		 
		}
	}

	/* mean */
	if(onlyk==0)
	{
		for(p1=0;p1<vector_dimen;++p1)
		{
			total_average[p1]/=row_id;
			fprintf(fid_mean,"%.9lf ",total_average[p1]);
		}
	}
	else
	{
		switch(onlyk)
		{
		case 1:   
			for(p1=0;p1<4;p1++)
			{
				total_average[p1]/=row_id;
				fprintf(fid_mean,"%.9lf ",total_average[p1]);
			}
			break;
		case 2:   
			for( p2=4,q2=0;q2<16;p2++,q2++)
			{
				total_average[p2]/=row_id;
				fprintf(fid_mean,"%.9lf ",total_average[p2]);
			}
			break;
		case 3:   
			for( p3=20,q3=0;q3<64;p3++,q3++)
			{
				total_average[p3]/=row_id;
				fprintf(fid_mean,"%.9lf ",total_average[p3]);
			}
			break;
		case 4:   
			for( p4=84,q4=0;q4<256;p4++,q4++)
			{
				total_average[p4]/=row_id;
				fprintf(fid_mean,"%.9lf ",total_average[p4]);
			}
			break;
		case 5:   
			for( p5=340,q5=0;q5<1024;p5++,q5++)
			{
				total_average[p5]/=row_id;
				fprintf(fid_mean,"%.9lf ",total_average[p5]);
			}
			break;
		case 6:   
			for( p6=1364,q6=0;q6<4096;p6++,q6++)
			{
				total_average[p6]/=row_id;
				fprintf(fid_mean,"%.9lf ",total_average[p6]);
			}
			break;
		
		}
	}

	fputc('\n',fid_mean);
	fclose(fid_mean);
	fclose(fid_vector);
	fclose(fid_sequence_file); 
 	free(buff_data);/* free */
	free(vector_data);/* free */
	free(total_average);/* free */

	if (mer_count>=2) free(count_mer2);	/* free */
	if (mer_count>=3) free(count_mer3);	/* free */
	if (mer_count>=4) free(count_mer4);	/* free */	
	if (mer_count>=5) free(count_mer5);	/* free */
	if (mer_count>=6) free(count_mer6);	/* free */
 

	if (REMOVE_TEMP_FILES) remove(sequence_only_file);

  } /* end of frequency with class label */ 


  void file_to_libsvm_format(char old_file[],char new_file[],long classlabel,long mer_count,long onlyk)
  {

	  FILE *fid_origin, *file_new; long i;
	  //long rowid=0;
	  double x;
	  long vector_dimen;
	  //const char "%ld:%.9lf "[]="%ld:%.9lf "; 

	  printf("trans %s to svm format ... \n",old_file);
	  
	  
	  file_new=fopen(new_file,"w");
	  if(file_new==NULL )
	  {
		  printf("Cannot open %s\n",new_file);
		  return;
	  }
	  
	  if (onlyk==0)
	  {
		  if(mer_count==1) vector_dimen=4;
		  if(mer_count==2) vector_dimen=20;
		  if(mer_count==3) vector_dimen=84;
		  if(mer_count==4) vector_dimen=340;
		  if(mer_count==5) vector_dimen=1364;
		  if(mer_count==6) vector_dimen=5460;
		 
	  }
	  else
	  {
		  if(mer_count==1) vector_dimen=4;
		  if(mer_count==2) vector_dimen=16;
		  if(mer_count==3) vector_dimen=64;
		  if(mer_count==4) vector_dimen=256;
		  if(mer_count==5) vector_dimen=1024;
		  if(mer_count==6) vector_dimen=4096;
	 
	  }
	  
	  fid_origin=fopen(old_file,"r"); 
	  if(fid_origin==NULL)
	  {
		  printf("Cannot open %s\n",old_file);
		  return;
	  }
	  
	  while(!feof(fid_origin))
	  {
		  fprintf(file_new,"%ld ",classlabel);  
		  for(i=0;i<vector_dimen;++i)
		  {
			  fscanf(fid_origin,"%lf",&x); 		 
			  fprintf(file_new,"%ld:%.9lf ",i+1,x);  
		  }
		  fputc('\n',file_new);   
		  
	  }
	  
	  fclose(fid_origin);
	  fclose(file_new); 

  } //file_to_libsvm_format



  void file_to_libsvm_format2(char old_file[],char new_file[],long mer_count,long onlyk)
  {

	  FILE *fid_origin, *file_new; long i;
	  long rowid=0;
	  double x;
	  long vector_dimen,classlabel;
	  //const char "%ld:%.9lf "[]="%ld:%.9lf ";
	  printf("trans %s to svm format ... \n",old_file);
	  
	  
	  file_new=fopen(new_file,"w");
	  if(file_new==NULL )
	  {
		  printf("Cannot open %s\n",new_file);
		  return;
	  }
	  
	  if (onlyk==0)
	  {
		  if(mer_count==1) vector_dimen=4;
		  if(mer_count==2) vector_dimen=20;
		  if(mer_count==3) vector_dimen=84;
		  if(mer_count==4) vector_dimen=340;
		  if(mer_count==5) vector_dimen=1364;
		  if(mer_count==6) vector_dimen=5460;
	 
	  }
	  else
	  {
		  if(mer_count==1) vector_dimen=4;
		  if(mer_count==2) vector_dimen=16;
		  if(mer_count==3) vector_dimen=64;
		  if(mer_count==4) vector_dimen=256;
		  if(mer_count==5) vector_dimen=1024;
		  if(mer_count==6) vector_dimen=4096;
		  
	  }
	  
	  fid_origin=fopen(old_file,"r"); 
	  if(fid_origin==NULL)
	  {
		  printf("Cannot open %s\n",old_file);
		  return;
	  }
	  
	  while(!feof(fid_origin))
	  {
		  fscanf(fid_origin,"%ld",&classlabel);
		  fprintf(file_new,"%ld ",classlabel);  
		  for(i=0;i<vector_dimen;++i)
		  {
			  fscanf(fid_origin,"%lf",&x); 		 
			  fprintf(file_new,"%ld:%.9lf ",i+1,x);  
		  }
		  fputc('\n',file_new);   
		  rowid++;
		  if(rowid%10000==0)
			printf("to svm format, %s, row index: %ld\n", old_file,rowid);
		  
	  }
	  
	  fclose(fid_origin);
	  fclose(file_new); 

  } //file_to_libsvm_format



  void file_to_willows_format(char old_file[],char new_file[],long classlabel,long mer_count,long onlyk)
  {

	  FILE *fid_origin, *file_new; long i;
	  //long rowid=0;
	  double x;
	  long vector_dimen;
	  
	  printf("to willows format ... %s \n",old_file);
	  
	  
	  file_new=fopen(new_file,"w");
	  if(file_new==NULL )
	  {
		  printf("Cannot open %s\n",new_file);
		  return;
	  }
	  
	  if (onlyk==0)
	  {
		  if(mer_count==1) vector_dimen=4;
		  if(mer_count==2) vector_dimen=20;
		  if(mer_count==3) vector_dimen=84;
		  if(mer_count==4) vector_dimen=340;
		  if(mer_count==5) vector_dimen=1364;
		  if(mer_count==6) vector_dimen=5460;
	 
	  }
	  else
	  {
		  if(mer_count==1) vector_dimen=4;
		  if(mer_count==2) vector_dimen=16;
		  if(mer_count==3) vector_dimen=64;
		  if(mer_count==4) vector_dimen=256;
		  if(mer_count==5) vector_dimen=1024;
		  if(mer_count==6) vector_dimen=4096;
		   
	  }
	  
	  fid_origin=fopen(old_file,"r"); 
	  if(fid_origin==NULL)
	  {
		  printf("Cannot open %s\n",old_file);
		  return;
	  }

	  //r o o o ...
	  fprintf(file_new,"%c ",'r');
	  for (i=0;i<vector_dimen;i++)
	  {
		fprintf(file_new,"%c ",'o');
	  }
	  fprintf(file_new,"%c",'\n');

	  //a a1 a2 a3 ...
	  fprintf(file_new,"%c ",'a');
	  for (i=0;i<vector_dimen;i++)
	  {
		  fprintf(file_new,"a%ld ",i);
	  }
	  fprintf(file_new,"%c",'\n');

	  
	  while(!feof(fid_origin))
	  {
		  fprintf(file_new,"%ld ",classlabel);  
		  for(i=0;i<vector_dimen;++i)
		  {
			  fscanf(fid_origin,"%lf",&x); 		 
			  fprintf(file_new,"%.9lf ",x);  
		  }
		  fputc('\n',file_new);   
		  
	  }
	  
	  fclose(fid_origin);
	  fclose(file_new); 
}

  void file_to_willows_format_no_header(char old_file[],char new_file[],long classlabel,long mer_count,long onlyk)
  {
	  
	  FILE *fid_origin, *file_new; long i;
	  //long rowid=0;
	  double x;
	  long vector_dimen;
	  
	  printf("trans %s to willows format ... \n",old_file);
	  
	  
	  file_new=fopen(new_file,"w");
	  if(file_new==NULL )
	  {
		  printf("Cannot open %s\n",new_file);
		  return;
	  }
	  
	  if (onlyk==0)
	  {
		  if(mer_count==1) vector_dimen=4;
		  if(mer_count==2) vector_dimen=20;
		  if(mer_count==3) vector_dimen=84;
		  if(mer_count==4) vector_dimen=340;
		  if(mer_count==5) vector_dimen=1364;
		  if(mer_count==6) vector_dimen=5460;
	 
	  }
	  else
	  {
		  if(mer_count==1) vector_dimen=4;
		  if(mer_count==2) vector_dimen=16;
		  if(mer_count==3) vector_dimen=64;
		  if(mer_count==4) vector_dimen=256;
		  if(mer_count==5) vector_dimen=1024;
		  if(mer_count==6) vector_dimen=4096;
	 
	  }
	  
	  fid_origin=fopen(old_file,"r"); 
	  if(fid_origin==NULL)
	  {
		  printf("Cannot open %s\n",old_file);
		  return;
	  }
	  
 
	  
	  while(!feof(fid_origin))
	  {
		  fprintf(file_new,"%ld ",classlabel);  
		  for(i=0;i<vector_dimen;++i)
		  {
			  fscanf(fid_origin,"%lf",&x); 		 
			  fprintf(file_new,"%.9lf ",x);  
		  }
		  fputc('\n',file_new);   
		  
	  }
	  
	  fclose(fid_origin);
	  fclose(file_new); 

  } // file to willows


 
 void frequency_label_svm(char sequence_only_file[],char vector_file_name[], long kmer_count, long onlyk,
					  long class_label,long balanced_count,
					  long steplength_is_1,char seq_desc[],char allsvm[],long denomi_type,
					  char out_file[],FILE*logs,long mustbebalanced,char allsvmdesc[]
					  , long outmsg, long rmtempfile) //long groups[],

{
 
	double denominator ; 

	double *vector_data; /*  frequency vector */ 
	double *total_average; /* positive average, negative average */
	char *buff_data; 

	long i,len; /* length of sequence */
	long row_id=0; /* row counter */
	long vector_dimen;
 
	double count_mer1[4]={0.000};  /* the count of A,T,C,G; 1-mer */
	double *count_mer2; /* the count of AA,AT,AC,AG...; 2-mer */
	double *count_mer3;
	double *count_mer4;
	double *count_mer5;
	double *count_mer6;

	long i1,i2,i3,i4,i5,i6 ;
	long begin_mer1=0,begin_mer2=0,begin_mer3=0,begin_mer4=0,begin_mer5=0,begin_mer6=0 ;			
	long a1,a2,a3,a4,a5,a6 ;
	long p1,p2,p3,p4,p5,p6 ;
	long q2,q3,q4,q5,q6 ;
	FILE *fid_sequence_file,*fid_vector,*fid_allsvm;
	FILE  *fid_input_desc,*fid_allsvmdesc; 	

	char *my_buff; /* store all the data, then write them to file for only one time */
	long my_buff_length;
	char temp_str[200];

	if(kmer_count==1) vector_dimen=4;
	if(kmer_count==2) vector_dimen=20;
	if(kmer_count==3) vector_dimen=84;
	if(kmer_count==4) vector_dimen=340;
	if(kmer_count==5) vector_dimen=1364;
	if(kmer_count==6) vector_dimen=5460;
 

	buff_data=(char*)malloc(MAX_LINE_LENGTH*sizeof(char)); 

	vector_data=(double*)malloc(vector_dimen*sizeof(double));
	for(i=0;i<vector_dimen;i++)
		vector_data[i]=0.0;
	total_average=(double*)malloc(vector_dimen*sizeof(double));
	for(i=0;i<vector_dimen;i++)
		total_average[i]=0.0;

	if (kmer_count>=2)
	{
		count_mer2=(double*)malloc(16*sizeof(double));
	}
	if (kmer_count>=3)
	{
		count_mer3=(double*)malloc(64*sizeof(double));
	}	
	if (kmer_count>=4)
	{
		count_mer4=(double*)malloc(256*sizeof(double));
	}	
	if (kmer_count>=5)
	{
		count_mer5=(double*)malloc(1024*sizeof(double));
	}	
	if (kmer_count>=6)
	{
		count_mer6=(double*)malloc(4096*sizeof(double));
	}	
 

	fid_sequence_file=fopen(sequence_only_file,"r");
	fid_vector=fopen(vector_file_name,"w");
	 
	fid_allsvm=fopen(allsvm,"a"); 
	fid_allsvmdesc=fopen(allsvmdesc,"a"); 
	if(fid_vector==NULL  )
	{
		printf("Cannot open %s\n",vector_file_name);
		fprintf(logs,"Cannot open %s\n",vector_file_name);
		fflush(logs);
		return;
	}
	if(fid_sequence_file==NULL )
	{
		printf("Cannot open %s\n",sequence_only_file);
		fprintf(logs, "Cannot open %s\n",sequence_only_file);
		fflush(logs);
		return;
	}
 
	if( fid_allsvm==NULL)
	{
		printf("Cannot open %s\n",allsvm);
		fprintf(logs, "Cannot open %s\n",allsvm);
		fflush(logs);
		return;
	}
	
	if( fid_allsvmdesc==NULL)
	{
		printf("Cannot open %s\n",allsvmdesc);
		fprintf(logs, "Cannot open %s\n",allsvmdesc);
		fflush(logs);
		return;
	}
	if(outmsg) printf("\nCalculating k-mer usage frequencies ... \n");
	fprintf(logs,"\nCalculating k-mer usage frequencies ... \n");
	fflush(logs);
 
	my_buff_length=(vector_dimen+300)*40;
 
	my_buff=(char*)malloc(sizeof(char)*my_buff_length);


	fid_input_desc=fopen(seq_desc,"r");
	if (fid_input_desc ==NULL)
	{
		printf("Cannot open file.\n");
		fprintf(logs, "Cannot open file.\n");
		fflush(logs);
		return ;
	}	


	while(!feof(fid_sequence_file))
	{
		if(fgets(buff_data,MAX_LINE_LENGTH,fid_sequence_file))
		{
		
			/* initial to zero */
			if (kmer_count>=1)
			{
				for (i=0;i<4;i++)
				{
					count_mer1[i]=0.0;
				}
			}

			if (kmer_count>=2)
			{
				for (i=0;i<16;i++)
				{
					count_mer2[i]=0.0;
				}
			}
			
			if (kmer_count>=3)
			{
				for (i=0;i<64;i++)
				{
					count_mer3[i]=0.0;
				}
			}
			
			if (kmer_count>=4)
			{
				for (i=0;i<256;i++)
				{
					count_mer4[i]=0.0;
				}
			}				
			
			if (kmer_count>=5)
			{
				for (i=0;i<1024;i++)
				{
					count_mer5[i]=0.0;
				}
			}				
			
			if (kmer_count>=6)
			{
				for (i=0;i<4096;i++)
				{
					count_mer6[i]=0.0;
				}
			}				
			
		 
			if(row_id%1000==0)
			{
				if(outmsg) printf("   %s, row index: %ld\n", vector_file_name,row_id);
				fprintf(logs, "   %s, row index: %ld\n", vector_file_name,row_id);
				fflush(logs);
			}
			

			/* computation */

			len=strlen(buff_data);
			if(buff_data[len-1]==10) /* carriage return */
			{
				len=len-1;
			}
			if(buff_data[len-1]==13) /* line break */
			{
				len=len-1;
			}


		 
			if(kmer_count>=1 && (onlyk==0 || onlyk==1))
			{
				//1mer 
				long kth=1;
				begin_mer1=0;
				while(begin_mer1<len)
				{
					for(i1=0;i1<4;++i1)
					{
						if(buff_data[begin_mer1]==mer_1[i1])
						{
							count_mer1[i1]=count_mer1[i1]+1;
						}
					}
					begin_mer1++;
				}
				for(a1=0;a1<4;++a1)
				{				 
					denominator=denominator_value(kth,len,denomi_type,kmer_count);
					count_mer1[a1]=count_mer1[a1]/denominator;
				}
			}

			if(kmer_count>=2 && (onlyk==0 || onlyk==2))
			{
				
				//2mer 
				long kth=2;
				begin_mer2=0;
				while(begin_mer2<len-1)
				{
					for(i2=0;i2<16;++i2)
					{			
						if(buff_data[begin_mer2]==mer_2[i2][0] &&
							buff_data[begin_mer2+1]==mer_2[i2][1])
						{
							count_mer2[i2]=count_mer2[i2]+1;					
							break; 
						}			
					}
					begin_mer2=steplength_is_1?begin_mer2+1:begin_mer2+2;
				}
				for( a2=0;a2<16;++a2)
				{
					denominator=denominator_value(kth,len,denomi_type,kmer_count);
					count_mer2[a2]=count_mer2[a2]/denominator;
				}
			}

			if(kmer_count>=3 && (onlyk==0 || onlyk==3))
			{
				//3mer
				long kth=3;
				begin_mer3=0;
				while(begin_mer3<len-2)
				{
					for(i3=0;i3<64;++i3)
					{				
						if(buff_data[begin_mer3]==mer_3[i3][0] &&
							buff_data[begin_mer3+1]==mer_3[i3][1] &&
							buff_data[begin_mer3+2]==mer_3[i3][2])
						{
							count_mer3[i3]=count_mer3[i3]+1;
							break;
						}
					}
					begin_mer3=steplength_is_1?begin_mer3+1:begin_mer3+3;
				}
				for( a3=0;a3<64;++a3)
				{		
					denominator=denominator_value(kth,len,denomi_type,kmer_count);
					count_mer3[a3]=count_mer3[a3]/denominator;
				}	
			}

			if(kmer_count>=4 && (onlyk==0 || onlyk==4))
			{
				//4mer
				long kth=4;
				begin_mer4=0;
				while(begin_mer4<len-3)
				{
					for(i4=0;i4<256;++i4)
					{			
						if(buff_data[begin_mer4]==mer_4[i4][0] &&
							buff_data[begin_mer4+1]==mer_4[i4][1] &&
							buff_data[begin_mer4+2]==mer_4[i4][2] &&
							buff_data[begin_mer4+3]==mer_4[i4][3] )
						{
							count_mer4[i4]+=1;								
							break;
						}
								
					}
					begin_mer4=steplength_is_1?begin_mer4+1:begin_mer4+4;
				}

				for( a4=0;a4<256;++a4)
				{
					denominator=denominator_value(kth,len,denomi_type,kmer_count);
					count_mer4[a4]=count_mer4[a4]/denominator;				
				}
			}


			if(kmer_count>=5 && (onlyk==0 || onlyk==5))
			{
				//5mer 
				long kth=5;
				begin_mer5=0;
				while(begin_mer5<len-4)
				{
					for(i5=0;i5<1024;++i5)
					{			
						if(buff_data[begin_mer5]==mer_5[i5][0] &&
							buff_data[begin_mer5+1]==mer_5[i5][1] &&
							buff_data[begin_mer5+2]==mer_5[i5][2] &&
							buff_data[begin_mer5+3]==mer_5[i5][3] &&
							buff_data[begin_mer5+4]==mer_5[i5][4] )
						{
							count_mer5[i5]+=1;								
							break;
						}
								
					}
					begin_mer5=steplength_is_1?begin_mer5+1:begin_mer5+5;
				}
				for( a5=0;a5<1024;++a5)
				{
					denominator=denominator_value(kth,len,denomi_type,kmer_count);
					count_mer5[a5]=count_mer5[a5]/denominator;				
				}
			}

			if(kmer_count>=6 && (onlyk==0 || onlyk==6))
			{
				//6mer 
				long kth=6;
				begin_mer6=0;
				while(begin_mer6<len-5)
				{
					for(i6=0;i6<4096;++i6)
					{			
						if(buff_data[begin_mer6]==mer_6[i6][0] &&
							buff_data[begin_mer6+1]==mer_6[i6][1] &&
							buff_data[begin_mer6+2]==mer_6[i6][2] &&
							buff_data[begin_mer6+3]==mer_6[i6][3] &&
							buff_data[begin_mer6+4]==mer_6[i6][4] &&
							buff_data[begin_mer6+5]==mer_6[i6][5])
						{
							count_mer6[i6]+=1;								
							break;
						}
						
					}
					begin_mer6=steplength_is_1?begin_mer6+1:begin_mer6+6;
				}
				for( a6=0;a6<4096;++a6)
				{
					denominator=denominator_value(kth,len,denomi_type,kmer_count);
					count_mer6[a6]=count_mer6[a6]/denominator;				
				}
			}

		

		
			/* combine all frequency together */
			sprintf(my_buff,"%ld ",class_label);

			if(kmer_count>=1 && (onlyk==0 || onlyk==1))
			{
				for(p1=0;p1<4;p1++)
				{
					vector_data[p1]=count_mer1[p1];
					total_average[p1]+=vector_data[p1];
					if (onlyk)
					{
						sprintf(temp_str,"%ld:%.9lf ",p1+1,vector_data[p1]);
					}
					else
					{
						sprintf(temp_str,"%ld:%.9lf ",p1+1,vector_data[p1]);
					}
					
					strcat(my_buff,temp_str);
					
				}
			}

			if(kmer_count>=2 && (onlyk==0 || onlyk==2))
			{
				for( p2=4,q2=0;q2<16;p2++,q2++)
				{
					vector_data[p2]=count_mer2[q2];
					total_average[p2]+=vector_data[p2];  
					if (onlyk)
					{
						sprintf(temp_str,"%ld:%.9lf ",q2+1,vector_data[p2]);
					}
					else
					{
						sprintf(temp_str,"%ld:%.9lf ",p2+1,vector_data[p2]);
					}
					strcat(my_buff,temp_str);
					
				}
			}

			if(kmer_count>=3 && (onlyk==0 || onlyk==3))
			{
				for( p3=20,q3=0;q3<64;p3++,q3++)
				{
					vector_data[p3]=count_mer3[q3];
					total_average[p3]+=vector_data[p3];  
					if (onlyk)
					{
						sprintf(temp_str,"%ld:%.9lf ",q3+1,vector_data[p3]);
					}
					else
					{
						sprintf(temp_str,"%ld:%.9lf ",p3+1,vector_data[p3]);
					}
					strcat(my_buff,temp_str);
					
				}
			}
		
			if(kmer_count>=4 && (onlyk==0 || onlyk==4))
			{
				for( p4=84,q4=0;q4<256;p4++,q4++)
				{
					vector_data[p4]=count_mer4[q4];
					total_average[p4]+=vector_data[p4];  
					if (onlyk)
					{
						sprintf(temp_str,"%ld:%.9lf ",q4+1,vector_data[p4]);
					}
					else
					{
						sprintf(temp_str,"%ld:%.9lf ",p4+1,vector_data[p4]);
					}
					strcat(my_buff,temp_str);
					
				}			
			}

			if(kmer_count>=5 && (onlyk==0 || onlyk==5))
			{ 
				for( p5=340,q5=0;q5<1024;p5++,q5++)
				{
					vector_data[p5]=count_mer5[q5];
					total_average[p5]+=vector_data[p5];  
					if (onlyk)
					{
						sprintf(temp_str,"%ld:%.9lf ",q5+1,vector_data[p5]);
					}
					else
					{
						sprintf(temp_str,"%ld:%.9lf ",p5+1,vector_data[p5]);
					}
					strcat(my_buff,temp_str);
				
				}
			}

			if(kmer_count>=6 && (onlyk==0 || onlyk==6))
			{ 
				for( p6=1364,q6=0;q6<4096;p6++,q6++)
				{
					vector_data[p6]=count_mer6[q6];
					total_average[p6]+=vector_data[p6]; 
					if (onlyk)
					{
						sprintf(temp_str,"%ld:%.9lf ",q6+1,vector_data[p6]);
					}
					else
					{
						sprintf(temp_str,"%ld:%.9lf ",p6+1,vector_data[p6]);
					}
					strcat(my_buff,temp_str);
					
				}
			}			

			//change		 
			//write

			if (mustbebalanced)
			{	

			 	fprintf(fid_allsvm,"%s\n",my_buff); //	=================================

				fgets(buff_data,MAX_LINE_LENGTH,fid_input_desc);
				fprintf(fid_allsvmdesc,"%s",buff_data); //
			}
			else
			{
				fprintf(fid_allsvm,"%s\n",my_buff); //	
				
				fgets(buff_data,MAX_LINE_LENGTH,fid_input_desc);
				fprintf(fid_allsvmdesc,"%s",buff_data); //

				fprintf(fid_vector,"%s\n",my_buff); //
			}			
		
			row_id++;	//==========================================================
		 
		 
		}
	}

	if(outmsg) printf("   Total lines: %ld\n\n",row_id);
 
	fclose(fid_vector);
	fclose(fid_sequence_file); //
	fclose(fid_input_desc);
	fclose(fid_allsvm);	
	fclose(fid_allsvmdesc);
 	free(buff_data);  /* free */
	free(vector_data);  /* free */
	free(total_average);  /* free */
	free(my_buff);  /* free */

	if (kmer_count>=2) free(count_mer2);	/* free */
	if (kmer_count>=3) free(count_mer3);	/* free */
	if (kmer_count>=4) free(count_mer4);	/* free */	
	if (kmer_count>=5) free(count_mer5);	/* free */
	if (kmer_count>=6) free(count_mer6);	/* free */
 

	if (rmtempfile) 
	{
		remove(sequence_only_file);	 
		remove(vector_file_name);				 
		remove(seq_desc);	
		//remove(allsvm);
	}

  } /* end of frequency_label_svm */


