///////////////////////////////////////////////////////////////////
//                                                               
//  PLEK - predictor of lncRNAs and mRNAs based on k-mer scheme  
//  Authors: Aimin Li, Junying Zhang                             
//  Contacts: LiAiminMail@gmail.com, jyzhang@mail.xidian.edu.cn  
//  Webcite: https://sourceforge.net/projects/plek/                      
//  Updated on: May 25, 2014                                                   
//                                                               
///////////////////////////////////////////////////////////////////

/* LANG=C  gcc -g -Wall PLEK_spsn.c -o PLEK_spsn -lm */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "time.h"  
 
#ifndef MAX_SAMPLE_COUNT  
#define MAX_SAMPLE_COUNT 2000000  
#endif 
#ifndef MAX_LINE_LENGTH  
#define MAX_LINE_LENGTH  500000  
#endif 


int main(int argc,char *argv[])
{	
	time_t rawtime;
	struct tm * timeinfo;
	char test_file[200]="", predict_file[200]="";
	char desc_file[200]="", desc_file_with_class[200]="";
	char output_file[200]="lnc_logs.txt", *buff_data;
 	FILE *f_test, *f_predict, *f_desc, *f_desc_with_class, *logs;
 	long i,j, TP=0,FP=0,TN=0,FN=0, outmsg=0, rmtempfile=1 ;
	long *testlabel, *predictlabel;
	double *predict_dec_value;
	long totalsamples=0 , pos_class_label=1, neg_class_label=0, input_type=1; 


	if (argc>=4)
	{
		for (i=0;i<argc;i++)
		{
			if (strcmp(argv[i],"-svm")==0)
			{
				strcpy(test_file,argv[i+1]); //input: get the pre-defined class label from the first column; svm-file			 
			}
			if (strcmp(argv[i],"-predict")==0)
			{
				strcpy(predict_file,argv[i+1]);	//input:  svm-predict results		 
			}
			if (strcmp(argv[i],"-desc")==0)
			{
				strcpy(desc_file,argv[i+1]); //input: desc 			 
			}
			if (strcmp(argv[i],"-descclass")==0)
			{
				strcpy(desc_file_with_class,argv[i+1]);	 //output:  output_predicted
			}
			if (strcmp(argv[i],"-output")==0)
			{
				strcpy(output_file,argv[i+1]);	//output: logs, 		 
			}
			if (strcmp(argv[i],"-isoutmsg")==0)
			{
				outmsg=atol(argv[i+1]);
			}
			if (strcmp(argv[i],"-isrmtempfile")==0)
			{
				rmtempfile=atol(argv[i+1]);
			}
			if (strcmp(argv[i],"-input_type")==0)
			{
				input_type=atol(argv[i+1]);
			}
			//  for modeling: 	_options.input_type=1
			//  for known Pos & Neg:   _options.input_type=2
			// 	for known Pos:	_options.input_type=3 			 
			//  for unkown: _options.input_type=4
			//  for known Neg:	_options.input_type=5
			
		}
	} 
	if(argc<4)
	{		
		printf("usage on linux: ./PLEK_spsn -predict predictX -desc desc_file -descclass desc_file_with_class -output output\n");
		return -1;
	}
					 
	logs=fopen(output_file,"a"); // output
	if(logs==NULL)
	{
		printf("WARNING: Cannot open %s\n","lnc_logs.txt");
		return -1;
	}
	f_test=fopen(test_file,"r"); // input
	if(f_test==NULL)
	{
		printf("WARNING: Cannot open %s\n",test_file);
		return -1;
	}	
	f_predict=fopen(predict_file,"r"); // input
	if(f_predict==NULL)
	{
		printf("WARNING: Cannot open %s\n",predict_file);
		return -1;
	}
	if(argc>5)
	{
		f_desc=fopen(desc_file,"r"); // input
		if(f_desc==NULL)
		{
			printf("WARNING: Cannot open %s\n",desc_file);
			return -1;
		}
		f_desc_with_class=fopen(desc_file_with_class,"w"); // output
		if(f_desc_with_class==NULL)
		{
			printf("WARNING: Cannot open %s\n",desc_file_with_class);
			return -1;
		}

	}

	time(&rawtime);timeinfo=localtime(&rawtime);
	fprintf(logs, "	The current date/time is: %s\n", asctime (timeinfo) ); 
	//if (outmsg)	printf("\n   Statistic %s, %s ...\n",test_file,predict_file);
	

	testlabel=(long *)malloc(sizeof(long)*MAX_SAMPLE_COUNT);
	predictlabel=(long *)malloc(sizeof(long)*MAX_SAMPLE_COUNT);
	predict_dec_value=(double *)malloc(sizeof(double)*MAX_SAMPLE_COUNT);
	
	for(i=0; i<MAX_SAMPLE_COUNT; i++)
		predict_dec_value[i]=0;

	i=0;
	buff_data=(char*)malloc(MAX_LINE_LENGTH*sizeof(char));
	// read svm file
	if (input_type==2 || input_type==3 || input_type==5)
	{
		while(!feof(f_test))
		{
			if(fgets(buff_data, MAX_LINE_LENGTH, f_test) != NULL)/* read the first number then next line */
			{
				if(sscanf(buff_data,"%ld",&testlabel[i])==1) 
					i++;
			}
			else
				break	;
			//blank line at the end of the file.
			
		}
	}

	// read predict file
	if (input_type==2 || input_type==3 || input_type==5 || input_type==4)
	{
		j=0;
		while(!feof(f_predict))
		{
			if(fscanf(f_predict,"%ld%lf",&predictlabel[j], &predict_dec_value[j])==2)
				j++;
			else
				break;
			//blank line at the end of the file.
		}
		totalsamples=j;
	}
	

	//known class label
	if (input_type==2 || input_type==3 || input_type==5)
	{
		for (i=0;i<totalsamples;i++)
		{
			if(argc>5) fgets(buff_data, MAX_LINE_LENGTH, f_desc);

			if (testlabel[i]==pos_class_label && predictlabel[i]==pos_class_label)
			{
				TP++;
				if(argc>5) fprintf(f_desc_with_class,"%s ","TP");
			}
			if (testlabel[i]==pos_class_label && predictlabel[i]==neg_class_label)
			{
				FN++;
				if(argc>5) fprintf(f_desc_with_class,"%s ","FN");
			}
			if (testlabel[i]==neg_class_label && predictlabel[i]==pos_class_label)
			{
				FP++;
				if(argc>5) fprintf(f_desc_with_class,"%s ","FP");
			}
			if (testlabel[i]==neg_class_label && predictlabel[i]==neg_class_label)
			{
				TN++;
				if(argc>5) fprintf(f_desc_with_class,"%s ","TN");
			}

			if(argc>5) fprintf(f_desc_with_class,"%s",buff_data);
		}
	}

	//unknown
	if (input_type==4)
	{
		for (i=0;i<totalsamples;i++)
		{
			if(argc>5) fgets(buff_data, MAX_LINE_LENGTH, f_desc);
			if (predictlabel[i]==pos_class_label)
			{
				 fprintf(f_desc_with_class,"%s\t%lf\t", "Coding", predict_dec_value[i]);
			}
			else
			{
				fprintf(f_desc_with_class,"%s\t%lf\t","Non-coding", predict_dec_value[i]);
			}			
			
			if(argc>5) fprintf(f_desc_with_class,"%s",buff_data);
		}
	}
	

	fclose(f_test);
	fclose(f_predict);
	if(argc>5) fclose(f_desc);
	if(argc>5) fclose(f_desc_with_class);
	//remove 

	free(buff_data);
	free(testlabel);
	free(predictlabel);
	
 if (input_type==2 || input_type==3 || input_type==5)
 {
	 if (outmsg)
	 {
		 printf("Calculate Precision, Sensitivity, Specificity, Accuracy ...\n");
		 printf("  TP=%ld\n",TP);
		 printf("  FN=%ld\n",FN);
		 printf("  FP=%ld\n",FP);
		 printf("  TN=%ld\n",TN);
		 printf("  Precision=TP/(TP+FP)=%lf\n",TP*1.0/(TP+FP));
		 printf("  Sensitivity=TP/(TP+FN)=%lf\n",TP*1.0/(TP+FN)); 
		 printf("  Specificity=TN/(TN+FP)=%lf\n",TN*1.0/(TN+FP)); 
		 printf("  Accuracy=(TP+TN)/(TP+TN+FP+FN)=%lf\n",(TP+TN)*1.0/(TP+TN+FP+FN)); 
	 }
	
	 fprintf(logs,"Calculate Precision, Sensitivity, Specificity, Accuracy ...\n");
	 fprintf(logs,"  TP=%ld\n",TP);
	 fprintf(logs,"  FN=%ld\n",FN);
	 fprintf(logs,"  FP=%ld\n",FP);
	 fprintf(logs,"  TN=%ld\n",TN);
	 fprintf(logs,"  Precision=TP/(TP+FP)=%lf\n",TP*1.0/(TP+FP));
	 fprintf(logs,"  Sensitivity=TP/(TP+FN)=%lf\n",TP*1.0/(TP+FN)); 
	 fprintf(logs,"  Specificity=TN/(TN+FP)=%lf\n",TN*1.0/(TN+FP)); 
	 fprintf(logs,"  Accuracy=(TP+TN)/(TP+TN+FP+FN)=%lf\n",(TP+TN)*1.0/(TP+TN+FP+FN));
	 fprintf(logs,"  MCC=(TP*TN-FP*FN)/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))=%lf\n",
		 1.0*(TP*TN-FP*FN)/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))); 
	 fprintf(logs,"  %s\n","Precision   Sensitivity   Specificity   Accuracy   MCC"); 
	 fprintf(logs,"  %lf %lf %lf %lf %lf\n",TP*1.0/(TP+FP),TP*1.0/(TP+FN),TN*1.0/(TN+FP),
		(TP+TN)*1.0/(TP+TN+FP+FN),1.0*(TP*TN-FP*FN)/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))); 
 }


	fprintf(logs,"  test_file: %s\n",test_file);
	fprintf(logs,"  predict_file: %s\n",predict_file);
	
	 
	time(&rawtime);timeinfo=localtime(&rawtime);
	fprintf(logs, "The current date/time is: %s\n", asctime (timeinfo) ); 
	fprintf(logs, "%s\n", "==========================================================="); 
	fclose(logs);	
	
	if (rmtempfile)
	{
		remove(test_file);
		remove(predict_file);
		remove(desc_file);
	}

	return -1;

}
