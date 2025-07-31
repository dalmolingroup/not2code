#!/usr/bin/env python
#################################################################
#                                                               
#  PLEK - predictor of lncRNAs and mRNAs based on k-mer scheme  
#  Authors: Aimin Li, Junying Zhang   
#  Contacts: LiAiminMail@gmail.com, jyzhang@mail.xidian.edu.cn  
#  Webcite: https://sourceforge.net/projects/plek/              
#  Updated on: June 19, 2014                                                 
#                                                               
#################################################################

__all__ = ['']
import os, sys, traceback, getpass, time


if __name__ == '__main__':
    def exit_with_help():
        print("""\
usage: python PLEK_setup.py
           """)
        sys.exit(1)      
    
    def compile_c(): # re-compile source
        
        os.system("g++ -c svm.cpp")
        os.system("gcc -g -Wall svm-train.c svm.o -o svm-train -lstdc++ -lm ")
        os.system("gcc -g -Wall svm-predict.c svm.o -o svm-predict  -lstdc++ -lm ")
        os.system("gcc -g -Wall svm-scale.c svm.o -o svm-scale  -lstdc++ -lm ")
        
        os.system("gcc -g -Wall PLEK_main.c -o PLEK -lm")
        os.system("gcc -g -Wall PLEK_spsn.c -o PLEK_spsn -lm ")

        dirname = os.path.dirname(__file__)
        model_pathname = os.path.join(dirname, './PLEK.model')
        if not os.path.exists(model_pathname):
            print("Warning: PLEK.model does not exist.")
            #os.system("cat PLEK.model0 PLEK.model1 PLEK.model2 > PLEK.model")

    try:
        print("Compiling source ...")
        compile_c()       
        print("""\
Compilation complete.
           """)
        sys.exit(1)
        
    except (IOError,ValueError) as e:
        sys.stderr.write(str(e) + '\n')
        
        sys.exit(1)
