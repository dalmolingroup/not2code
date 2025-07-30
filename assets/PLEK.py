#!/usr/bin/env python
#######################################################
#                                                                                      
#  PLEK - predictor of lncRNAs and mRNAs based on k-mer scheme
#  Authors: Aimin Li, JunYing Zhang                                                    
#  Contacts: LiAiminMail@gmail.com, jyzhang@mail.xidian.edu.cn
#  Webcite: https://sourceforge.net/projects/plek/
#  Version: 1.2
#  Updated on: June 19, 2014
#                                                                                      
#######################################################

#usage: python PLEK.py -fasta fasta_file -out prefix_of_output_file -thread number_of_threads

#example: python PLEK.py -fasta PLEK_test.fa -out myfile -thread 10
#output: myfile

__all__ = ['']
import os, sys, traceback, getpass, time, re
from threading import Thread
from subprocess import *

class SVMScaleThread(Thread):  # class for svm-scale
        def __init__(self, range_file, svm_origin, svm_scaled, script_dir):
                Thread.__init__(self)
                self.range_file = range_file
                self.svm_origin = svm_origin
                self.svm_scaled = svm_scaled
                self.script_dir = script_dir

        def run(self):
                cmdline= (self.script_dir +'svm-scale -r {0} {1} > {2} ').format\
                                (self.range_file, self.svm_origin, self.svm_scaled)
                os.system(cmdline)

class SVMPredictThread(Thread):   # class for svm-predict
        def __init__(self, svm_scaled, svm_model, predicted, script_dir):
                Thread.__init__(self)
                self.svm_scaled = svm_scaled
                self.svm_model = svm_model
                self.predicted = predicted
                self.script_dir = script_dir

        def run(self):
                cmdline= (self.script_dir + 'svm-predict {0} {1} {2} ').format\
                                (self.svm_scaled,self.svm_model,self.predicted)
                os.system(cmdline)


class GridOption: # set or get input parameters
        def __init__(self, options):
                dirname = os.path.dirname(__file__)
                self.svmtrain_pathname = os.path.join(dirname, './svm-train')
                self.svmpredict_pathname = os.path.join(dirname, './svm-predict')
                self.svmscale_pathname = os.path.join(dirname, './svm-scale')
                self.pos_file = ""
                self.neg_file = ""
                self.prefix = "plek_output_"
                self.svmrangefile = "PLEK.range"
                self.is_posneg_balanced = 0
                self.thread_count = 5
                self.modelfile="PLEK.model"
                self.kmer=5
                self.min_seq_length=200
                self.unkown=0
                self.is_recompile=0
                self.input_type=1
                self.isoutmsg=0
                self.isrmtempfile=1
                self.script_dir="./"
                self.parse_options(options)

        def parse_options(self, options):
                args=options # save options to args
                if type(options) == str:
                        options = options.split()
                i = 0
                pass_through_options = []
                # get python script file's name, args[0]
                # determine it is with dir or only file_name
                script_name=args[0]
                #print strHello
                if script_name.rfind('/')>=0:
                        self.script_dir=script_name[:script_name.rfind('/')+1]

                self.svmtrain_pathname = self.script_dir + 'svm-train'
                self.svmpredict_pathname = self.script_dir + 'svm-predict'
                self.svmscale_pathname = self.script_dir + 'svm-scale'
                self.svmrangefile = self.script_dir +"PLEK.range"
                self.modelfile = self.script_dir +"PLEK.model"

                while i < len(options):
                        if options[i] == '-fasta': # The name of a fasta file, its sequences are to be predicted.
                                i = i + 1
                                self.pos_file = options[i]
                                self.unkown=1
                        elif options[i] == '-thread': # The number of threads for running the PLEK program. The bigger this number is, the faster PLEK runs.
                                i = i + 1
                                self.thread_count = options[i]
                        elif options[i] == '-out': # The file name for the results of prediction. Predicted positive samples are labeled as "Coding", and negative as "Non-coding".
                                i = i + 1
                                self.prefix = options[i]
                        elif options[i] == '-minlength': # The minimum length of sequences. The sequences whose lengths are more than minlength will be processed.
                                i = i + 1
                                self.min_seq_length = options[i]
                        elif options[i] == '-isoutmsg': # Output messages to stdout(screen) or not. "0" means that PLEK be run quietly. "1" means that PLEK outputs the details of processing.
                                i = i + 1
                                self.isoutmsg = options[i]
                        elif options[i] == '-isrmtempfile': # Remove temporary files or not. "0" means that PLEK retains temporary files. "1" means that PLEK remove temporary files.
                                i = i + 1
                                self.isrmtempfile = options[i]
                        elif options[i] == '-pos': # positive class.
                                i = i + 1
                                self.pos_file = options[i]
                        elif options[i] == '-neg': # negative class.
                                i = i + 1
                                self.neg_file = options[i]
                        elif options[i] == '-range': # svm range file.
                                i = i + 1
                                self.svmrangefile = options[i]
                        elif options[i] == '-k': # range of k.
                                i = i + 1
                                self.kmer = options[i]
                        elif options[i] == '-model': # svm model file.
                                i = i + 1
                                self.modelfile = options[i]
                        elif options[i] == '-balance':
                                self.is_posneg_balanced = 1 #  NOTE: -b , need to balance; NO -b, not balance.
                        elif options[i] == '-isrecompile': # re-compile source.
                                i = i + 1
                                self.is_recompile = options[i]
                        else:
                                pass_through_options.append(options[i])
                        i = i + 1

                self.pass_through_string = ' '.join(pass_through_options)


if __name__ == '__main__':
        def exit_with_help():
                print("""\
=====================
  USAGE AND EXAMPLES
=====================
Usage:
python PLEK.py -fasta fasta_file -out output_file -thread number_of_threads
   -minlength min_length_of_sequence -isoutmsg 0_or_1 -isrmtempfile 0_or_1

   -fasta               The name of a fasta file, its sequences are to be predicted.

   -out           The file name for the results of prediction. Predicted positive
                                 samples are labeled as "Coding", and negative as "Non-coding".

   -thread         (Optional) The number of threads for running the PLEK program.
                                  The bigger this number is, the faster PLEK runs. Default value: 5.

   -minlength   (Optional) The minimum length of sequences. The sequences whose
                                  lengths are more than minlength will be processed. Default
                                 value: 200.

   -isoutmsg     (Optional) Output messages to stdout(screen) or not. "0" means
                                 that PLEK be run quietly. "1" means that PLEK outputs the details
                                 of processing. Default value: 0.

   -isrmtempfile (Optional) Remove temporary files or not. "0" means that PLEK
                                  retains temporary files. "1" means that PLEK remove temporary
                                  files. Default value: 1.



Examples:
1. $ python PLEK.py -fasta PLEK_test.fa -out predicted -thread 10

   NOTE: To predict the sequences in the 'PLEK_test.fa' file, run the PLEK program with
   10 threads. The program outputs the predicted sequences  in the file 'predicted'.


2. $ python PLEK.py -fasta PLEK_test.fa -out predicted -thread 10 -minlength 300

   NOTE: To predict the sequences in the 'PLEK_test.fa' file, run the PLEK program with
   10 threads. The program outputs the predicted sequences  in the file 'predicted'.
   The sequences with the length of >300nt will be processed (remained).


3. $ python PLEK.py -fasta PLEK_test.fa -out predicted -thread 10 -isrmtempfile 0

   NOTE: To predict the sequences in the 'PLEK_test.fa' file, run the PLEK program with
   10 threads. The program outputs the predicted sequences  in the file 'predicted'.
   The details of PLEK run will output to the files with "predicted" as prefix.


4. $ python PLEK.py -fasta PLEK_test.fa -out predicted -thread 10 -isoutmsg 1

   NOTE: To predict the sequences in the 'PLEK_test.fa' file, run the PLEK program with
   10 threads. The program outputs the predicted sequences  in the file 'predicted'.
   The details of PLEK run will output to user's screen(stdio).

=====================
  CONTACTS
=====================
Aimin Li: LiAiminMail@gmail.com
Junying Zhang: jyzhang@mail.xidian.edu.cn

======================
  WEBSITE
=====================
https://sourceforge.net/projects/plek/
                   """)
                sys.exit(1)

        def file_id_by_lineid(totalcount,filecount,n): # for splitting input files.
                # totalcount,  total line count
                #filecount,  the number of files
                #n , current row id
                countperfile=totalcount/filecount
                for i in range(1,filecount):
                        if n>(i-1)*countperfile and n<=i*countperfile:
                                return i
                return filecount

        def compile_c(_options): # re-compile source
                print('[{0}] Compiling svm, svm-train, svm-predict, svm-scale'.format(time.strftime('%Y-%m-%d %H:%M:%S')))
                os.system("g++ -c " + _options.script_dir + "svm.cpp -o " + _options.script_dir + "svm.o")
                os.system("LNAG=C gcc -g -Wall " + _options.script_dir + "svm-train.c " + _options.script_dir + "svm.o -o " + _options.script_dir + "svm-train -lstdc++ -lm")
                os.system("LNAG=C gcc -g -Wall " + _options.script_dir + "svm-predict.c " + _options.script_dir + "svm.o -o " + _options.script_dir + "svm-predict  -lstdc++ -lm")
                os.system("LNAG=C gcc -g -Wall " + _options.script_dir + "svm-scale.c " + _options.script_dir + "svm.o -o " + _options.script_dir + "svm-scale  -lstdc++ -lm")

                print('[{0}] Compiling PLEK_main, PLEK_spsn'.format(time.strftime('%Y-%m-%d %H:%M:%S')))
                os.system("LNAG=C gcc -g -Wall " + _options.script_dir + "PLEK_main.c -o " + _options.script_dir + "PLEK -lm")
                os.system("LNAG=C gcc -g -Wall " + _options.script_dir + "PLEK_spsn.c -o " + _options.script_dir + "PLEK_spsn -lm ")


        if len(sys.argv) < 2:
                exit_with_help()
        options = sys.argv
        try:
                print('[{0}] Beginning PLEK run (Version 1.2) '.format(time.strftime('%Y-%m-%d %H:%M:%S')))
                # get input options
                _options = GridOption(options);
                cmdline=None;

                # recompile source
                if _options.is_recompile or (not os.path.isfile(_options.script_dir + 'PLEK')):
                        compile_c(_options)

                # check if the model file exists
                if not os.path.isfile( _options.modelfile):
                        print('Building model')
                        os.system("cat " + _options.script_dir + "PLEK.model0 " + _options.script_dir + "PLEK.model1 " + _options.script_dir + "PLEK.model2 > "  +  _options.modelfile )
                        if not os.path.isfile( _options.modelfile):
                                print("ERROR: No such file '" + _options.modelfile + "'")
                                sys.exit(1)

                svm_file=str(_options.prefix)+"_allsvm";
                print('[{0}] PLEK is running'.format(time.strftime('%Y-%m-%d %H:%M:%S')))

                # calculate k-mer usage frequencies
                print('[{0}] Calculating k-mer usage'.format(time.strftime('%Y-%m-%d %H:%M:%S')))
                if _options.pos_file!="" and _options.neg_file!="" and _options.is_posneg_balanced==1:  # for modeling
                        _options.input_type=1
                        if not os.path.isfile(_options.pos_file):
                                print("ERROR: No such file '" + _options.pos_file + "'")
                                sys.exit(1)
                        if not os.path.isfile(_options.neg_file):
                                print("ERROR: No such file '" + _options.neg_file + "'")
                                sys.exit(1)
                        cmdline = (_options.script_dir + 'PLEK -s 1 -d 5 -p {0} -n {1} -o {2} -k {3} -l {4}  -b -isoutmsg {5} -isrmtempfile {6}').format\
                                (_options.pos_file, _options.neg_file, _options.prefix, _options.kmer, _options.min_seq_length, _options.isoutmsg, _options.isrmtempfile)
                if _options.pos_file!="" and _options.neg_file!="" and _options.is_posneg_balanced==0: # for known Pos & Neg
                        _options.input_type=2
                        if not os.path.isfile(_options.pos_file):
                                print("ERROR: No such file '" + _options.pos_file + "'")
                                sys.exit(1)
                        if not os.path.isfile(_options.neg_file):
                                print("ERROR: No such file '" + _options.neg_file + "'")
                                sys.exit(1)
                        cmdline = (_options.script_dir + 'PLEK   -s 1 -d 5 -p {0} -n {1} -o {2} -k {3} -l {4}  -b -isoutmsg {5} -isrmtempfile {6}').format\
                                (_options.pos_file, _options.neg_file, _options.prefix, _options.kmer, _options.min_seq_length, _options.isoutmsg, _options.isrmtempfile)
                if _options.pos_file!="" and _options.neg_file=="" :  # for known Pos; for unkown Neg
                        _options.input_type=3
                        if _options.unkown==1:
                                _options.input_type=4
                        if not os.path.isfile(_options.pos_file):
                                print("ERROR: No such file '" + _options.pos_file + "'")
                                sys.exit(1)
                        cmdline = (_options.script_dir + 'PLEK   -s 1 -d 5 -p {0} -o {1}  -k {2} -l {3}  -isoutmsg {4} -isrmtempfile {5}').format\
                                (_options.pos_file, _options.prefix, _options.kmer, _options.min_seq_length, _options.isoutmsg, _options.isrmtempfile)
                if _options.pos_file=="" and _options.neg_file!="" :  # for known Neg
                        _options.input_type=5
                        if not os.path.isfile(_options.neg_file):
                                print("ERROR: No such file '" + _options.neg_file + "'")
                                sys.exit(1)
                        cmdline = (_options.script_dir + 'PLEK   -s 1 -d 5  -n {0} -o {1} -k {2} -l {3}  -isoutmsg {4} -isrmtempfile {5}').format\
                                (_options.neg_file, _options.prefix, _options.kmer, _options.min_seq_length, _options.isoutmsg, _options.isrmtempfile)

                os.system(cmdline)

                # split file (input, file of k-mer usage frequencies )

                #total number of rows
                count = -1
                for count, line in enumerate(open(svm_file, 'rU')):
                        pass
                        count += 1
                #
                file_count=int(_options.thread_count)
                #print('          Number of sequence: {0}'.format(count))

                file_array=[]
                for fn in range(1,file_count+1):
                        fp=open(_options.prefix+'_temp_'+str(fn),'w')
                        file_array.append(fp)
                #
                if not os.path.isfile(svm_file):
                        print("ERROR: No such file '" + svm_file + "'")
                        sys.exit(1)
                fv=open(svm_file,'r')
                n=1
                for line in fv:
                        file_id=file_id_by_lineid(count,file_count,n)
                        file_p=file_array[file_id-1]
                        file_p.write(line)
                        n+=1

                for fn in range(1,file_count+1):
                        file_p=file_array[fn-1]
                        file_p.close()
                fv.close()

                # svm-scale
                print('[{0}] Scaling data'.format(time.strftime('%Y-%m-%d %H:%M:%S')))
                scale_array=[]
                for fn in range(1,file_count+1):
                        scale_obj=SVMScaleThread(_options.svmrangefile, _options.prefix+'_temp_'+str(fn), _options.prefix+'_temp_'+str(fn)+'_scaled', _options.script_dir)
                        scale_array.append(scale_obj);
                        scale_obj.start()

                # svm-predict
                for fn in range(1,file_count+1):
                        scale_obj=scale_array[fn-1]
                        scale_obj.join()

                # remove temporary files
                if _options.isrmtempfile==1:
                        for fn in range(1,file_count+1):
                                os.remove(_options.prefix+'_temp_'+str(fn))

                print('[{0}] Predicting'.format(time.strftime('%Y-%m-%d %H:%M:%S')))
                predict_array=[]
                for fn in range(1,file_count+1):
                        predict_obj=SVMPredictThread(_options.prefix+'_temp_'+str(fn)+'_scaled', _options.modelfile, _options.prefix+'_temp_'+str(fn)+'_predicted', _options.script_dir)
                        predict_array.append(predict_obj);
                        predict_obj.start()

                # merge:  predicted file
                for fn in range(1,file_count+1):
                        predict_obj=predict_array[fn-1]
                        predict_obj.join()

                # remove temporary files
                if _options.isrmtempfile==1:
                        for fn in range(1,file_count+1):
                                os.remove(_options.prefix+'_temp_'+str(fn)+'_scaled')

                f_final=open(_options.prefix+"_predicted",'w')

                for fn in range(1,file_count+1):
                        f_predicted=open(_options.prefix+'_temp_'+str(fn)+'_predicted')
                        for line in f_predicted:
                                f_final.write(line)
                        f_predicted.close()
                        f_final.flush()

                f_final.close()

                # compare svm_file & predicted file
                if _options.input_type==2 or _options.input_type==3 or _options.input_type==5 or _options.input_type==4:
                        cmdline=(_options.script_dir + 'PLEK_spsn  -svm {0} -predict {1} -desc {2} -descclass {3} -output {4} -input_type {5}  -isoutmsg {6} -isrmtempfile {7}').format\
                                (svm_file, _options.prefix+"_predicted", _options.prefix+"_allsvmdesc", _options.prefix+"_result", _options.prefix+"_logs", _options.input_type, _options.isoutmsg, _options.isrmtempfile)
                        os.system(cmdline);

                if _options.isrmtempfile==1:
                        os.remove(_options.prefix+"_logs")
                os.rename(_options.prefix+"_result", _options.prefix)

                # remove temporary files
                if _options.isrmtempfile==1:
                        for fn in range(1,file_count+1):
                                os.remove(_options.prefix+'_temp_'+str(fn)+'_predicted')

                print('[{0}] Run complete'.format(time.strftime('%Y-%m-%d %H:%M:%S')))
                print(' Result file: {0}'.format(_options.prefix))

                # statistics
                Total_count=0
                Noncoding_count=0
                for line in open(_options.prefix, 'r'):
                        Total_count=Total_count+1
                        #if line.find("Non-coding")>=0:
                        if re.compile(r'^Non-coding').match(line):
                                Noncoding_count=Noncoding_count+1
                print(' Coding: {0}/{1}={2}%, Non-coding: {3}/{4}={5}%'.format(
                        Total_count-Noncoding_count, Total_count, 1.0*(Total_count-Noncoding_count)/Total_count*100,
                        Noncoding_count, Total_count, 1.0*(Noncoding_count)/Total_count*100))


                sys.exit(1)

        except (IOError,ValueError) as e:
                sys.stderr.write(str(e) + '\n')
                sys.stderr.write('Try "python {0}" for more information.\n'.format(sys.arg[0]))
                sys.exit(1)