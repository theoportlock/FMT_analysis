# FMT_analysis
A slowly increasing collection of scripts for the analysis of data pertaining to the effectiveness of FMT for the treatment of liver disease

Download data from ena with the script using:

'''
wget https://github.com/enasequence/enaBrowserTools/archive/refs/tags/v1.6.zip
unzip v1.6.zip
'''

then, make sure that you have the python requests library installed and run:

'''
python enaBrowserTools-1.6/python3/enaGroupGet.py -f fastq PRJEB6337
'''
