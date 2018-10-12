#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

It downloads the lastest available sequences of the chromosomes from the Ensembl database.



Author: Daniel Nicorici, Daniel.Nicorici@gmail.com

Copyright (c) 2009-2018 Daniel Nicorici

This file is part of FusionCatcher.

FusionCatcher is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

FusionCatcher is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with FusionCatcher (see file 'COPYING.txt').  If not, see
<http://www.gnu.org/licenses/>.

By default, FusionCatcher is running BLAT aligner
<http://users.soe.ucsc.edu/~kent/src/> but it offers also the option to disable
all its scripts which make use of BLAT aligner if you choose explicitly to do so.
BLAT's license does not allow to be used for commercial activities. If BLAT
license does not allow to be used in your case then you may still use
FusionCatcher by forcing not use the BLAT aligner by specifying the option
'--skip-blat'. Fore more information regarding BLAT please see its license.

Please, note that FusionCatcher does not require BLAT in order to find
candidate fusion genes!

This file is not running/executing/using BLAT.
"""



import os
import sys
import ftplib
import gzip
import socket
import optparse
from numpy import concatenate
import shutil


if __name__ == '__main__':

    #command line parsing

    usage = "%prog [options]"
    description = """It downloads the lastest available sequences of the chromosomes from the Ensembl database."""
    version = "%prog 0.10 beta"

    parser = optparse.OptionParser(usage=usage,description=description,version=version)

    parser.add_option("--organism",
                      action = "store",
                      type = "string",
                      dest = "organism",
                      default = "homo_sapiens",
                      help="""The name of the organism for which the chromosomes are downloaded, e.g. homo_sapiens, mus_musculus, etc. Default is '%default'.""")

    parser.add_option("--output",
                      action="store",
                      type="string",
                      dest="output_directory",
                      default = '.',
                      help="""The output directory where the chromosomes are stored. Default is '%default'.""")

    parser.add_option("--server",
                      action="store",
                      type="string",
                      dest="server",
                      default = "ftp.ensembl.org",
                      help="""The Ensembl server from where the chromosomes are downloaded. Default is '%default'.""")

    parser.add_option("--server-path",
                      action="store",
                      type="string",
                      dest="server_path",
                      default = "/pub/current_fasta/", #ftp://ftp.ensembl.org/pub/release-75/fasta/
                      help="""The path of Ensembl server from where the chromosomes are downloaded. Default is '%default'.""")

    (options,args) = parser.parse_args()

    # validate options
    if not (options.output_directory
            ):
        parser.print_help()
        sys.exit(1)


    # timeout in seconds
    timeout = 1800
    socket.setdefaulttimeout(timeout)

    #ftp.ensembl.org/pub/current_fasta/homo_sapiens/dna
    #
    #Files:
    #
    #Homo_sapiens.NCBI36.*.dna.chromosome.**.fa.gz
    #
    #where:
    #* is 54 (ensembl version)
    #** is chromosome name e.g. 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y,MT

    #ftp://ftp.ensembl.org/pub/current_fasta/mus_musculus/dna/
    #ftp://ftp.ensembl.org/pub/current_fasta/homo_sapiens/dna/
    #ftp://ftp.ensembl.org/pub/current_fasta/rattus_norvegicus/dna/

    url = '%s/%s/dna' % (options.server_path,options.organism.lower())
    ftp_path = ''
    print "Downloading the genome of organism '%s' from Ensembl!" % (options.organism.lower())
    try:
        ftp = ftplib.FTP(options.server)
        print ftp.login()
        ftp.cwd(url)
        ftp_path = ftp.pwd()
        list_files = ftp.nlst()

        fasta_name = '.dna.primary_assembly.fa'
        if options.organism.lower() != 'homo_sapiens':
            fasta_name = '.dna.toplevel.fa'

        list_files = [el for el in list_files if el.lower().startswith(options.organism.lower()) and (el.lower().find('.dna.chromosome.mt.fa.gz')>-1 or el.lower().endswith('.dna.chromosome.mito.fa.gz') or el.lower().find('%s.gz' %(fasta_name,))>-1)]

        if len(list_files) == 2:
            new_files = []
            for filename in list_files:
                #print "Downloading: %s/pub/current_fasta/%s/dna/%s" % (options.server,options.organism.lower(),filename)
                print "Downloading: "+url+"/"+filename
                nf = os.path.join(options.output_directory,filename)
                new_files.append(nf)
                fid = open(nf,'wb')
                ftp.retrbinary("RETR " + filename, fid.write)
                fid.close()
        else:
            print "ERROR: Cannot find the genome and mitochondria files!"
            print list_files
            sys.exit(1)
        ftp.close()
    except ftplib.all_errors, e:
        print 'FTP Error = ' + str(e)
        sys.exit(1)
    except Exception, e:
        print "Error: Generic exception!",str(e)
        sys.exit(1)

    new_list_files=[]
    for filename in new_files:
#        f = gzip.open(filename, 'rb')
#        file_content = f.read()
#        f.close()
        fn = filename.replace('.fa.gz','.fa')
        new_list_files.append(fn)
#        fod = open(fn,'wb')
#        fod.write(file_content)
#        fod.close()
        t = "gzip -d -f -c '%s' > '%s'" % (filename,fn)
        print "Executing:  ",t
        r = os.system(t)
        os.remove(filename)

    new_list_files = sorted(new_list_files)

    (first_dir,first_file) = os.path.split(new_list_files[0])
    info = first_file.split('.dna.chromosome.')[0].split('.')
    gv = ""
    if info:
        ni = len(info)
        if ni == 2:
            gv = info[1]
        elif ni > 2:
            gv = '.'.join(info[1:])
    txt = [
    "Ensembl database version: "+ftp_path.split('/')[2].replace('release-','')+"\n",
    "Organism: "+info[0].replace('_',' ')+"\n",
    "Genome version: "+gv+"\n"
    ]
    file(os.path.join(options.output_directory,'version.txt'),'w').writelines(txt)
    # file: version_ensembl.txt
    #
    # Ensembl version 60
    # Homo sapiens GRCh37/hg19 (February 2009 reference sequence)
    # example of filename: Homo_sapiens.GRCh37.60.dna.chromosome.4.fa

    #
    # Ensembl version 76
    # Homo sapiens GRCh38/hg38 (December 2013 reference sequence)
    # example of filename: Homo_sapiens.GRCh38.dna.chromosome.4.fa


    # find the FASTA file containing the MT and save it as mt.fa
    mt_filename = [el for el in new_list_files if el.lower().endswith('.dna.chromosome.mt.fa') or el.lower().endswith('.dna.chromosome.mito.fa')]
    mt_filename = mt_filename.pop(0)
    (mt_dir,mt_new_filename) = os.path.split(mt_filename)
    mt_new_filename = os.path.join(mt_dir,'mt.fa')
    shutil.move(mt_filename,mt_new_filename)

    #
    # find the FASTA file containing the genome and save it as genome.fa
    #
    genome_filename = [el for el in new_list_files if el.lower().endswith(fasta_name)]
    genome_filename = genome_filename.pop(0)
    (genome_dir,genome_new_filename) = os.path.split(genome_filename)
    file(os.path.join(options.output_directory,'genome_information.txt'),'w').write("%s/%s/%s" % (options.server,ftp_path,genome_new_filename))
    genome_new_filename = os.path.join(genome_dir,'genome.fa')
    shutil.move(genome_filename,genome_new_filename)
    #
