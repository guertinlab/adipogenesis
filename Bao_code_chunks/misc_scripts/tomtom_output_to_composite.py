#! /usr/bin/python
import sys
import getopt
import os
import itertools
def tomtom_output(infilename):
    outfilename = infilename.split('.tomtom_output/tomtom.xml')[0].split('/')[-1]
    outfile = open(outfilename+'_test_index_pswm.txt', 'w')
    dic = {}
    order_dic ={}
    exp_name ='NULL'
    count = 0
    for line in file(infilename):
        splitline=line.split()
        if splitline[0]=="<queries>":
            count = 'query'
        if splitline[0]=="<targets>":
            count = -1
            print count
        if splitline[0]=="<motif":
            if count != 'query':
                count = count + 1
                print line
        if splitline[0]=="<pos":
            print splitline
            afreq = float(line.split('<pos A="')[1].split('" C=')[0])
            cfreq = float(line.split('" C="')[1].split('" G=')[0])
            gfreq = float(line.split('" G="')[1].split('" T=')[0])
            tfreq = float(line.split('" T="')[1].split('"/>')[0])
            outfile.write('%s\t%s\t%s\t%s\t%s\n'%(str(count), str(afreq), str(cfreq), str(gfreq), str(tfreq)))
    outfile.close()
    
    outfile2 = open(outfilename+'_test_index_rc_offset.txt', 'w')
    for line in file(infilename):
        splitline=line.split()
        if splitline[0]=="<target":
            print splitline
            idx = line.split('idx="')[1].split('" rc="')[0]
            rc = line.split('" rc="')[1].split('" off="')[0]
            offset = line.split('" off="')[1].split('" pv="')[0]
            outfile2.write('%s\t%s\t%s\n'%(str(idx), str(rc), str(offset)))
    outfile2.close()
    return



def main(argv):
    try:
        opts, args = getopt.getopt(argv, "hi:", ["help", "infile="])
    except getopt.GetoptError, err:
        print str(err)
        sys.exit(2)
    inf = False
    outf = False
    for opt, arg in opts:
        if opt in ('-i', '--infile'):
            inf = arg
        elif opt in ('-h', '--help'):
            print '-i\ninput file\n'
            sys.exit()
    if inf:
        tomtom_output(inf)
if __name__ == "__main__":
    main(sys.argv[1:])    
