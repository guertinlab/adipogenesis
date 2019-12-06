#! /usr/bin/python
import sys
import getopt
def matrix(filename):
    infile=open(filename, 'r')
    while 1:
        line = infile.readline()
        if not line: break
        splitline = line.split()
        cluster = filename.split('.meme_output/meme.txt')[0].split('/')[-1]
        count = 0
        if line.startswith('MEME version'):
            meme_line = line
        if line.startswith('ALPHABET'):
            alphabet_line = line
        if line.startswith('strands:'):
            strands_line = line
        if line.startswith('A 0.'):
            bkg_freq = line
        if line.startswith('MOTIF'):
            if len(splitline) == 3:
                outfilename_prefix = splitline[2]
            else:
                outfilename_prefix = splitline[1]
            outfile=open(str(outfilename_prefix) + '_'+ str(cluster) +'_meme.txt', 'w')
            outfile.write(meme_line)
            outfile.write('\n')
            outfile.write(alphabet_line)
            outfile.write('\n')
            outfile.write(strands_line)
            outfile.write('\n')
            outfile.write(bkg_freq)
            outfile.write('\n')
            #outfile.write(line)
            outfile.write('%s\t%s\n'%(splitline[0], splitline[1]))
        if line.startswith('letter-probability matrix:'):# or 'log-odds matrix'):
            outfile.write(line)
            while 1:
                next = infile.readline()
                if len(next.split()) != 4: break
                outfile.write(next)
            outfile.write('\n')
            outfile.close()
    return

def main(argv):
    try:
        opts, args = getopt.getopt(argv, "hi:", ["help", "input="])
    except getopt.GetoptError, err:
        print str(err)
        sys.exit(2)
    name = False
    for opt, arg in opts:
        if opt in ('-i', '--input'):
            name = arg
        elif opt in ('-h', '--help'):
            print 'python ~/pyscripts/MEME_individual_from_db.py -i /Users/guertinlab/Desktop/tomtom_db/combined_motif_db.txt'
            sys.exit()
    if name:
        matrix(name)
    else:
        print 'python ~/pyscriptsMEME_individual_from_db.py -i /Users/guertinlab/Desktop/tomtom_db/combined_motif_db.txt'
if __name__ == "__main__":
    main(sys.argv[1:])


