#! /usr/bin/python
import sys
import getopt
def matrix(filename, outfilename_prefix):
    infile=open(filename, 'r')
    nmotifs = 0
    outfile=open(str(outfilename_prefix) + '_meme.txt', 'w')
    outfile.write('%s\n'%("MEME version 4"))
    outfile.write('\n')
    outfile.write('%s\n'%("ALPHABET= ACGT"))
    outfile.write("\n")
    outfile.write('%s\n'%("strands: + -"))
    outfile.write("\n")
    outfile.write('%s\n'%("A 0.288 C 0.212 G 0.212 T 0.288"))
    count = 0 
    while 1:
        line = infile.readline()
        if not line: break
        splitline = line.split()
        directory = filename.split('meme.txt')[0]
        if line.startswith('>'):
            count = count + 1
            outfile.write("\n")
            namemotif = str.split(splitline[1], "(")[0]
            print namemotif
            outfile.write('%s%s%s%s\n'%("MOTIF ", splitline[1], ' ', namemotif))
            str(len(splitline[0]) - 1)
            outfile.write('%s\n'%("letter-probability matrix: alength= 4"))
        else:
            outfile.write(line)
    print count
    outfile.close()
    return

def main(argv):
    try:
        opts, args = getopt.getopt(argv, "hi:o:", ["help", "input=", "out="])
    except getopt.GetoptError, err:
        print str(err)
        sys.exit(2)
    name = False
    outname = False
    for opt, arg in opts:
        if opt in ('-i', '--input'):
            name = arg
        if opt in ('-o', '--out'):
            outname = arg
        elif opt in ('-h', '--help'):
            print 'python /Users/guertinlab/pyscripts/HOMER_MEME_conversion.py -i /Users/guertinlab/Desktop/custom.motifs -o custom.motifs'
            sys.exit()
    if name and outname:
        matrix(name, outname)
if __name__ == "__main__":
    main(sys.argv[1:])


