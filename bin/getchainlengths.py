#!/usr/bin/env python3

import re
import sys, os, io
import tarfile, gzip
from Bio import PDB
from Bio.PDB import PDBParser, MMCIFParser, PDBIO
from Bio.PDB.PDBIO import Select
from optparse import OptionParser

description = (
  "List all chains along with their lengths encountered in a single, "
  "optionally gzipped, structure file (.pdb, .cif) or across "
  "all files in a .tar archive.")

def ParseArguments():
    """Parse command-line options.
    """
    parser = OptionParser(description=description)

    parser.add_option("-i", "--infile", dest="input",
              help="input, optionally gzipped, structure file (.pdb, .cif) or .tar "
                "archive containing structure files", metavar="FILE")

    (options, args) = parser.parse_args()

    if not (options.input):
        sys.stderr.write("ERROR: Input file is not provided.\n")
        sys.exit()

    if options.input and not os.path.isfile(options.input):
        sys.stderr.write("ERROR: Input file does not exist: " + options.input + "\n")
        sys.exit()

    return options



class ChainSelect(Select):
    def __init__(self, cid):
        self.cid = cid

    def accept_chain(self, chain):
        if chain.get_id() == self.cid:
            return 1
        else:
            return 0



def GetStructure(tfo, filename, extn, filestring):
    """read a structure from filehandle
    """
    erc = 1; structure = None
    try:
        if tfo is not None:
            try:
                with tfo.extractfile(filename) as tfefh:
                    tfefh = gzip.open(tfefh,'rt') if extn == '.gz' else io.TextIOWrapper(tfefh)
                    structure = PDBParser().get_structure('?', tfefh)
                    if len([_ for _ in structure.get_chains()]) < 1: raise ValueError
            except:# Exception as e:
                #print(e)
                with tfo.extractfile(filename) as tfefh:
                    tfefh = gzip.open(tfefh,'rt') if extn == '.gz' else io.TextIOWrapper(tfefh)
                    structure = MMCIFParser().get_structure('?', tfefh)
        else:
            try:
                fo = gzip.open(filename,'rt') if extn == '.gz' else open(filename)
                structure = PDBParser().get_structure('?', fo)
                if len([_ for _ in structure.get_chains()]) < 1:
                    fo.close()
                    raise ValueError
            except:
                fo = gzip.open(filename,'rt') if extn == '.gz' else open(filename)
                structure = MMCIFParser().get_structure('?', fo)
            finally:
                fo.close()
    except Exception as e:
        #print(e)
        sys.stderr.write("W: File ignored: " + filestring + '\n')
        erc = 0
    return [erc, structure]



if __name__ == "__main__":
    options = ParseArguments()

    basename = os.path.basename(options.input)
    dirname = os.path.dirname(options.input)
    name, extension = os.path.splitext(basename)

    tfentrydir = os.path.join(dirname, basename + '_ext')

    code = 0

    if extension == '.tar':
        with tarfile.open(options.input,'r') as tfo:
            for tfentry in tfo.getmembers():
                tfentn = tfentry.name
                tfentnmod = re.sub(r'[^\sa-zA-Z0-9_+-]', '_', tfentn)
                tfebname, tfeextn0 = os.path.splitext(os.path.basename(tfentn))
                tfebnam2, tfeextn1 = os.path.splitext(os.path.basename(tfebname))
                code, structure = GetStructure(tfo, tfentry, tfeextn0, options.input + ':' + tfentn)
                if not code: continue
                pdbcode = structure.header.get('idcode')
                if pdbcode == '': pdbcode = '_'
                for model in structure.get_models():
                    for chain in model.get_chains():
                        length = len([_ for _ in chain.get_residues() if PDB.is_aa(_,standard=False)])
                        print(f'D:  "{tfentnmod}"  {pdbcode}  {chain.id}  {length}  {model.serial_num}')

    else:
        code, structure = GetStructure(None, options.input, extension, options.input)

    if code and not extension == '.tar':
        pdbcode = structure.header.get('idcode')
        basenamemod = re.sub(r'[^\sa-zA-Z0-9_+-]', '_', basename)
        if pdbcode == '': pdbcode = '_'
        for model in structure.get_models():
            for chain in model.get_chains():
                length = len([_ for _ in chain.get_residues() if PDB.is_aa(_,standard=False)])
                print(f'D:  "{basenamemod}"  {pdbcode}  {chain.id}  {length}  {model.serial_num}')


