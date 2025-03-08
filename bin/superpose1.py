#!/usr/bin/env python3

import sys, os, io
import tarfile, gzip
import numpy as np
from Bio import PDB
from Bio.PDB import PDBParser, MMCIFParser, PDBIO
from Bio.PDB.PDBIO import Select
from optparse import OptionParser

class MyException(Exception):
    pass

description = (
  "Superimpose a pair of 3D structure (model) chains")

def get_archive_file(option, opt, value, parser):
    setattr(parser.values, option.dest, value.rsplit(':',1))

def get_list(option, opt, value, parser):
    setattr(parser.values, option.dest, value.split(','))

def ParseArguments():
    """Parse command-line options.
    """
    parser = OptionParser(description=description)

    parser.add_option("--i1", "--file1", type=str,
              action='callback', callback=get_archive_file, dest='file1',
              help="First input, optionally gzipped, structure file (.pdb, .cif) or " \
                ".tar archive containing structure files", metavar="FILE")
    parser.add_option("--i2", "--file2", type=str,
              action='callback', callback=get_archive_file, dest='file2',
              help="Second input, optionally gzipped, structure file (.pdb, .cif) or " \
                ".tar archive containing structure files", metavar="FILE")
    parser.add_option("--c1", "--chain1", dest="cid1", default="",
              help="First structure's chain id; " \
                "by default, the first chain is selected", metavar="ID")
    parser.add_option("--c2", "--chain2", dest="cid2", default="",
              help="Second structure's chain id; " \
                "by default, the first chain is selected", metavar="ID")
    parser.add_option("--m1", "--model1", dest="mod1", type=int, default=-999999,
              help="First structure's model number (optional)", metavar="NUMBER")
    parser.add_option("--m2", "--model2", dest="mod2", type=int, default=-999999,
              help="Second structure's model number (optional)", metavar="NUMBER")
    parser.add_option("-r", "--rotation", type=str,
              action='callback', callback=get_list, dest='rot',
              help="Rotation matrix in row-major order with elements separated by a " \
                "comma", metavar="MATRIX")
    parser.add_option("-t", "--translation", type=str,
              action='callback', callback=get_list, dest='trl',
              help="Translation vector with elements separated by a comma",
              metavar="VECTOR")
    parser.add_option("-2", "--referenced",
              action="store_true", dest="referenced", default=False,
              help="Apply transformation to the second structure")
    parser.add_option("-s", "--save1",
              action="store_true", dest="save1", default=False,
              help="Include the first structure in the output")
    parser.add_option("-o", "--outfile", dest="output",
              help="Output file in PDB format", metavar="FILE")

    (options, args) = parser.parse_args()

    if not options.file1[0] or not os.path.isfile(options.file1[0]):
        sys.stderr.write("E: First input file not found: " + options.file1[0] + "\n")
        sys.exit(1)

    if not options.file2[0] or not os.path.isfile(options.file2[0]):
        sys.stderr.write("E: Second input file not found: " + options.file2[0] + "\n")
        sys.exit(1)

    if not options.output:
        sys.stderr.write("E: Output file is not provided.\n")
        sys.exit(1)

    if options.rot and len(options.rot) != 9:
        sys.stderr.write("E: Rotation matrix should consist of 9 elements.\n")
        sys.exit(1)

    if options.trl and len(options.trl) != 3:
        sys.stderr.write("E: Translation vector should consist of 3 elements.\n")
        sys.exit(1)

    if (options.rot or options.trl) and not (options.rot and options.trl):
        sys.stderr.write("E: Both rotation and translation should be provided.\n")
        sys.exit(1)

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
        sys.stderr.write("E: File ignored: " + filestring + '\n')
        erc = 0
    return [erc, structure]


def GetStructureModelChain(fnamelst, chainid, modelnum):
    """
    read a structure file and extract the required chain and
    model from it
    """
    basename = os.path.basename(fnamelst[0])
    dirname = os.path.dirname(fnamelst[0])
    name, extension = os.path.splitext(basename)
    erc = 1; model = None; chain = None

    try:
        if extension == '.tar':
            if len(fnamelst) < 2:
                raise MyException('E: Member not specified (should follow :) for .tar archive ' + fnamelst[0])
            with tarfile.open(fnamelst[0],'r') as tfo:
                tfentry = tfo.getmember(fnamelst[1])
                tfentn = tfentry.name
                tfebname, tfeextn0 = os.path.splitext(os.path.basename(tfentn))
                tfebnam2, tfeextn1 = os.path.splitext(os.path.basename(tfebname))
                code, structure = GetStructure(tfo, tfentry, tfeextn0, fnamelst[0] + ':' + tfentn)

        else:
            code, structure = GetStructure(None, fnamelst[0], extension, fnamelst[0])

        if not code: raise MyException('E: Reading structure failed.')

        if modelnum > -99999:
            for mod in structure.get_models():
                if mod.serial_num == modelnum: model = mod; break
        if model is None:
            model = structure.get_list()[0] #1st model
        if chainid:
            for chn in model.get_chains():
                if chn.get_id() == chainid: chain = chn; break
        else: chain = model.get_list()[0] #1st chain

        if chain is None: raise MyException('E: Chain unidentified: ' + fnamelst[0])

        chain = chain.copy()

    except Exception as e:
        ##if hasattr(e,'message'):
        sys.stderr.write(str(e) + '\n')
        erc = 0

    return [erc, chain]


if __name__ == "__main__":
    options = ParseArguments()

    code = 0

    code, chain1 = GetStructureModelChain(options.file1, options.cid1, options.mod1)

    if not code: sys.exit(1)

    code, chain2 = GetStructureModelChain(options.file2, options.cid2, options.mod2)

    if not code: sys.exit(1)

    if options.rot and options.trl:
        trl = np.array(options.trl)
        rot = np.array(options.rot)
        rot = rot.reshape(3, 3, order='F')
        ##print(rot.astype(np.float32));print(trl.astype(np.float32))
        if options.referenced:
            for atom in chain2.get_atoms(): atom.transform(rot.astype(np.float32), trl.astype(np.float32))
        else:
            for atom in chain1.get_atoms(): atom.transform(rot.astype(np.float32), trl.astype(np.float32))

    builder = PDB.StructureBuilder.StructureBuilder()
    builder.init_structure('Sup')
    builder.init_model(0,0)

    outstr = builder.get_structure()

    chain1.id = 'A'
    chain2.id = 'B'

    if options.save1: outstr[0].add(chain1)
    outstr[0].add(chain2)

    io = PDBIO()
    io.set_structure(outstr)
    io.save(options.output)

