#!/usr/bin/env python3

import datetime
import sys, os, io
import tarfile, gzip
import numpy as np
from Bio import PDB
from Bio.PDB import PDBParser, MMCIFParser, PDBIO
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
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

    if options.save1 and not (options.file1[0] and os.path.isfile(options.file1[0])):
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
    erc = 1; structure = None; header = {}
    try:
        if tfo is not None:
            try:
                with tfo.extractfile(filename) as tfefh:
                    tfefh = gzip.open(tfefh,'rt') if extn == '.gz' else io.TextIOWrapper(tfefh)
                    structure = PDBParser().get_structure('?', tfefh)
                    header = structure.header
                    if len([_ for _ in structure.get_chains()]) < 1: raise ValueError
            except:# Exception as e:
                #print(e)
                with tfo.extractfile(filename) as tfefh:
                    tfefh = gzip.open(tfefh,'rt') if extn == '.gz' else io.TextIOWrapper(tfefh)
                    structure = MMCIFParser().get_structure('?', tfefh)
                    tfefh.seek(0)
                    header = MMCIF2Dict(tfefh)
        else:
            try:
                fo = gzip.open(filename,'rt') if extn == '.gz' else open(filename)
                structure = PDBParser().get_structure('?', fo)
                header = structure.header
                if len([_ for _ in structure.get_chains()]) < 1:
                    fo.close()
                    raise ValueError
            except:
                fo = gzip.open(filename,'rt') if extn == '.gz' else open(filename)
                structure = MMCIFParser().get_structure('?', fo)
                fo.seek(0)
                header = MMCIF2Dict(fo)
            finally:
                fo.close()
    except Exception as e:
        #print(e)
        sys.stderr.write("E: File ignored: " + filestring + '\n')
        erc = 0
    return [erc, structure, header]


def GetStructureModelChain(fnamelst, chainid, modelnum):
    """
    read a structure file and extract the required chain and
    model from it
    """
    basename = os.path.basename(fnamelst[0])
    dirname = os.path.dirname(fnamelst[0])
    name, extension = os.path.splitext(basename)
    erc = 1; model = None; chain = None; header = {}

    try:
        if extension == '.tar':
            if len(fnamelst) < 2:
                raise MyException('E: Member not specified (should follow :) for .tar archive ' + fnamelst[0])
            with tarfile.open(fnamelst[0],'r') as tfo:
                tfentry = tfo.getmember(fnamelst[1])
                tfentn = tfentry.name
                tfebname, tfeextn0 = os.path.splitext(os.path.basename(tfentn))
                tfebnam2, tfeextn1 = os.path.splitext(os.path.basename(tfebname))
                code, structure, header = GetStructure(tfo, tfentry, tfeextn0, fnamelst[0] + ':' + tfentn)

        else:
            code, structure, header = GetStructure(None, fnamelst[0], extension, fnamelst[0])

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

    return [erc, chain, header]


def PrintHeader(fh, header):
    """
    write header to filehandle
    """
    if fh is not None:
        code = ''; name = ''; head = ''; date = ''
        cmol = ''; scie = ''; comm = ''; taxi = ''
        cmolalt = ''

        try:
            if 'idcode' in header: code = header['idcode']
            elif '_entry.id' in header: code = header['_entry.id'][0]
        except Exception as e: code = ''

        try:
            if 'head' in header: head = header['head']
            elif '_struct_keywords.pdbx_keywords' in header: head = header['_struct_keywords.pdbx_keywords'][0]
        except Exception as e: head = ''

        try:
            if 'name' in header: name = header['name']
            elif '_struct.title' in header: name = header['_struct.title'][0]
            elif '_entry.id' in header: name = header['_entry.id'][0]
        except Exception as e: name = ''

        try:
            if 'compound' in header:
                for nd_key, nd_value in header['compound'].items():
                    if 'molecule' in nd_value: cmol = nd_value['molecule']; break
                    if 'misc' in nd_value and nd_value['misc']: cmolalt = nd_value['misc']
            elif '_entity.pdbx_description' in header: cmol = header['_entity.pdbx_description'][0]
            if cmolalt and not cmol: cmol = cmolalt
        except Exception as e: cmol = ''

        try:
            if 'source' in header:
                for nd_key, nd_value in header['source'].items():
                    if 'organism_scientific' in nd_value: scie = nd_value['organism_scientific']; break
            elif '_entity_src_nat.pdbx_organism_scientific' in header: scie = header['_entity_src_nat.pdbx_organism_scientific'][0]
            elif '_entity_src_gen.pdbx_gene_src_scientific_name' in header: scie = header['_entity_src_gen.pdbx_gene_src_scientific_name'][0]
            elif '_pdbx_entity_src_syn.organism_scientific' in header: scie = header['_pdbx_entity_src_syn.organism_scientific'][0]
            elif '_ma_target_ref_db_details.organism_scientific' in header: scie = header['_ma_target_ref_db_details.organism_scientific'][0]
        except Exception as e: scie = ''

        try:
            if 'source' in header:
                for nd_key, nd_value in header['source'].items():
                    if 'organism_common' in nd_value: comm = nd_value['organism_common']; break
            elif '_entity_src_nat.common_name' in header: comm = header['_entity_src_nat.common_name'][0]
            elif '_entity_src_gen.gene_src_common_name' in header: comm = header['_entity_src_gen.gene_src_common_name'][0]
            elif '_pdbx_entity_src_syn.organism_common_name' in header: comm = header['_pdbx_entity_src_syn.organism_common_name'][0]
        except Exception as e: comm = ''

        try:
            if 'source' in header:
                for nd_key, nd_value in header['source'].items():
                    if 'organism_taxid' in nd_value: taxi = nd_value['organism_taxid']; break
            elif '_entity_src_gen.pdbx_gene_src_ncbi_taxonomy_id' in header: taxi = header['_entity_src_gen.pdbx_gene_src_ncbi_taxonomy_id'][0]
            elif '_pdbx_entity_src_syn.ncbi_taxonomy_id' in header: taxi = header['_pdbx_entity_src_syn.ncbi_taxonomy_id'][0]
            elif '_ma_target_ref_db_details.ncbi_taxonomy_id' in header: taxi = header['_ma_target_ref_db_details.ncbi_taxonomy_id'][0]
        except Exception as e: taxi = ''

        try:
            if 'deposition_date' in header: date = header['deposition_date']
            elif '_database_pdb_rev.date_original' in header: date = header['_database_pdb_rev.date_original'][0]
            elif '_pdbx_database_status.recvd_initial_deposition_date' in header: date = header['_pdbx_database_status.recvd_initial_deposition_date'][0]
            do = datetime.datetime.strptime(date, '%Y-%m-%d')
            date = datetime.date.strftime(do, '%d-%b-%y')
        except Exception as e: date = ''

        code = code.upper(); name = name.upper(); head = head.upper(); date = date.upper()
        cmol = cmol.upper(); scie = scie.upper(); comm = comm.upper(); taxi = taxi.upper()

        s = "HEADER    {myhead:<39} {mydate:<9} {mycode:<17}\n"
        fh.write(s.format(myhead = head[:39], mydate = date[:9], mycode = code[:17]))
        s1 = "TITLE     {myname:<68}\n"
        s2 = "TITLE    {myi} {myname:<67}\n"; i = 2
        fh.write(s1.format(myname = name[:68])); name = name[68:];
        while 0 < len(name):
            fh.write(s2.format(myi = i, myname = name[:67]));
            name = name[67:]; i = min(9, i + 1)
        s = "COMPND    MOL_ID: 1;" + (" "*58) + "\n"
        fh.write(s)
        s1 = "COMPND   {myi} MOLECULE: {mycmpd:<57}\n"; i = 2
        s2 = "COMPND   {myi} {mycmpd:<67}\n";
        if cmol:
            fh.write(s1.format(myi = i, mycmpd = cmol[:57])); cmol = cmol[57:]; i = i + 1
            while 0 < len(cmol):
                fh.write(s2.format(myi = i, mycmpd = cmol[:67]));
                cmol = cmol[67:]; i = min(9, i + 1)
        ##s = "COMPND   {myi} CHAIN: B;" + (" "*58) + "\n"
        ##fh.write(s.format(myi = i))
        s = "SOURCE    MOL_ID: 1;" + (" "*58) + "\n"
        fh.write(s)
        s1 = "SOURCE   {myi} ORGANISM_SCIENTIFIC: {myscie:<46}\n"; i = 2
        s2 = "SOURCE   {myi} {myscie:<67}\n";
        if scie:
            fh.write(s1.format(myi = i, myscie = scie[:46])); scie = scie[46:]; i = i + 1
            while 0 < len(scie):
                fh.write(s2.format(myi = i, myscie = scie[:67]));
                scie = scie[67:]; i = min(9, i + 1)
        s1 = "SOURCE   {myi} ORGANISM_COMMON: {mycomm:<50}\n";
        s2 = "SOURCE   {myi} {mycomm:<67}\n";
        if comm:
            fh.write(s1.format(myi = i, mycomm = comm[:50])); comm = comm[50:]; i = min(9, i + 1)
            while 0 < len(comm):
                fh.write(s2.format(myi = i, mycomm = comm[:67]));
                comm = comm[67:]; i = min(9, i + 1)
        s = "SOURCE   {myi} ORGANISM_TAXID: {mytaxi:<51}\n";
        if taxi: fh.write(s.format(myi = i, mytaxi = taxi[:51])); i = min(9, i + 1)
        s = "REMARK   1 SUPERPOSITION BY GTALIGN" + (" "*44) + "\n"
        fh.write(s)


if __name__ == "__main__":
    options = ParseArguments()

    code = 1

    if options.save1:
        code, chain1, header1 = GetStructureModelChain(options.file1, options.cid1, options.mod1)

    if not code: sys.exit(1)

    code, chain2, header2 = GetStructureModelChain(options.file2, options.cid2, options.mod2)

    if not code: sys.exit(1)

    if options.rot and options.trl:
        trl = np.array(options.trl)
        rot = np.array(options.rot)
        rot = rot.reshape(3, 3, order='F')
        ##print(rot.astype(np.float32));print(trl.astype(np.float32))
        if options.referenced:
            for atom in chain2.get_atoms(): atom.transform(rot.astype(np.float32), trl.astype(np.float32))
        elif options.save1:
            for atom in chain1.get_atoms(): atom.transform(rot.astype(np.float32), trl.astype(np.float32))

    builder = PDB.StructureBuilder.StructureBuilder()
    builder.init_structure('Sup')
    builder.init_model(0,0)

    outstr = builder.get_structure()

    if options.save1: chain1.id = 'A'
    chain2.id = 'B'

    if options.save1: outstr[0].add(chain1)
    outstr[0].add(chain2)

    stream = io.StringIO()

    pio = PDBIO()
    pio.set_structure(outstr)
    pio.save(stream)
    ##pio.save(options.output)

    ##print(header2)

    with open(options.output,'w') as of:
        PrintHeader(of, header2)
        of.write(stream.getvalue())

    stream.close()

