# Copyright 2022 Informatics Matters Ltd.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

from argparse import ArgumentError
import gzip
from rdkit import Chem
import utils

class SdfWriter:

    def __init__(self, outfile):
        if outfile.endswith('.gz'):
            self.gzip = gzip.open(outfile, 'wt')
            self.writer = Chem.SDWriter(self.gzip)
        else:
            self.gzip = None
            self.writer = Chem.SDWriter(outfile)

    def write(self, smiles=None, mol=None, mol_id=None, existing_props=None, prop_names=None, new_props=None, smiles_prop_name=None):
        if not mol:
            mol = Chem.MolFromSmiles(smiles)
        if mol_id is not None:
            mol.SetProp('_Name', mol_id)
        if smiles_prop_name is not None:
            mol.SetProp(smiles_prop_name, smiles)
        if not prop_names:
            prop_names = []
        if not new_props:
            new_props = []
            
        for i, prop_name in enumerate(prop_names):
            value = new_props[i]
            if prop_name is not None:
                mol.SetProp(prop_name, str(value))

        self.writer.write(mol)

    def write_header(self, values):
        utils.log("INFO: asked to write header for an SDF. No action will be taken.")

    def close(self):
        self.writer.close()
        if self.gzip:
            self.gzip.close()


class SmilesWriter:

    def __init__(self, outfile, sep):
        self.writer = open(outfile, 'w')
        if sep is None:
            self.sep = ' '
        else:
            self.sep = sep

    def write_header(self, values):
        line = self.sep.join(values)
        self.writer.write(line + "\n")

    def write(self, smiles=None, mol=None, mol_id=None, existing_props=None, prop_names=None, new_props=None, smiles_prop_name=None):
        values = [smiles]
        if not existing_props:
            existing_props = []
        if not prop_names:
            prop_names = []
        if not new_props:
            new_props = []
            
        for prop in existing_props:
            if prop is not None:
                values.append(prop)
            else:
                values.append('')

        for prop in new_props:
            if prop is not None:
                values.append(str(prop))
            else:
                values.append('')
        line = self.sep.join(values)
        self.writer.write(line + "\n")

    def close(self):
        self.writer.close()


class SdfReader:

    def __init__(self, input_file, id_col, recs_to_read):

        self.field_names = []
        # read a number of records to determine the field names
        if recs_to_read:
            r = self.create_reader(input_file)
            for i in range(recs_to_read):
                if r.atEnd():
                    break
                else:
                    try:
                        mol = next(r)
                        if mol:
                            names = mol.GetPropNames()
                            for name in names:
                                if name not in self.field_names:
                                    self.field_names.append(name)
                    except StopIteration:
                        break

        # now create the real reader
        self.reader = self.create_reader(input_file)
        self.id_col = id_col


    @staticmethod
    def create_reader(input_file):
        if input_file.endswith('.gz'):
            reader = Chem.ForwardSDMolSupplier(gzip.open(input_file))
        else:
            reader = Chem.ForwardSDMolSupplier(input_file)
        return reader

    def read(self):
        try:
            props = []
            mol = next(self.reader)
            if not mol:
                raise TypeError(f'{self}: Error parsing molecule')
            smi = Chem.MolToSmiles(mol)
                
            if self.id_col:
                id = mol.GetProp(self.id_col)
            else:
                id = None

            for name in self.field_names:
                if mol.HasProp(name):
                    props.append(mol.GetProp(name))
                else:
                    props.append(None)
            t = (mol, smi, id, props)
            return t

        except StopIteration as ex:
            raise StopIteration from ex

    def get_extra_field_names(self):
        return self.field_names

    def close(self):
        pass

    def __str__(self):
        return str(self.__class__)


class SmilesReader:

    def __init__(self, input_file, read_header, delimiter, id_col):
        if input_file.endswith('.gz'):
            self.reader = gzip.open(input_file, 'rt')
        else:
            self.reader = open(input_file, 'rt')
        self.delimiter = delimiter
        if id_col is None:
            self.id_col = None
        else:
            self.id_col = int(id_col)
        self.field_names = None
        # skip header lines
        if read_header:
            line = self.reader.readline()
            tokens = self.tokenize(line)
            self.field_names = []
            for token in tokens:
                self.field_names.append(token.strip())

    def tokenize(self, line):
        line = line.strip()
        if self.delimiter is None:
            tokens = line.split()
        else:
            tokens = line.split(self.delimiter)
        stripped = []
        for token in tokens:
            stripped.append(token.strip())
        return stripped

    def read(self):
        line = self.reader.readline()
        if line:
            tokens = self.tokenize(line)
            smi = tokens[0]
            if self.id_col:
                id = tokens[self.id_col]
            else:
                id = None

            mol = Chem.MolFromSmiles(smi)
            if not mol:
                raise TypeError(f'{self}: Error parsing molecule')            
            props = []

            for i, token in enumerate(tokens):
                token = token.strip()
                if i != 0:
                    props.append(token)
                    if mol:
                        if self.field_names:
                            mol.SetProp(self.field_names[i], token)
                        else:
                            mol.SetProp('field' + str(i), token)


            t = (mol, smi, id, props)
            return t
        else:
            raise StopIteration

    def get_extra_field_names(self):
        if self.field_names:
            results = []
            for i, name in enumerate(self.field_names):
                if i != 0 and i != self.id_col:
                    results.append(name)
            return results
        else:
            return None

    def close(self):
        self.reader.close()

    def __str__(self):
        return str(self.__class__)        


def generate_header_values(extra_field_names, num_orig_props, calc_prop_names):
    headers = ['smiles']
    if extra_field_names:
        headers.extend(extra_field_names)
    else:
        for i in range(num_orig_props):
            headers.append('field' + str(i + 2))
    headers.extend(calc_prop_names)
    return headers


def create_reader(input_file, type=None, id_column=None, sdf_read_records=100, read_header=False, delimiter='\t'):
    if type is None:
        if input_file.endswith('.sdf') or input_file.endswith('.sdf.gz') or input_file.endswith('.sd') or input_file.endswith('.sd.gz'):
            type = 'sdf'
        else:
            type = 'smi'

    if type == 'sdf':
        return SdfReader(input_file, id_column, sdf_read_records)
    elif type == 'smi':
        return SmilesReader(input_file, read_header, delimiter, id_column)
    else:
        raise ValueError('Unexpected file type', type)


def create_writer(outfile, delimiter='\t'):
    if outfile.endswith('.sdf') or outfile.endswith('sd'):
        return SdfWriter(outfile)
    else:
        return SmilesWriter(outfile, delimiter)


def updateChargeFlagInAtomBlock(mb):
    """
    Add data for the charges to the atom block. This data is now deprecated and should be specified using "M  CHG" lines
    but some old software such as rDock only handle the old syntax. RDKit only supports the new syntax so this method
    handles adding the old syntax as well.
    This function is based on work by Jose Manuel Gally that can be found here:
    See https://sourceforge.net/p/rdkit/mailman/message/36425493/
    """
    f="{:>10s}"*3+"{:>2}{:>4s}"+"{:>3s}"*11
    chgs = []    # list of charges
    lines = mb.split("\n")
    #if mb[0] == '' or mb[0] == "\n":
    #    del lines[0]
    CTAB = lines[3]
    atomCount = int(CTAB.split()[0])
    # parse mb line per line
    for l in lines:
        # look for M CHG property
        if l[0:6] == "M  CHG":
            records = l.split()[3:]    # M  CHG X is not needed for parsing, the info we want comes afterwards
            # record each charge into a list
            for i in range(0,len(records),2):
                idx = records[i]
                chg = records[i+1]
                chgs.append((int(idx), int(chg)))    # sort tuples by first element?
            break    # stop iterating

    # sort by idx in order to parse the molblock only once more
    chgs = sorted(chgs, key=lambda x: x[0])

    # that we have a list for the current molblock, attribute each charges
    for chg in chgs:
        i=4
        while i < 4+atomCount:    # do not read from beginning each time, rather continue parsing mb!
            # when finding the idx of the atom we want to update, extract all fields and rewrite whole sequence
            if i-3 == chg[0]:    # -4 to take into account the CTAB headers, +1 because idx begin at 1 and not 0
                fields = lines[i].split()
                x=fields[0]
                y=fields[1]
                z=fields[2]
                symb=fields[3]
                massDiff=fields[4]
                charge=fields[5]
                sp=fields[6]
                hc=fields[7]
                scb=fields[8]
                v=fields[9]
                hd=fields[10]
                nu1=fields[11]
                nu2=fields[12]
                aamn=fields[13]
                irf=fields[14]
                ecf=fields[15]
                # update charge flag
                if chg[1] == -1:
                    charge = '5'
                elif chg[1] == -2:
                    charge = '6'
                elif chg[1] == -3:
                    charge = '7'
                elif chg[1] == 1:
                    charge = '3'
                elif chg[1] == 2:
                    charge = '2'
                elif chg[1] == 3:
                    charge = '1'
                else:
                    utils.log("ERROR! " + str(lines[0]) + "unknown charge flag: " + str(chg[1]))
                    break
                # update modatom block line
                lines[i] = f.format(x,y,z,symb,massDiff,charge,sp,hc,scb,v,hd,nu1,nu2,aamn,irf,ecf)
            i+=1
    upmb = "\n".join(lines)
    return(upmb)


def get_num_chiral_centers(mol):
    chiral_centers = Chem.FindMolChiralCenters(mol, force=True, includeUnassigned=True)
    undef_cc = 0
    for cc in chiral_centers:
        if cc[1] == '?':
            undef_cc += 1
    return len(chiral_centers), undef_cc


def get_num_sp3_centres(mol):
    return sum((x.GetHybridization() == Chem.HybridizationType.SP3) for x in mol.GetAtoms())


def sdf_read_mol(input_file):
    txt = sdf_read_block(input_file)
    if txt == None:
        return None
    mol = Chem.MolFromMolBlock(txt)
    return (txt, mol)


def sdf_record_gen(hnd):
    """A generator for text records fom a SD file
    """
    mol_text_tmp = ""
    while 1:
        line = hnd.readline()
        if not line:
            return
        line = line.decode("utf-8")
        mol_text_tmp += line
        if line.startswith("$$$$"):
            mol_text = mol_text_tmp
            mol_text_tmp = ""
            yield mol_text


def rdk_read_single_mol(input_file):
    # read the molecule. Can be SDF or Mol format
    # if SDF then the first molecule is used.
    if input_file.endswith('.mol'):
        mol = Chem.MolFromMolFile(input_file)
    elif input_file.endswith('.sdf'):
        suppl = Chem.SDMolSupplier(input_file)
        mol = next(suppl)
    elif input_file.endswith('.sdf.gz'):
        with gzip.open(input_file, 'rb') as gz:
            suppl = Chem.ForwardSDMolSupplier(gz)
            mol = next(suppl)
    else:
        raise ValueError('Unsupported file type. Must be .mol .sdf or .sdf.gz. Found ' + input_file)
    return mol


def rdk_read_mols(input_file):
    """
    Read molecules from a single file (.mol or .sdf)
    """
    if input_file.endswith('.mol'):
        mol = Chem.MolFromMolFile(input_file)
        return [mol]
    else:
        supplr = rdk_mol_supplier(input_file)
        mols = [m for m in supplr]
        return mols


def rdk_read_molecule_files(inputs):
    """
    Read input molecules.
    A list of inputs is specified, each one containing either a single filename or a comma separated list of filenames
    (no spaces). The files can either be .mol files with a single molecule or .sdf files with multiple molecules.
    The molecules in all the files specified are read and returned as a list of molecules.
    Examples:
    ['mols.sdf']                          - single SDF with one or more molecules
    ['mol1.mol', 'mol2.mol', 'mol3.mol']  - 3 molfiles as separate elements of the list
    ['mol1.mol,mol2.mol,mol3.mol']        - 3 molfiles as a single comma separated element of the list
    ['mol1.mol,mol2.mol', 'mols.sdf']     - 3 molfiles plus 1 SDF

    :param input_files: Input_File filenames
    :return: List of molecules that have been read
    """
    mols = []
    for input_file in inputs:
        tokens = input_file.split(',')
        for token in tokens:
            if token.endswith('.mol'):
                m = Chem.MolFromMolFile(token)
                if m:
                    mols.append(m)
                else:
                    DmLog.emit_event('WARNING: could not process', token)
            else:
                ms = rdk_read_mols(token)
                for i, m in enumerate(ms):
                    if m:
                        mols.append(m)
                    else:
                        DmLog.emit_event('WARNING: could not process molecule', i, 'from', token)
    return mols


def rdk_merge_mols(inputs):
    """
    Merge multiple molecules into a single molecule
    :param input_file:
    :return:
    """
    merged_mol = Chem.RWMol()
    mols = rdk_read_molecule_files(inputs)
    for i, mol in enumerate(mols):
        if mol:
            merged_mol.InsertMol(mol)
        else:
            DmLog.emit_event('WARNING: could not process molecule', i)
    Chem.SanitizeMol(merged_mol)
    return merged_mol, i + 1


def rdk_mol_supplier(input_file):
    """
    Generate a ForwardSDMolSupplier from a .sdf or .sdf.gz file
    """
    if input_file.endswith('.sdf'):
        suppl = Chem.ForwardSDMolSupplier(input_file)
    elif input_file.endswith('.sdf.gz'):
        gz = gzip.open(input_file, 'rb')
        suppl = Chem.ForwardSDMolSupplier(gz)
    else:
        raise ValueError('Unsupported file type. Must be .sdf or .sdf.gz. Found ' + input_file)
    return suppl


def fragment(mol, mode):
    frags = Chem.GetMolFrags(mol, asMols=True)

    if len(frags) == 1:
        return mol
    else:
        # TODO - handle ties
        biggest_index = -1
        i = 0
        if mode == 'hac':
            biggest_count = 0
            for frag in frags:
                hac = frag.GetNumHeavyAtoms()
                if hac > biggest_count:
                    biggest_count = hac
                    biggest_mol = frag
                    biggest_index = i
                i+=1
            utils.log("Chose fragment", biggest_index, "from", len(frags), "based on HAC")
        elif mode == 'mw':
            biggest_mw = 0
            for frag in frags:
                mw = Descriptors.MolWt(frag)
                if mw > biggest_mw:
                    biggest_mw = mw
                    biggest_mol = frag
                    biggest_index = i
                i+=1
            utils.log("Chose fragment", biggest_index, "from", len(frags), "based on MW")
        else:
            raise ValueError('Invalid fragment mode:',mode)

        # copy the properties across
        for name in mol.GetPropNames():
            biggest_mol.SetProp(name, mol.GetProp(name))

        # _Name is a magical property that is not in the ones returned by GetPropNames
        if '_Name' in mol.GetPropNames():
            biggest_mol.SetProp("_Name", mol.GetProp("_Name"))

        return biggest_mol


def fragmentAndFingerprint(reader, mols, data, fps, descriptor, fragmentMethod='hac', outputFragment=False):
    """
    Fragment the molecule if it has multiple fragments and generate fingerprints on the fragment.

    :param reader: MolSupplier from which to read the molecules
    :param mols: List to which the molecules are added
    :param fps: List to which the fingerprints are added
    :param descriptor: Function to generate the fingerprints from the molecule
    :param fragmentMethod: The fragmentation method to use when there are multiple fragments (hac or mw)
    :param outputFragment: Boolean that specifies whether to write the fragment or the original molecule to the mols list
    :return: The number of errors encountered
    """
    errors = 0
    count = 0
    while True:
        count += 1
        t = reader.read()
        # break if no more data to read
        if not t:
            break
        mol, smi, id, props = t
        if not mol:
            utils.log('Failed to read molecule', count)
            errors += 1
            continue

        frag = fragment(mol, fragmentMethod)
        d = descriptor(frag)
        if d:
            if outputFragment:
                mols.append(frag)
                data.append((id, Chem.MolToSmiles(frag), props))
            else:
                mols.append(mol)
                data.append((id, Chem.MolToSmiles(mol), props))
            fps.append(d)
            continue
    return errors
