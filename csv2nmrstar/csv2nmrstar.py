import csv
import ntpath
import pynmrstar
from datetime import date
import os
import json
import logging
import sys



def read_csv_file(file_name,path):
    file_name='{}/{}'.format(path,file_name)
    fpath,fname = ntpath.split(file_name)
    meta_data = {}
    meta_info = fname.split("_")
    #print (meta_info)
    meta_data['software_name'] = meta_info[0]
    meta_data['proteome_id'] = meta_info[2].split(".")[0]
    meta_data['uniprot_id'] = meta_info[1].split("-")[1]
    #print (meta_data)
    data=[]
    with open(file_name) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter = ',')
        for row in csv_reader:
            data.append(row)
    return meta_data,data[1:]

def create_entry_info_sf(meta_data):
    entry_id = '{}-{}'.format(meta_data['proteome_id'], meta_data['uniprot_id'])
    sf = pynmrstar.Saveframe.from_scratch(sf_name='entry_information', tag_prefix='Entry')
    sf.add_tag('Sf_category', 'entry_information')
    sf.add_tag('Sf_framecode', 'entry_information')
    sf.add_tag('ID', entry_id)
    sf.add_tag('NMR_STAR_version', '3.2')
    lp = pynmrstar.Loop.from_scratch('Related_entries')
    lp.add_tag(['Database_name', 'Database_accession_code', 'Entry_ID'])
    lp.add_data(['UniProt', meta_data['uniprot_id'], entry_id])
    sf.add_loop(lp)
    lp = pynmrstar.Loop.from_scratch('Release')
    lp.add_tag(['Release_number', 'Date', 'Type', 'Detail', 'Entry_ID'])
    today = date.today()
    lp.add_data([1, today.strftime("%Y-%m-%d"), '{} prediction added'.format(meta_data['software_name']), 'AF chemical shift prediction', entry_id])
    sf.add_loop(lp)
    return sf

def create_software_sf(meta_data,id):
    entry_id = '{}-{}'.format(meta_data['proteome_id'], meta_data['uniprot_id'])
    sf = pynmrstar.Saveframe.from_scratch(meta_data['software_name'], tag_prefix='Software')
    sf.add_tag('Sf_category', 'software')
    sf.add_tag('Sf_framecode', meta_data['software_name'])
    sf.add_tag('ID', id)
    sf.add_tag('Entry_ID', entry_id)
    sf.add_tag('Name', meta_data['software_name'])
    sf.add_tag('Version', '.')
    return sf

def create_chem_shift_sf(meta_data,data,id):
    flg,msg,atom_dict = load_json_data('../atomDict.json')
    entry_id = '{}-{}'.format(meta_data['proteome_id'], meta_data['uniprot_id'])
    sf = pynmrstar.Saveframe.from_scratch('predicted_chem_shift_{}'.format(meta_data['software_name']),
                                          tag_prefix='Assigned_chem_shift_list')
    sf.add_tag('Sf_category', 'assigned_chemical_shifts')
    sf.add_tag('Sf_framecode', 'predicted_chem_shift_{}'.format(meta_data['software_name']))
    sf.add_tag('ID', id)
    sf.add_tag('Entry_ID', entry_id)
    lp = pynmrstar.Loop.from_scratch('Chem_shift_software')
    lp.add_tag(['Software_ID', 'Entry_ID', 'Assigned_chem_shift_list_ID'])
    lp.add_data([id, entry_id, id])
    sf.add_loop(lp)
    lp = pynmrstar.Loop.from_scratch('Atom_chem_shift')
    lp.add_tag(['ID', 'Entity_assembly_ID', 'Entity_ID', 'Comp_index_ID', 'Comp_ID', 'Atom_ID',
                'Atom_type', 'Atom_isotope_number', 'Val', 'Val_err', 'Ambiguity_code', 'Entry_ID',
                'Assigned_chem_shift_list_ID'])
    i = 1
    for dat in data:
        if 'C' in dat[2]:
            atype = 'C'
            isono = 13
        elif 'N' in dat[2]:
            atype = 'N'
            isono = '15'
        else:
            atype = 'H'
            isono = 1
        if dat[2] not in atom_dict[dat[1]]:
            if dat[1]=='MET' and dat[2]=='HE':
                lp.add_data([i, 1, 1, dat[0], dat[1], 'HE1', atype, isono, dat[3], '.', '1', entry_id, id])
                lp.add_data([i, 1, 1, dat[0], dat[1], 'HE2', atype, isono, dat[3], '.', '1', entry_id, id])
                lp.add_data([i, 1, 1, dat[0], dat[1], 'HE3', atype, isono, dat[3], '.', '1', entry_id, id])
            elif dat[1]=='LYS' and dat[2]=='HZ':
                lp.add_data([i, 1, 1, dat[0], dat[1], 'HZ1', atype, isono, dat[3], '.', '1', entry_id, id])
                lp.add_data([i, 1, 1, dat[0], dat[1], 'HZ2', atype, isono, dat[3], '.', '1', entry_id, id])
                lp.add_data([i, 1, 1, dat[0], dat[1], 'HZ3', atype, isono, dat[3], '.', '1', entry_id, id])
            elif dat[1]=='LEU' and dat[2]=='HD1':
                lp.add_data([i, 1, 1, dat[0], dat[1], 'HD11', atype, isono, dat[3], '.', '1', entry_id, id])
                lp.add_data([i, 1, 1, dat[0], dat[1], 'HD12', atype, isono, dat[3], '.', '1', entry_id, id])
                lp.add_data([i, 1, 1, dat[0], dat[1], 'HD13', atype, isono, dat[3], '.', '1', entry_id, id])
            elif dat[1]=='LEU' and dat[2]=='HD2':
                lp.add_data([i, 1, 1, dat[0], dat[1], 'HD21', atype, isono, dat[3], '.', '1', entry_id, id])
                lp.add_data([i, 1, 1, dat[0], dat[1], 'HD22', atype, isono, dat[3], '.', '1', entry_id, id])
                lp.add_data([i, 1, 1, dat[0], dat[1], 'HD23', atype, isono, dat[3], '.', '1', entry_id, id])
            elif dat[1]=='GLY' and dat[2]=='HA':
                lp.add_data([i, 1, 1, dat[0], dat[1], 'HA2', atype, isono, dat[3], '.', '2', entry_id, id])
                lp.add_data([i, 1, 1, dat[0], dat[1], 'HA3', atype, isono, dat[3], '.', '2', entry_id, id])
            elif dat[1]=='VAL' and dat[2]=='HG1':
                lp.add_data([i, 1, 1, dat[0], dat[1], 'HG11', atype, isono, dat[3], '.', '1', entry_id, id])
                lp.add_data([i, 1, 1, dat[0], dat[1], 'HG12', atype, isono, dat[3], '.', '1', entry_id, id])
                lp.add_data([i, 1, 1, dat[0], dat[1], 'HG13', atype, isono, dat[3], '.', '1', entry_id, id])
            elif dat[1]=='VAL' and dat[2]=='HG2':
                lp.add_data([i, 1, 1, dat[0], dat[1], 'HG21', atype, isono, dat[3], '.', '1', entry_id, id])
                lp.add_data([i, 1, 1, dat[0], dat[1], 'HG22', atype, isono, dat[3], '.', '1', entry_id, id])
                lp.add_data([i, 1, 1, dat[0], dat[1], 'HG23', atype, isono, dat[3], '.', '1', entry_id, id])
            elif dat[1]=='MET' and dat[2]=='HE':
                lp.add_data([i, 1, 1, dat[0], dat[1], 'HE1', atype, isono, dat[3], '.', '1', entry_id, id])
                lp.add_data([i, 1, 1, dat[0], dat[1], 'HE2', atype, isono, dat[3], '.', '1', entry_id, id])
                lp.add_data([i, 1, 1, dat[0], dat[1], 'HE3', atype, isono, dat[3], '.', '1', entry_id, id])
            elif dat[1]=='ALA' and dat[2]=='HB':
                lp.add_data([i, 1, 1, dat[0], dat[1], 'HB1', atype, isono, dat[3], '.', '1', entry_id, id])
                lp.add_data([i, 1, 1, dat[0], dat[1], 'HB2', atype, isono, dat[3], '.', '1', entry_id, id])
                lp.add_data([i, 1, 1, dat[0], dat[1], 'HB3', atype, isono, dat[3], '.', '1', entry_id, id])
            elif dat[1]=='THR' and dat[2]=='HG2':
                lp.add_data([i, 1, 1, dat[0], dat[1], 'HG21', atype, isono, dat[3], '.', '1', entry_id, id])
                lp.add_data([i, 1, 1, dat[0], dat[1], 'HG22', atype, isono, dat[3], '.', '1', entry_id, id])
                lp.add_data([i, 1, 1, dat[0], dat[1], 'HG23', atype, isono, dat[3], '.', '1', entry_id, id])
            elif dat[1]=='ILE' and dat[2]=='HG2':
                lp.add_data([i, 1, 1, dat[0], dat[1], 'HG21', atype, isono, dat[3], '.', '1', entry_id, id])
                lp.add_data([i, 1, 1, dat[0], dat[1], 'HG22', atype, isono, dat[3], '.', '1', entry_id, id])
                lp.add_data([i, 1, 1, dat[0], dat[1], 'HG23', atype, isono, dat[3], '.', '1', entry_id, id])
            elif dat[1]=='ILE' and dat[2]=='HD1':
                lp.add_data([i, 1, 1, dat[0], dat[1], 'HD11', atype, isono, dat[3], '.', '1', entry_id, id])
                lp.add_data([i, 1, 1, dat[0], dat[1], 'HD12', atype, isono, dat[3], '.', '1', entry_id, id])
                lp.add_data([i, 1, 1, dat[0], dat[1], 'HD13', atype, isono, dat[3], '.', '1', entry_id, id])
            else:
                print (dat)
        lp.add_data([i, 1, 1, dat[0], dat[1], dat[2], atype, isono, dat[3], '.', '1', entry_id, id])
        i += 1
    sf.add_loop(lp)
    return sf




def write_nmrstar(meta_data,data,out_path=None):
    flg,msg,atom_dict = load_json_data('../atomDict.json')
    entry_id = '{}-{}'.format(meta_data['proteome_id'],meta_data['uniprot_id'])
    str_file='{}/{}.str'.format(out_path,entry_id)
    #print (str_file)
    if not os.path.isfile(str_file):
        logging.info('Creating  {} file using {}'.format(str_file,meta_data['software_name']))
        ent = pynmrstar.Entry.from_scratch(entry_id)
        ent_info_sf = create_entry_info_sf(meta_data)
        ent.add_saveframe(ent_info_sf)
        software_sf = create_software_sf(meta_data,1)
        ent.add_saveframe(software_sf)
        cs_sf = create_chem_shift_sf(meta_data,data,1)
        ent.add_saveframe(cs_sf)
        ent.normalize()
        with open(str_file, 'w') as wstarfile:
            wstarfile.write(str(ent))
    else:
        ent=pynmrstar.Entry.from_file(str_file)
        software_name = ent.get_tag('_Software.Name')
        if  meta_data['software_name'] in software_name:
            logging.info('{} prediction already exists; NMR-STAR file not updated'.format(meta_data['software_name']))
        else:
            logging.info('Updating  {} file using {}'.format(str_file, meta_data['software_name']))
            id = max([int(i) for i in ent.get_tag('_Software.ID')]) + 1
            for sf in ent:
                for lp in sf:
                    if lp.category == '_Release':
                        relese_num = max([int(i) for i in ent.get_tag('_Release.Release_number')]) + 1
                        today = date.today()
                        lp.add_data([relese_num, today.strftime("%Y-%m-%d"),
                                     '{} prediction added'.format(meta_data['software_name']),
                                     'AF chemical shift prediction', entry_id])
            sofrware_sf = create_software_sf(meta_data,id)
            ent.add_saveframe(sofrware_sf)
            chem_shift_sf = create_chem_shift_sf(meta_data,data,id)
            ent.add_saveframe(chem_shift_sf)
            ent.normalize()
            with open(str_file, 'w') as wstarfile:
                wstarfile.write(str(ent))




def load_json_data(json_file):
    """
    Loads json data files from lib folder
    :param json_file: json file
    :return: dictionay
    """
    try:
        with open(json_file, 'r') as jsonF:
            data_dict = json.loads(jsonF.read())
        is_ok = True
        msg = "{} file is read!".format(json_file)
    except IOError:
        msg = "{} file is missing!".format(json_file)
        is_ok = False
        data_dict = []
    return is_ok, msg, data_dict

def generate_nmrstar(filename,inpath,outpath):
    meta_data,data=read_csv_file(filename,inpath)
    write_nmrstar(meta_data,data,outpath)
if __name__ == "__main__":
    flist=sys.argv[1]
    inpath=sys.argv[2]
    outpath=sys.argv[3]
    f=open(flist).read().split("\n")
    for l in f:
        generate_nmrstar(l,inpath,outpath)
