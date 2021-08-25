import csv
import ntpath
import pynmrstar
from datetime import date

def read_csv_file(file_name):
    fpath,fname = ntpath.split(file_name)
    meta_data = {}
    meta_info = fname.split("_")
    meta_data['software_name'] = meta_info[0]
    meta_data['proteome_id'] = meta_info[1]
    meta_data['uniprot_id'] = meta_info[2].split("-")[1]
    #print (meta_data)
    data=[]
    with open(file_name) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter = ',')
        for row in csv_reader:
            data.append(row)
    return meta_data,data[1:],fpath

def write_nmrstar(meta_data,data,out_path=None):
    entry_id = '{}-{}'.format(meta_data['proteome_id'],meta_data['uniprot_id'])
    ent = pynmrstar.Entry.from_scratch(entry_id)
    sf = pynmrstar.Saveframe.from_scratch(sf_name='entry_information',tag_prefix='Entry')
    sf.add_tag('Sf_category','entry_information')
    sf.add_tag('Sf_framecode','entry_information')
    sf.add_tag('ID',entry_id)
    sf.add_tag('NMR_STAR_version','3.2')
    lp = pynmrstar.Loop.from_scratch('Related_entries')
    lp.add_tag(['Database_name','Database_accession_code','Entry_ID'])
    lp.add_data(['UniProt',meta_data['uniprot_id'],entry_id])
    sf.add_loop(lp)
    lp = pynmrstar.Loop.from_scratch('Release')
    lp.add_tag(['Release_number','Date','Type','Detail','Entry_ID'])
    today = date.today()
    lp.add_data([1,today.strftime("%Y-%m-%d"),'Initial','AF prediction',entry_id])
    sf.add_loop(lp)
    ent.add_saveframe(sf)
    sf = pynmrstar.Saveframe.from_scratch(meta_data['software_name'],tag_prefix='Software')
    sf.add_tag('Sf_category', 'software')
    sf.add_tag('Sf_framecode',meta_data['software_name'] )
    sf.add_tag('ID',1)
    sf.add_tag('Entry_ID',entry_id)
    sf.add_tag('Name',meta_data['software_name'])
    sf.add_tag('Version','.')
    ent.add_saveframe(sf)
    sf = pynmrstar.Saveframe.from_scratch('predicted_chem_shift_{}'.format(meta_data['software_name']),
                                          tag_prefix='Assigned_chem_shift_list')
    sf.add_tag('Sf_category','assigned_chemical_shifts')
    sf.add_tag('Sf_framecode','predicted_chem_shift_{}'.format(meta_data['software_name']))
    sf.add_tag('ID',1)
    lp= pynmrstar.Loop.from_scratch('Chem_shift_software')
    lp.add_tag(['Software_ID','Entry_ID','Assigned_chem_shift_list_ID'])
    lp.add_data([1,entry_id,1])
    sf.add_loop(lp)
    lp = pynmrstar.Loop.from_scratch('Atom_chem_shift')
    lp.add_tag(['ID','Entity_assembly_ID','Entity_ID','Comp_index_ID','Comp_ID','Atom_ID',
                'Atom_type','Atom_isotope_number','Val','Val_err','Ambiguity_code','Entry_ID',
                'Assigned_chem_shift_list_ID'])
    i=1
    for dat in data:
        if 'C' in dat[2]:
            atype='C'
            isono=13
        elif 'N' in dat[2]:
            atype='N'
            isono='15'
        else:
            atype='H'
            isono = 1
        lp.add_data([i,1,1,dat[0],dat[1],dat[2],atype,isono,dat[3],'.','.',entry_id,1])
        i+=1
    sf.add_loop(lp)
    ent.add_saveframe(sf)
    out_file = '{}/{}.str'.format(out_path,entry_id)
    with open(out_file, 'w') as wstarfile:
        wstarfile.write(str(ent))

if __name__ == "__main__":
    m,d,p=read_csv_file('/Users/kumaran/Downloads/sparta-plus_UP000002485_AF-O94312-F1-model_v1.csv')
    write_nmrstar(m,d,p)