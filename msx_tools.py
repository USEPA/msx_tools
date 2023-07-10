# -*- coding: utf-8 -*-
"""
TOOLS ASSOCIATED WITH RUNNING AND USING EPANET-MSX

@author: JBurkhar
"""

import wntr
import numpy as np
import pandas as pd
import os
import subprocess
import re ## text processing

from msx_library import msx_dict


std_options = ['area_units', 'rate_units', 'solver', 'timestep', 'coupling', 
               'rtol', 'atol', 'compiler', 'segments', 'peclet']

file_loc = os.path.abspath(__file__)
split_file = file_loc.split("\\")
msx_location = '/'.join(split_file[:-1])


class MSXobj():
    def __init__(self, **kw):
        '''
        

        Parameters
        ----------
        **kw : TYPE
            DESCRIPTION.
            
        Defaults in brackets []
        title      = title of MSX file
        OPTIONS:
        area_units = [M2] or FT2     : surface area units
        rate_units = SEC, [MIN], HR  : rate units
        solver     = [RK5], ROS2, EUL: solver to use
        timestep   = number of seconds for the timestep [INT, default = 1 minute, 60 seconds]
        Returns
        -------
        None.

        '''
        self.title = kw.get('title', 'MSX Builder - Python Generated File')
        ### Options Definitions
        self.options = {}
        self.area_units = kw.get('area_units', 'M2')
        self.options['area_units'] = {'val': self.area_units, 
                                      'note': 'Surface concentration in mass/area_units'}
        self.rate_units = kw.get('rate_units', 'MIN')
        self.options['rate_units'] = {'val': self.rate_units, 
                                     'note': 'Reaction rates are concentration/rate_units'}
        self.solver = kw.get('solver', 'RK5')
        self.options['solver'] = {'val': self.solver,
                                  'note': 'RK5, 5th order Runge-Kutta Integrator or  ROS2, 2nd order Rosenbrock or EUL, Euler'}
        self.timestep = kw.get('timestep', 60)
        self.options['timestep'] = {'val': int(self.timestep),
                                    'note': 'in seconds, Default 60 seconds'}
        self.coupling = kw.get('coupling', 'NONE')
        self.options['coupling'] = {'val': self.coupling,
                                    'note': 'Default NONE?'}
        self.rtol = kw.get('rtol', 1e-4)
        self.options['rtol'] = {'val': self.rtol,
                                'note': 'Relative concentration tolerance, 1e-4 default'}
        self.atol = kw.get('atol', 1e-4)
        self.options['atol'] = {'val': self.atol,
                                'note': "Absolute concentration tolerance, 1e-4 default"}
        self.compiler = kw.get('compiler', '')
        if self.compiler != '':
            self.options['compiler'] = {'val': self.compiler,
                                        'note': ''}
        self.segments = kw.get('segments', 5000)
        self.options['segments'] = {'val': self.segments,
                                    'note': 'Maximum number of segments per pipe (MSX 2.0 or newer only)'}
        self.peclet = kw.get('peclet', 1000)
        self.options['peclet'] = {'val': self.peclet,
                                  'note': 'Peclet threshold for applying dispersion, 1000 default (MSX 2.0 or newer only)'}
        
        ## create additional storage places
        ### containers that may have duplication of species are lists, if unique species, species acts as key for dictionary
        self.species = {} 
        self.coefficients = {}
        self.terms = {}
        self.pipes = {}
        self.tanks = {}
        self.sources = [] 
        self.parameters = []
        self.quality = [] 
        self.patterns = {}
        self.report = []
        self.dispersion = {}
        self.msx_info_dict = 'None'
        
        ## read in preconfigured MSX reactions. Set up msx_helper for user's reference
        self.msx_dict = msx_dict
        self.msx_helper = {}
        for key in self.msx_dict.keys():
            self.msx_helper[key] = self.msx_dict[key]['description']
        
        ## Is user supplying a file, or reaction type??
        self.read_file = kw.get('file_name', '')
        self.rxn_type = kw.get('rxn_type', '')
        
        if self.rxn_type != '' and self.read_file != '':
            print('WARNING: User provided both a reaction type and MSX file, using MSX file.')
            print('IF you wish to only supply a reaction type, remove "file_name" keyword')
            self.rxn_type = ''
        
        ###################################
        ## set up where to find file information
        
        if self.read_file != '':
            print(f"INFO: Reading data from {self.read_file}\n")
            try:
                self.msx_info_dict = _read_msx_file_to_dict(self.read_file)
            except Exception as e:
                print(f'ERROR READING IN {self.read_file}: {e}')
        else:
            pass
            # print('No File Provided')  ## will delete
        
        
        
        
# =============================================================================
#       ### SUPPLIED REACTION NAME ############ USES MSX_LIBRARY ##############  
# =============================================================================
        ### works with preloaded reaction dictionaries
        if self.rxn_type != '':
            print('SELECTED:    ', self.rxn_type.lower())
            print("Description: ", self.msx_helper[self.rxn_type.lower()],'\n')
            if self.rxn_type.lower() in self.msx_dict.keys():
                self.msx_info_dict = self.msx_dict[self.rxn_type.lower()]
            else:
                print(f"ERROR: Dictionary for {self.rxn_type} is not available")
        
# =============================================================================
#       #### BUILD FROM DICTIONARY INTO CLASS OBJECT###########################
#       USED by either actual file, or reaction specification
# =============================================================================
        if self.msx_info_dict != 'None':
            key_list = list(self.msx_info_dict.keys())
            if 'title' in key_list:
                self.title = self.msx_info_dict['title']
            
            if 'options' in key_list:
                for option in std_options:
                    if option in self.msx_info_dict['options'].keys(): 
                        self.options[option] = self.msx_info_dict['options'][option]
                    elif option.upper() in self.msx_info_dict['options'].keys():
                        self.options[option] = self.msx_info_dict['options'][option.upper()]
                        
            if 'species' in key_list:
                for key in self.msx_info_dict['species'].keys():
                    self.species[key] = self.msx_info_dict['species'][key]
                    
            if 'coefficients' in key_list:
                for key in self.msx_info_dict['coefficients'].keys():
                    self.coefficients[key] = self.msx_info_dict['coefficients'][key]
            
            if 'terms' in key_list:
                for key in self.msx_info_dict['terms'].keys():
                    self.terms[key] = self.msx_info_dict['terms'][key]
            
            if 'pipes' in key_list:
                for key in self.msx_info_dict['pipes'].keys():
                    self.pipes[key] = self.msx_info_dict['pipes'][key]
            
            if 'tanks' in key_list:
                for key in self.msx_info_dict['tanks'].keys():
                    self.tanks[key] = self.msx_info_dict['tanks'][key]
            
            if 'sources' in key_list:
                for item in self.msx_info_dict['sources']:
                    if len(item) == 4: ## adds empty for pattern and note
                        item.append('')
                        item.append('')
                    elif len(item) == 5: ## adds empty for note
                        item.append('')
                    self.sources.append(item)
            
            if 'dispersion' in key_list:
                for key in self.msx_info_dict['dispersion'].keys():
                    self.dispersion[key] = self.msx_info_dict['dispersion'][key]
            
            if 'parameters' in key_list:
                for item in self.msx_info_dict['parameters']:
                    self.parameters.append(item)
            
            if 'quality' in key_list:
                #### TODO: NEEDS TO HANDLE type NODE/LINK SPECIES VALUE ; NOTE STYLE

                for item in self.msx_info_dict['quality']:
                    if len(item) == 3: ## add empty note if missing
                        item.append('')
                    self.quality.append(item)
            
            if 'patterns' in key_list:
                for key in self.msx_info_dict['patterns'].keys():
                    self.patterns[key] = self.msx_info_dict['patterns'][key]
            
            if 'report' in key_list:
                for item in self.msx_info_dict['report']:
                    self.report.append(item)
                        
# =============================================================================
#   #################################### FUNCTIONS ############################    
# =============================================================================
        
    def add_species(self, species_name, units, note='', loc_type='BULK' ):
        if species_name in self.species.keys():
            print(f"WARNING: SPECIES: Adding {species_name} again, new data will overwite previous values.")
        self.species[species_name] = {'val': units,
                                      'type': loc_type,
                                      'note': note}
        
    def add_coefficient(self, coeff_name, value, note='', coeff_type='CONSTANT'):
        if coeff_name in self.coefficients.keys():
            print(f"WARNING: COEFFICIENT: {coeff_name} already defined, new data will overwrite previous values.")
        self.coefficients[coeff_name] = {'val': value,
                                         'type': coeff_type,
                                         'note': note}
    def update_coefficient(self, coeff_name, value, note='', coeff_type='CONSTANT'):
        if coeff_name not in self.coefficients.keys():
            print(f"WARNING: COEFFICIENTS: Attempting to update {coeff_name} that does not exist: Adding instead.")
            self.add_coefficient(coeff_name, value, note, coeff_type)
        else:
            if note == '' and self.coefficients[coeff_name]['note'] != '':
                ## preserves the original notes if present
                note = self.coefficients[coeff_name]['note']
            if coeff_type == 'CONSTANT' and coeff_type != self.coefficients[coeff_name]['type']:
                # print(f"Warning: Coefficient Type mismatch for updating {coeff_name}, using original")
                coeff_type = self.coefficients[coeff_name]['type']
            self.coefficients[coeff_name]['val'] = value * 1
            self.coefficients[coeff_name]['type'] = coeff_type
            self.coefficients[coeff_name]['note'] = note
    
    def add_term(self, term_name, eq, note=''):
        if term_name in self.terms.keys():
            print(f"WARNING: TERM: {term_name} already defined, new data will overwite previous values")
        self.terms[term_name] = {'val': eq,
                                 'note': note}
    
    def add_pipe_rxn(self, rxn_type, species, rate, note=''):
        ### currently assumes you can only have rate or equil per species... might need to tweak if not limited
        self.pipes[species] = {'type': rxn_type,
                               'val': rate,
                               'note': note}
    
    def add_tank_rxn(self, rxn_type, species, rate, note=''):
        ### currently assumes you can only have rate or equil per species... might need to tweak if not limited
        self.tanks[species] = {'type': rxn_type,
                               'val': rate,
                               'note': note}
        
    def add_source(self, source_type, node_ID, species, value, pattern='', note=''):
        self.sources.append([source_type, node_ID, species, value, pattern, note])
        
    def add_dispersion (self, species, value, extra='', note=''):
        self.dispersion[species] = {'val': value,
                                    'extra': extra,
                                    'note': note}
    
    def add_parameter(self, loc_type, ID, coeff_name, value, note=''):
        self.parameters.append([loc_type, ID, coeff_name, value, note])
    
    def add_quality(self, loc_type, species, value, note=''):
        self.quality.append([loc_type, species, value, note])
        
    def add_pattern(self, pattern_ID, pattern):
        if pattern_ID in self.patterns.keys():
            print(f"WARNING: PATTERNS: {pattern_ID} already defined, new data will overwrite previous values")
        self.patterns[pattern_ID] = pattern
    
    
    def add_report_type(self, loc, value):
        self.report.append([loc, value])
    
    def add_report_item(self, report_type, species, report='YES', sig_figs=5):
        self.report.append([report_type, species, report, sig_figs])
    
    def _available_species(self):
        print(list(self.species.keys()))
    
    def build_msx_file(self, file_name='temp.msx', style='MSX2'):
        '''
        

        Parameters
        ----------
        file_name : string, filename, optional
            DESCRIPTION. The default is 'temp.msx'.
        style : string, [MSX2 or MSX11], optional
            DESCRIPTION. The default is 'MSX2'.

        Returns
        -------
        None.

        '''
        file_string = '[TITLE]\n'
        file_string += self.title + '\n'*2
        
        file_string += '[OPTIONS]\n'
        for key in self.options.keys():
            skip = False
            if style != 'MSX2':
                if key in ['peclet', 'segments', ]:
                    skip = True
            
            if not skip:
                lval = len(key)
                file_string += f"  {key.upper()}{' '*(11-lval)}{self.options[key]['val']} \t\t\t;{self.options[key]['note']}\n"
        
        file_string += '\n[SPECIES]\n'
        for key in self.species.keys():
            lval = len(key)
            file_string += f"  {self.species[key]['type'].upper()} {key}{' '*(11-lval)}{self.species[key]['val']}\t\t\t;{self.species[key]['note']}\n"
        
        file_string += '\n[COEFFICIENTS]\n'
        for key in self.coefficients.keys():
            lval = len(key)
            file_string += f"  {self.coefficients[key]['type'].upper()} {key}{' '*(11-lval)}{self.coefficients[key]['val']}\t\t\t;{self.coefficients[key]['note']}\n"
        
        file_string += '\n[TERMS]\n'
        for key in self.terms.keys():
            file_string += f"  {key}\t\t{self.terms[key]['val']}\t\t\t;{self.terms[key]['note']}\n"
        
        file_string += '\n[PIPES]\n'
        for key in self.pipes.keys():
            file_string += f"  {self.pipes[key]['type'].upper()}\t{key}\t\t{self.pipes[key]['val']}\t\t\t;{self.pipes[key]['note']}\n"
        
        
        file_string += '\n[TANKS]\n'
        for key in self.tanks.keys():
            file_string += f"  {self.tanks[key]['type'].upper()}\t{key}\t\t{self.tanks[key]['val']}\t\t\t;{self.tanks[key]['note']}\n"
        
        
        file_string += '\n[SOURCES]\n'
        for item in self.sources:
            file_string += f"  {item[0].upper()}\t{item[1]}\t\t\t\t{item[2]}\t\t{item[3]}\t\t\t{item[4]}\t\t\t;{item[5]}\n"
        
        if style == 'MSX2':
            file_string += '\n[DIFFUSIVITY]\n'
            for key in self.species.keys():
                if key not in self.dispersion.keys():
                    file_string += f"  {key}\t1.0\n"
                else:
                    file_string += f"  {key}\t{self.dispersion[key]['val']}\t\t{self.dispersion[key]['extra']}\t\t\t;{self.dispersion[key]['note']}\n"
        
        file_string += '\n[PARAMETERS]\n'
        for item in self.parameters:
            file_string += f"  {item[0].upper()}\t{item[1]}\t\t\t\t{item[2]}\t\t\t{item[3]}\t\t\t;{item[4]}\n"
        
        
        file_string += '\n[QUALITY]\n'
        for item in self.quality:
            #### TODO: NEEDS TO HANDLE type NODE/LINK SPECIES VALUE ; NOTE STYLE
            file_string += f"  {item[0].upper()}\t{item[1]}\t\t\t\t{item[2]}\t\t\t;{item[3]}\n"

        
        file_string += '\n[PATTERNS]'
        for key in self.patterns.keys():
            patt_vals = self.patterns[key]
            for i in range(len(patt_vals)):
                if i % 10 == 0:
                    file_string += f"\n  {key}\t{patt_vals[i]}"
                else:
                    file_string += f"\t{patt_vals[i]}"
        
        file_string += '\n\n[REPORT]\n'
        for item in self.report:
            if len(item) > 1:
                for sub_item in item:
                    if item[0] == sub_item:
                        file_string += f"  {sub_item.upper()}\t\t"
                    else:
                        file_string += f"{sub_item}\t"
                file_string += "\n"
            # if len(item) == 2:
            #     file_string += f"  {item[0].upper()}\t\t{item[1]}\n"
            # else:
            #     file_string += f"  {item[0].upper()}\t\t{item[1]}\t\t\t{item[2]}\t{item[3]}\n"
        
        f = open(file_name, 'w')
        f.write(file_string)
        f.close()
        print(f"INFO: MSX FILE {file_name} generated.")
        
        # print(file_string) ## will be replaced with output step later
        

# =============================================================================
# #################### STAND ALONE FUNCTIONS #################################
# =============================================================================

def MSXBinReader(filename, epanetinpfile):
    wn = wntr.network.WaterNetworkModel(epanetinpfile)
    duration = int(wn.options.time.duration)
    
    with open(filename, 'rb') as fin:
          ftype = '=f4'
          idlen = 32
          prolog = np.fromfile(fin, dtype = np.int32, count=6)
          magic1 = prolog[0]
          version = prolog[1]
          nnodes = prolog[2] 
          nlinks = prolog[3]
          nspecies = prolog[4]
          reportstep = prolog[5]
          
          node_list = wn.node_name_list
          link_list = wn.link_name_list
          
          species_list = []
          species_mass = []
          
          if version >= 200000:
              for i in range(nspecies):
                  species_len = int(np.fromfile(fin, dtype=np.int32, count=1))
                  species_name = ''.join(chr(f) for f in np.fromfile(fin, dtype=np.uint8, count=species_len) if f!=0)
                  species_list.append(species_name)
       
              # for i in range(nspecies):
                  species_mass.append(''.join(chr(f) for f in np.fromfile(fin, dtype = np.uint8, count=16) if f != 0))
          else:
          # elif version == 100000:
              # older version of MSX had names then units, rather than name-unit/name-unit as in 2.0 or newer
              for i in range(nspecies):
                  species_len = int(np.fromfile(fin, dtype=np.int32, count=1))
                  species_name = ''.join(chr(f) for f in np.fromfile(fin, dtype=np.uint8, count=species_len) if f!=0)
                  species_list.append(species_name)
       
              for i in range(nspecies):
                  species_mass.append(''.join(chr(f) for f in np.fromfile(fin, dtype = np.uint8, count=16) if f != 0))
                      
                
          timerange = range(0, duration+1, reportstep)
          tr = len(timerange)
          
          # print(species_list, species_mass)
          
          row1 = ['node']*nnodes*len(species_list)+['link']*nlinks*len(species_list)
          row2 = []
          for i in [nnodes,nlinks]:
                  for j in species_list:
                        row2.append([j]*i)
          row2 = [y for x in row2 for y in x]
          row3 = [node_list for i in species_list] + [link_list for i in species_list]
          row3 = [y for x in row3 for y in x]    
          
          tuples = list(zip(row1, row2, row3))
          index = pd.MultiIndex.from_tuples(tuples, names = ['type','species','name'])
          
          try:
                  data = np.fromfile(fin, dtype = np.dtype(ftype), count = tr*(len(species_list*(nnodes + nlinks))))
                  data = np.reshape(data, (tr, len(species_list*(nnodes + nlinks))))
          except Exception as e:
              print(e)
              print ("oops")
             
          postlog = np.fromfile(fin, dtype = np.int32, count=4)
          offset = postlog[0]
          numreport = postlog[1]
          errorcode = postlog[2]
          magicnew = postlog[3] 
          if errorcode !=0:
              print(f'ERROR CODE: {errorcode}')
              print(offset, numreport)
              
          df_fin = pd.DataFrame(index=index, columns=timerange).transpose()
          if magic1 == magicnew:
              # print("Magic# Match")
              df_fin = pd.DataFrame(data.transpose(), index=index, columns=timerange)
              df_fin = df_fin.transpose()
              
          else:
              print("Magic#s do not match!")
    return df_fin


def line_parser(line):
    # print(line)
    note = ''
    if ';' in line:
        head, note = line.split(';')
    else:
        head = line
    
    if '\t' in head:
        a = re.split(r"[\t ]", head)
    else:
        a = head.split(' ')
        
    parsed_list = [i for i in a if i != ''] # removes empty blanks
    
    return parsed_list, note

def line_parser_terms(line):
    note = ''
    parsed_list = []
    if ';' in line:
        head, note = line.split(';')
    else:
        head = line
    
    if len(head) > 5: # must have enough information to actually have data
        if '\t' in head:
            a = re.split(r"[ \t]", head)
        else:
            a = head.split(' ')
            
        int_parsed = [i for i in a if i != '']
        
        parsed_list = ['', '']
        parsed_list[0] = int_parsed[0] ### should be term name
        
        value = "".join(int_parsed[1:]) ## join remaining into 2nd term equation
        parsed_list[1] = value
        
    return parsed_list, note

def line_parser_eq(line):
    note = ''
    parsed_list = []
    if ';' in line:
        head, note = line.split(';')
    else:
        head = line
    
    if len(head) > 5: # must have enough information to actually have data
        if '\t' in head:
            a = re.split(r"[ \t]", head)
        else:
            a = head.split(' ')
            
        int_parsed = [i for i in a if i != '']
        
        parsed_list = ['', '','']
        parsed_list[:2] = int_parsed[:2]
        
        value = "".join(int_parsed[2:])
        parsed_list[2] = value
    
    return parsed_list, note

def _read_msx_file_to_dict(filename):
    #### TODO: NEEDS TO HANDLE type NODE/LINK SPECIES VALUE ; NOTE STYLE

    f = open(filename, 'r')
    file = f.read()
    
    read_in_dict = {}
    key_val = None
    split_file = file.split('\n')
    
    for line in split_file:
        if len(line) > 3: ## only bothers doing something if there is real data in the line
            if '[' in line and ']' in line:
                st_idx = line.index('[')
                en_idx = line.index(']')
                key_val = line[st_idx+1: en_idx].lower()
                
                read_in_dict[key_val] = {}
                if key_val in ['parameters', 'sources', 'quality', 'report']:
                    read_in_dict[key_val] = []
                
            else: # do something with the data
                if key_val == 'title':
                    read_in_dict[key_val] = line
                elif key_val == 'options':
                    parsed, note = line_parser(line)
                    key, val = parsed
                    read_in_dict[key_val][key] = {'val': val,
                                                  'note': note}
                elif key_val in ['species', 'coefficients']:
                    parsed, note = line_parser(line)
                    type_name, key, val = parsed
                    read_in_dict[key_val][key] = {'type': type_name, 
                                                  'val': val,
                                                  'note': note}
                elif key_val == 'terms':
                    parsed, note = line_parser_terms(line)    
                    
                    key, val = parsed
                    read_in_dict[key_val][key] = {'val': val,
                                                  'note': note}
                elif key_val in ['pipes', 'tanks']:
                    parsed, note = line_parser_eq(line)
                    type_name, key, val = parsed
                    read_in_dict[key_val][key] = {'type': type_name,
                                                  'val': val,
                                                  'note': note}
                elif key_val in ['sources', 'quality', 'parameters']:
                    parsed, note = line_parser(line)
                    read_in_dict[key_val].append(parsed + [note])
                elif key_val == 'report':
                    parsed, note = line_parser(line)
                    read_in_dict[key_val].append(parsed)
                elif key_val == 'patterns':
                    parsed, note = line_parser(line)
                    if len(parsed) > 0:
                        pattern = parsed[0]
                        if pattern not in read_in_dict[key_val].keys():
                            read_in_dict[key_val][pattern] = []
                        for item in parsed[1:]:
                            read_in_dict[key_val][pattern].append(item)
    
    return read_in_dict



def run_msx(inpfile, msxfile, engine='32', version='2'):
    cwd = os.getcwd()
    if version == '2' or version == 2:
        if engine == '32' or engine == '64':
            os.chdir(f'{msx_location}/EPANET_MSX_2/Release{engine}')
        elif engine == 32 or engine == 64:
            os.chdir(f'{msx_location}/EPANET_MSX_2/Release{engine}')
    elif version == '11' or version == '1.1' or version == 1.1:
        os.chdir(f"{msx_location}/EPANETMSX")
    
    job_base, _ = msxfile.split('.')
    try:
        subprocess.Popen(['runepanetmsx.exe', 
                          cwd+'/'+inpfile, 
                          cwd+'/'+msxfile,
                          cwd+'/'+job_base+'.rpt', 
                          cwd+'/'+job_base+'.bin']).communicate()
    except Exception as e:
        print('ERROR: Failure to Run EPANET-MSX ')
        print(e)
        # os.chdir(cwd)
    os.chdir(cwd)
    
    return cwd+'/'+job_base+'.bin', cwd+'/'+job_base+'.rpt'

