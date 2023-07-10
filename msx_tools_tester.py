# -*- coding: utf-8 -*-
"""
Created on Wed Jun 29 11:19:51 2022

@author: JBurkhar
"""

import msx_tools
import numpy as np


test1 = True ## Build from scratch
test2 = True ## Read from msx file
test3 = True ## build from preconfigured dictionary

test4 = True ## run msx file with msx_tools.run_msx()

test5 = True ## Lead model



if test1:
    ### Test 1: Build from Scratch
    print('Test #1 - Building MSX file from scratch')
    msxobj = msx_tools.MSXobj()
    
    
    msxobj.add_species('Nx', 'MG', 'Nicotine', 'BULK')
    msxobj.add_species('HOCL', 'MG', 'Chlorine', 'BULK')
    
    msxobj.add_coefficient('kd', 3e-5, 'Decay rate (min^-1)', 'CONSTANT')
    msxobj.add_coefficient('K1', 0.056, 'decay constant for chlorine: L (mg Nic)^-1 min^-1', 'CONSTANT')
    msxobj.add_coefficient('K2', 0.157, 'decay constant for nicotine: L (mg Cl)^-1 min^-1', 'CONSTANT')
    
    msxobj.add_term('RXCL', 'kd*HOCL + K1*Nx*HOCL')
    msxobj.add_term('RXN', 'K2*Nx*HOCL', 'Nicotine Reaction Rate')
    
    msxobj.add_pipe_rxn('RATE', 'Nx', '-RXN')
    msxobj.add_pipe_rxn('RATE', 'HOCL', '-RXCL')
    
    msxobj.add_tank_rxn('RATE', 'Nx', 0)
    msxobj.add_tank_rxn('RATE', 'HOCL', 0)
    
    msxobj.add_source('CONCEN', '17', 'HOCL', 0.88)
    msxobj.add_source('CONCEN', 'InjPort', 'Nx', 66.6162, 'P2')
    
    msxobj.add_quality('GLOBAL', 'HOCL', 0.88)
    
    # msxobj.add_parameter('PIPE', 2, 'F', 1.0, 'LSL')
    
    pattern = np.zeros(60*60) ## 60 minute simulation
    msxobj.add_pattern('P2', pattern)
    
    msxobj.add_dispersion('Nx', 100)
    
    msxobj.add_report_type('NODES', 'ALL')
    msxobj.add_report_item('SPECIE', 'HOCL')
    msxobj.add_report_item('SPECIE','Nx', sig_figs=6)
    
    
    msxobj.build_msx_file()
# 
    # msxobj.build_msx_file(style='MSX11')


if test2:
    print('Test #2: Read from existing MSX file')
    ## Not finished
    ## Read from MSX file
    msxobj = msx_tools.MSXobj(file_name='temp.msx')
    
    msxobj.build_msx_file(file_name='new.msx')
    

if test3:
    print('Test #3: Build from a preconfigured dictionary')
    ## build from preconfigured dictionary 
    

    msxobj = msx_tools.MSXobj(rxn_type='nicotine')
    # print(msxobj.msx_info_dict)
    
    msxobj.build_msx_file()
    
if test4:
    print('Test #4: Run MSX')
    inp_fn = 'input_files/updated.inp'  ## name of input file
    bin_name, rpt_name = msx_tools.run_msx(inp_fn, 'new.msx')
    
    rpt_file = open(rpt_name, 'r').readlines()
    no_errors = True
    for line in rpt_file:
        if 'Error' in line:
            print('ERROR found in report file')
            no_errors = False
    
    if no_errors:
        # msxobj = msx_tools.MSXobj(file_name='new.msx')
        # msxobj.species.keys()
        results = msx_tools.MSXBinReader(bin_name, inp_fn) 
        ### MSXBinReader must be provided with bin_file_name and associated
        ### EPANET INP file
        #### -- Saved as a results object ('node', 'link') ('Species') ('nodeID' or 'linkID')
        
        species = list(set(results.columns.get_level_values('species').values))
        for specie in species:
            results.node[specie].plot.line()

if test5:
    print('Test #5: Read in and change item')
    msxobj = msx_tools.MSXobj(rxn_type='lead_ppm')
    msxobj.update_coefficient('M', 0.12)
    msxobj.build_msx_file(file_name='lead_test.msx')
        
