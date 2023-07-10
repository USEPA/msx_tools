# -*- coding: utf-8 -*-
"""
Created on Tue Jun 15 10:02:56 2021

MSX Builder - Python Script for building and running MSX Files


NOTE: Values in libarary are for test purposes ONLY. These highlight the reaction form/relationship but not necessarily the specific system kinetics. Citations are provided where appropriate for sources of values or models. In most cases updates will be necessary to use for different models -- That is, sources, pipes, nodes are for a given network, and must be updated for different models.



@author: JBurkhar
"""
import numpy as np

mins = 1 * 60
injL = 10 # 10 minute long injection
P2 = np.zeros(60*mins)
P2[:injL*60] = 1.

injectConc = np.round(0.4 * 22/0.1321,7)


msx_nicotine_1 = {'title': 'MSX Builder - Python Generated File',
                  'description': 'Nicotine - Chlorine reaction\nSources for Test Model, please update sources.', 
                  'options': {'area_units': {'val': 'M2',   'note': 'Surface concentration is mass/m2'},
                              'rate_units': {'val': 'MIN',  'note': 'Reaction rates are concen/sec'},
                              'solver':     {'val': 'RK5',  'note': '5-th order Runge-Kutta integrator'},
                              'timestep':   {'val': 1,      'note': 'seconds'},
                              'COUPLING':   {'val': 'NONE', 'note': 'Default NONE (???)'},
                              'RTOL':       {'val': 1e-10,  'note': 'Relative concentration tolerance'},
                              'ATOL':       {'val': 1e-10,  'note': 'Absolute concentration tolerance'},
                              #'COMPILER':   {'val': 'VC',   'note': ''},
                              'SEGMENTS':   {'val': 5000,   'note': 'Max # of segments per pipe: Default value is 5000'},
                              'PECLET':     {'val': 1000.,  'note': 'default value is 1000 for Peclet Threshold for applying dispersion'}
                              },
                  'species': {'Nx':     {'type': 'BULK',        'val': 'MG',        'note': 'Nicotine'},
                              'HOCL':   {'type': 'BULK',        'val': 'MG',        'note': 'Free Chlorine'}},
                  'coefficients': {'kd': {'type': 'constant',   'val': 2.33e-3,     'note': 'Decay rate : min$^-$$^1$)'},
                                   'K1': {'type': 'constant',   'val': 5.92e-2,     'note': 'decay constant for chlorine: L (mg Nic)$^-$$^1$ min$^-$$^1$'},
                                   'K2': {'type': 'constant',   'val': 1.84e-1,     'note': 'decay constant for nicotine: L (mg Cl)$^-$$^1$ min$^-$$^1$'}},
                  'pipes': {'Nx':        {'type': 'RATE',       'val': '-RXN',      'note': ''},
                            'HOCL':      {'type': 'RATE',       'val': '-RXCL',     'note': ''}},
                  'tanks': {'Nx':        {'type': 'RATE',       'val': '0',         'note': ''},
                            'HOCL':      {'type': 'RATE',       'val': '0',         'note': ''}},
                  'terms': {'RXCL':      {'val': 'kd*HOCL + K1*Nx*HOCL',            'note': ''},
                            'RXN':       {'val': 'K2*Nx*HOCL',                      'note': ''}},
                  'sources': [['CONCEN', '17',      'HOCL', 1.01],
                              ['CONCEN', 'InjPort', 'Nx',   injectConc, 'P2']], ### might want to delete this, and have it handled seperately. 
                  'dispersion': {'HOCL': {'val': 1.0, 'extra': '', 'note': ''},
                                'Nx':    {'val': 1.0, 'extra': '', 'note': ''}},  # relative value compared to chlorine
                  'quality': [['GLOBAL', 'HOCL', 1.01]],
                  'patterns': {'P2': P2},
                  'report': [['NODES', 'ALL'],
                             ['SPECIE','HOCL','YES','5'],
                             ['SPECIE','Nx', 'YES','5']]
                  }


injL = 10 # 10 minute long injection
P2 = np.zeros(60*mins)
P2[:injL*60] = 1.

msx_nicotine_2 = {'title': 'MSX Builder - Python Generated File',
                  'description': 'Nicotine - Chlorine reaction with reactive intermediate\nSources for test using EPANET Net3 Model, please update sources',
              'options': {'area_units': {'val': 'M2',   'note': 'Surface concentration is mass/m2'},
                          'rate_units': {'val': 'MIN',  'note': 'Reaction rates are concen/sec'},
                          'solver':     {'val': 'RK5',  'note': '5-th order Runge-Kutta integrator'},
                          'timestep':   {'val': 1,      'note': 'seconds'},
                          'RTOL':       {'val': 1e-10,  'note': 'Relative concentration tolerance'},
                          'ATOL':       {'val': 1e-10,  'note': 'Absolute concentration tolerance'},
                          #'COMPILER':   {'val': 'VC',   'note': ''},
                          'SEGMENTS':   {'val': 5000,   'note': 'default value is 5000'}},
              'species': {'Nx':     {'type': 'BULK',        'val': 'MG',    'note': 'Nicotine'},
                          'HOCL':   {'type': 'BULK',        'val': 'MG',    'note': 'Free Chlorine'},
                          'NX2':    {'type': 'BULK',        'val': 'MG',    'note': 'Int. Nicotine Reactive'}},
              'coefficients': {'kd': {'type': 'constant',   'val': 3.0e-5,  'note': 'Decay rate : (min$^-$$^1$)'},
                               'K1': {'type': 'constant',   'val': 9.75e-2, 'note': 'decay constant for chlorine: L (mg Nic)$^-$$^1$ min$^-$$^1$'},
                               'K2': {'type': 'constant',   'val': 5.73e-1, 'note': 'decay constant for nicotine: L (mg Cl)$^-$$^1$ min$^-$$^1$'},
                               'K3': {'type': 'constant',   'val': 1.34e-2, 'note': 'decay constant for nicotine: L (mg N2)$^-$$^1$ min$^-$$^1$'},
                               'K4': {'type': 'constant',   'val': 2.19e-2, 'note': 'decay constant for nicotine: L (mg Cl)$^-$$^1$ min$^-$$^1$'},
                               },
              'pipes': {'Nx':   {'type': 'RATE',        'val': '-RXN',  'note': ''},
                        'HOCL': {'type': 'RATE',        'val': '-RXCL', 'note': ''},
                        'NX2':  {'type': 'RATE',        'val': 'RXNX2', 'note': ''}},
              'tanks': {'Nx':   {'type': 'RATE',        'val': '-RXN',  'note': ''},
                        'HOCL': {'type': 'RATE',        'val': '-RXCL', 'note': ''},
                        'NX2':  {'type': 'RATE',        'val': 'RXNX2', 'note': ''}},
              'terms': {'RXCL':     {'val': 'kd*HOCL + K1*Nx*HOCL + K3*NX2*HOCL', 'note': ''},
                        'RXN':      {'val': 'K2*Nx*HOCL', 'note': '', 'note': ''},
                        'RXNX2':    {'val': 'K2*Nx*HOCL - K4*NX2*HOCL', 'note': ''}},
              'sources': [['CONCEN','River','HOCL', 1.01],
                          ['CONCEN','Lake', 'HOCL', 1.01],
                          ['CONCEN', '120', 'Nx', 3e5, 'P2']],
              'dispersion': {'HOCL': {'val': 1.0, 'extra': '', 'note': ''},
                             'Nx':   {'val': 1, 'extra': '', 'note': ''},
                             'NX2':  {'val': 1, 'extra': '', 'note': ''}},  # relative value compared to chlorine
              'quality': [['GLOBAL', 'HOCL', 1.01]],
              'patterns': {'P2': P2},
              'report': [['NODES', 'ALL'],
                         ['SPECIE','HOCL','YES','5'],
                         ['SPECIE','Nx', 'YES','5'],
                         ['SPECIE', 'NX2', 'YES', '5']]
                        
              }

msx_lead = {'title': 'MSX Builder - Python Generated File',
            'description': 'Lead Plumbosolvency Model (from Burkhardt et al 2020)\nParameters for EPA HPS Simulator Model',
              'options': {'area_units': {'val': 'M2',   'note': 'Surface concentration is mass/m2'},
                          'rate_units': {'val': 'SEC',  'note': 'Reaction rates are concen/sec'},
                           'solver':     {'val': 'RK5',  'note': '5-th order Runge-Kutta integrator'},
                          # 'solver':     {'val': 'ROS2',  'note': 'a second order Rosenbrock 2(1) method'},
                          #'solver': {'val': 'EUL', 'note': 'Eulerian Solver?'},
                          'timestep':   {'val': 1,      'note': 'seconds'},
                          'COUPLING': {'val': 'NONE', 'note': 'Default NONE (???)'},
                          'RTOL':       {'val': 1e-8,  'note': 'Relative concentration tolerance'},
                          'ATOL':       {'val': 1e-8,  'note': 'Absolute concentration tolerance'},
                          # 'COMPILER':   {'val': 'VC',   'note': ''},
                          'SEGMENTS':   {'val': 5000,   'note': 'Max # of segments per pipe: Default value is 5000'},
                          'PECLET': {'val': 1000., 'note': 'default value is 1000 for Peclet Threshold for applying dispersion'}
                          },
              'species': {'PB2':     {'type': 'BULK', 'val': 'UG', 'note': 'Dissolved Lead (Pb)'},
                          },
              'coefficients': {'M': {'type': 'constant', 'val': 0.117, 'note': 'Desorption Rate (ug/m^2/s)'},
                               'E': {'type': 'constant', 'val': 140.0, 'note': 'Saturation/Plumbosolvency Level (ug/L)'},
                               'F': {'type': 'PARAMETER', 'val': 0.0, 'note': 'Determines which pipes have reactions'}},
              'pipes': {'PB2': {'type': 'RATE', 'val': 'F*Av*M*(E - PB2)/E', 'note':''},
                        },
              'tanks': {'PB2': {'type': 'RATE', 'val': 0, 'note': ''}},
              'terms': {},
              'sources': [],
              'dispersion': {'PB2': {'val': 1.0, 'extra':'', 'note':''}},
              'parameters': [['PIPE', '1', 'F', 1.0, 'Sets Desorption from Lead Pipe']],
              'quality': [],
              'patterns': {},
              'report': [[]]
                           # [['NODES','ALL'],
                           # ['LINKS','ALL'],
                           # ['SPECIE','PB2','YES','5']]
              }


msx_dict = {'nicotine': msx_nicotine_1,
            'nicotine_ri': msx_nicotine_2,
            'lead_ppm': msx_lead}
