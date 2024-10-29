import sys
import os
from drormd import wrappers
from conditions import conditions, generic

metrics = {
    # these are single atoms that we are selecting
    'distances': {
        # TM6 Intracellular Position
        '6.31-4.42': ('{6x31} and name CA', '{4x42} and name CA'),
        '6.31-3.50': ('{6x31} and name CA', '{3x50} and name CA')
    },

    'dihedrals': {
        '7.53-chi1': ('{7x53} and name N',
                      '{7x53} and name CA',
                      '{7x53} and name CB',
                      '{7x53} and name CG'),
        '3.50-chi2': ('{3x50} and name CA',
                      '{3x50} and name CB',
                      '{3x50} and name CG',
                      '{3x50} and name CD')},
    
    'angles':
        {'6.50-6.43-6.35': ('{6x50} and name CA',
                            '{6x43} and name CA',
                            '{6x35} and name CA'),
        '3.50-ang': ('{3x50} and name CA',
                            '{3x50} and name CG',
                            '{3x50} and name CZ'),
        '3.50-7.56-7.54': ('{3x50} and name CA',
                            '{7x56} and name CA',
                            '{7x54} and name CA')}
    
    # the first argument is the thing you want to compare the RMSD
    # of, the second argument is how you are aligning the simulation
    # and the reference, and the third is what molecule is the
    # reference molecule. 'self' refers the the same molecule.

    'rmsds': {
        'UK5099': ('{UK} and noh, {prot} and noh', 'self'), 
    },  

    # you must keep an 'align' dictionary if you are using avg_vector_angles
    # avg_vector_angles value is a nested list, each of which contains 4 atoms.
    # the first two atoms make up the first vector on the simulation molecule (self)
    # and the second two atoms make up the first vector on the reference molecule 
    # eventually the two vectors will be computed and the angle between them
    # will be computed. Then it will average the angle across all quadruples.
    'avg_vector_angles': {'TM7 twist':[['{7x49} and name CA','{7x50} and name CA','{7x49_ref} and name CA','{7x50_ref} and name CA' ],
                                        ['{7x50} and name CA','{7x51} and name CA','{7x50_ref} and name CA','{7x51_ref} and name CA' ],
                                        ['{7x51} and name CA','{7x52} and name CA','{7x51_ref} and name CA','{7x52_ref} and name CA' ],
                                        ['{7x52} and name CA','{7x53} and name CA','{7x52_ref} and name CA','{7x53_ref} and name CA' ]]
                        }
    # the key of the align dictionary MUST be 'ref', but everything else can change
    # the first entry is the selection you will align on in the simulation molecule,
    # the second entry is the selection you will align on in the reference molecule,
    # and the third entry is the name of the entry in structs that you will be comparing to
    # (it can be self, then it will just compare to the first frame)
    'align':{'ref':['({sim_selection}) and name CA', '({ref_selection}) and name CA','ref'] }
    },  
    # if there's a metric you need but no example, I would go and see what the function needs in wrappers.py
    # in add_{metric}, and you can see what each input is in analysis.py
}

structs = {'ref': '/oak/stanford/groups/rondror/users/evejfine/simulation_analysis/oprm_mouse_6DDF.pdb'}


wrappers.main(sys.argv[1:], os.path.abspath(__file__), conditions, metrics,
              generic=generic, structs=structs)

            