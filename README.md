# drormd_v1.2
drormd but with my edits and arguably more documentation

Documentation for Joe's code (analysis.py and wrappers.py)

Basically, the architecture is that you will provide a metrics
and conditions file, wrappers.py will substitute all of the {} in
your metrics files with the dictionary in conditions, and then it
will compute the relevant metrics using functions from the analysis.py file.

I'm goin to go through each function that is not just a simple
readable helper function so that everything is documented and if
changes need to be made the script can be understood better.

I'm going to go through the workflow so it's easily understood.
You start with resolve_selections in wrappers.py. 
k1 and v1 just refer to the key adn value in metrics, so this
would be somthing like 'distances' and {'dist1': ['{sel1}','{sel2}], 'dit2':['{sel3}','{sel4}]}
It loops through these different types of metrics, and then loops through the different
measurements of this type that you are taking and resolves the item by formating
with the generic and selections dictionaries provided. k2 adn v2 refer to 'dist1' and 
['{sel1}','{sel2}] for example. The script will fix v2 as long as it's either a string or a 
list of strings (or list of lists). It will return these in the same format as they 
came, but now just withthe dictionary replacements.

Next we go to compute in wrappers. It will load in each struct provided. It will have a 
molid (molecule ID) associated with it in a dictionary with the key being the name you
gave to the struct in the input. Basically this function will go through and see what
metrics you have (dihedrals, angles, distances, rmsds, avg_vector_angles, centroid_distances)
and call the appropriate helper function. In this, molids is the dictionary of molids (so of
the simulation would be molids['self']). metrics['key'] will return the corresponding dictionary
in your metrics file (with the resolved selections). data refers to what you will get out
and eventually write to a csv. molids['key'] will give a number that you will put into
the molids dictionary to get out the coordinates from the molecule you want.

Each of the helper functions in wrappers will call a function in analysis that will
actually compute your desired metric. In each of these, they unpack each metric by
name (sel1, sel2, ...) in measurement.items(), which in my previous example would be
'dist1' and ['selection1', 'selection2'] for one iteration. It will then pass this information
onto a function in analysis.

Each function in analaysis is fairly well documented/can be figured out especially
with this information. The only caveat is that for add_angle_btwn_vectors, to compare
the angle between two vectors between a simulation and some reference structure, you 
need to provide an align dictionary, which is explained in metrics_example.py.
This is a MUST! Can be coded to be different but I don't feel like it right now.
