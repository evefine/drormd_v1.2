"""
Methods for specifying a set of simple metrhics, running them,
and saving them to files.
"""

from drormd import analysis
import pandas
import numpy as np

def write_to_csv(data, fname):
	pandas.DataFrame.from_dict(data).to_csv(fname, mode='w')

def add_distances(data, distances, molid):
	for name, (sel1, sel2) in distances.items():
		data[name] = analysis.min_distance(sel1, sel2, molid)

def add_angles(data, angles, molid):
	for name, (sel1, sel2, sel3) in angles.items():
		data[name] = analysis.angle(sel1, sel2, sel3, molid)

def add_dihedrals(data, dihedrals, molid):
	for name, (sel1, sel2, sel3, sel4) in dihedrals.items():
		data[name] = analysis.dihedral(sel1, sel2, sel3, sel4, molid)

def add_rmsd_to_avg(data, rmsd_to_avg, molid):
	for name, (rmsd_sel, align_sel) in rmsd_to_avg.items():
		analysis.align(align_sel, molid)
		data[name] = analysis.rmsd_to_average(rmsd_sel, molid)


def add_rmsds(data, rmsds, molids_dict, molid):
	#if len(rmsds.items()[1]) == 4:
		#sel_ref = rmsds.items()[3]
	#else:
		#sel_ref = rmsds.items()[1]
	for name, (rmsd_sel, align_sel, name_ref) in rmsds.items():
		analysis.align(align_sel, molid, molid_ref = molids_dict[name_ref])
		data[name] = analysis.rmsd(rmsd_sel, molid, molid_ref = molids_dict[name_ref])

# EJ added
def add_centroids(data, distances, molid):
	for name, (sel1, sel2) in distances.items():
		data[name] = analysis.centroid_distance(sel1, sel2, molid)

def add_centroids_noz(data, distances, molid):
        for name, (sel1, sel2) in distances.items():
                data[name] = analysis.centroid_distance_noz(sel1, sel2, molid)

def add_centroids_z(data, distances, molid):
        for name, (sel1, sel2) in distances.items():
                data[name] = analysis.centroid_distance_z(sel1, sel2, molid)


def add_angle_btwn_vectors(data, vector_grps, molids_dict, align_info):
        for name, vectors in vector_grps.items():
                angle_difs = []
                ref_mol = align_info['ref'][2]
                for sel1, sel2, refsel1, refsel2 in vectors:
                        analysis.align(align_info['ref'][0], molids_dict['self'],sel_ref=align_info['ref'][1], molid_ref=molids_dict[ref_mol])
                        angle_difs.append(analysis.angle_btwn_vector(sel1, sel2, refsel1, refsel2, molids_dict[ref_mol], molids_dict['self']))
                angle_difs = np.stack(angle_difs,axis=0)
                avg_angle_difs = np.mean(angle_difs,axis=0)
                data[name] = avg_angle_difs



def resolve_selections(metrics, condition, generic, conditions):
    # Create merged dictionary of selections.
    selections = {**generic}  # Merge dictionaries with generic selections
    if 'selections' in conditions[condition]:
        for name, sel in conditions[condition]['selections'].items():
            selections[name] = sel.format(**generic)

    # Helper function to recursively resolve selections
    def resolve_item(item):
        #print('resolving ', item)
        if isinstance(item, str):
            return item.format(**selections)
        elif isinstance(item, (list, tuple)):  # Check for both lists AND tuples
            # Convert result to the same type as input (list or tuple)
            return type(item)(resolve_item(sub_item) for sub_item in item)
        else:
            return item

    # Create a new version of metrics where all `str`s have been formatted
    resolved = {}
    for k1, v1 in metrics.items():
        #print('step1',k1,v1)
        resolved[k1] = {}
        for k2, v2 in v1.items():
            #print('step2',k2,v2)
            resolved[k1][k2] = resolve_item(v2)

            # Remove entries that contain 'N/A'
            if isinstance(resolved[k1][k2], (list, tuple)):
                resolved[k1][k2] = type(resolved[k1][k2])(
                    x for x in resolved[k1][k2] if 'N/A' not in str(x)
                )

    return resolved











def compute(fname, topology, trajectory, metrics, structs):
    data, molids = {}, {}
    
    # Load trajectory and any structures to be used in rmsd computations.
    molids['self'] = analysis.load_traj(topology, trajectory)
    for struct, path in structs.items():
        molids[struct] = analysis.load_struct(path) 
    # Add all metrics that are specified.
    if 'dihedrals' in metrics:
        add_dihedrals(data, metrics['dihedrals'], molids['self'])
    if 'angles' in metrics:
        add_angles(data, metrics['angles'], molids['self'])
    if 'distances' in metrics:
    	add_distances(data, metrics['distances'], molids['self'])
    if 'rmsds' in metrics:
       	add_rmsds(data, metrics['rmsds'], molids, molids['self'])
    # EJ added this
    if 'avg_rmsds' in metrics:
        add_rmsd_to_avg(data, metrics['avg_rmsds'], molids['self'])
    if 'avg_vector_angles' in metrics:
        add_angle_btwn_vectors(data,metrics['avg_vector_angles'],molids,metrics['align'])    
    if 'centroid_distances' in metrics:
        add_centroids(data, metrics['centroid_distances'], molids['self'])
    if 'centroid_distances_noz' in metrics:
        add_centroids_noz(data, metrics['centroid_distances_noz'], molids['self'])
    if 'centroid_distances_z' in metrics:
        add_centroids_z(data, metrics['centroid_distances_z'], molids['self'])
    write_to_csv(data, fname)

def setup(script_path, conditions, metrics, generic):
	"""
	Use this command to preview the resolved metrics and
	print the commands that need to be run.

	Pipe the output to a file when you're ready to run.
	"""
	# Resolved metrics for each condition
	# This is just for double checking purposes
	for condition in conditions:
		print('# {}'.format(condition))
		resolved = resolve_selections(metrics, condition, generic, conditions)
		for k1, v1 in metrics.items():
			print('# {}'.format(k1))
			for k2 in v1:
				print('# Original: {}'.format(metrics[k1][k2]))
				if k2 in resolved[k1]:
					print('# Resolved: {}'.format(resolved[k1][k2]))
				else:
					print('# Resolved: {}'.format('N/A'))

	# Print commands that need to be run.
	for condition, info in conditions.items():
		for replicate in info['reps']:
			print('python {} {} {}'.format(script_path, condition, replicate))

	print('# setup_jobarray.sh run.sh drormd 00:30:00 12GB 10 1  > run.sbatch')

def main(args, script_path, conditions, metrics, generic={}, structs={}, base=''):
	if len(args) == 1:
		assert args[0] == 'setup'
		setup(script_path, conditions, metrics, generic)
	else:
		condition, replicate = args

		fname = '{}{}_r{}.csv'.format(base, condition, replicate)
		top = conditions[condition]['topology']
		print(conditions[condition]['trajectory'])
		print(replicate)
		try:
			traj = conditions[condition]['trajectory'].format(replicate)
		except:
			traj = conditions[condition]['trajectory'].format(replicate,replicate)		
		selections = resolve_selections(metrics, condition, generic, conditions)

		compute(fname, top, traj, selections, structs)
