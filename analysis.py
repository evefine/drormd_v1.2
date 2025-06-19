"""
Functions that handle direct interfacing with the MD trajectories through
vmd-python.
"""

import numpy as np
from vmd import molecule, atomsel, vmdnumpy

# I/O.
def load_traj(topology, traj, start=0, stop=-1, stride=1):
	"""
	Loads a PSF structure file and an NC trajectory path.
	By default, loads in all frames from start to finish. Use `start`, `stop`,
	and `stride` to adjust which frames to load. By default, this function
	will not return until the molecule is fully loaded.

	Args:
        topology (str): Path to a PSF file. Must have extension '.psf'.
        traj (str): Path to a netcdf file. Must have extension '.nc'.
        start (int, optional): First frame to load. Defaults to 0.
        stop (int, optional): Last frame to load. Defaults to -1, the last frame.
        stride (int, optional): Load every `stride` frames. Defaults to 1.

    Returns:
        int: The molecule ID of the newly-imported molecule.
	"""
	if traj[-3:] == '.nc':
		traj_filetype = 'netcdf'
	elif traj[-4:] == '.dcd':
		traj_filetype = 'dcd'
	else:
		assert false
	assert topology[-4:] in ['.psf','.pdb']
	topology_filetype = topology.split('.')[-1]
	
	trajid = molecule.load(topology_filetype, topology)
	molecule.read(trajid, traj_filetype, traj,
	              beg=start, end=stop, skip=stride, waitfor=-1)
	return trajid

def load_struct(struct):
	"""
	Loads an MAE or PDB formatted structure file.

	Args:
        struct (str): The first parameter.
        trajectory (str): Path to a netcdf file. Must have extension '.nc'.
    Returns:
        int: The molecule ID of the newly-imported molecule.
	"""
	filetype = struct.split('.')[-1]
	assert filetype in ['pdb', 'mae'], \
		'{} is not a pdb or mae file extension'.format(filetype)
	
	trajid = molecule.load(filetype, struct)
	return trajid

def write_pdb(outpath, molid, frame):
	"""
	Writes frame `frame` of molecule `molid` to `outpath`.

	In the future, this could be updated to allow writing of
	multiple frames or only a subset of the atoms.

	Args:
		outpath (str): Path to write PDB file.
		molid (int): Molecule for which to write coordinates.
		frame (int): Frame for which to write coordinates.

	Returns:
		None.
	"""
	vmd.molecule.write(molid, "pdb", outpath, beg=frame, end=frame)

# Analysis.
def get_coords(sel, molid, frames=None):
	"""
	Returns an (num frames x num atoms x 3) numpy array of the
	xyz coordinates of atoms in sel for the molecule molid.

	Args:
		sel (str): VMD atom selection for which to extract coordinates.
		molid (int): molid of a loaded Molecule.
		frames (iterable yielding int, optional): Iterable of frames for which
			to compute coordinates. Default is to compute coordinated for all
			available frames (None).
			Ex: for first frame only use [0].
			Ex: for every tenth frame use range(0, molecule.numframes(molid), 10).
	Returns:
		np.array(float): shape is (num frames x num atoms x 3).
	"""
	# Use all the frames by default.
	print(sel,molid,range(molecule.numframes(molid)))
	if frames is None:
		frames = range(molecule.numframes(molid))
	
	mask = vmdnumpy.atomselect(molid=molid, frame=0, selection=sel)
	if not np.any(mask):
		raise ValueError('Atomsel "{}" references no atoms.'.format(sel))

	return np.stack([np.compress(mask, vmdnumpy.timestep(molid, frame), axis = 0)
	                 for frame in frames])

def align(sel, molid, sel_ref=None, molid_ref=None, frame_ref=0):
	"""
	Aligns molecule molid on atomsel sel.

	Args:
		sel (str): VMD atom selection to align on.
		molid (int): Specifies molecule to align.
		sel_ref (str, optional): VMD atom selection to use for the
			reference structure. Defaults to `sel`.
		molid_ref (str, optional): Specifies molecule to use as reference
			structure. Defaults to `molid`.
		frame_ref (int, optional): Frame of `molid_ref` to use as the
			reference structure. Defaults to 0.
	Returns:
		None.
	"""
	if molid_ref is None:
		molid_ref = molid
	if sel_ref is None:
		sel_ref = sel

	ref = atomsel(selection=sel_ref, molid=molid_ref, frame=frame_ref)
	for frame in range(molecule.numframes(molid)):
		sel_t = atomsel(selection=sel, molid=molid, frame=frame)
		M = sel_t.fit(ref)
		atomsel('all', molid=molid, frame=frame).move(M)

def rmsd(sel, molid, sel_ref=None, molid_ref=None, frame_ref=0):
	"""
	Computes the rmsd of sel for molecule molid.

	Args:
		sel (str): VMD atom selection to compute the RMSD for.
		molid (int): Specifies molecule to compute an RMSD for.
		sel_ref (str, optional): VMD atom selection to use for the
			reference structure. Defaults to `sel`.
		molid_ref (str, optional): Specifies molecule to use as reference
			structure. Defaults to `molid`.
		frame_ref (int, optional): Frame of `molid_ref` to use as the
			reference structure. Defaults to 0.
	Returns:
		np.array(float): shape is (num frames,)
	"""
	if molid_ref is None: molid_ref = molid
	if sel_ref is None: sel_ref = sel

	ref = atomsel(sel_ref, molid=molid_ref, frame = frame_ref)
	return np.array([ref.rmsd(atomsel(sel, molid=molid, frame=frame))
	                 for frame in range(molecule.numframes(molid))])

def rmsd_to_average(sel, molid):
	"""
	Computes the rmsd of sel for molecule molid to the average position
	of these atoms.

	The average over time of this is the RMSF.

	** Note that this implementation does not take into account symmetry,
		i.e. this will give too high of a value for a side chain like leucine. **

	Args:
		sel (str): VMD atom selection to compute the RMSD for.
		molid (int): Specifies molecule to compute an RMSD for.
	Returns:
		np.array(float): shape is (num frames,)
	"""
	coords = get_coords(sel, molid)
	# axis = 0 gives time average
	# keepdims = True gives shape of (1, N, 3)
	avg_coords = coords.mean(axis=0, keepdims = True)
	
	deviation = (coords - avg_coords)**2
	deviation = np.sum((coords - avg_coords)**2, axis=-1) # sum over x, y, z
	deviation = np.mean(deviation, axis=1) # mean over atoms
	return np.sqrt(deviation)

def centroid_position(sel, molid):
	"""
	Computes the centroid at each timepoint.

	Args:
		sel (str): VMD atom selection.
		molid (int): molecule ID.
	Returns:
		np.array(float): shape is (num frames, 3)
	"""
	coords = get_coords(sel, molid)
	return np.mean(coords, axis=1)

def centroid_distance(sel1, sel2, molid):
	"""
	Computes the distance between the centroid of sel1 and sel2
	at each timepoint.

	Args:
		sel1, sel2 (str): VMD atom selection.
		molid (int): molecule ID.
	Returns:
		np.array(float): shape is (num frames,)
	"""
	center1 = centroid_position(sel1, molid)
	center2 = centroid_position(sel2, molid)
	return np.linalg.norm(center2 - center1, axis=1)



def centroid_position_z(sel, molid):
    """
    Computes the Z centroid at each timepoint.

    Args:
        sel (str): VMD atom selection.
        molid (int): molecule ID.

    Returns:
        np.array(float): shape is (num frames,)
    """
    coords = get_coords(sel, molid)  # shape: T x N x 3
    z_coords = coords[:, :, 2]       # shape: T x N
    return np.mean(z_coords, axis=1) # shape: T


def centroid_distance_z(sel1, sel2, molid):
    """
    Computes the Z-axis distance between the centroids of sel1 and sel2
    at each timepoint.

    Args:
        sel1, sel2 (str): VMD atom selection.
        molid (int): molecule ID.

    Returns:
        np.array(float): shape is (num frames,)
    """
    z1 = centroid_position_z(sel1, molid)  # shape: T
    z2 = centroid_position_z(sel2, molid)  # shape: T
    return np.abs(z2 - z1)                 # shape: T



def centroid_position_noz(sel, molid):
    """
    Computes the XY centroid at each timepoint (ignores Z).

    Args:
        sel (str): VMD atom selection.
        molid (int): molecule ID.
    Returns:
        np.array(float): shape is (num frames, 2)
    """
    coords = get_coords(sel, molid)
    xy_coords = coords[:, :, :2]  # Keep only x and y
    return np.mean(xy_coords, axis=1)

def centroid_distance_noz(sel1, sel2, molid):
    """
    Computes the 2D distance (ignoring Z) between the centroid of sel1 and sel2
    at each timepoint.

    Args:
        sel1, sel2 (str): VMD atom selection.
        molid (int): molecule ID.
    Returns:
        np.array(float): shape is (num frames,)
    """
    center1 = centroid_position_noz(sel1, molid)
    center2 = centroid_position_noz(sel2, molid)
    return np.linalg.norm(center2 - center1, axis=1)


def min_distance(sel1, sel2, molid):
	"""
	Computes the minimum distance at each timepoint between atoms
	in sel1 and sel2.

	Args:
		sel1, sel2 (str): VMD atom selection.
		molid (int): molecule ID.
	Returns:
		np.array(float): shape is (num frames,)
	"""
	coords1 = get_coords(sel1, molid)
	coords2 = get_coords(sel2, molid)
	return _min_distance(coords1, coords2)


def min_z_distance(sel1, sel2, molid):
        """
        Computes the minimum distance at each timepoint between atoms
        in sel1 and sel2.

        Args:
                sel1, sel2 (str): VMD atom selection.
                molid (int): molecule ID.
        Returns:
                np.array(float): shape is (num frames,)
        """
        coords1 = get_coords(sel1, molid)
        coords2 = get_coords(sel2, molid)
        return _min_z_distance(coords1, coords2)


def angle(sel1, sel2, sel3, molid):
	"""
	Computes the dihedral angle defined by sel1, sel2, sel3, sel4 for
	all available frames of molecule `molid`.

	Args:
		sel1, sel2, sel3, sel4 (str): VMD atomsels containing that evaluate
			to one atom each.
		molid (int): Specifies molecule to compute dihedrals for.

	Returns:
		np.array(float): shape is (num frames,)
	"""
	coords = []
	for sel in [sel1, sel2, sel3]:
		coords += [get_coords(sel, molid)]

		if coords[-1].shape[1] != 1:
			msg = 'sel: "{}"" in angle defined by\n'.format(sel)
			msg += '\tsel1: "{}"\n'.format(sel1)
			msg += '\tsel2: "{}"\n'.format(sel2)
			msg += '\tsel3: "{}"\n'.format(sel3)
			msg += 'has wrong number of atoms: {}'.format(coords[-1].shape[1])
			raise ValueError(msg)

	coords = np.stack(coords, axis=1)
	coords = np.squeeze(coords, axis=2)
	return _angle(coords)

def dihedral(sel1, sel2, sel3, sel4, molid):
	"""
	Computes the dihedral angle defined by sel1, sel2, sel3, sel4 for
	all available frames of molecule `molid`.

	Args:
		sel1, sel2, sel3, sel4 (str): VMD atomsels containing that evaluate
			to one atom each.
		molid (int): Specifies molecule to compute dihedrals for.

	Returns:
		np.array(float): shape is (num frames,)
	"""
	coords = []
	for sel in [sel1, sel2, sel3, sel4]:
		coords += [get_coords(sel, molid)]

		if coords[-1].shape[1] != 1:
			msg = 'sel: "{}"" in dihedral defined by\n'.format(sel)
			msg += '\tsel1: "{}"\n'.format(sel1)
			msg += '\tsel2: "{}"\n'.format(sel2)
			msg += '\tsel3: "{}"\n'.format(sel3)
			msg += '\tsel4: "{}"\n'.format(sel4)
			msg += 'has wrong number of atoms: {}'.format(coords[-1].shape[1])
			raise ValueError(msg)

	coords = np.stack(coords, axis=1)
	coords = np.squeeze(coords, axis = 2)
	return _dihedral(coords)

def vector_between_points(point1, point2):
    """
    Calculates the vector between two points in 2D space.
    
    Args:
        point1 (array-like): Coordinates of the first point [x1, y1].
        point2 (array-like): Coordinates of the second point [x2, y2].
    
    Returns:
        np.array: The vector from point1 to point2 as [vx, vy].
    """
    return np.array(point2) - np.array(point1)

def angle_between_vectors(vector1, vector2):
    """
    Calculates the angle between two vectors in 2D space.
    
    Args:
        vector1, vector2 (array-like): Vectors [vx, vy] between points.
    
    Returns:
        float: Angle in radians between vector1 and vector2.
    """
    dot_product = np.dot(vector1, vector2)
    norm_product = np.linalg.norm(vector1) * np.linalg.norm(vector2)
    return np.arccos(dot_product / norm_product)




def signed_angle_between_vectors_2d(vector1, vector2):
    """
    Calculates the signed angle from vector1 to vector2 in 2D space.
    Positive angle means counterclockwise rotation from vector1 to vector2.
    Negative angle means clockwise rotation from vector1 to vector2.
    
    Args:
        vector1, vector2 (array-like): 2D vectors [x, y]
    
    Returns:
        float: Signed angle in degrees from vector1 to vector2, range [-180, 180]
    """
    # For 2D vectors, the cross product is just a single number
    # representing the z-component
    cross_product = vector1[0] * vector2[1] - vector1[1] * vector2[0]
    dot_product = np.dot(vector1, vector2)
    
    # atan2 returns angle in range [-π, π]
    angle = np.arctan2(cross_product, dot_product)
    return np.degrees(angle)






def scalar_projection(vector1,vector2):
	"""
	Computes the magnitude of vector1 projected onto vector2.

	Args:
		vector1 (np.array float)
		vector2 (np.array float)

	Returns:
		np.array float
	"""
	return np.dot(vector1, vector2) / np.linalg.norm(vector2)



def angle_btwn_vector(sel1, sel2, refsel1, refsel2, molidref, molid):
        """
        Computes the angle between two vectors.

        Args:
                sel1, sel2, sel3, sel4 (str): VMD atomsels containing that evaluate
                        to one atom each.
                molid (int): Specifies molecule to compute the vector for sel1 and sel2.
				molidref: molecule to compute vector for for refsel1 and refsel2

        Returns:
                np.array(float): shape is (num frames,)
        """
        coords_2d = []
        for sel in [sel1, sel2]:
                coords_2d += [get_coords(sel, molid)[:, :, :2]]

                if coords_2d[-1].shape[1] != 1:
                        msg = 'sel: "{}"" in angle defined by\n'.format(sel)
                        msg += '\tsel1: "{}"\n'.format(sel1)
                        msg += '\tsel2: "{}"\n'.format(sel2)
                        msg += 'has wrong number of atoms: {}'.format(coords[-1].shape[1])
                        raise ValueError(msg)
        coords_2d_ref = []
        for sel in [refsel1, refsel2]:
                coords_2d_ref += [get_coords(sel, molidref)[:, :, :2]]

                if coords_2d_ref[-1].shape[1] != 1:
                        msg = 'sel: "{}"" in angle defined by\n'.format(sel)
                        msg += '\tsel1: "{}"\n'.format(sel1)
                        msg += '\tsel2: "{}"\n'.format(sel2)
                        msg += 'has wrong number of atoms: {}'.format(coords[-1].shape[1])
                        raise ValueError(msg)

        coords_2d = np.stack(coords_2d, axis=1)
        coords_2d = np.squeeze(coords_2d, axis=2)

        coords_2d_ref = np.stack(coords_2d_ref, axis=1)
        coords_2d_ref = np.squeeze(coords_2d_ref, axis=2)
        ind = 133
        print('index ', ind)
        print('coords ', coords_2d[ind], 'coords ref ', coords_2d_ref[0])
            # Calculate vectors between points for each frame
        vectors_2d = coords_2d[:, 1] - coords_2d[:, 0]           # Vector sel2 - sel1 for each frame
        vectors_2d_ref = coords_2d_ref[:, 1] - coords_2d_ref[:, 0]  # Vector refsel2 - refsel1 for each frame
        if vectors_2d_ref.shape[0] ==1:
                vectors_2d_ref = np.tile(vectors_2d_ref, (vectors_2d.shape[0], 1)) 
        print('vector ',vectors_2d[ind])
        print('ref vector ',vectors_2d_ref[ind])
        # Calculate angles between the vectors for each frame
        angle_diffs = []
        for i in range(vectors_2d.shape[0]):
                angle = signed_angle_between_vectors_2d(vectors_2d[i], vectors_2d_ref[i])
                angle_diffs.append(angle)  # No need to negate
        return np.array(angle_diffs)





# Helper functions extracted from above to facilitate testing.
def _min_distance(coords1, coords2):
	"""
	Computes the minimum distance between coodinates in `coords1`
	and `coords2` for each timepoint.

	Args: coords1, coords2 (np.array float): Arrays should have dimensions
		N timepoints (T) x N atoms (N1, N2) x 3.

	Returns:
		np.array(float): 1-D array of length T.
	"""
	# coords are T x N1 x 3 and T x N2 x 3
	# expand_dims gives T x 1 x N1 x 3 and T x N2 x 1 x 3
	coords1 = np.expand_dims(coords1, 1)
	coords2 = np.expand_dims(coords2, 2)
	
	# Since above made coords1 and coords2 orthoganol, subtracting them
	# results in all pairwise combionsinations being computed.
	# subtraction gives T x N2 x N1 x 3
	displacements = coords1 - coords2

	# norm gives T x N2 x N1
	distances = np.linalg.norm(displacements, axis=3)
	
	# min gives T
	return distances.min(axis=(1, 2))



def _min_z_distance(coords1, coords2):
    """
    Computes the minimum absolute distance between z-coordinates 
    of atoms in `coords1` and `coords2` for each timepoint.

    Args:
        coords1, coords2 (np.array float): Arrays with dimensions
            T timepoints x N1 atoms x 3
            T timepoints x N2 atoms x 3

    Returns:
        np.array(float): 1-D array of minimum z-distances per timepoint (length T).
    """
    # Extract z-coordinates (last dimension index 2)
    z1 = coords1[:, :, 2]  # shape: T x N1
    z2 = coords2[:, :, 2]  # shape: T x N2

    # Expand dims to compute pairwise z-differences: T x N1 x 1 and T x 1 x N2
    z1 = np.expand_dims(z1, 2)  # T x N1 x 1
    z2 = np.expand_dims(z2, 1)  # T x 1 x N2

    # Compute absolute z-coordinate differences: T x N1 x N2
    z_diff = np.abs(z1 - z2)

    # Take minimum over all atom pairs: T
    return z_diff.min(axis=(1, 2))


def _angle(coords):
	v1 = coords[:, 0] - coords[:, 1]
	v2 = coords[:, 2] - coords[:, 1]

	x = (v1*v2).sum(axis=1)
	x /= np.linalg.norm(v1, axis=1)
	x /= np.linalg.norm(v2, axis=1)
	return np.arccos(x)*180.0/np.pi

def _dihedral(coords):
	"""
	Computes dihedral angles for np.array of coordinates for 4 points.

	Formula was adapted from stackoverflow:
    https://stackoverflow.com/questions/20305272/
    dihedral-torsion-angle-from-four-points-in-cartesian-coordinates-in-python

	** Is there an efficient built-in method to get dihedrals? **

	Args:
		coords (np.array of floats of shape T x 4 x 3)

	Returns:
		length T np.array: Dihedral angles in degrees for all timepoints.
	"""
	b0 = -(coords[:, 1] - coords[:, 0])
	b1 = coords[:, 2] - coords[:, 1]
	b2 = coords[:, 3] - coords[:, 2]
	b1 /= np.linalg.norm(b1, axis = 1, keepdims = True)
	# v = projection of b0 onto plane perpendicular to b1
	#   = b0 minus component that aligns with b1
	# w = projection of b2 onto plane perpendicular to b1
	#   = b2 minus component that aligns with b1
	v = b0 - np.sum(b0*b1, axis = 1, keepdims = True)*b1
	w = b2 - np.sum(b2*b1, axis = 1, keepdims = True)*b1
	# angle between v and w in a plane is the torsion angle
	# v and w may not be normalized but that's fine since tan is y/x
	x = np.sum(v*w, axis = 1)
	y = np.sum(np.cross(b1, v, axis = 1)*w, axis = 1)
	return np.degrees(np.arctan2(y, x))
