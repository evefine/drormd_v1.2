generic = {

    '6x31': 'protein and resid 319',
    '4x42': 'protein and resid 192',
    '6x48': 'protein and resid 336',
    '6x33': 'protein and resid 321',
    '6x35': 'protein and resid 323',
    '6x37': 'protein and resid 325',
    '5x58': 'protein and resid 254',

    '7x54': 'protein and resid 381',
    '2x43': 'protein and resid 113',
    '7x50': 'protein and resid 377',
    '2x47': 'protein and resid 117',
    '5x54': 'protein and resid 189',

    '1x50': 'protein and resid 92',
    '7x46': 'protein and resid 373',
    '7x47': 'protein and resid 374',

    '8x47': 'protein and resid 384',
    '2x39': 'protein and resid 109',

    '5x57': 'protein and resid 253',
    '4x45': 'protein and resid 195',

    '3x50': 'protein and resid 173',
    '7x53': 'protein and resid 380',

    '6x50': 'protein and resid 338',
    '6x43': 'protein and resid 331',

    '7x54': 'protein and resid 381',
    '7x56': 'protein and resid 383'
}

root = '/oak/stanford/groups/rondror/projects/MD_simulations/amber/5HT2AR'

conditions = {
    'apo': {
        'topology': '{}/apo/v2/prep/dabble/system_dabbled.psf'.format(root),
        'trajectory': '{}/apo/v2/run_{{}}/summary_traj_w_eq_stride5.nc'.format(root),
        'reps': range(1, 13),
    },
    'BOL': {
        'topology': '{}/BOL/v3/prep/dabble/system_dabbled.psf'.format(root),
        'trajectory': '{}/BOL/v3/run_{{}}/summary_traj_w_eq_stride5.nc'.format(root),
        'reps': range(1, 7),
    },
    'DMT': {
        'topology': '{}/DMT/v1/prep/dabble/system_dabbled.psf'.format(root),
        'trajectory': '{}/DMT/v1/run_{{}}/summary_traj_w_eq_stride5.nc'.format(root),
        'reps': range(1, 13),
    },
    'DOPR': {
        'topology': '{}/DOPR/v1/prep/dabble/system_dabbled.psf'.format(root),
        'trajectory': '{}/DOPR/v1/run_{{}}/summary_traj_w_eq_stride5.nc'.format(root),
        'reps': range(1, 7),
    },
    'lisuride': {
        'topology': '{}/lisuride/v3/prep/dabble/system_dabbled.psf'.format(root),
        'trajectory': '{}/lisuride/v3/run_{{}}/summary_traj_w_eq_stride5.nc'.format(root),
        'reps': range(1, 13),
    },
    'LSD': {
        'topology': '{}/LSD/v3/prep/dabble/system_dabbled.psf'.format(root),
        'trajectory': '{}/LSD/v3/run_{{}}/summary_traj_w_eq_stride5.nc'.format(root),
        'reps': range(1, 13),
    },
    'mescaline': {
        'topology': '{}/mescaline/v1/prep/dabble/system_dabbled.psf'.format(root),
        'trajectory': '{}/mescaline/v1/run_{{}}/summary_traj_w_eq_stride5.nc'.format(root),
        'reps': range(1, 13),
    },
    'psilocin': {
        'topology': '{}/psilocin/v4/prep/dabble/system_dabbled.psf'.format(root),
        'trajectory': '{}/psilocin/v4/run_{{}}/summary_traj_w_eq_stride5.nc'.format(root),
        'reps': range(1, 13),
    },
    'serotonin': {
        'topology': '{}/serotonin/v3/prep/dabble/system_dabbled.psf'.format(root),
        'trajectory': '{}/serotonin/v3/run_{{}}/summary_traj_w_eq_stride5.nc'.format(root),
        'reps': range(1, 13),
    },
}
