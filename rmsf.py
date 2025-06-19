import MDAnalysis as mda
from MDAnalysis.analysis.rms import RMSF
from MDAnalysis.analysis.align import alignto
import numpy as np

def get_rmsf(top, ncs, align_sel, rmsf_sel):

        u = mda.Universe(top, ncs)

        # Select protein backbone or alpha-carbons
        selection = u.select_atoms(rmsf_sel)

        # Align to first frame or an external reference
        # (optional but improves RMSF interpretation)
        alignto(u, u, select=align_sel)  # align all to first frame

        # Run RMSF
        rmsf = RMSF(selection).run()
        return np.mean(rmsf.results.rmsf)

lig_dirs = ['SK_tao1', 'PCP_tao']

reps = [6, 6]

align_sels = ['name CA and (protein and ((resid 58 to 69) or (resid 186 to 193) or (resid 286 to 294) or (resid 107 to 116) or (resid 135 to 144) or (resid 313 to 319) or (resid 227 to 231)))',
'name CA and (protein and ((resid 58 to 72) or (resid 108 to 119) or (resid 135 to 142) or (resid 186 to 193) or (resid 227 to 235) or (resid 286 to 293) or (resid 312 to 320)))']

rmsf_sels = ['resname SKET and not type H',
'resname PCP and not type H']

root = '/oak/stanford/groups/rondror/projects/MD_simulations/amber/KOR/KOR_SK_PCP/v1/'

files = {}

for j, lig_dir in enumerate(lig_dirs):
        files[lig_dir] = {}
        files[lig_dir]['top'] = f'{root}/{lig_dir}/prep/dabble/system_dabbled.psf'
        files[lig_dir]['trajs'] = [f'{root}/{lig_dir}/run_{x}/summary_traj_w_eq_stride5.nc' for x in range(1,reps[j]+1)]

rmsfs = {}

for i, lig_dir in enumerate(lig_dirs):
        rmsfs[lig_dir] = get_rmsf(files[lig_dir]['top'], files[lig_dir]['trajs'], align_sels[i], rmsf_sels[i])
        print(f'{lig_dir} RMSF: {round(rmsfs[lig_dir],2)}')
