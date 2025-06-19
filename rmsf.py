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

root = '/oak/stanford/groups/rondror/projects/MD_simulations/amber/KOR/KOR_SK_PCP/v1/SK_tao1'

SK_top = '/oak/stanford/groups/rondror/projects/MD_simulations/amber/KOR/KOR_SK_PCP/v1/SK_tao1/prep/dabble/system_dabbled.psf'
PCP_top = '/oak/stanford/groups/rondror/projects/MD_simulations/amber/KOR/KOR_SK_PCP/v1/PCP_tao/prep/dabble/system_dabbled.psf'

SK_ncs = []
PCP_ncs = []


for i in range(1, 7):
        SK_ncs.append(f'/oak/stanford/groups/rondror/projects/MD_simulations/amber/KOR/KOR_SK_PCP/v1/SK_tao1/run_{i}/summary_traj_w_eq_stride5.nc')
        PCP_ncs.append(f'/oak/stanford/groups/rondror/projects/MD_simulations/amber/KOR/KOR_SK_PCP/v1/PCP_tao/run_{i}/summary_traj_w_eq_stride5.nc')

align_SK = 'name CA and (protein and ((resid 58 to 69) or (resid 186 to 193) or (resid 286 to 294) or (resid 107 to 116) or (resid 135 to 144) or (resid 313 to 319) or (resid 227 to 231)))'
align_PCP = 'name CA and (protein and ((resid 58 to 72) or (resid 108 to 119) or (resid 135 to 142) or (resid 186 to 193) or (resid 227 to 235) or (resid 286 to 293) or (resid 312 to 320)))'

rmsf_sel_SK = 'resname SKET and not type H'
rmsf_sel_PCP = 'resname PCP and not type H'

rmsf_SK = get_rmsf(SK_top, SK_ncs, align_SK, rmsf_sel_SK)
rmsf_PCP = get_rmsf(PCP_top, PCP_ncs, align_PCP, rmsf_sel_PCP)

print(f"""
RMSF SK: {round(rmsf_SK,3)}
RMSF PCP: {round(rmsf_PCP,3)}
""")
