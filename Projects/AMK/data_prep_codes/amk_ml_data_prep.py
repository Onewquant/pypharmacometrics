from tools import *
from pynca.tools import *
from datetime import datetime, timedelta
from scipy.stats import spearmanr

# result_type = 'Phoenix'
# result_type = 'R'

prj_name = 'AMK'
prj_dir = f'./Projects/{prj_name}'
resource_dir = f'{prj_dir}/resource'
output_dir = f"{prj_dir}/results"
nonmem_dir = f'C:/Users/ilma0/NONMEMProjects/{prj_name}'
if not os.path.exists(output_dir):
    os.mkdir(output_dir)

aki_df = pd.read_csv(f"{output_dir}/amk_aki.csv")
# aki_df.columns
# modeling_df.columns
modeling_df = pd.read_csv(f"{nonmem_dir}/amk_modeling_df_covar.csv")
# simres_df = pd.read_csv(f"{nonmem_dir}/amk_modeling_df_covar.csv")
# modeling_df['UID']

totlab_df = pd.read_csv(f"{output_dir}/totlab_df.csv")
comed_df = pd.read_csv(f"{output_dir}/final_comed_df.csv")
cm_df = pd.read_csv(f"{output_dir}/final_comorbidity_df.csv")
vs_df = pd.read_csv(f"{output_dir}/final_vs_data.csv")
proc_df = pd.read_csv(f"{output_dir}/final_procedure_df.csv")
micbio_df = pd.read_csv(f"{output_dir}/final_microbiology_df.csv")
losmot_df = pd.read_csv(f"{output_dir}/final_locmotality_df.csv")