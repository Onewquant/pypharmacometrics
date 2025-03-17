from tools import *


project_dict = {'KSCPTSPRWS25MULTI':['SGLT2INH'], 'KSCPTSPRWS25SINGLE':['SGLT2INH']}
modeling_dir_path = "C:/Users/ilma0/PycharmProjects/pypharmacometrics/resource/KSCPTSPRWS25"
prepconc_dir_path = "C:/Users/ilma0/PycharmProjects/pypharmacometrics/resource/KSCPTSPRWS25/prep_data"
resource_dir_path = "C:/Users/ilma0/PycharmProjects/pypharmacometrics/resource/KSCPTSPRWS25/resources"


bpd_df = pd.read_excel(f"{resource_dir_path}/PD prep_Final_0707_BASE group.xlsx")
upd_df = pd.read_excel(f"{resource_dir_path}/PD urine_prep_Final_0707_prep_BASE group.xlsx")
