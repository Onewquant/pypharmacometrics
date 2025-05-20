from emr_scraper.tools import *

# df = pd.read_csv("C:/Users/snubh/PycharmProjects/snubhprac/emr_scraper/anti_TNFa_patients.csv")

## Lab 수집

existing_files = glob.glob("C:/Users/snubh/PycharmProjects/snubhprac/emr_scraper/lab/IBD_PGx_lab(*).xlsx")
# existing_id_list = [int(f.split('IBD_PGx_lab(')[-1].split('_')[0]) for f in existing_files]
for finx, f in enumerate(existing_files):
    df = pd.read_excel(f)
    print(f"({finx}) {f.split('IBD_PGx_lab(')[-1].split(').xlsx')[0]} / {list(df.iloc[-1])} / {list(df.iloc[0])}")
    # raise ValueError
