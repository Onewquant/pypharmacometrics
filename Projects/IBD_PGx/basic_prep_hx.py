from tools import *
from pynca.tools import *

# result_type = 'Phoenix'
# result_type = 'R'

prj_name = 'IBD_PGX'
prj_dir = 'C:/Users/ilma0/PycharmProjects/pypharmacometrics/Projects/IBD_PGx'
resource_dir = f'{prj_dir}/resource'
output_dir = f"{prj_dir}/results"
if not os.path.exists(output_dir):
    os.mkdir(output_dir)

## Orders
drug_order_set = set()

hx_files = glob.glob(f'{resource_dir}/hx/IBD_PGx_hx(*).txt')
# dose_result_df = list()
# wierd_result_df = list()
# result_cols = ['ID','NAME','DATETIME','DRUG','DOSE','ROUTE','ACTING','PERIOD','ETC_INFO','PLACE']
# no_dup_cols = [c for c in result_cols if c!='NAME']

isornot_dict = {'유':1.0, '무':0.0}

pms_df = list()
ms_df = list()
cdai_df = list()
cdsurvey_df = list()
for finx, fpath in enumerate(hx_files): #break

    pid = fpath.split('(')[-1].split('_')[0]
    pname = fpath.split('_')[-1].split(')')[0]

    print(f"({finx}) {pname} / {pid}")
    # if pname=='이학준':
    #     raise ValueError

    # if pid in ("15322168", "19739357", "34835292", "37366865", "21618097", "36898756", "36975211", "37858047"):       # lab, order 파일 다시 수집 필요
    #     continue

    with open(fpath, "r", encoding="utf-8") as f:
        f_content = f.read()
        # print(content)
    for fct_inx, fct_str in enumerate([fct.split('\n작성자')[0] for fct in f_content.split('작성과: ')][1:]): #break

        if 'Partial Mayo Score for Ulcerative Colitis' in fct_str:
            print(f'Partial Mayo Score for Ulcerative Colitis')
            pms_dict = {'ID':pid, 'NAME':pname}
            pms_dict['DATE'] = re.findall(r'[\d][\d][\d][\d]-[\d][\d]-[\d][\d]',fct_str)[0]
            try: pms_dict['PMS_TOTALSCORE'] = float(re.findall(r'Total\s*:\s*\d+',fct_str)[0].split(':')[-1].strip())
            except: pms_dict['PMS_TOTALSCORE'] = np.nan
            try: pms_dict['PMS_STLCNT'] = float(re.findall(r'대변 횟수\s*\d+',fct_str)[0].split('대변 횟수')[-1].strip())
            except: pms_dict['PMS_STLCNT'] = np.nan
            try: pms_dict['PMS_NORMSTLCNT'] = float(re.findall(r'\(정상 횟수:\s*\d+\s*회',fct_str)[0].split(':')[-1].split('회')[0].strip())
            except: pms_dict['PMS_NORMSTLCNT'] = np.nan
            try: pms_dict['PMS_HEMATOCHEZIA'] = float(re.findall(r'혈변\s*\d+\s*',fct_str)[0].split('혈변')[-1].strip())
            except: pms_dict['PMS_HEMATOCHEZIA'] = np.nan
            try: pms_dict['PMS_TOTALEVAL'] = float(re.findall(r'전체적인 평가\s*\d+\s*',fct_str)[0].split('전체적인 평가')[-1].strip())
            except: pms_dict['PMS_TOTALEVA'] = np.nan
            try: pms_dict['PMS_URGENCY'] = isornot_dict[re.findall(r'급박감\s*[유|무]\s*',fct_str)[0].split('급박감')[-1].strip()]
            except: pms_dict['PMS_URGENCY'] = np.nan
            try: pms_dict['PMS_TENESMUS'] = isornot_dict[re.findall(r'잔변감\s*[유|무]\s*',fct_str)[0].split('잔변감')[-1].strip()]
            except: pms_dict['PMS_TENESMU'] = np.nan
            try: pms_dict['PMS_MUCOUSSTL'] = isornot_dict[re.findall(r'점액변\s*[유|무]\s*', fct_str)[0].split('점액변')[-1].strip()]
            except: pms_dict['PMS_MUCOUSSTL'] = np.nan
            try: pms_dict['PMS_NOCTSTL'] = isornot_dict[re.findall(r'야간배변\s*[유|무]\s*',fct_str)[0].split('야간배변')[-1].strip()]
            except: pms_dict['PMS_NOCTSTL'] = np.nan
            try:
                adherencepct = re.findall(r'약제 순응도\s*[처방받은 약제 없음|\d]+\s*',fct_str.replace('(%)',''))[0].split('약제 순응도')[-1].strip()
                if adherencepct in ['처방받은 약제 없음','']:
                    adherencepct = np.nan
                else:
                    try: adherencepct = float(adherencepct)
                    except:
                        print("ADHERENCE // ", adherencepct)
                        raise ValueError
            except:
                adherencepct = np.nan

            pms_dict['PMS_ADHERENCEPCT'] = adherencepct

            # try:
            suppository_usepct = re.findall(r'관장액 사용도\(%\)\s*[처방받은 약제 없음|\d]*\n', fct_str)
            if len(suppository_usepct)==0:
                suppository_usepct = np.nan
            else:
                # raise ValueError
                suppository_usepct = suppository_usepct[0].split('\n')[0].split(')')[-1].strip()
                if suppository_usepct in ['처방받은 약제 없음','']:
                    suppository_usepct = np.nan
                else:
                    try: suppository_usepct = float(suppository_usepct)
                    except:
                        print(suppository_usepct)
                        raise ValueError
            # except:
            #     suppository_usepct = np.nan

            pms_dict['PMS_SUPPOSITORY_USEPCT'] = suppository_usepct
            try: pms_dict['PMS_HEIGHT'] = float(re.findall(r'[신장|키]\s*:?\s*\d+[\.]?\d*\s*',fct_str.replace('(cm)',''))[0].split(':')[-1].split('키')[-1].split('cm')[0].strip())
            except: pms_dict['PMS_HEIGHT'] = np.nan
            try: pms_dict['PMS_WEIGHT'] = float(re.findall(r'체중\s*:?\s*\d+[\.]?\d*\s*',fct_str.replace('(kg)',''))[0].split(':')[-1].split('체중')[-1].split('kg')[0].strip())
            except: pms_dict['PMS_WEIGHT'] = np.nan
            try: pms_dict['PMS_BMI'] = float(re.findall(r'BMI\s*:?\s*\d+[\.]?\d*\s*',fct_str.replace('(kg/㎡)',''))[0].split(':')[-1].split('BMI')[-1].split('kg')[0].strip())
            except: pms_dict['PMS_BMI'] = np.nan

            pms_df.append(pms_dict)
        elif 'Mayo Score' in fct_str:
            print('Mayo Score for UC')
            continue
            raise ValueError
        elif "Crohn's Disease Activity index" in fct_str:
            print("Crohn's Disease Activity index")
            continue
            raise ValueError
        elif "크론병 질병 활성도 평가를 위한 설문" in fct_str:
            print('크론병 질병 활성도')
            continue
            raise ValueError
        else:
            continue
            raise ValueError

        # break
pms_df = pd.DataFrame(pms_df)
pms_df.to_csv(f"{output_dir}/pdmarker_pms_df.csv", encoding='utf-8-sig', index=False)
print('COMPLETED')