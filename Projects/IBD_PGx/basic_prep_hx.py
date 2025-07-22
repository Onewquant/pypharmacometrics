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
mss_df = list()
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
        elif 'Mayo Scoring system' in fct_str:
            print('Mayo Score for UC')
            mss_dict = {'ID': pid, 'NAME': pname}
            mss_dict['DATE'] = re.findall(r'[\d][\d][\d][\d]-[\d][\d]-[\d][\d]', fct_str)[0]

            mss_total = fct_str.split('합계')[-1].strip()
            try:
                mss_dict['MSS_TOTALSCORE'] = float(mss_total)
            except:
                mss_dict['MSS_TOTALSCORE'] = np.nan
                print('MSS_TOTALSCORE ERROR')
                # raise ValueError

            mss_stlcnt = fct_str.split('5회 이상 증가')[-1].split('혈변')[0].strip()
            try: mss_dict['MSS_STLCNT'] = float(mss_stlcnt)
            except:
                mss_dict['MSS_STLCNT'] = np.nan
                print('MSS_STLCNT ERROR')
                # raise ValueError

            mss_hematochezia = fct_str.split('피만 나옴')[-1].split('내시경 관찰')[0].strip()
            try: mss_dict['MSS_HEMATOCHEZIA'] = float(mss_hematochezia)
            except:
                mss_dict['MSS_HEMATOCHEZIA'] = np.nan
                print('MSS_HEMATOCHEZIA ERROR')
                # raise ValueError

            mss_endoscopy = fct_str.split('중증(지속적 출혈, 궤양)')[-1].split('전체적인 평가')[0].strip()
            try: mss_dict['MSS_ENDOSCOPY'] = float(mss_endoscopy)
            except:
                mss_dict['MSS_ENDOSCOPY'] = np.nan
                print('MSS_ENDOSCOPY ERROR')
                # raise ValueError

            mss_ovreval = fct_str.split('전체적인 평가')[-1].split('중증')[-1].split('(마지막칸을')[0].strip()
            try:
                mss_dict['MSS_OVREVAL'] = float(mss_ovreval)
            except:
                mss_dict['MSS_OVREVAL'] = np.nan
                print('MSS_OVREVAL ERROR')
                # raise ValueError
            # fct_str['MSS_TOTALSCORE']


            # continue
            mss_df.append(mss_dict)
        elif "Crohn's Disease Activity Index" in fct_str:
            print("Crohn's Disease Activity Index")
            cdai_dict = {'ID': pid, 'NAME': pname}
            cdai_dict['DATE'] = re.findall(r'[\d][\d][\d][\d]-[\d][\d]-[\d][\d]', fct_str)[0]
            # fct_str_ori = fct_str
            fct_str = fct_str.split("Crohn's Disease Activity Index\n\n")[-1]

            try:cdai_dict['CDAI_TOTALSCORE'] = float(re.findall(r'Total\s*:\s*\d+', fct_str)[0].split(':')[-1].strip())
            except:cdai_dict['CDAI_TOTALSCORE'] = np.nan

            try:cdai_dict['CDAI_DIARHEACNT'] = float(re.findall(r'설사 :   하루평균\s*\d+', fct_str)[0].split('설사 :   하루평균')[-1].strip())
            except:cdai_dict['CDAI_DIARHEACNT'] = np.nan
            try:cdai_dict['CDAI_DIARHEASCORE'] = float(re.findall(r'x 7일 x 2 =\s*\d+', fct_str)[0].split('x 7일 x 2 =')[-1].strip())
            except:cdai_dict['CDAI_DIARHEASCORE'] = np.nan

            try:cdai_dict['CDAI_ABDPAINCNT'] = float(re.findall(r'복통 :\s*\d+', fct_str)[0].split('복통 :')[-1].strip())
            except:cdai_dict['CDAI_ABDPAINCNT'] = np.nan
            try:cdai_dict['CDAI_ABDPAINSCORE'] = float(re.findall(r'x 7일 x 5 =\s*\d+', fct_str)[0].split('x 7일 x 5 =')[-1].strip())
            except:cdai_dict['CDAI_ABDPAINSCORE'] = np.nan

            try:cdai_dict['CDAI_GENWELLBEING'] = float(re.findall(r'전신 안녕감 :\s*\d+', fct_str)[0].split('전신 안녕감 :')[-1].strip())
            except:cdai_dict['CDAI_GENWELLBEING'] = np.nan
            try:cdai_dict['CDAI_GENWBSCORE'] = float(re.findall(r'x 7일 x 7 =\s*\d+', fct_str)[0].split('x 7일 x 7 =')[-1].strip())
            except:cdai_dict['CDAI_GENWBSCORE'] = np.nan

            try:cdai_dict['CDAI_ARTHSX'] = float(re.findall(r'관절염\/ 관절통 :\s*\d+', fct_str)[0].split(':')[-1].strip())
            except:cdai_dict['CDAI_ARTHSX'] = np.nan
            try:cdai_dict['CDAI_ARTHSXSCORE'] = float(fct_str.split('B. 홍채염/ 포도막염')[0].split('x 20 =')[-1].strip())
            except:cdai_dict['CDAI_ARTHSXSCORE'] = np.nan

            try:cdai_dict['CDAI_EYESX'] = float(re.findall(r'홍채염\/ 포도막염 :\s*\d+', fct_str)[0].split(':')[-1].strip())
            except:cdai_dict['CDAI_EYESX'] = np.nan
            try:cdai_dict['CDAI_EYESXSCORE'] = float(fct_str.split('C. 결절홍반/ 괴저농피증')[0].split('x 20 =')[-1].strip())
            except:cdai_dict['CDAI_EYESXSCORE'] = np.nan

            try:cdai_dict['CDAI_SKINSX'] = float(re.findall(r'결절홍반\/ 괴저농피증\/ 아프타구내염 :\s*\d+', fct_str)[0].split(':')[-1].strip())
            except:cdai_dict['CDAI_SKINSX'] = np.nan
            try:cdai_dict['CDAI_SKINSXSCORE'] = float(fct_str.split('D. 항문열창, 치루 또는 농양')[0].split('x 20 =')[-1].strip())
            except:cdai_dict['CDAI_SKINSXSCORE'] = np.nan

            try:cdai_dict['CDAI_ANALSX'] = float(re.findall(r'항문열창, 치루 또는 농양 :\s*\d+', fct_str)[0].split(':')[-1].strip())
            except:cdai_dict['CDAI_ANALSX'] = np.nan
            try:cdai_dict['CDAI_ANALSXSCORE'] = float(fct_str.split('E. 기타 누공')[0].split('x 20 =')[-1].strip())
            except:cdai_dict['CDAI_ANALSXSCORE'] = np.nan

            try:cdai_dict['CDAI_OTHERFIST'] = float(re.findall(r'기타 누공 :\s*\d+', fct_str)[0].split(':')[-1].strip())
            except:cdai_dict['CDAI_OTHERFIST'] = np.nan
            try:cdai_dict['CDAI_OTHERFISTSCORE'] = float(fct_str.split('F. 최근 7일동안 37.8')[0].split('x 20 =')[-1].strip())
            except:cdai_dict['CDAI_OTHERFISTSCORE'] = np.nan

            try:cdai_dict['CDAI_FEVER'] = float(re.findall(r'38\.3도 \(항문\) 이상의 열 :\s*\d+', fct_str)[0].split(':')[-1].strip())
            except:cdai_dict['CDAI_FEVER'] = np.nan
            try:cdai_dict['CDAI_FEVERSCORE'] = float(fct_str.split('최근 7일동안 지사제 치료')[0].split('x 20 =')[-1].strip())
            except:cdai_dict['CDAI_FEVERSCORE'] = np.nan

            try:cdai_dict['CDAI_ANTIDIARHMEDI'] = float(re.findall(r'지사제 치료를 받은적이 있는 경우 :\s*\d+', fct_str)[0].split(':')[-1].strip())
            except:cdai_dict['CDAI_ANTIDIARHMEDI'] = np.nan
            try:cdai_dict['CDAI_ANTIDIARHMEDISCORE'] = float(fct_str.split('복부 종괴 :')[0].split('x 30 =')[-1].strip())
            except:cdai_dict['CDAI_ANTIDIARHMEDISCORE'] = np.nan

            try:cdai_dict['CDAI_ABDMASS'] = float(re.findall(r'복부 종괴 :\s*\d+', fct_str)[0].split(':')[-1].strip())
            except:cdai_dict['CDAI_ABDMASS'] = np.nan
            try:cdai_dict['CDAI_ABDMASSSCORE'] = float(fct_str.split('헤마토크릿  남성   :')[0].split('x 10 =')[-1].strip())
            except:cdai_dict['CDAI_ABDMASSSCORE'] = np.nan


            if ('헤마토크릿\n\n남성' in fct_str) or ('헤마토크릿  남성' in fct_str):
                try:cdai_dict['CDAI_HEMATOCRIT'] = float(re.findall(r'헤마토크릿[\s|\n][\s|\n]남성   \: \( 47 \-\s*\d+', fct_str)[0].split('-')[-1].strip())
                except:cdai_dict['CDAI_HEMATOCRIT'] = np.nan
                try:cdai_dict['CDAI_HCTSCORE'] = float(fct_str.split('x 6 =')[-1].split('조정값')[0].strip())
                except:cdai_dict['CDAI_HCTSCORE'] = np.nan
                try:cdai_dict['CDAI_ADJHCTSCORE'] = float(fct_str.split('\n\n신장 :')[0].split('조정값')[-1].strip())
                except:cdai_dict['CDAI_ADJHCTSCORE'] = np.nan
            elif ('헤마토크릿\n\n여성' in fct_str) or ('헤마토크릿  여성' in fct_str):
                try:cdai_dict['CDAI_HEMATOCRIT'] = float(re.findall(r'헤마토크릿[\s|\n][\s|\n]여성   \: \( 47 \-\s*\d+', fct_str)[0].split('-')[-1].strip())
                except:cdai_dict['CDAI_HEMATOCRIT'] = np.nan
                try:cdai_dict['CDAI_HCTSCORE'] = float(fct_str.split('x 6 =')[-1].split('조정값')[0].strip())
                except:cdai_dict['CDAI_HCTSCORE'] = np.nan
                try:cdai_dict['CDAI_ADJHCTSCORE'] = float(fct_str.split('\n\n신장 :')[0].split('조정값')[-1].strip())
                except:cdai_dict['CDAI_ADJHCTSCORE'] = np.nan
            else:
                cdai_dict['CDAI_HEMATOCRIT'] = np.nan
                cdai_dict['CDAI_HCTSCORE'] = np.nan
                cdai_dict['CDAI_ADJHCTSCORE'] = np.nan


                # try:cdai_dict['CDAI_HEMATOCRIT'] = float(re.findall(r'헤마토크릿  남성   \: \( 47 \-\s*\d+', fct_str)[0].split('-')[-1].strip())
                # except:cdai_dict['CDAI_HEMATOCRIT'] = np.nan


            try:cdai_dict['CDAI_HCTSCORE'] = float(fct_str.split('x 6 =')[-1].split('조정값')[0].strip())
            except:cdai_dict['CDAI_HCTSCORE'] = np.nan
            try:cdai_dict['CDAI_ADJHCTSCORE'] = float(fct_str.split('\n\n신장 :')[0].split('조정값')[-1].strip())
            except:cdai_dict['CDAI_ADJHCTSCORE'] = np.nan

            try:cdai_dict['CDAI_STDWEIGHT'] = float(re.findall(r'표준체중\s*:?\s*\d+[\.]?\d*\s*', fct_str)[0].split('표준체중')[-1].strip())
            except:cdai_dict['CDAI_STDWEIGHT'] = np.nan
            try:cdai_dict['CDAI_WEIGHT'] = float(re.findall(r'체중\s*:?\s*\d+[\.]?\d*\s*', fct_str.replace('(kg)', ''))[0].split(':')[-1].split('체중')[-1].split('kg')[0].strip())
            except:cdai_dict['CDAI_WEIGHT'] = np.nan
            try:cdai_dict['CDAI_WTSCORE'] = float(fct_str.split('x 100 =')[-1].split('조정값')[0].strip())
            except:cdai_dict['CDAI_WTSCORE'] = np.nan
            try:cdai_dict['CDAI_ADJWTSCORE'] = float(fct_str.split("Crohn’s Disease Obstructive Score")[0].split("조정값")[-1].strip())
            except:cdai_dict['CDAI_ADJWTSCORE'] = np.nan

            try:cdai_dict['CDAI_HEIGHT'] = float(re.findall(r'[신장|키]\s*:?\s*\d+[\.]?\d*\s*', fct_str.replace('(cm)', ''))[0].split(':')[-1].split('키')[-1].split('cm')[0].strip())
            except:cdai_dict['CDAI_HEIGHT'] = np.nan
            try:cdai_dict['CDAI_BMI'] = float(re.findall(r'BMI\s*:?\s*\d+[\.]?\d*\s*', fct_str.replace('(kg/㎡)', ''))[0].split(':')[-1].split('BMI')[-1].split('kg')[0].strip())
            except:cdai_dict['CDAI_BMI'] = np.nan
            try:cdai_dict['CDAI_STDWEIGHT'] = float(re.findall(r'실제체중\s*:?\s*\d+[\.]?\d*\s*', fct_str)[0].split('표준체중')[-1].strip())
            except:cdai_dict['CDAI_STDWEIGHT'] = np.nan

            # raise ValueError
            fct_cdos_split = fct_str.split('Crohn’s Disease Obstructive Score')
            cdos_str = fct_cdos_split[-1].split('약제 순응도')[0].strip()
            # try:cdai_dict['CDAI_CDOSTOTAL'] =
            # except:cdai_dict['CDAI_CDOSTOTAL'] = np.nan
            if cdos_str.split('통증 강도 :')[-1].split('연관 증상 :')[0].split('통증 기간 :')[0].strip()=='없음':
                cdai_dict['CDAI_CDOS'] = 0
                cdai_dict['CDAI_CDOSPAININT'] = np.nan
                cdai_dict['CDAI_CDOSPAINPERIOD'] = np.nan
                cdai_dict['CDAI_CDOSASSOSX'] = np.nan
                cdai_dict['CDAI_CDOSASSOSXPERIOD'] = np.nan
                cdai_dict['CDAI_CDOSEMERGENCY'] = np.nan

            elif len(fct_cdos_split)>1:
                cdai_dict['CDAI_CDOS'] = 1

                try:cdai_dict['CDAI_CDOSPAININT'] = re.findall(r'통증\s*강도\s*:\s*(\S+)', cdos_str)[0].strip()
                except:cdai_dict['CDAI_CDOSPAININT'] = np.nan

                try:cdai_dict['CDAI_CDOSPAINPERIOD'] = re.findall(r'통증\s*기간\s*:\s*(\S+)', cdos_str)[0].strip()
                except:cdai_dict['CDAI_CDOSPAINPERIOD'] = np.nan

                try:cdai_dict['CDAI_CDOSASSOSX'] = re.findall(r'연관\s*증상\s*:\s*(\S+)', cdos_str)[0].strip()
                except:cdai_dict['CDAI_CDOSASSOSX'] = np.nan

                try:cdai_dict['CDAI_CDOSASSOSXPERIOD'] = re.findall(r'연관\s*증상기간\s*:\s*(\S+)', cdos_str)[0].strip()
                except:cdai_dict['CDAI_CDOSASSOSXPERIOD'] = np.nan

                try:cdai_dict['CDAI_CDOSEMERGENCY'] = re.findall(r'응급실\s*방문\s*혹은\s*입원\s*:\s*(\S+)', cdos_str)[0].strip()
                except:cdai_dict['CDAI_CDOSEMERGENCY'] = np.nan
            else:
                cdai_dict['CDAI_CDOS'] = 0
                cdai_dict['CDAI_CDOSPAININT'] = np.nan
                cdai_dict['CDAI_CDOSPAINPERIOD'] = np.nan
                cdai_dict['CDAI_CDOSASSOSX'] = np.nan
                cdai_dict['CDAI_CDOSASSOSXPERIOD'] = np.nan
                cdai_dict['CDAI_CDOSEMERGENCY'] = np.nan
                # raise ValueError

            try:cdai_dict['CDAI_ADHERENCEPCT'] = float(fct_str.split('약제 순응도 (%)')[-1].split('Objective findings')[0].strip())
            except:cdai_dict['CDAI_ADHERENCEPCT'] = np.nan

            cdai_df.append(cdai_dict)

        elif "크론병 질병 활성도 평가를 위한 설문" in fct_str:
            print('크론병 질병 활성도 설문')
            cdsurvey_dict = {'ID': pid, 'NAME': pname}
            cdsurvey_dict['DATE'] = re.findall(r'[\d][\d][\d][\d]-[\d][\d]-[\d][\d]', fct_str)[0]
            fct_str = fct_str.split("< 크론병 질병 활성도 평가를 위한 설문 >")[-1]
            # if '복통 :' in fct_str:
            #     raise ValueError

            try:cdsurvey_dict['CDAI_TOTALSCORE'] = float(re.findall(r'Total\s*:\s*\d+', fct_str)[0].split(':')[-1].strip())
            except:
                try:cdsurvey_dict['CDAI_TOTALSCORE'] = float(re.findall(r'Total\s*\d+', fct_str)[0].split('Total')[-1].strip())
                except:cdsurvey_dict['CDAI_TOTALSCORE'] = np.nan
            # raise ValueError
            try:cdsurvey_dict['CDAI_DIARHEACNT'] = float(re.findall(r'1주일간의\s*설사횟수\s*하루평균\s*\d+', fct_str)[0].split('하루평균')[-1].strip())
            except:
                try:cdsurvey_dict['CDAI_DIARHEACNT'] = float(re.findall(r'설사\s*:\s*하루평균\s*\d+', fct_str)[0].split('하루평균')[-1].strip())
                except:cdsurvey_dict['CDAI_DIARHEACNT'] = np.nan
            try:cdsurvey_dict['CDAI_DIARHEASCORE'] = float(re.findall(r'x 7일  x 2 \s*\d+', fct_str)[0].split('7일  x 2')[-1].strip())
            except:
                try:cdsurvey_dict['CDAI_DIARHEASCORE'] = float(re.findall(r'회\s*x\s*7일\s*x\s*2\s*=\s*\d+', fct_str)[0].split('x 2 =')[-1].strip())
                except:cdsurvey_dict['CDAI_DIARHEASCORE'] = np.nan


            try:cdsurvey_dict['CDAI_ABDPAINDEG'] = float(re.findall(r'복통의 정도\s*0=없음\s*1=경증\s*\d+', fct_str)[0].split('경증')[-1].strip())
            except:cdsurvey_dict['CDAI_ABDPAINDEG'] = np.nan
            try:cdsurvey_dict['CDAI_ABDPAINSCORE'] = float(re.findall(r'\s*7일\s*x\s*5\s*\d+', fct_str)[0].split('x 5')[-1].strip())
            except:
                try:cdsurvey_dict['CDAI_ABDPAINSCORE'] = float(re.findall(r'x\s*7일\s*x\s*5\s*=\s*\d+', fct_str)[0].split('x 5 =')[-1].strip())
                except:cdsurvey_dict['CDAI_ABDPAINSCORE'] = np.nan
            if np.isnan(cdsurvey_dict['CDAI_ABDPAINDEG']) and not np.isnan(cdsurvey_dict['CDAI_ABDPAINSCORE']):
                cdsurvey_dict['CDAI_ABDPAINDEG'] = cdsurvey_dict['CDAI_ABDPAINSCORE']/(7*5)

            try:cdsurvey_dict['CDAI_GENWELLBEING'] = float(re.findall(r'일반적으로 전신 안녕감\s*0\s*=\s*나쁘지 않음\s*\d+', fct_str)[0].split('나쁘지 않음')[-1].strip())
            except:cdsurvey_dict['CDAI_GENWELLBEING'] = np.nan
            try:cdsurvey_dict['CDAI_GENWBSCORE'] = float(re.findall(r'\s*7일\s*x\s*7\s*\d+', fct_str)[0].split('x 7')[-1].strip())
            except:
                try:cdsurvey_dict['CDAI_GENWBSCORE'] = float(re.findall(r'x\s*7일\s*x\s*7\s*=\s*\d+', fct_str)[0].split('x 7 =')[-1].strip())
                except:cdsurvey_dict['CDAI_GENWBSCORE'] = np.nan
            if np.isnan(cdsurvey_dict['CDAI_GENWELLBEING']) and not np.isnan(cdsurvey_dict['CDAI_GENWBSCORE']):
                cdsurvey_dict['CDAI_GENWELLBEING'] = cdsurvey_dict['CDAI_GENWBSCORE']/(7*7)

            try:cdsurvey_dict['CDAI_ARTHSX'] = float(re.findall(r'관절염\/ 관절통\s*0\s*=\s*아니오\s*1\s*=\s*예\s*\d+', fct_str)[0].split('예')[-1].strip())
            except:cdsurvey_dict['CDAI_ARTHSX'] = np.nan
            try:cdsurvey_dict['CDAI_ARTHSXSCORE'] = float(re.findall('x\s*20\s*\d+\s*B\. 홍채염/ 포도막염',fct_str)[0].split('x 20')[-1].split('B. 홍채염/')[0].strip())
            except:
                try:cdsurvey_dict['CDAI_ARTHSXSCORE'] = float(re.findall('관절염\/\s*관절통\s*:\s*0\s*1\s*x\s*20\s*=\s*\d+', fct_str)[0].split('x 20 =')[-1].strip())
                except:cdsurvey_dict['CDAI_ARTHSXSCORE'] = np.nan
            if np.isnan(cdsurvey_dict['CDAI_ARTHSX']) and not np.isnan(cdsurvey_dict['CDAI_ARTHSXSCORE']):
                cdsurvey_dict['CDAI_ARTHSX'] = cdsurvey_dict['CDAI_ARTHSXSCORE']/20

            try:cdsurvey_dict['CDAI_EYESX'] = float(re.findall(r'홍채염\/ 포도막염\s*0\s*=\s*아니오\s*1\s*=\s*예\s*\d+', fct_str)[0].split('예')[-1].strip())
            except:cdsurvey_dict['CDAI_EYESX'] = np.nan
            try:cdsurvey_dict['CDAI_EYESXSCORE'] = float(re.findall('x\s*20\s*\d+\s*C\. 결절홍반/ 괴저농피증/ 아프타구내염',fct_str)[0].split('x 20')[-1].split('C. 결절홍반/')[0].strip())
            except:
                try:cdsurvey_dict['CDAI_EYESXSCORE'] = float(re.findall('홍채염\/ 포도막염\s*:\s*0\s*1\s*x\s*20\s*=\s*\d+', fct_str)[0].split('x 20 =')[-1].strip())
                except:cdsurvey_dict['CDAI_EYESXSCORE'] = np.nan
            if np.isnan(cdsurvey_dict['CDAI_EYESX']) and not np.isnan(cdsurvey_dict['CDAI_EYESXSCORE']):
                cdsurvey_dict['CDAI_EYESX'] = cdsurvey_dict['CDAI_EYESXSCORE']/20

            try:cdsurvey_dict['CDAI_SKINSX'] = float(re.findall(r'결절홍반\/ 괴저농피증\/ 아프타구내염\s*0\s*=\s*아니오\s*1\s*=\s*예\s*\d+', fct_str)[0].split('예')[-1].strip())
            except:cdsurvey_dict['CDAI_SKINSX'] = np.nan
            try:cdsurvey_dict['CDAI_SKINSXSCORE'] = float(re.findall('x\s*20\s*\d+\s*D\. 항문열창, 치루 또는 농양',fct_str)[0].split('x 20')[-1].split('D. 항문열창')[0].strip())
            except:
                try:cdsurvey_dict['CDAI_SKINSXSCORE'] = float(re.findall('결절홍반\/ 괴저농피증\/ 아프타구내염\s*:\s*0\s*1\s*x\s*20\s*=\s*\d+', fct_str)[0].split('x 20 =')[-1].strip())
                except:cdsurvey_dict['CDAI_SKINSXSCORE'] = np.nan
            if np.isnan(cdsurvey_dict['CDAI_SKINSX']) and not np.isnan(cdsurvey_dict['CDAI_SKINSXSCORE']):
                cdsurvey_dict['CDAI_SKINSX'] = cdsurvey_dict['CDAI_SKINSXSCORE']/20

            try:cdsurvey_dict['CDAI_ANALSX'] = float(re.findall(r'항문열창, 치루 또는 농양\s*0\s*=\s*아니오\s*1\s*=\s*예\s*\d+', fct_str)[0].split('예')[-1].strip())
            except:cdsurvey_dict['CDAI_ANALSX'] = np.nan
            try:cdsurvey_dict['CDAI_ANALSXSCORE'] = float(re.findall('x\s*20\s*\d+\s*E\. 기타 누공',fct_str)[0].split('x 20')[-1].split('E. 기타 누공')[0].strip())
            except:
                try:cdsurvey_dict['CDAI_ANALSXSCORE'] = float(re.findall('항문열창, 치루 또는 농양\s*:\s*0\s*1\s*x\s*20\s*=\s*\d+', fct_str)[0].split('x 20 =')[-1].strip())
                except:cdsurvey_dict['CDAI_ANALSXSCORE'] = np.nan
            if np.isnan(cdsurvey_dict['CDAI_ANALSX']) and not np.isnan(cdsurvey_dict['CDAI_ANALSXSCORE']):
                cdsurvey_dict['CDAI_ANALSX'] = cdsurvey_dict['CDAI_ANALSXSCORE']/20

            try:cdsurvey_dict['CDAI_OTHERFIST'] = float(re.findall(r'기타 누공\s*0\s*=\s*아니오\s*1\s*=\s*예\s*\d+', fct_str)[0].split('예')[-1].strip())
            except:cdsurvey_dict['CDAI_OTHERFIST'] = np.nan
            try:cdsurvey_dict['CDAI_OTHERFISTSCORE'] = float(re.findall('x\s*20\s*\d+\s*F\. 최근 7일동안 37',fct_str)[0].split('x 20')[-1].split('F. 최근 7일')[0].strip())
            except:
                try:cdsurvey_dict['CDAI_OTHERFISTSCORE'] = float(re.findall('기타 누공\s*:\s*0\s*1\s*x\s*20\s*=\s*\d+', fct_str)[0].split('x 20 =')[-1].strip())
                except:cdsurvey_dict['CDAI_OTHERFISTSCORE'] = np.nan
            if np.isnan(cdsurvey_dict['CDAI_OTHERFIST']) and not np.isnan(cdsurvey_dict['CDAI_OTHERFISTSCORE']):
                cdsurvey_dict['CDAI_OTHERFIST'] = cdsurvey_dict['CDAI_OTHERFISTSCORE']/20

            try:cdsurvey_dict['CDAI_FEVER'] = float(re.findall(r'\(구강\)\/\s*0\s*=\s*아니오\s*1\s*=\s*예\s*\d+', fct_str)[0].split('예')[-1].strip())
            except:cdsurvey_dict['CDAI_FEVER'] = np.nan
            try:cdsurvey_dict['CDAI_FEVERSCORE'] = float(re.findall('x\s*20\s*\d+\s*38.3 ℃ \(항문\) 이상의 열',fct_str)[0].split('x 20')[-1].split('38.3 ℃')[0].strip())
            except:
                try:cdsurvey_dict['CDAI_FEVERSCORE'] = float(re.findall('\(항문\) 이상의 열\s*:\s*0\s*1\s*x\s*20\s*=\s*\d+', fct_str)[0].split('x 20 =')[-1].strip())
                except:cdsurvey_dict['CDAI_FEVERSCORE'] = np.nan
            if np.isnan(cdsurvey_dict['CDAI_FEVER']) and not np.isnan(cdsurvey_dict['CDAI_FEVERSCORE']):
                cdsurvey_dict['CDAI_FEVER'] = cdsurvey_dict['CDAI_FEVERSCORE']/20

            try:cdsurvey_dict['CDAI_ANTIDIARHMEDI'] = float(re.findall(r'최근 7일동안 지사제 치료를\s*0\s*=\s*아니오\s*1\s*=\s*예\s*\d+', fct_str)[0].split('예')[-1].strip())
            except:cdsurvey_dict['CDAI_ANTIDIARHMEDI'] = np.nan
            try:cdsurvey_dict['CDAI_ANTIDIARHMEDISCORE'] = float(re.findall('x\s*30\s*\d+\s*받은 적이 있는 경우',fct_str)[0].split('x 30')[-1].split('받은 적이 있는 경우')[0].strip())
            except:
                try:cdsurvey_dict['CDAI_ANTIDIARHMEDISCORE'] = float(re.findall('지사제 치료를 받은적이 있는 경우\s*:\s*0\s*1\s*x\s*30\s*=\s*\d+', fct_str)[0].split('x 30 =')[-1].strip())
                except: cdsurvey_dict['CDAI_ANTIDIARHMEDISCORE'] = np.nan
            if np.isnan(cdsurvey_dict['CDAI_ANTIDIARHMEDI']) and not np.isnan(cdsurvey_dict['CDAI_ANTIDIARHMEDISCORE']):
                cdsurvey_dict['CDAI_ANTIDIARHMEDI'] = cdsurvey_dict['CDAI_ANTIDIARHMEDISCORE']/30

            try:cdsurvey_dict['CDAI_ABDMASS'] = float(re.findall(r'복부 종괴\s*0\s*=\s*없음\s*\d+', fct_str)[0].split('없음')[-1].strip())
            except:cdsurvey_dict['CDAI_ABDMASS'] = np.nan
            try:cdsurvey_dict['CDAI_ABDMASSSCORE'] = float(re.findall('x\s*10\s*\d+\s*2\s*=\s*의심됨',fct_str)[0].split('x 10')[-1].split('2 \t \t= 의심됨')[0].strip())
            except:
                try:cdsurvey_dict['CDAI_ABDMASSSCORE'] = float(re.findall('복부 종괴\s*:\s*0\s*2\s*5\s*x\s*10\s*=\s*\d+', fct_str)[0].split('x 10 =')[-1].strip())
                except:cdsurvey_dict['CDAI_ABDMASSSCORE'] = np.nan
            if np.isnan(cdsurvey_dict['CDAI_ABDMASS']) and not np.isnan(cdsurvey_dict['CDAI_ABDMASSSCORE']):
                cdsurvey_dict['CDAI_ABDMASS'] = cdsurvey_dict['CDAI_ABDMASSSCORE']/10


            try:cdsurvey_dict['CDAI_HEMATOCRIT_MALE'] = float(re.findall(r'헤마토크릿\s*남성\s*\:\s*[\(]*\s*47\s*\-\s*\d+', fct_str)[0].split('-')[-1].strip())
            except:cdsurvey_dict['CDAI_HEMATOCRIT_MALE'] = np.nan
            try: cdsurvey_dict['CDAI_HEMATOCRIT_FEMALE'] = float(re.findall(r'여성\s*\:\s*[\(]*\s*42\s*\-\s*\d+', fct_str)[0].split('-')[-1].strip())
            except:cdsurvey_dict['CDAI_HEMATOCRIT_FEMALE'] = np.nan
            if not np.isnan(cdsurvey_dict['CDAI_HEMATOCRIT_FEMALE']):
                cdsurvey_dict['CDAI_HEMATOCRIT'] = cdsurvey_dict['CDAI_HEMATOCRIT_FEMALE']
            elif not np.isnan(cdsurvey_dict['CDAI_HEMATOCRIT_MALE']):
                cdsurvey_dict['CDAI_HEMATOCRIT'] = cdsurvey_dict['CDAI_HEMATOCRIT_MALE']
            else:
                cdsurvey_dict['CDAI_HEMATOCRIT'] = np.nan

            try:
                cdsurvey_dict['CDAI_HCTSCORE'] = float(re.findall(r'x\s*6\s*\d+', fct_str)[0].split('x 6')[-1].strip())
                cdsurvey_dict['CDAI_ADJHCTSCORE'] = cdsurvey_dict['CDAI_HCTSCORE']
            except:
                try:
                    cdai_hctscore_pattern = re.findall('x\s*6\s*=\s*\d+\s*조정값\s*\d+', fct_str)[0]
                    cdsurvey_dict['CDAI_HCTSCORE'] = float(cdai_hctscore_pattern.split('x 6 =')[-1].split('조정값')[0].strip())
                    cdsurvey_dict['CDAI_ADJHCTSCORE'] = float(cdai_hctscore_pattern.split('조정값')[-1].strip())
                except:
                    cdsurvey_dict['CDAI_HCTSCORE'] = np.nan
                    cdsurvey_dict['CDAI_ADJHCTSCORE'] = np.nan

            try:cdsurvey_dict['CDAI_STDWEIGHT'] = float(re.findall(r'표준체중\s*:?\s*\d+[\.]?\d*\s*', fct_str)[0].split('표준체중')[-1].strip())
            except:cdsurvey_dict['CDAI_STDWEIGHT'] = np.nan
            try:cdsurvey_dict['CDAI_WEIGHT'] = float(re.findall(r'실제체중\s*:?\s*\d+[\.]?\d*\s*', fct_str.replace('(kg)', ''))[0].split(':')[-1].split('체중')[-1].split('kg')[0].strip())
            except:cdsurvey_dict['CDAI_WEIGHT'] = np.nan

            try:
                cdsurvey_dict['CDAI_WTSCORE'] = float(re.findall(r'x100\s*-?\d+', fct_str)[0].split('100')[-1].strip())
                cdsurvey_dict['CDAI_ADJWTSCORE'] = float(fct_str.split("Crohn’s Disease Obstructive Score")[0].split("조정값")[-1].strip())
            except:
                try:
                    cdai_wtscore_pattern = re.findall('x 100\s*=\s*-?\d+\s*조정값\s*-?\d+', fct_str)[0]
                    cdsurvey_dict['CDAI_WTSCORE'] = float(cdai_wtscore_pattern.split('x 100 =')[-1].split("조정값")[0].strip())
                    cdsurvey_dict['CDAI_ADJWTSCORE'] =  float(cdai_wtscore_pattern.split("조정값")[-1].strip())
                except:
                    cdsurvey_dict['CDAI_WTSCORE'] = np.nan
                    cdsurvey_dict['CDAI_ADJWTSCORE'] = np.nan
            # raise ValueError

            try:cdsurvey_dict['CDAI_HEIGHT'] = float(re.findall(r'환자신장\s*:?\s*\d+[\.]?\d*\s*', fct_str.replace('(cm)', ''))[0].split(':')[-1].split('신장')[-1].split('cm')[0].strip())
            except:
                try:cdsurvey_dict['CDAI_HEIGHT'] = float(re.findall(r'신장\s*:?\s*\d+[\.]?\d*\s*', fct_str.replace('(cm)', ''))[0].split(':')[-1].split('신장')[-1].split('cm')[0].strip())
                except:cdsurvey_dict['CDAI_HEIGHT'] = np.nan

            try:cdsurvey_dict['CDAI_BMI'] = float(re.findall(r'BMI\s*:?\s*\d+[\.]?\d*\s*', fct_str.replace('(kg/㎡)', ''))[0].split(':')[-1].split('BMI')[-1].split('kg')[0].strip())
            except:
                try:cdsurvey_dict['CDAI_BMI'] = float(re.findall(r'BMI\s*:?\s*\d+[\.]?\d*\s*', fct_str.replace('(kg/㎡)', ''))[0].split(':')[-1].split('BMI')[-1].split('kg')[0].strip())
                except:cdsurvey_dict['CDAI_BMI'] = np.nan


            # fct_cdos_split = fct_str.split('Crohn’s Disease Obstructive Score')
            # cdos_str = fct_cdos_split[-1].split('약제 순응도')[0].strip()
            # if cdos_str.split('통증 강도 :')[-1].split('연관 증상 :')[0].split('통증 기간 :')[0].strip() == '없음':
            #     cdsurvey_dict['CDAI_CDOS'] = 0
            #     cdsurvey_dict['CDAI_CDOSPAININT'] = np.nan
            #     cdsurvey_dict['CDAI_CDOSPAINPERIOD'] = np.nan
            #     cdsurvey_dict['CDAI_CDOSASSOSX'] = np.nan
            #     cdsurvey_dict['CDAI_CDOSASSOSXPERIOD'] = np.nan
            #     cdsurvey_dict['CDAI_CDOSEMERGENCY'] = np.nan
            #
            # elif len(fct_cdos_split) > 1:
            #     cdsurvey_dict['CDAI_CDOS'] = 1
            #
            #     try:cdsurvey_dict['CDAI_CDOSPAININT'] = re.findall(r'통증\s*강도\s*:\s*(\S+)', cdos_str)[0].strip()
            #     except:cdsurvey_dict['CDAI_CDOSPAININT'] = np.nan
            #
            #     try:cdsurvey_dict['CDAI_CDOSPAINPERIOD'] = re.findall(r'통증\s*기간\s*:\s*(\S+)', cdos_str)[0].strip()
            #     except:cdsurvey_dict['CDAI_CDOSPAINPERIOD'] = np.nan
            #
            #     try:cdsurvey_dict['CDAI_CDOSASSOSX'] = re.findall(r'연관\s*증상\s*:\s*(\S+)', cdos_str)[0].strip()
            #     except:cdsurvey_dict['CDAI_CDOSASSOSX'] = np.nan
            #
            #     try:cdsurvey_dict['CDAI_CDOSASSOSXPERIOD'] = re.findall(r'연관\s*증상기간\s*:\s*(\S+)', cdos_str)[0].strip()
            #     except:cdsurvey_dict['CDAI_CDOSASSOSXPERIOD'] = np.nan
            #
            #     try:cdsurvey_dict['CDAI_CDOSEMERGENCY'] = re.findall(r'응급실\s*방문\s*혹은\s*입원\s*:\s*(\S+)', cdos_str)[0].strip()
            #     except:cdsurvey_dict['CDAI_CDOSEMERGENCY'] = np.nan
            # else:
            #     cdsurvey_dict['CDAI_CDOS'] = 0
            #     cdsurvey_dict['CDAI_CDOSPAININT'] = np.nan
            #     cdsurvey_dict['CDAI_CDOSPAINPERIOD'] = np.nan
            #     cdsurvey_dict['CDAI_CDOSASSOSX'] = np.nan
            #     cdsurvey_dict['CDAI_CDOSASSOSXPERIOD'] = np.nan
            #     cdsurvey_dict['CDAI_CDOSEMERGENCY'] = np.nan
            #     # raise ValueError
            #
            # try:cdsurvey_dict['CDAI_ADHERENCEPCT'] = float(fct_str.split('약제 순응도 (%)')[-1].split('Objective findings')[0].strip())
            # except:cdsurvey_dict['CDAI_ADHERENCEPCT'] = np.nan

            cdsurvey_df.append(cdsurvey_dict)
        else:
            continue
            raise ValueError

        # break
pms_df = pd.DataFrame(pms_df)
pms_df.to_csv(f"{output_dir}/pdmarker_pms_df.csv", encoding='utf-8-sig', index=False)
mss_df = pd.DataFrame(mss_df)
mss_df.to_csv(f"{output_dir}/pdmarker_mss_df.csv", encoding='utf-8-sig', index=False)
cdai_df = pd.DataFrame(cdai_df)
cdai_df.to_csv(f"{output_dir}/pdmarker_cdai_df.csv", encoding='utf-8-sig', index=False)
cdsurvey_df = pd.DataFrame(cdsurvey_df)
cdsurvey_df.to_csv(f"{output_dir}/pdmarker_cdsurvey_df.csv", encoding='utf-8-sig', index=False)

print('COMPLETED')