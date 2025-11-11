import pandas as pd
import numpy as np
import glob
import os
from scipy.stats import mannwhitneyu, fisher_exact
import matplotlib.pyplot as plt
import seaborn as sns
# 성별 차이 Subgroup analysis

# 데이터 로드

output_dir = 'C:/Users/ilma0/PycharmProjects/pypharmacometrics/Projects/LINEZOLID/results'

subgroup_files = glob.glob(f"{output_dir}/b1da/mvlreg_output/datasubset/b1da_lnz_mvlreg_datasubset(Total_Adult)(Lactate).csv")
for file_path in subgroup_files: #break
    # raise ValueError
    subset_group = file_path.split(')(')[0].split('(')[-1]
    pd_endpoint = file_path.split(')(')[-1].split(')')[0]

    df = pd.read_csv(file_path)

    for group_col in ['SEX', 'EV']:
        numeric_cols = [c for c in df.columns if len(df[c].unique())>2]

        rows = []
        rows.append({'Covar':'N',
                     f'{group_col} = 0':len(df[df[group_col]==0]),
                     f'{group_col} = 1':len(df[df[group_col]==1]),
                     'p_value': 0,
                     })

        ###################################################
        # Subgroup간 Covariates 값 비교
        ###################################################

        # 1️⃣ 연속형 변수: mean (SD)와 Mann-Whitney U p-value
        for col in numeric_cols:
            g0 = df.loc[df[group_col] == 0, col].dropna()
            g1 = df.loc[df[group_col] == 1, col].dropna()
            pval = mannwhitneyu(g0, g1, alternative='two-sided')[1] if len(g0) and len(g1) else np.nan
            rows.append({
                'Covar': col,
                f'{group_col} = 0': f"{g0.mean():.1f} ({g0.std():.1f})",
                f'{group_col} = 1': f"{g1.mean():.1f} ({g1.std():.1f})",
                'p_value': round(pval,3)
            })

        # 2️⃣ EV 컬럼: N (%)와 Fisher의 정확 검정
        g0_total = (df[group_col] == 0).sum()
        g1_total = (df[group_col] == 1).sum()
        binary_cols = [bc for bc in list(set(df.columns)-set(numeric_cols)) if bc!=group_col]
        for binary_col in binary_cols:
            g0_ev1 = df.loc[df[group_col] == 0, binary_col].sum()
            g1_ev1 = df.loc[df[group_col] == 1, binary_col].sum()

            contingency = [[g0_ev1, g0_total - g0_ev1],
                           [g1_ev1, g1_total - g1_ev1]]
            pval_ev = fisher_exact(contingency)[1]

            rows.append({
                'Covar': binary_col,
                f'{group_col} = 0': f"{g0_ev1} ({g0_ev1 / g0_total * 100:.1f})",
                f'{group_col} = 1': f"{g1_ev1} ({g1_ev1 / g1_total * 100:.1f})",
                'p_value': round(pval_ev,3)
            })

        # raise ValueError

        # 3️⃣ 결과 테이블 생성
        demographic_table = pd.DataFrame(rows).sort_values('p_value').reset_index(drop=True)
        demographic_table_saving = demographic_table.copy()
        demographic_table_saving['p_value'] = demographic_table_saving['p_value'].replace(0,'<0.001')
        if not os.path.exists(f'{output_dir}/b1da/mvlreg_output/subgroup_analysis'):
            os.mkdir(f'{output_dir}/b1da/mvlreg_output/subgroup_analysis')
        demographic_table_saving.to_csv(f"{output_dir}/b1da/mvlreg_output/subgroup_analysis/b1da_lnz_subgrpa({subset_group})({pd_endpoint})_({group_col} subgroup).csv", index=False, encoding='utf-8-sig')
        # print(demographic_table)

    ###################################################
    # (Event 발생 그룹 vs 미발생 그룹) 에서의 WT
    ###################################################

    # 색상 설정 (고급스러운 tone)
    colors = ['#1f77b4', '#ff7f0e']  # EV=0 → blue, EV=1 → orange
    ggdf = df.copy()
    ggdf['SEX'] = ggdf['SEX'].map({0:'Male',1:'Female'})
    # for y2type, unit_str in {'WT':('(kg)','Weight'),'DOSE24PERWT':('(mg/kg/day)','Weight-normalized mean daily dose')}.items():
    #
    #     plt.figure(figsize=(10, 8))
    #     sns.boxplot(
    #         data=ggdf,
    #         x='EV',
    #         y=y2type,
    #         hue='SEX',
    #         hue_order=['Male','Female'],
    #         palette=colors,
    #         width=0.5,
    #         showfliers=False
    #     )
    #
    #     # 스타일링
    #     # plt.xticks([0, 1], ['EV = 0', 'EV = 1'])
    #     plt.yticks(fontsize=12)
    #     plt.xticks([0, 1], ['Event = 0', 'Event = 1'],fontsize=12)
    #     plt.xlabel('Lactic acidosis event group',fontsize=12)
    #     # plt.ylabel(f'{y2type} {unit_str[0]}')
    #     plt.ylabel(f'{unit_str[1]} {unit_str[0]}',fontsize=12)
    #     # plt.title(f'Comparison of {y2type} by Event Occurrence')
    #     plt.title(f'Comparison of {unit_str[1]} by event occurrence',fontsize=12)
    #     plt.grid(alpha=0.2)
    #     plt.tight_layout()
    #
    #     plt.savefig(f"{output_dir}/b1da/mvlreg_output/subgroup_analysis/EV_to_{y2type}_boxplot({pd_endpoint}).png")  # PNG 파일로 저장
    #
    #     plt.cla()
    #     plt.clf()
    #     plt.close()

    fig, axes = plt.subplots(2, 1, figsize=(12, 15), sharex=True)

    plot_items = {
        'WT': ('(kg)', 'weight'),
        'DOSE24PERWT': ('(mg/kg/day)', 'weight-normalized mean daily dose')
    }
    x_order = [0, 1]
    ggfontsize = 18

    for i, (y2type, unit_str) in enumerate(plot_items.items()):
        i_alph = 'A' if i == 0 else 'B'
        ax = axes[i]

        sns.boxplot(
            data=ggdf,
            x='EV', y=y2type,
            hue='SEX', hue_order=['Male', 'Female'],
            palette=colors, width=0.5, showfliers=False,
            order=x_order, ax=ax
        )
        ax.set_xticklabels(['Event = 0', 'Event = 1'], fontsize=ggfontsize)
        ax.set_xlabel('Lactic acidosis event group', fontsize=ggfontsize)
        ax.set_ylabel(f'{unit_str[1]} {unit_str[0]}', fontsize=ggfontsize)
        ax.set_title(f'({i_alph}) Comparison of {unit_str[1]} by event occurrence', fontsize=ggfontsize)
        ax.grid(alpha=0.2)

    # ⬇⬇⬇ 추가: 위/아래 subplot 모두 x축 라벨 보이게
    axes[0].tick_params(axis='x', which='both', labelbottom=True)
    axes[1].tick_params(axis='x', which='both', labelbottom=True)

    # 통합 legend
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, title="Sex", loc='upper right', fontsize=ggfontsize, title_fontsize=ggfontsize)
    axes[0].legend_.remove()
    axes[1].legend_.remove()

    plt.tight_layout(rect=[0, 0, 0.95, 1])
    plt.savefig(f"{output_dir}/b1da/mvlreg_output/subgroup_analysis/EV_to_WT_DOSE24PERWT_boxplot({pd_endpoint}).png")
    plt.cla()
    plt.clf()
    plt.close()

    ###################################################
    # Subgroup간 Covariates 값 비교
    ###################################################

    # 1️⃣ 체중을 10 단위 구간으로 나누기 (10~150)
    df['WT_BIN'] = pd.cut(df['WT'], bins=range(10, 160, 10), right=False)

    # 2️⃣ SEX별, 체중구간별 Event 발생률 계산
    rate_df = (
        df.groupby(['SEX', 'WT_BIN'])['EV'].agg(['sum', 'count']).reset_index()
    )

    rate_df['EV_count'] = rate_df['sum']
    rate_df['EV_rate'] = rate_df['sum']/rate_df['count']

    # 3️⃣ 색상 정의 (고급스러운 톤)
    colors = {
        0: '#1f77b4',  # navy blue
        1: '#ff7f0e'  # muted orange
    }

    for ytype in ['EV_count','EV_rate']:
        # 4️⃣ subplot 설정
        fig, axes = plt.subplots(1, 2, figsize=(12, 5), sharey=True)

        for i, sex in enumerate([0, 1]):
            subset = rate_df[rate_df['SEX'] == sex]

            axes[i].bar(
                subset['WT_BIN'].astype(str),
                subset[ytype],
                color=colors[sex],
                alpha=0.8,
                width=0.7
            )

            # dashed line y=0, y=1
            axes[i].axhline(0, color='gray', linestyle='--', linewidth=1)
            axes[i].axhline(1, color='gray', linestyle='--', linewidth=1)

            ylimval = rate_df[ytype].max() * 1.1
            axes[i].set_ylim(0, ylimval)
            axes[i].set_xlabel('Weight category (kg)')
            axes[i].set_ylabel(f'{ytype.replace("_"," ")} (EV=1 ratio)' if i == 0 else '')
            axes[i].set_title(f'SEX = {sex} ({"Male" if sex == 0 else "Female"})')
            axes[i].tick_params(axis='x', rotation=45)
            axes[i].grid(alpha=0.2)

        plt.tight_layout()
        # plt.tight_layout()
        plt.savefig(f"{output_dir}/b1da/mvlreg_output/subgroup_analysis/SubGrp(SEX)_WTbin_to_{ytype.replace('_','')}_plot({pd_endpoint}).png")  # PNG 파일로 저장

        plt.cla()
        plt.clf()
        plt.close()

    ###################################################
    # SEX subgroup 내에서 WT vs EV
    ###################################################

    # SEX subgroup별로 나누기
    df_0 = df[df['SEX'] == 0].copy()
    df_1 = df[df['SEX'] == 1].copy()

    # SEX=0: 남성(blue계열), SEX=1: 여성(orange계열)
    colors = {
        0: '#1f77b4',  # 고급스러운 navy blue
        1: '#ff7f0e'  # muted orange
    }

    # Subplot 생성 (y축 공유)
    fig, axes = plt.subplots(1, 2, figsize=(12, 6), sharey=True)

    # SEX=0 (left subplot)
    axes[0].scatter(df_0['WT'], df_0['EV'], color=colors[0], alpha=0.7, edgecolor='none')
    axes[0].axhline(0, color='gray', linestyle='--', linewidth=1)
    axes[0].axhline(1, color='gray', linestyle='--', linewidth=1)
    axes[0].set_xlim(10, 150)
    axes[0].set_ylim(0, 1.2)
    axes[0].set_xlabel('Weight (WT)')
    axes[0].set_ylabel('Event (EV)')
    axes[0].set_title('SEX = 0 (Male)')
    axes[0].grid(alpha=0.2)

    # SEX=1 (right subplot)
    axes[1].scatter(df_1['WT'], df_1['EV'], color=colors[1], alpha=0.7, edgecolor='none')
    axes[1].axhline(0, color='gray', linestyle='--', linewidth=1)
    axes[1].axhline(1, color='gray', linestyle='--', linewidth=1)
    axes[1].set_xlim(10, 150)
    axes[1].set_ylim(0, 1.2)
    axes[1].set_xlabel('Weight (WT)')
    axes[1].set_title('SEX = 1 (Female)')
    axes[1].grid(alpha=0.2)

    # 전체 스타일 정리
    plt.tight_layout()
    # plt.tight_layout()
    plt.savefig(f"{output_dir}/b1da/mvlreg_output/subgroup_analysis/SubGrp(SEX)_WT_to_EV_plot({pd_endpoint}).png")  # PNG 파일로 저장

    plt.cla()
    plt.clf()
    plt.close()

    # for sgval, sg_df in df.groupby('SEX'):
    #     df['WT'].max()