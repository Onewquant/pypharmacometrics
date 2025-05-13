from tools import *
from pynca.tools import *


def load_data_dict(drug_list, filename_format, input_file_dir_path):
    drug_prep_df_dict = dict()
    for drug in drug_list:
        result_file_path = f"{input_file_dir_path}/" + filename_format.replace('[drug]',drug)
        if filename_format.split('.')[-1]=='csv':
            drug_prep_df_dict[drug] = pd.read_csv(result_file_path)
        if filename_format.split('.')[-1] == 'xls':
            drug_prep_df_dict[drug] = pd.read_excel(result_file_path)
        # drug_prep_df_dict[drug]['FEEDING'] = drug_prep_df_dict[drug]['FEEDING'].replace('FASTING','FASTED')
        # drug_prep_df_dict[drug]['Subject'] = drug_prep_df_dict[drug].apply(lambda row:f'{row["ID"]}|{row["FEEDING"]}',axis=1)
    return drug_prep_df_dict

result_type = 'Phoenix'
result_type = 'R'

drug_list = ['W2406']
input_file_dir_path = 'C:/Users/ilma0/PycharmProjects/pypharmacometrics/Projects/glpharma/resource/prep_data'
result_file_dir_path = 'C:/Users/ilma0/PycharmProjects/pypharmacometrics/Projects/glpharma/resource/Figures'

drug_prep_df_dict = load_data_dict(drug_list=drug_list, filename_format=f"GLPHARMA_ConcPrep_[drug]_{result_type}.csv", input_file_dir_path=input_file_dir_path)
# drug_prep_df_dict['W2406']['CONC'].max()
# drug_prep_df_dict['W2406']['CONC'].std()
# drug_prep_df_dict['W2406']['CONC'].mean()
# drug_prep_df_dict['W2406']['ATIME'].max()

def time_to_conc_graph_ckd(gdf, sid_list, drug, hue, result_file_dir_path, hue_order=None, file_format='png', dpi=300, estimator=np.mean, yscale='linear', save_fig=True):

    g_palette = 'Dark2'
    g_palette_colors = sns.color_palette('Dark2')
    sns.set_style("whitegrid", {'grid.linestyle': ':',
                                })

    mode = 'Individual' if len(sid_list)==1 else 'Population'
    if mode=='Individual':
        title_str = sid_list[0]
        last_tag = '('+sid_list[0]+')'
        time_col = 'ATIME'
    else:
        # if errorbar[0]=='sd':
        #     if errorbar[1]==1: errorbar_str = f' ({errorbar[0].upper()})'
        #     else: errorbar_str = f' ({errorbar[1]} {errorbar[0].upper()})'
        # elif errorbar[0]=='ci':
        #     errorbar_str = f' ({errorbar[1]}% {errorbar[0].upper()})'
        # else:
        #     errorbar_str = ''
        # title_str = f'Sample Mean{errorbar_str}'
        # last_tag = 'sample'+str(tuple(sid_list)).replace(",)",")").replace("'","")
        last_tag = 'sample'
        time_col = 'NTIME'
    filename = f'{mode}_{drug}_{last_tag}'

    act_gdf = gdf[gdf['ID'].isin(sid_list)].copy()

    # marker_list = ['o', '^', 'v','<', '>', 's', 'p', '*', 'h', 'H', 'D', 'd', '+', 'x', '|', '_']
    marker_list = ['o']*36
    # g = sns.relplot(data=act_gdf, x=time_col,y='CONC', palette=g_palette, marker='o',hue=hue, hue_order=hue_order, markersize=7, markeredgecolor='white', markeredgewidth=1, kind='line', linewidth=1.5, linestyle='--', errorbar=errorbar, estimator=estimator, err_style=err_style)
    g = sns.relplot(data=act_gdf, x=time_col, y='CONC', palette=g_palette, markers=marker_list[:len(hue_order)], hue=hue, hue_order=hue_order, style=hue, style_order=hue_order, markersize=10, markeredgecolor='white', markeredgewidth=1, kind='line', linewidth=2, estimator=estimator, errorbar=None)
    # g = sns.relplot(data=act_gdf, x=time_col, y='CONC', palette=g_palette, marker='o', hue=hue, hue_order=hue_order, markersize=7, markeredgecolor='white', markeredgewidth=1, kind='line', linewidth=1.5, linestyle='--', estimator=estimator, errorbar=None)
    # errorbar = ("sd", 2), err_style = 'band',
    # plt.setp(plt.gca().get_lines()[1], fillstyle='none')

    if mode=='Population':

        ## 에러바 데이터 생성

        eb_df_dict = dict()
        for hue_inx, hue_act_gdf in act_gdf.groupby(hue):
            for_eb_df = hue_act_gdf.groupby('NTIME')['CONC'].agg([np.nanmean, np.nanstd]).reset_index(drop=False)
            eb_x = tuple(for_eb_df['NTIME'])
            eb_y = tuple(for_eb_df['nanmean'])
            eb_y_errbar = tuple(for_eb_df['nanstd'])

            eb_df_dict[hue_inx] = {'for_eb_df':for_eb_df,'eb_x':eb_x, 'eb_y':eb_y, 'eb_y_errbar':eb_y_errbar}

        hue_order_dict = dict([(ho,i) for i, ho in enumerate(hue_order)])
        for hue_eb_key, hue_eb_val in eb_df_dict.items():
            g.ax.errorbar(hue_eb_val['eb_x'], hue_eb_val['eb_y'], yerr=[tuple(np.zeros(len(eb_y))), hue_eb_val['eb_y_errbar']], fmt='o', ecolor=g_palette_colors[hue_order_dict[hue_eb_key]], capsize=2.5, capthick=2,barsabove=True)

    else:
        plt.title(sid_list[0], fontsize=18)

    # eb.get_children()[3].set_linestyle('--')  ## 에러 바 라인 스타일
    # eb.get_children()[1].set_marker('v') ## 에러 바 아래쪽 마커 스타일
    # eb.get_children()[2].set_marker('^') ## 에러 바 위쪽 마커 스타일
    # palette = sns.color_palette('Dark2')


    # g.fig.set_size_inches(15,11)
    g.fig.set_size_inches(15, 11)
    # g.fig.subplots_adjust(left=0.1, right=0.1)

    if yscale=="log": g.set(yscale="log")
    else: pass
    # g.set_axis_labels('Time (hr)', 'Concentration (mg/L)')
    # sns.move_legend(g, 'upper right', frameon=True)
    # g.fig.subplots_adjust(top=0.85)

    sns.move_legend(g, 'center right', title=None, frameon=False, fontsize=18)
    # if mode == 'Population':
    #     sns.move_legend(g, 'upper right', title=None, frameon=False, fontsize=18)
    # else:
    #     sns.move_legend(g, 'center right', title=None, frameon=False, fontsize=18)
    # sns.move_legend(g, 'upper center', ncol=2, title=None, frameon=False, fontsize=15)
    # g.fig.suptitle("A001", fontsize=20, fontweight='bold')

    plt.tight_layout(pad=3.5)

    plt.xlabel('Time (h)', fontsize=20, labelpad=8)
    plt.ylabel(f'{drug} plasma concentration (μg/L)', fontsize=20, labelpad=8)

    # plt.xticks(np.arange(-6,30, step=6), fontsize=18)
    # plt.xlim(-1,30)

    plt.xticks(np.arange(-1, 4, step=1), fontsize=18)
    plt.xlim(-1, 4)

    if drug=='W2406':
        if yscale=='linear':
            # plt.yticks(np.linspace(0, 2500, 11, endpoint=True), fontsize=18)
            # plt.ylim(-50, 2500)
            plt.yticks(np.linspace(0, 10, 11, endpoint=True), fontsize=18)
            plt.ylim(-1, 10)
        elif yscale=='log':
            plt.yticks([0,1,10], fontsize=18)
            plt.ylim(1, 10)


    if save_fig:
        if not os.path.exists(f"{result_file_dir_path}"): os.mkdir(f"{result_file_dir_path}")
        if not os.path.exists(f"{result_file_dir_path}/{yscale}"): os.mkdir(f"{result_file_dir_path}/{yscale}")
        if not os.path.exists(f"{result_file_dir_path}/{yscale}/{mode}"): os.mkdir(f"{result_file_dir_path}/{yscale}/{mode}")
        if not os.path.exists(f"{result_file_dir_path}/{yscale}/{mode}/{drug}"): os.mkdir(f"{result_file_dir_path}/{yscale}/{mode}/{drug}")
        plt.savefig(f"{result_file_dir_path}/{yscale}/{mode}/{drug}/{filename}.{file_format}", dpi=dpi)



############################

hue = 'TRT'
hue_order = ['T','R']
estimator=np.mean

for yscale in ['linear','log']:
    for drug in drug_list:
        # if (drug!='Metformin'):
        #     continue

        gdf = drug_prep_df_dict[drug]

        ## Population

        time_to_conc_graph_ckd(gdf=gdf, sid_list=list(gdf['ID'].unique()), drug=drug, hue=hue, result_file_dir_path=result_file_dir_path, hue_order=hue_order, estimator=estimator, yscale=yscale, save_fig=True)

        plt.cla()
        plt.clf()
        plt.close()

        ## Individual

        for sid in gdf['ID'].unique():

            time_to_conc_graph_ckd(gdf=gdf, sid_list=[sid,], drug=drug, hue=hue, result_file_dir_path=result_file_dir_path, hue_order=hue_order, estimator=estimator, yscale=yscale, save_fig=True)

            plt.cla()
            plt.clf()
            plt.close()

############################


fpp_df = pd.read_excel(f"{input_file_dir_path}/Final Parameters Pivoted - BSH5.xls")
for drug in drug_list:
    fpp_df[fpp_df['DRUG']==drug].reset_index(drop=True).to_csv(f"{input_file_dir_path}/Final Parameters Pivoted - ({drug}).csv", index=False)


fpp_df = fpp_df.iloc[1:].copy()
comp_criteria = 'TRT'
params = ['AUClast','Cmax']
# ylabel_dict = {'AUClast':r'$AUC_t$ (ng·h/mL)','Cmax':r'$C_{max}$ (ng/mL)'}
ylabel_dict = {'AUClast':r'$AUC_t$ (μg·h/L)','Cmax':r'$C_{max}$ (μg/L)'}
save_fig=True
# hue='id', style=None, palette=["blue"]
comp_order = list(fpp_df[comp_criteria].unique())
if len(comp_order)==2:
    # x_position_dict = dict(zip(comp_order,[1/5, 4/5]))
    x_position_dict = dict(zip(comp_order,[1/3, 2/3]))
else:
    raise ValueError

# yscale = "Linear"
# drug = 'Lobeglitazone'
# p = 'Cmax'

# for yscale in ["linear","log"]:
for yscale in ["linear", ]:
    for drug in drug_list:
        for p in params:

            # fpp_df[fpp_df['ID'].isin(['A001','A002','A003']) & (fpp_df['DRUG']=='Empagliflozin')][['ID','FEEDING','Cmax']].to_csv(f"{input_file_dir_path}/ind_plot_demo.csv")
            dfpp_df = fpp_df[fpp_df['DRUG']==drug].copy()
            dfpp_df[f'{comp_criteria}_NUM'] = dfpp_df[comp_criteria].map(x_position_dict)

            # if (drug=='Empagliflozin') and (p=='AUClast'):
            #     raise ValueError

            plt.figure(figsize=(10, 10))

            # Draw the line plot with hollow markers and dark grey color for lines and markers
            g = sns.lineplot(
                x=f'{comp_criteria}_NUM',
                y=p,
                # hue='ID',
                data=dfpp_df,
                marker='o',
                markersize=12,  # Increase the marker size
                linestyle='-',  # Line segments
                # style='ID',
                hue='ID',
                # palette=["darkgrey"]*len(fpp_df['ID'].unique()),
                # palette=["black"] * len(fpp_df['ID'].unique()),
                palette=["dimgrey"] * len(fpp_df['ID'].unique()),
                # color='darkgrey',  # Set the line and marker edge color to dark grey
                markerfacecolor=(1, 1, 1, 0),  # Keep the marker hollow
                # markeredgecolor='black',  # Dark grey marker edge
                markeredgecolor='dimgrey',  # Dark grey marker edge
                # markeredgecolor='darkgrey',  # Dark grey marker edge
                legend=False,
                errorbar=('ci', False)
            )

            if yscale=="log":
                g.set(yscale="log")

            # Customize the x-axis to show FASTED and FED at custom positions


            plt.xticks(list(x_position_dict.values()), list(x_position_dict.keys()), fontsize=18)

            plt.tight_layout(pad=5)

            plt.xlim(0, 1)
            # ds = dfpp_df[p]
            ymaxlim_value = get_lim_num_for_graph(dfpp_df[p])
            ymin_digit_num = 10**(get_digit_count(dfpp_df[p].min())-1)
            yminlim_value = ymin_digit_num*1.5 if ymin_digit_num*1.5 < dfpp_df[p].min() else ymin_digit_num
            # (ymaxlim_value-yminlim_value)/25
            new_max = 2 * dfpp_df[p].max() + dfpp_df[p].min() - ymaxlim_value
            new_ymaxlim_value = get_lim_num_for_graph(pd.Series([new_max, new_max]))
            if yscale == 'linear':
                # plt.yticks(np.linspace(0, 2500, 11, endpoint=True), fontsize=18)
                # plt.ylim(-50, 2500)

                plt.yticks(np.linspace(0, new_ymaxlim_value, 11, endpoint=True), fontsize=18)
                plt.ylim(-10, new_ymaxlim_value)
                # plt.yticks(np.linspace(0, ymaxlim_value, 11, endpoint=True), fontsize=18)
                # plt.ylim(-10, ymaxlim_value)
                # plt.ylim(-10, ymaxlim_value)
                # plt.ylim(yminlim_value, ymaxlim_value)

            elif yscale == 'log':
                ylim_log_range = [0,]
                for exp_n in range(10):
                    add_cand = 10**exp_n
                    if add_cand < dfpp_df[p].max():
                        ylim_log_range.append(add_cand)
                    else:
                        ylim_log_range += [ymaxlim_value]
                        break

                # ylim_log_range = [0, 1, 10, 100, 1000]
                # if drug == 'Metformin':
                #     ylim_log_range += [ylim_value]
                plt.yticks(ylim_log_range, fontsize=18)
                plt.ylim(yminlim_value, ylim_log_range[-1])

            # ds = dfpp_df[p].copy()
            # get_lim_num_for_graph(ds)
            # if drug == 'Metformin':
            #     if yscale == 'linear':
            #         # plt.yticks(np.linspace(0, 2500, 11, endpoint=True), fontsize=18)
            #         # plt.ylim(-50, 2500)
            #         pmax_value = dfpp_df[p].max()
            #
            #         plt.yticks(np.linspace(0, 3100, 11, endpoint=True), fontsize=18)
            #         plt.ylim(-50, 2500)
            #     elif yscale == 'log':
            #         plt.yticks([0, 1, 10, 100, 1000, 3500], fontsize=18)
            #         plt.ylim(1, 3500)
            # elif drug == 'Empagliflozin':
            #     if yscale == 'linear':
            #         # plt.yticks(np.linspace(0, 650, 11, endpoint=True), fontsize=18)
            #         # plt.ylim(-10,650)
            #         plt.yticks(np.linspace(0, 550, 11, endpoint=True), fontsize=18)
            #         plt.ylim(-10, 650)
            #     elif yscale == 'log':
            #         plt.yticks([0, 1, 10, 100, 1000], fontsize=18)
            #         plt.ylim(1, 1000)
            # elif drug == 'Sitagliptin':
            #     if yscale == 'linear':
            #         plt.yticks(np.linspace(0, 650, 11, endpoint=True), fontsize=18)
            #         plt.ylim(-10, 650)
            #     elif yscale == 'log':
            #         plt.yticks([0, 1, 10, 100, 1000], fontsize=18)
            #         plt.ylim(1, 1000)
            # elif drug == 'Lobeglitazone':
            #     if yscale == 'linear':
            #         plt.yticks(np.linspace(0, 60, 11, endpoint=True), fontsize=18)
            #         plt.ylim(-10, 65)
            #     elif yscale == 'log':
            #         plt.yticks([0, 1, 10, 100, 1000], fontsize=18)
            #         plt.ylim(1, 1000)


            # Add titles and labels
            # plt.title("Cmax by FEEDING status for different IDs")
            # plt.xlabel("FEEDING", fontsize=20, labelpad=8)
            plt.xlabel(" ", fontsize=20, labelpad=8)
            plt.ylabel(f"{drug} {ylabel_dict[p]}", fontsize=20, labelpad=8)

            dpi=300
            filename = f'Individual_Plot_({drug}_{p})({yscale})'
            file_format='png'
            if save_fig:
                if not os.path.exists(f"{result_file_dir_path}"): os.mkdir(f"{result_file_dir_path}")
                # if not os.path.exists(f"{result_file_dir_path}/{yscale}"): os.mkdir(f"{result_file_dir_path}/{yscale}")
                # if not os.path.exists(f"{result_file_dir_path}/{yscale}/{drug}"): os.mkdir(f"{result_file_dir_path}/{yscale}/{mode}/{drug}")
                plt.savefig(f"{result_file_dir_path}/{filename}.{file_format}", dpi=dpi)

            plt.cla()
            plt.clf()
            plt.close()

