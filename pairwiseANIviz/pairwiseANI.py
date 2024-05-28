# the classification infor is optional
# the ANI and GTDBTK analysis is optional


import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import seaborn as sns
import scipy.cluster.hierarchy as hc
import os
import sys   

# read the ani result and classification result
def readIput(anifile, classifile=None):
    anidf = pd.read_csv(anifile, sep='\t', header=None).rename(columns={0:'q',1:'r',2:'ANI'}).iloc[:,:3]
    anidf = anidf.pivot(index='q', columns='r', values='ANI')

    # test if the numbers of cols and rows are the same
    cols = set(anidf.columns.tolist())
    rows = set(anidf.index.tolist())
    if cols==rows:
        print('ANI result is read and checked.')
    else:
        print('The query and reference genomes are different. Please check the ANI result or rerun the ANI analysis.')
        sys.exit()

    if classifile:
        classidf = pd.read_csv(classifile, sep='\t')
        # test if all genomes have the classification info
        classi_gs = set(classidf.iloc[:,0].tolist())
        if cols.issubset(classi_gs):
            print('Classification result is read and checked.')
        else:
            absent_gs = list(cols-classi_gs)
            print(f'Genomes {absent_gs} used in ANI analysis have no classification information in {classifile}.')
            print('Please check the classification result or rerun the classification analysis')
            sys.exit()

        taxaLevels = 'domain phylum class order family genus species'.split()
        for idx in range(7):
            classidf[taxaLevels[idx]] = classidf.iloc[:,1].apply(lambda x: x.split(';')[idx])

        for i in taxaLevels:
            classidf[i] = classidf[i].apply(lambda x: f'Unknown_{i}' if not x.split('__')[1] else x)
            dftmp = classidf[classidf[i].apply(lambda x: True if 'Unknown' in x else False)].copy()
            unknownNums = len(dftmp)
            dftmp['idx'] = [idx for idx in range(unknownNums)]
            dftmp['tmp'] = dftmp[i]+dftmp['idx'].astype(str)
            classidf.loc[dftmp.index.tolist(),i] = dftmp['tmp']
        return anidf, classidf
    return anidf
    

# define row colors
def rowColors(classidf, cols, taxaLevel='phylum', palette='hls'):
    classidf.set_index(classidf.columns.tolist()[0], inplace=True)
    classidf = classidf.loc[list(cols),:].copy()
    uniqueTaxa = classidf[taxaLevel].unique().tolist()
    colorNums = len(uniqueTaxa)
    rgb_list = [i for i in sns.color_palette(palette,n_colors=colorNums)]
    # print(colorNums)
    # print(uniqueTaxa)
    lut = dict(zip(classidf[taxaLevel].unique(), rgb_list))
    row_colors_rgb = classidf[taxaLevel].map(lut)

    return (lut, row_colors_rgb)

# draw clustermap
def plot(anidf, 
         outdir='./pairwiseANIviz', 
         method="average", 
         metric="euclidean", 
         cmap='Blues', 
         figsize=(15,15),
         linewidths=0.5, 
         linecolor='grey', 
         row_cluster=False,
         col_cluster=False,
         annotation=False,
         vmax=100,
         classidf=None,
         taxaLevel='phylum',
         palette='hls',
         ):

    plt.rcParams['svg.fonttype'] = 'none'
    plt.rcParams["figure.dpi"] = 300
    # clustering
    linkage = hc.linkage(anidf.fillna(0), method=method, metric=metric)

    # custome cmap
    custom_cmap=plt.cm.__dict__[cmap]
    custom_cmap.set_under(color='lightgrey')
    custom_cmap.set_over(color='red')

    # see if there is classification info, calculate row colors
    if classidf is not None:
        cols=set(anidf.columns.tolist())
        lut_colors = rowColors(classidf, cols, taxaLevel=taxaLevel, palette=palette)
        lut = lut_colors[0]
        handles = [Patch(facecolor=lut[name]) for name in lut]

        ax=sns.clustermap(
            anidf.fillna(0), 
            cmap=custom_cmap, 
            vmin=round(anidf.describe().loc['min',:].describe()['min']), 
            vmax=vmax, 
            linewidths=linewidths, 
            linecolor=linecolor, 
            figsize=figsize,
            row_linkage=linkage, 
            col_linkage=linkage, 
            row_cluster=row_cluster,
            col_cluster=col_cluster,
            annot=annotation,
            row_colors=lut_colors[1],
            xticklabels = False,
            yticklabels = True, 
            cbar_kws={
                "label": "ANI (%)",
                "orientation": "vertical",
                "spacing": "proportional"
            },
            tree_kws={"linewidths": 1.5},
            fmt=".3g"
            )
        
        plt.legend(handles, lut, title=taxaLevel,
                bbox_to_anchor=[0.9, 0.9], bbox_transform=plt.gcf().transFigure, loc='upper right')
    else:
        ax=sns.clustermap(
            anidf.fillna(0), 
            cmap=custom_cmap, 
            vmin=anidf.describe().loc['min',:].describe()['min'], 
            vmax=vmax, 
            linewidths=linewidths, 
            linecolor=linecolor, 
            figsize=figsize,
            row_linkage=linkage, 
            col_linkage=linkage, 
            row_cluster=row_cluster,
            col_cluster=col_cluster,
            annot=annotation,
            cbar = True,
            xticklabels = False,
            yticklabels = True, 
            cbar_kws={
                "label": "ANI (%)",
                "orientation": "vertical",
                "spacing": "proportional"
            },
            tree_kws={"linewidths": 1.5},
            fmt=".3g"
            )
    
    # save figs
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    
    suffix=['svg','png','jpg','pdf','tiff', 'eps']
    outfiles=['pairwiseANIviz.'+i for i in suffix]
    for file in outfiles:
        figure=os.path.join(outdir, file)
        ax.figure.savefig(figure, bbox_inches='tight')
