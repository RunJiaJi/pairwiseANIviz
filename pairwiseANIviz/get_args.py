import argparse

def get_args():
    parser = argparse.ArgumentParser(
        prog="pairwiseANIviz",
        description="Pairwise ANI(Average Nucleotide Identity) visulization tool.",
        usage='%(prog)s [options]',
        formatter_class=argparse.RawTextHelpFormatter,
        epilog=
        "General usage\n----------------\n"
        "ANI result visulization without classification info:\n"
        "   $ pairwiseANIviz ani_result.txt\n\n"
        "ANI result visulization with classification info:\n"
        "   $ pairwiseANIviz ani_result.txt --classificationFile classification_result.tsv\n\n"
        "Runjia Ji, 2023"
    )

    parser.add_argument('-v', '--version', action='version', version='%(prog)s 1.0', help='Show pairwiseANIviz version number and exit.')
    parser.add_argument('anifile', type=str, help='File containing pairwise ANI analysis result.')
    parser.add_argument('-o', '--outdir', type=str, default='./pairwiseANIviz', help="Directory to save the output figures (default 'pairwiseANIviz').")
    parser.add_argument('--method', type=str, choices=['single', 'complete', 'average', 'weighted', 'centroid', 'median', 'ward'], 
                        default='average',
                        help="Linkage method to use for calculating clusters (default 'average').\nSee https://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.linkage.html#scipy.cluster.hierarchy.linkage")
    parser.add_argument('--metric', type=str, 
                        choices=[
                        'braycurtis', 'canberra', 'chebyshev', 'cityblock', 'correlation', 'cosine', 'dice', 'euclidean', 'hamming', 'jaccard', 'jensenshannon', 'kulczynski1', 'mahalanobis', 'matching', 'minkowski', 'rogerstanimoto', 'russellrao', 'seuclidean', 'sokalmichener', 'sokalsneath', 'sqeuclidean', 'yule'
                        ],
                        default='euclidean',
                        help="The distance metric to use (default 'euclidean').\nSee https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.distance.pdist.html#scipy.spatial.distance.pdist"
                    )
    parser.add_argument('-cmap', '--colormap', type=str, 
                        default='Blues',
                        help="Matplotlib colormap used when drawing the heatmap of ANI values (default 'Blues').\nSee https://matplotlib.org/stable/users/explain/colors/colormaps.html")
    parser.add_argument('--figWidth', type=int, default=15, help="Figure width (default '15').")
    parser.add_argument('--figHeight', type=int, default=15, help="Figure height (default '15').")
    parser.add_argument('--linewidth', type=float, default=0.5, help="Line width of the main heatmap (default 0.5)")
    parser.add_argument('--linecolor', type=str, default='grey', help="Line color of the main heatmap (default 'grey').")
    parser.add_argument('--rowCluster', help="Draw the row cluster.", action='store_true')
    parser.add_argument('--colCluster', help="Draw the column cluster.", action='store_true')
    parser.add_argument('--annotation', help="Show ANI values on the plot.", action='store_true')
    parser.add_argument('--illustrateOverValue',type=float, 
                        default=100,
                        help="Cells have ANI values over specific threshold set to red (eg. cells have ANI value >=0.95 set to red).")
    

    parser.add_argument('-c', '--classificationFile', type=str, 
                        help='File containing classification result generated by GTDBTk(https://github.com/Ecogenomics/GTDBTk).')
    parser.add_argument('-t', '--taxaLevel', type=str, 
                        default='phylum',
                        choices=['domain', 'phylum', 'class', 'order', 'family', 'genus', 'species'], help='Taxa level illustrated on the plot.\nNote that this parameter only works if classification result was input.')
    parser.add_argument('--colorPalette', type=str, 
                        default='hls',
                        help="Color palette used to return a specified number of evenly spaced hues which are then used to illustrate different taxa (default 'hls').\n Note that this parameter only works if classification result was input.")
    

    args = parser.parse_args()
    return args
