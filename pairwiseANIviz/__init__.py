from .get_args import get_args
from .pairwiseANI import readIput, plot

def main():
    args = get_args()
    if args.classificationFile:
        anidf, classidf = readIput(args.anifile, args.classificationFile)
    else:
        anidf = readIput(args.anifile)
        classidf = None
    plot(anidf, 
        outdir=args.outdir, 
        method=args.method, 
        metric=args.metric, 
        cmap=args.colormap, 
        figsize=(args.figWidth,args.figHeight),
        linewidths=args.linewidth, 
        linecolor=args.linecolor, 
        row_cluster=args.rowCluster,
        col_cluster=args.colCluster,
        annotation=args.annotation,
        vmax=args.outrangeValue,
        classidf=classidf,
        taxaLevel=args.taxaLevel,
        palette=args.colorPalette,
        )

#############################################
if __name__ == '__main__':
    main()