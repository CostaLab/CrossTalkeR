## ----vignete0-----------------------------------------------------------------
suppressPackageStartupMessages({require(CrossTalkeR)})
suppressPackageStartupMessages({require(igraph)})
suppressPackageStartupMessages({require(ggraph)})
suppressPackageStartupMessages({require(ggplot2)})


## ----vignete1, results="hide", warning=FALSE,message=F, output=FALSE----------
paths <- c('CTR' = system.file("extdata",
                               "ctr_nils_bm_human.csv",
                               package = "CrossTalkeR"),
           'EXP' = system.file("extdata",
                               "exp_nils_bm_human.csv",
                               package = "CrossTalkeR"))
genes <- c('TGFB1','PF4','PPBP')
data <- generate_report(paths,
                        genes,
                        out_path='~/Documents/',
                        threshold=0,
                        out_file = 'vignettes_example.html',
                        output_fmt = "html_document",
                        report = FALSE)

## ----vignete2, results="hide", warning=FALSE,message=F, output=FALSE,fig.width=8,fig.height=8----
plot_cci(graph = data@graphs$CTR,
        colors = data@colors,
        plt_name = 'Example 1',
        coords = data@coords[V(data@graphs$CTR)$name,],
        emax = NULL,
        leg = FALSE,
        low = 0,
        high = 0,
        ignore_alpha = FALSE,
        log = FALSE,
        efactor = 8,
        vfactor = 12)


## ----vignete3, results="hide", warning=FALSE,message=F, output=FALSE,fig.width=15,fig.height=15----
plot_ggi(graph = data@graphs_ggi$EXP_x_CTR,
         color = data@colors,name="EXP_x_CTR")


## ----vignete4, results="hide", warning=FALSE,message=F, output=FALSE,fig.width=15,fig.height=10----
plot_sankey(lrobj_tbl = data@tables$EXP_x_CTR,
            target = c("TGFB1"),
            ligand_cluster = NULL,
            receptor_cluster = NULL,
            plt_name = "TGFB1")


## ----sessionInfo, echo=FALSE--------------------------------------------------
sessionInfo()

