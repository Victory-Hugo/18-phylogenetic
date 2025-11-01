      seqfile = /mnt/c/Users/Administrator/Desktop/BIG.fasta
     treefile = /mnt/c/Users/Administrator/Desktop/BIG.treefile

      outfile = /mnt/f/OneDrive/文档（科研）/脚本/Download/12-PAML/output/mlb       * main result file
        noisy = 3   * 屏幕输出冗余度（0=安静, 3=详细）
      verbose = 0   * 输出结果详细程度（0=简洁, 1=详细）

        model = 4   * 4: HKY85核苷酸替换模型 
                    * 0:JC69, 1:K80, 2:F81, 3:F84, 4:HKY85
                    * 5:T92, 6:TN93, 7:REV, 8:UNREST, 9:REVu; 10:UNRESTu
        Mgene = 0   * 不区分基因，整个对齐用同一模型
                    * 0:rates, 1:separate; 2:diff pi, 3:diff kapa, 4:all diff

        clock = 1   * 1: 全局分子钟模式 * 0:no clock, 1:clock; 2:local clock; 3:CombinedAnalysis
    fix_kappa = 0   * 0: estimate kappa; 1: fix kappa at value below
        kappa = 2.0 * initial or fixed kappa
* fix_blength = 2

    fix_alpha = 1   * 0: estimate alpha; 1: fix alpha at value below
        alpha = 0 * initial or fixed alpha, 0:infinity (constant rate)
       Malpha = 0   * 1: different alpha's for genes, 0: one alpha
        ncatG = 4   * # of categories in the dG, AdG, or nparK models of rates

         getSE = 2   * 0: don't want them, 1: want S.E.s of estimates
* RateAncestor = 0   * (0,1,2): rates (alpha>0) or ancestral states
        method = 0   * 0: simultaneous; 1: one branch at a time
    Small_Diff = 1e-6 * “变化很小，可以结束”的标准
     cleandata = 1   * remove sites with ambiguity data (1:yes, 0:no)?
