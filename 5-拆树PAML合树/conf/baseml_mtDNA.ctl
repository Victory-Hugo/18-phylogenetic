	seqfile	  =	.fasta    * [输入序列文件，此处为模板，请勿修改]
	treefile  =	.treefile * [输入树文件，此处为模板，请勿修改]

	outfile	  =	mlb		  * [主结果文件，此处为模板，请勿修改]
	noisy	  =	3			* 屏幕输出冗余度（0=安静, 3=详细）
	verbose	= 0   		  * 输出结果详细程度（0=简洁, 1=详细）
	runmode = 0         * 0: 使用 treefile 中给定拓扑

	model 	= 7   		    * 7: REV（通常对应 GTR）核苷酸替换模型
									* 0:JC69, 1:K80, 2:F81, 3:F84, 4:HKY85
									* 5:T92, 6:TN93, 7:REV, 8:UNREST, 9:REVu; 10:UNRESTu
	Mgene  = 0   			* 不区分基因，整个对齐用同一模型
                    				* 0:共享速率, 1:分开估计; 2:碱基频率不同, 3:kappa不同, 4:全部不同

	clock 		    = 1  		* 1: 全局分子钟模式
	fix_kappa    = 0   		* 0:估计 kappa; 1:将 kappa 固定为下面的值
	kappa 	   	  = 2.0 	* kappa 的初始值或固定值
	fix_blength = 0  	  * 0:估计枝长; 1:固定为初始值; 2:固定为树文件中的值

	fix_alpha  = 0   		* 0:估计 alpha; 1:将 alpha 固定为下面的值
	alpha 		= 0.5 		* alpha 的初始值或固定值，0 表示无穷大（恒定速率）
	Malpha	  = 0  		   * 1:不同基因使用不同 alpha，0:共用一个 alpha
	ncatG 		= 4   		 * dG、AdG 或 nparK 速率模型中的类别数

	getSE = 0   			   * 不应把getSE当成树支长置信区间来用;0:不输出标准误, 1:输出参数估计的标准误
	* RateAncestor = 0 	 * (0,1,2): 速率（alpha>0）或祖先状态
	method = 0  			* 0:同时优化; 1:每次优化一条分支
	Small_Diff = 1e-6 	  * “变化很小，可以结束”的标准
	cleandata = 1   		* 删除含模糊字符的数据位点吗（1:是, 0:否）？
