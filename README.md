# GATK_pre_pipeline
数据的预处理流程
## 最近做的数据处理，从sra数据得到mutation calling 以及indels，的过程，同样包括HLA-typing
得到的突变以及indels的信息在VCF文件中。从最初的测序数据sra得到VCF文件。
> 流程参考的GATK官网推荐的流程(https://software.broadinstitute.org/gatk/best-practices/workflow?id=11165),计算所用平台是上海科技大学计算平台。PBS脚本皆使用@ShixiangWang师兄写的工具批量创建提交(https://github.com/ShixiangWang/sync-deploy).

每一个样本生成一个PBS文件。

### 环境创建
首先我们使用的是conda管理环境，创建一个名为wes的环境.
### 1.sra数据处理
  我们得到的数据是加密过的sra数据，在经过解密之后，需要解压成为fastq文件才能进行以后的操作.
  srar文件批量解压为faastq文件：
  用到的工具是**sratools**中的**fastq-dump**，安装采用的是conda.具体命令：
> fastq-dump -outdir /path/to/dir --split-3 --skip-technical --clip --gzip /path/to/file/* .sra

fastq-dump的具体参数含义参照(https://www.jianshu.com/p/43680bdd42ae)

