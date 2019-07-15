# GATK_pre_pipeline
数据的预处理流程
## 最近做的数据处理，从sra数据得到mutation calling 以及indels，的过程，同样包括HLA-typing
得到的突变以及indels的信息在VCF文件中。从最初的测序数据sra得到VCF文件。
> 流程参考的[GATK官](https://software.broadinstitute.org/gatk/best-practices/workflow?id=11165)网推荐的流程
  计算所用平台是上海科技大学计算平台。PBS脚本皆使用[@ShixiangWang](https://github.com/ShixiangWang/sync-deploy)师兄写的工具批量创建提交

每一个样本生成一个PBS文件。

### 环境创建
首先我们使用的是conda管理环境，创建一个名为wes的环境.
### 1.sra数据处理
  我们得到的数据是加密过的sra数据，在经过解密之后，需要解压成为fastq文件才能进行以后的操作.
  sra文件批量解压为faastq文件：
  用到的工具是**sratools**中的**fastq-dump**，安装采用的是conda.具体命令：
> fastq-dump -outdir /path/to/dir --split-3 --skip-technical --clip --gzip /path/to/file/* .sra

fastq-dump的具体参数含义参照(https://www.jianshu.com/p/43680bdd42ae)
#### 1.1 sra解压脚本
~~~
#PBS -N fastq_<sample>
#PBS -l nodes=1:ppn=1
#PBS -l walltime=20:00:00
#PBS -l mem=10gb
#PBS -S /bin/bash
#PBS -q normal_3
#PBS -j oe

source activate wes
workdir=/path/to/work
cd $workdir
fastq-dump -outdir fastq_gz --split-3 --skip-technical --clip --gzip $workdir/sra/<sample>
~~~
这样输出的是`fastq.gz`文件,是fastq文件压缩过的，占用空间小。
### 2.QC
这一步的意义是将价压出来的文件做一个质量控制，可以方便下一步去接头以及其他处理.
所用工具：fastqc
经conda安装后，具体命令：
> fastqc -o /path/to/dir -t 8 /path/to/work/<sample>.fastq.gz
  
之后可以用`multiqc`工具将QC结果整合成一个网页，方便总体查看。
### 3.去接头
测序的时候会加上接头，在这里要去除.
工具：[Trim galore](https://www.jianshu.com/p/7a3de6b8e503),
Trim galore 适用于所有的高通量测序，在这里我们的数据使用的是ILLUMINA,主要功能：
* 去除低质量的碱基，然后取出3'末端的接头，如果没有指定具体接头信息，程序会自动检测匹配。
具体脚本
~~~
#PBS -N <sample>_trim
#PBS -l nodes=1:ppn=1
#PBS -l walltime=15:00:00
#PBS -S /bin/bash
#PBS -q normal_3
#PBS -j oe
source activate wes
workdir=/path/to/work
trim_galore -q 20 --phred33 --stringency 3 --length 20 -e 0.1 --paired $workdir/fastq_gz/<sample>_1.fastq.gz $workdir \
/fastq_gz/<sample>_2.fastq.gz -o $workdir/fastq_clean/
~~~
具体参数含义参照前文网页。

### 4.bwa




