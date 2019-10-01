# GATK_pre_pipeline
数据的预处理流程

相关流程加本都在`demo`中
## 最近做的数据处理，从sra数据得到mutation calling 以及indels，的过程，同样包括HLA-typing
得到的突变以及indels的信息在VCF文件中。从最初的测序数据sra得到VCF文件。
> 流程参考的[GATK官](https://software.broadinstitute.org/gatk/best-practices/workflow?id=11165)网推荐的流程
  计算所用平台是上海科技大学计算平台。PBS脚本皆使用[@ShixiangWang](https://github.com/ShixiangWang/sync-deploy)师兄写的工具批量创建提交每一个样本生成一个PBS文件。


![流程图](https://github.com/limbo1996/GATK_pre_pipeline/blob/master/IMG/%E5%BE%AE%E4%BF%A1%E5%9B%BE%E7%89%87_20190716222447.png)
## 环境创建
首先我们使用的是conda管理环境，创建一个名为wes的环境.
## 1.sra数据处理
  我们得到的数据是加密过的sra数据，在经过解密之后，需要解压成为fastq文件才能进行以后的操作.
  sra文件批量解压为faastq文件：
  用到的工具是**sratools**中的**fastq-dump**，安装采用的是conda.具体命令：
> fastq-dump -outdir /path/to/dir --split-3 --skip-technical --clip --gzip /path/to/file/* .sra

这样输出的是`fastq.gz`文件,是fastq文件压缩过的，占用空间小。
fastq-dump的具体参数含义参照(https://www.jianshu.com/p/43680bdd42ae)

## 2.QC
这一步的意义是将价压出来的文件做一个质量控制，可以方便下一步去接头以及其他处理.
所用工具：fastqc
经conda安装后，具体命令：
> fastqc -o /path/to/dir -t 8 /path/to/work/<sample>.fastq.gz
  
之后可以用`multiqc`工具将QC结果整合成一个网页，方便总体查看，通过查看结果中的每条序列质量等信息，对接下来的操作期参考作用。
## 3.去接头
测序的时候会加上接头，在这里要去除.
工具：[Trim galore](https://www.jianshu.com/p/7a3de6b8e503),
Trim galore 适用于所有的高通量测序，在这里我们的数据使用的是ILLUMINA,主要功能：
* 去除低质量的碱基，然后取出3'末端的接头，如果没有指定具体接头信息，程序会自动检测匹配。
具体参数含义参照前文网页。

## 4.bwa比对
  这一步的意义是将测序数据与参考基因组mapping，因为我们测序结果是一个一个的reads,所以需要组装。因为是重复的文章的数据，所以参考基因组用的和文章相同即hg19。
  用到的工具是：BWA，安装依旧使用的conda
  [BWA流程](http://starsyi.github.io/2016/05/24/BWA-%E5%91%BD%E4%BB%A4%E8%AF%A6%E8%A7%A3/)：
  BWA的比对有三个不同的算法：
  * BWA-backtrack：是用来比对ILLUMINA的序列，reads长度能达到100bp
  * BWA-SW：用于比对long-reads的，支持的长度为70bp-1Mbp；同时支持剪切性比对
  * BWA-MEM：是更新的算法，也是这次我们要用的算法，支持较长的reads同时也支持剪接性比对。
 
 这一步十分重要，因为你比对的结果就是接下来分析的输入数据，参考基因组要选好，以后每个需要参考基因组的步骤都和这是选择的参考基因组相同。
 
 人类的参考基因组有好几个版本，这次选择的hg19， 之后推荐使用hg38
 hg19也有不同的版本，具体hg9与b37或者GRCH37之间的区别可以参照[参考基因组的区别](http://www.zxzyl.com/archives/829)
 
  ### 4.1构建索引
  进行比对之前，需要对fastq文件构建索引。算法选用的`bwtsw`
  ~~~
  index   Usage：bwa index [ –p prefix ] [ –a algoType ] <in.db.fasta>
        Index database sequence in the FASTA format.
        OPTIONS:
        -P STR  输出数据库的前缀；【默认和输入的文件名一致，输出的数据库在其输入文件所在的文件夹，并以该文件名为前缀。】
        -a [is|bwtsw]   构建index的算法，有两个算法：
                        is  是默认的算法，虽然相对较快，但是需要较大的内存，当构建的数据库大于2GB的时候就不能正常工作了。
                        bwtsw   对于短的参考序列式不工作的，必须要大于等于10MB, 但能用于较大的基因组数据，比如人的全基因组。
~~~

构建方式：
> bwa index -a bwtsw /path/to/ref.fa
  ### 4.2 mem比对
  使用方法：
  > mem Usage: bwa mem [options] ref.fa reads.fq [mates.fq]

具体参数含义：
  ~~~
1 -t INT  线程数，默认是1。
2 -M 将 shorter split hits 标记为次优，以兼容 Picard’s markDuplicates 软件。
3 -p 若无此参数：输入文件只有1个，则进行单端比对；若输入文件有2个，则作为paired reads进行比对。若加入此参数：则仅以第1个文件作为输入(输入的文件若有2个，则忽略之)，该文件必须是read1.fq和read2.fa进行reads交叉的数据。
4 -R STR  完整的read group的头部，可以用 '\t' 作为分隔符， 在输出的SAM文件中被解释为制表符TAB. read group 的ID，会被添加到输出文件的每一个read的头部。
5 -T INT  当比对的分值比 INT 小时，不输出该比对结果，这个参数只影响输出的结果，不影响比对的过程。
6 -a      将所有的比对结果都输出，包括 single-end 和 unpaired paired-end的 reads，但是这些比对的结果会被标记为次优。
~~~

比对之后就得到了sam文件其实这里可以直接将sam文件转换为bam文件，但是因为是第一次做就没有连起来而是分为两步进行了。

## 5.将sam文件转化成bam文件
工具`samtools`，conda安装
[samtools](http://www.chenlianfu.com/?p=1399)
~~~

samtools view -bS /path/to/sample/<sample>.sam> /path/to/target <sample>.bam
~~~
## 6.对bam文件排序(sort_bam)
对bam文件按照染色体对应的条目按照坐标顺序从小到大排序，之后GATK流程不排序会报错。
工具：`picard`中的`sortsam`。`picard`已经被集成在`GATK`中了。所以conda下载安装GATK4就好。
这次是单独安装的`picard`。
脚本：
~~~
gatk --java-options "-Xmx20G -Djava.io.tmpdir=/public/home/wangshx/wx/tmp" SortSam \
	INPUT=/path/to/sam/<sample>.bam \
	OUTPUT=/path/to/bam/<sample>.sort.bam \
	SORT_ORDER=coordinate
~~~
**注意在sort bam的时候会产生一个很大的临时文件，可以自定义临时文件存放位置.**
## 7.Mark Duplicates
标记PCR重复，在之前的质量控制中会看到有很多重复序列，也就是在测序是经过PCR产生的重复，需要将其标记出来(不用去除)后续GATK会识别这些标记。
工具：`picard`中的`MarkDuplicates`
~~~
gatk --java-options "-Xmx20G -Djava.io.tmpdir=/path/to/tmp" MarkDuplicates \
        I=/path/to/sam/<sample>.sort.bam \
        O=/path/tp/work/<sample>.rmdup.bam \
        VALIDATION_STRINGENCY=LENIENT \
        MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 \
        M=/path/to/sort_martics/<sample>.sort.addhead.rmdup.metric

~~~
标记之后生成索引`.bai`文件。
~~~
samtools index /path/to/bam/<sample>.rmdup.bam /path/to/index/<sample>.bam.bai
~~~
## 8.BQSR(Base Quality Score Recalibration)
工具：`GATK`
这一步是对bam文件里reads的碱基质量值进行重新校正，使最后输出的bam文件中reads中碱基的质量值能够更加接近真实的与参考基因组之间错配的概率。这一步适用于多种数据类型，包括illunima、solid、454、CG等数据格式。在GATK2.0以上版本中还可以对indel的质量值进行校正，这一步对indel calling非常有帮助。
具体信息[GATK](http://starsyi.github.io/2016/05/25/%E5%8F%98%E5%BC%82%E6%A3%80%E6%B5%8B%EF%BC%88BWA-SAMtools-picard-GATK%EF%BC%89/)
分为两步：
* BaseRecalibrator
* Applybam
### 8.1 BaseRecalibrator
~~~
ref=/path/to/hg19
bam=/path/to/bam
workdir=/path/to/work

dbsnp=/path/to/dbsnp_148.hg19.vcf.gz 

dbsnp1000G=/path/to/1000G_phase1.snps.high_confidence.hg19.vcf.gz

dbindel100G=/path/to/Mills_and_1000G_gold_standard.indels.hg19.vcf.gz
 gatk --java-options "-Xmx20G -Djava.io.tmpdir=/path/to/tmp"  BaseRecalibrator \
    -R $ref/Homo_sapiens_assembly38.fasta \
    -I $bam/<head>.sort.marked.bam \
    -O $workdir/<head>.recal_data.table \
    --known-sites $dbsnp --known-sites $dbsnp1000G --known-sites $dbindel100G
~~~
注意这里的`konwn-sites`文件的选择十分重要
### 8.2 Applybam
~~~
ref=/path/to/hg19
bam=/path/to/bam
workdir=/path/to/work

dbsnp=/path/to/dbsnp_148.hg19.vcf.gz 

dbsnp1000G=/path/to/1000G_phase1.snps.high_confidence.hg19.vcf.gz

dbindel100G=/path/to/Mills_and_1000G_gold_standard.indels.hg19.vcf.gz

 gatk --java-options "-Xmx20G -Djava.io.tmpdir=/path/to/tmp"  ApplyBQSR  \
    -R $ref/Homo_sapiens_assembly38.fasta \
    -I $bam/<head>.sort.marked.bam \
    --bqsr-recal-file $workdir/<head>.recal_data.table \
    -O $workdir/<head>.sort.marked.BQSR.bam
~~~
