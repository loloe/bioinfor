
# 1.安装下载

	从bioconda源下载，对选项yes，并且创建环境jcvi，下载的软件也是jcvi
==这里发现不创建新环境会报错，所以选择新环境==
	
```
conda create -y -c bioconda -n jcvi jcvi
```

	使用时需要激活环境
```bash
conda activate jcvi
```
# 2.数据准备

ncbi或ensembl数据库下载基因组fa文件、cds文件和gff文件

```
拟南芥  https://phytozome-next.jgi.doe.gov/info/Athaliana_Araport11
二月兰  http://www.bioinformaticslab.cn/pubs/OV_data/
白菜 Brassica napus   https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/020/379/485/GCA_020379485.1_Da-Ae/
```

## 2.1 基因组fa文件

## 2.2  cds文件

cds文件需要进行以下调整

1. 提取最长转录本

- 批量提取运行sh文件，sh文件如下，我们输入的是cds.fa文件

```
mkdir -p output  #创建文件夹
for i in *.fa ; do  #遍历以.fa结尾的文件
filename="${i%.*}"   ## 获取文件名（不包括扩展名）
python /public/home/bsun/miniconda3/envs/OF/bin/primary_transcript.py "$i" > "output/${filename}_.fa"  # 运行脚本并将输出保存到新文件夹
```

- ==报错信息1==

![[Pasted image 20240524142108.png]]
	解决方案：把windows换行符改成linux换行符，如下命令就行
	
```
dos2unix longest.sh
```
	这样得到的是最长转录本的cds信息，后续分析


	

脚本如下

```
import os
import re
import sys
from collections import Counter, defaultdict

# Use the 'all' version rather than ab initio

def CheckFile(fn):
    """
    Checks for:
    - Duplicated accession lines
    """
    accs = set()
    with open(fn, 'r') as infile:
        for l in infile:
            if l.startswith(">"):
                a = l.rstrip()[1:]
                if a in accs:
                    print("\nERROR: duplicated sequence accession:\n%s" % a)
                    print("\nPlease correct this and then rerun the script.\n")
                    return False
                accs.add(a)
    return True

def ScanTags(fn):
    """
    For ensembl genomes, look for tag:id and count repeated ids
    :param fn:
    :return:
    """
    tags = set()
    tokens = []
    with open(fn, 'r') as infile:
        for line in infile:
            if not line.startswith(">"): continue
            tokens.append([t.split(":", 1) for t in line.rstrip().split() if ":" in t])
            tags.update([t[0] for t in tokens[-1]])
    for this_tag in tags:
        print(this_tag)
        # print(tokens[-1])
        c = Counter([idd for acc in tokens for t, idd in acc if t == this_tag])
        print(c.most_common(5))
        print("")

def ScanTags_NCBI(fn):
    genes = []
    with open(fn, 'r') as infile:
        for line in infile:
            if not line.startswith(">"): continue
            genes.append(line[1:].split(".", 1)[0])
    print("%d sequences, %d genes" % (len(genes), len(set(genes))))

def ScanTags_with_fn(fn, gene_name_fn):
    genes = []
    with open(fn, 'r') as infile:
        for line in infile:
            if not line.startswith(">"): continue
            genes.append(gene_name_fn(line))
    print("%d sequences, %d genes" % (len(genes), len(set(genes))))
    # print(genes[0])
    # print(sorted(genes)[:10])

def GetGeneName_Ensembl(acc_line):
    tokens = [(t.split("=") if "=" in t else t.split(":"))[1] for t in acc_line.rstrip().split() if ("gene:" in t or "gene=" in t or "locus:" in t or "locus=" in t)]
    if len(tokens) != 1: return None
    return tokens[0]

def IsNCBI(fn):
    with open(fn, 'r') as infile:
        for l in infile:
            if l.startswith(">"):
                l = l.rstrip()
                if l.startswith(">NP_") and l.endswith("]"): return True
                elif l.startswith(">XP_") and l.endswith("]"): return True
                elif l.startswith(">YP_") and l.endswith("]"): return True
                elif l.startswith(">WP_") and l.endswith("]"): return True
                return False
    return False

def GetGeneName_NCBI(acc_line):
    acc_line = acc_line[1:]
    original = acc_line
    # look for "isoform X[:d]+" or "isoform [:d]+"
    acc_line = re.sub("isoform [0-9, A-Z]+ ", "", acc_line)
    acc_line = re.sub("isoform X[0-9, A-Z]+ ", "", acc_line)
    # This last step is nasty. These are the same gene:
    # >XP_024356342.1 pyruvate decarboxylase 2-like isoform X1 [Physcomitrella patens]
    # >XP_024356343.1 pyruvate decarboxylase 2-like isoform X1 [Physcomitrella patens]
    # as the name is the same and they same 'isoform ...'
    # Whereas these are not the same gene, even though the names are identical
    # because they don't say isoform:
    # >XP_024390255.1 40S ribosomal protein S12-like [Physcomitrella patens]
    # >XP_024399722.1 40S ribosomal protein S12-like [Physcomitrella patens]
    #
    # To deal with that, we remove the ID (e.g. XP_024356342.1) if it says 'isoform'
    # so that the lines are identical, but not when it doesn't say 'isoform'
    # so that the lines are different. If I were writting the script from scratch
    # for NCBI files I'd do it a different way, but this is a way to handle it so 
    # that it works with the existing logic in the file.
    if original != acc_line:
        acc_line = acc_line.split(None, 1)[-1]
    return acc_line

def CreatePrimaryTranscriptsFile(fn, dout, gene_name_fn, q_use_original_accession_line):
    # Get genes and lengths
    max_gene_lens = defaultdict(int)
    with open(fn, 'r') as infile:
        lines = [l.rstrip() for l in infile]
    N = len(lines) - 1
    nAcc = 0
    nGeneUnidentified = 0
    acc_to_use = defaultdict(str)
    iLine = -1
    while iLine < N:
        iLine += 1
        line = lines[iLine]
        if not line.startswith(">"): continue
        nAcc += 1
        iLineAcc = iLine
        gene = gene_name_fn(line)
        if gene == None:
            nGeneUnidentified += 1
            continue
        # get length
        l = 0
        while iLine < N:
            iLine += 1
            line = lines[iLine]
            if line.startswith(">"):
                iLine -= 1
                break
            l += len(line.rstrip())
        if l > max_gene_lens[gene]:
            max_gene_lens[gene] = l
            acc_to_use[gene] = iLineAcc
    print("Found %d accessions, %d genes, %d unidentified transcripts" % (nAcc, len(max_gene_lens), nGeneUnidentified))
    # print(gene)
    # print(sorted(max_gene_lens.keys())[:10])
    # print(len(set(max_gene_lens.keys())))

    # Get longest version for each gene
    # Parse file second time and only write out sequences that are longest variant
    nGenesWriten = 0
    outfn = dout + os.path.basename(fn)
    with open(outfn, 'w') as outfile:
        iLine = -1
        while iLine < N:
            iLine += 1
            line = lines[iLine]
            if not line.startswith(">"): continue
            gene = gene_name_fn(line)
            # transcripts not identifying the gene should be written
            if gene != None and iLine != acc_to_use[gene]: continue
            if q_use_original_accession_line or gene == None:
                acc_line_out = line + "\n"
            else:
                 acc_line_out = ">%s\n" % gene
            nGenesWriten += 1
            outfile.write(acc_line_out)
            while iLine < N:
                iLine += 1
                line = lines[iLine]
                if line.startswith(">"):
                    iLine -= 1
                    break
                outfile.write(line + "\n")
    print("Wrote %d genes" % nGenesWriten)
    if nGenesWriten != len(max_gene_lens) + nGeneUnidentified:
        print("ERROR")
        raise Exception
    print(outfn)


def last_dot(text):
    return text[1:].rstrip().rsplit(".", 1)[0]


def space(text):
    return text[1:].rstrip().split(None, 1)[0]


function_dict = {"last_dot":last_dot, "space":space}

def main(args=None):
    print("")
    if args is None:
        args = sys.argv[1:]
    fn = args[0]

    if not CheckFile(fn):
        return

    dout = os.path.dirname(os.path.abspath(fn)) + "/primary_transcripts/"
    if not os.path.exists(dout):
        os.mkdir(dout)

    if len(sys.argv) == 3:
        gene_name_function_name = function_dict[sys.argv[2]]
        ScanTags_with_fn(fn, gene_name_function_name)
        CreatePrimaryTranscriptsFile(fn, dout, gene_name_function_name)
    else:
        # ScanTags(fn)
        # ScanTags_NCBI(fn)
        # ScanTags_second_dot(fn)
        if IsNCBI(fn):
            print("Identified as NCBI file")
            gene_name_function_name = GetGeneName_NCBI
            q_use_original_accession_line = True
        else:
            # This is the default, unless we can determine that it is an NCBI
            # file in which case it needs special processing
            gene_name_function_name = GetGeneName_Ensembl
            q_use_original_accession_line = False
            print('Looking for "gene=" of "gene:" to identify isoforms of same gene')
        CreatePrimaryTranscriptsFile(fn, dout, gene_name_function_name, q_use_original_accession_line)

if __name__ == "__main__":
    args = sys.argv[1:]
    main(args)

```

- ==结果==

运行完脚本会出现两个文件夹，output文件里面是提取的基本信息，多少个基因被写入等等，primary文件夹里面就是我们需要的最长cds序列。

![[Pasted image 20240524142702.png]]



2. 根据cds和gff文件提取所需的bed文件，这部分--type=mRNA --key=ID需要两文件内容一一对应
	==这里的cds是我们上面提取过最长转录本的cds==
	本步骤结合cds文件的ID去选择type和key的类型 ，gff3文件转化bed文件时注意type和key类型对应gff中第三列和第九列信息。type一般为mRNA、cds等，但是key注意你的gff文件是取Name还是ID。

- 二月兰cds文件
我们可以观察到基因ID是带.1的

![[Pasted image 20240524144321.png]]

- 二月兰gff文件
![[Pasted image 20240524144741.png]]

==取--type=mRNA --key=ID

****
- 拟南芥最长cds
![[Pasted image 20240524143320.png]]
- 拟南芥gff
![[Pasted image 20240524143733.png]]

==取--type=gene --key=Name==
****
油菜是有亚基因组的，我们可以从基因组fa文件中看到有A和C两种
- 油菜cds
![[Pasted image 20240524145620.png]]
- 油菜gff
![[Pasted image 20240524150501.png]]

==这里我们遇到遇到gff中无论是mrna还是cds还是gene都无法找到cds文件中的基因id；我们可以自己利用gffread进行提取，然后在提取bed==

```
module load gffread/0.12.6
#用gffread提取cds序列,蛋白序列,转录本序列
gffread GCA_020379485.1_Da-Ae_genomic.gff -g GCA_020379485.1_Da-Ae_genomic.fna -x Branap.cds.fa
gffread GCA_020379485.1_Da-Ae_genomic.gff -g GCA_020379485.1_Da-Ae_genomic.fna -y Branap.protein.fa
gffread GCA_020379485.1_Da-Ae_genomic.gff -g GCA_020379485.1_Da-Ae_genomic.fna -w Branap.transcripts.fa
```

接着对油菜进行提取最长转录本，至此我们的最长转录本都已经被提取出来了，接着观察油菜的cds和gff对应关系

- 油菜cds
![[Pasted image 20240525095852.png]]

油菜gff

![[Pasted image 20240525095940.png]]

==取--type=mRNA --key=ID==

## 2.3  gff文件

就是从网上下载的gff文件



## 2.4 bed文件

### 2.4.1 提取bed文件

就是用gff转换成bed格式

1. 拟南芥是gene行的 Name和cds文件的id一致，所以我们选择gene和 Name去提取bed文件

```bash


python -m jcvi.formats.gff bed --type=gene --key=Name Athaliana_447_Araport11.gene.gff3 -o Athaliana_447_Araport11.gene.bed


#- `"python -m jcvi.formats.gff bed --type=mRNA --key=ID zr.gff -o zr.bed"`: 这是要在作业中运行的命令。它使用Python中的jcvi模块中的脚本，将输入的GFF文件转换为BED格式的文件。具体来说，`bed`是指定的jcvi格式转换函数，`--type=mRNA`指定要转换的GFF中的特征类型为mRNA，`--key=ID`指定要在BED文件中使用的键名，`zr.gff`是输入的GFF文件，`-o zr.bed`指定了输出的BED文件名为`zr.bed`。

```

2. 二月兰是--type=mRNA --key=ID

```bash

python -m jcvi.formats.gff bed --type=mRNA --key=ID OV_sca.gff -o OV_sca.bed

```


3. 油菜是type=mRNA --key=ID

```bash
python -m jcvi.formats.gff bed --type=mRNA --key=ID GCA_020379485.1_Da-Ae_genomic.gff -o GCA_020379485.1_Da-Ae_genomic.bed
```

### 2.4.2 去重bed文件


对bed文件进行去重复（jcvi自带python脚本）

```
1.拟南芥去重
python3 -m jcvi.formats.bed uniq Athaliana_447_Araport11.gene.bed
2.油菜去重
python3 -m jcvi.formats.bed uniq GCA_020379485.1_Da-Ae_genomic.bed
3.二月兰去重
python3 -m jcvi.formats.bed uniq OV_sca.bed

#是使用 Python 3 中的 `jcvi` 模块下的 `formats.bed` 子模块对 `spottgar.bed` 文件中的条目进行去重
```

==**重要的事说三遍：删除bed文件中多余的contig或scaffold对应的列，只保存染色体列即可，用大文本编辑器打开uniqbed文件自己删就行，也可以用命令（否则后续绘图报错****）**==



### 2.4.3 根据uniq.bed文件提取最终的cds文件

```bash
seqkit grep -f <(cut -f 4 mis.uniq.bed) mis.cds.fa | seqkit seq -i > mis.cds

#这个命令的作用是根据 `mis.uniq.bed` 文件中的第四列提取序列文件 `mis.cds.fa` 中对应的序列，并将结果输出到一个名为 `mis.cds` 的文件中
```


```
1.拟南芥
seqkit grep -f <(cut -f 4 Athaliana_447_Araport11.uniq.bed) primary_transcripts/Athaliana_447_rt11.cds.fa | seqkit seq -i > Athaliana_447_Araport11.cdslast.fa

2.二月兰
seqkit grep -f <(cut -f 4 OV_sca.uniq.bed) primary_transcripts/OV_sca.cds.fa | seqkit seq -i > OV_sca.cdslast.fa

3.油菜
seqkit grep -f <(cut -f 4 GCA_020379485.uniq.bed) primary_transcripts/Branap.cds.fa | seqkit seq -i > Branap.cdslast.fa
```

# 3.运行

## 3.1 **提取共线信息**


根据共线性图的展示顺序对三个物种的cds序列进行两两比对，并保留基因数量较多的block便于绘图展示，本文的展示顺序为二月兰（上）、拟南芥（中）、油菜（下），所以只需要将二月兰与拟南芥比较，拟南芥与油菜比较即可

==我们比较的是cds，所以要保证格式都是物种名+后缀==

```
python3 -m jcvi.compara.catalog ortholog --no_strip_names --cscore  0.99 species1  species2

#`python3 -m jcvi.compara.catalog ortholog`: 这部分使用 Python 3 解释器来运行 JCvi 工具包中的 `compara.catalog` 模块下的 `ortholog` 子命令。这个命令用于执行基因的同源性分析，查找两个物种之间的同源基因
#`--no_strip_names`: 这个选项告诉程序在进行同源基因分析时不要剥离基因的名称。默认情况下，程序会尝试剥离基因名称中的版本信息（例如 `_v1`、`-RA` 等）以便更好地进行匹配。使用这个选项会跳过这个步骤
#`--cscore 0.99`: 这个选项指定了一个阈值，即最小的同源基因匹配的信心得分。在这个例子中，设置为 0.99，表示只有信心得分大于等于 0.99 的同源基因匹配会被考虑
#`species1` 和 `species2`: 这两个参数指定了要进行比较的两个物种。你需要将它们替换为实际的物种名称或缩写，例如 `species1` 可以是 `human`，`species2` 可以是 `mouse`
```

创建一个文件夹，把我们需要比对的物种cds都给移动进去，统一最终要使用的cds、及bed文件的名字前缀，使JCVI能够识别。最终文件命名为：species.cds；species.bed。

![[Pasted image 20240525145834.png]]

我们最终画的是拟南芥在中间，所以只需要比两次，如果画一个三角形，那就要两两比较了，根据自己想要的结构来的

```
1.二月兰和拟南芥比较
python3 -m jcvi.compara.catalog ortholog --no_strip_names --cscore 0.99 OV_sca Athaliana
2.拟南芥和油菜比较
python3 -m jcvi.compara.catalog ortholog --no_strip_names --cscore 0.99 Athaliana  Branap


```








最后结果会得到有pdf和anchors的文件

![[Pasted image 20240529202038.png]]

![[Pasted image 20240529095047.png]]

.anchors
![[Pasted image 20240529095121.png]]


anchors文件中每个共线性区块以###分隔,第一和第二列分别是两基因组的基因ID，第三列BLAST的bit score，越大可靠性越高。该文件的内容需要用来搜索自己的基因在不在block里，如果在.new.anchors文件中；要找到目的基因在其block的开始和结束基因的id，然后到.smple中标记颜色。

-  .last  基于LAST的比对结果
- .last.filtered: LAST的比对结果过滤串联重复和低分比对
- lifted.anchors:增加了额外的锚点，形成最终的共线性区块

## 3.2 建立简单的anchor.simple文件

```
1.二月兰和拟南芥
python -m jcvi.compara.synteny screen --minspan=30 --simple OV_sca.Athaliana.anchors OV_sca.Athaliana.anchors.new

2.拟南芥和油菜
python -m jcvi.compara.synteny screen --minspan=30 --simple Athaliana.Branap.anchors Athaliana.Branap.anchors.new
#这条命令执行了Python模块`jcvi.compara.synteny`中的`screen`功能，参数包括：

#- `--minspan=30`：设置最小跨度为30。
#- `--simple`：使用简单模式。
#- `OV_sca.Athaliana.anchors ` 和 `OV_sca.Athaliana.anchors.new`：输入和输出文件。
```

**.simple: 从.anchors文件创建的更简化格式，对应的就是anchors文件里面一块一块的内容**
![[Pasted image 20240529100644.png]]

## 3.3 **配置画图文件-layout.txt和seqids.txt文件**

### 3.3.1**layout.txt**,用于设置绘制的一些选项。

第一列，二列，三列控制的是track的位置。rotation是方向，color是颜色，label是标签。va是vertical alignment。

```
vim layout
三个轨道和两个边缘（edges）的定义
# y, xstart, xend, rotation, color, label, va, bed 
.7, .2, .8, 0, , Qden, top, QD.bed 
.5, .2, .8, 0, , Qali, middle, QA.bed 
.3, .2, .8, 0, , Qvar, bottom, QV.bed 
# edges 
e, 0, 1, QA.QD.anchors.simple 
e, 1, 2, QA.QV.anchors.simple


- y位置为0.7
- 水平起始位置为0.2
- 水平结束位置为0.8
- 无旋转
- 颜色未指定
- 标签为`Qden`
- 垂直对齐方式为`top`
- 数据来自`QD.bed`文件

```

### 3.3.2 seqids.txt

==-一行代表一个物种的bed文件中的染色体名称，必须和layout文件中的物种先后一致（否则一直报错）**==


```
OV01，OV02，OV03，OV04，OV05，OV06，OV07，OV08，OV09，OV10，OV11，OV12
Chr1，Chr2，Chr3，Chr4，Chr5
CM035446.1，CM035447.1，CM035448.1，CM035449.1，CM035450.1，CM035451.1，CM035452.1，CM035453.1，CM035454.1，CM035455.1，CM035456.1，CM035457.1，CM035458.1，CM035459.1，CM035460.1，CM035461.1，CM035462.1，CM035463.1，CM035464.1
```

# 4 画多个基因组间共线性图

```
python -m jcvi.graphics.karyotype seqids layout

#- `seqids`：这是输入文件，包含序列标识符。这些标识符通常是染色体或其他大型序列的名称。
#- `layout`：这是布局文件，定义了如何在图中布置和显示这些序列。
```


```
python -m jcvi.graphics.karyotype seqids.txt layout.txt --nocircle
#--nocircle  是不加染色体编号
```

最后出来的图如下
![[Pasted image 20240529203513.png]]

如果想标注一些想要的基因，则在.new.anchors文件中找到目的基因在其block的开始和结束基因的id，然后到.smple中标记颜色（g、r、y、b*；颜色首字符加*表示）；如果找不到自己的基因说明同源性不够好，可以放宽cscore的参数尝试。（r是红色，g是绿色等）


![[Pasted image 20240529151329.png]]


然后会变成

![[Pasted image 20240529151426.png]]

# 5 利用jcvi绘制基因水平的局部共线性图

需要准备的三个输入文件

- 记录物种内或者物种间的共线性基因对
    
- 记录基因坐标的bed文件
    
- 布局文件



第一步，基于已有的共线性分析结果(WGDI, MCscan, MCscanX等软件的分析结果），整理出你需要展示的区间的基因对。注意分隔符是制表符，我们保存为blocks.txt。我们可以通过前面已整理好的文件来进行分析。

方法是两两基因组做共线性分析。还是以二月兰，拟南芥和油菜为例，方法如下:
```
python -m jcvi.compara.synteny mcscan QA.bed QA.QV.lifted.anchors --iter=1 -o QA.QV.blocks

- `QA.bed`：第一个输入文件，
    
- `QA.QV.lifted.anchors`：第二个输入文件，包含锚点信息。锚点是指在两个基因组之间对齐的基因对，通常是通过比对方法识别出的保守基因对。
    
- `--iter=1`：指定迭代次数。`mcscan` 工具可以通过迭代方法来优化共线性块的识别结果，这里指定进行一次迭代。
    
- `-o QA.QV.blocks`：指定输出文件。`mcscan` 工具会将识别出的共线性块结果输出到 `QA.QV.blocks` 文件中。
```

因为我们是两两比对，根据想要展出的格式，我们只需要让二月兰和拟南芥比较，拟南芥和油菜比较，拟南芥是中心，所以命令如下
```
1.二月兰和拟南芥
python -m jcvi.compara.synteny mcscan Athaliana.bed OV_sca.Athaliana.lifted.anchors --iter=1 -o OV_sca.Athaliana.blocks
2.拟南芥和油菜
python -m jcvi.compara.synteny mcscan Athaliana.bed Athaliana.Branap.lifted.anchors --iter=1 -o Athaliana.Branap.blocks
```

blocks结果里面放的是基因对

![[Pasted image 20240529164101.png]]

合并文件

```
python -m jcvi.formats.base join grape.peach.i1.blocks grape.cacao.i1.blocks --noheader | cut -f1,2,4,6 > grape.blocks

- **`grape.peach.i1.blocks` 和 `grape.cacao.i1.blocks`**：
    
    - 这两个是输入文件，包含了需要合并的数据。
- **`--noheader`**：
    
    - 这个选项表示输入文件没有表头，因此在合并过程中不考虑表头行。
- **`| cut -f1,2,4,6`**：
    
    - 使用管道操作符 `|` 将 `join` 工具的输出传递给 `cut` 命令。
    - `cut` 命令用于提取特定的字段。`-f1,2,4,6` 表示提取第1、2、4和6列的内容。
- **`> grape.blocks`**：
    
    - 将最终的结果重定向到 `grape.blocks` 文件中。
```


```
python -m jcvi.formats.base join OV_sca.Athaliana.blocks Athaliana.Branap.blocks --noheader | cut -f1,2,4,6 > ath.blocks
```
格式如下，注意列的顺序需要跟期望的排列顺序一致,如果不是可以用命令交换列
`awk '{print $2 "\t" $1 "\t" $3}' ath.blocks > ath1.blocks
![[Pasted image 20240529193540.png]]

选取感兴趣的基因,这里展示Ath1，OV6，Bra10,自己通过编辑软件筛选






配置文件如下
```
# x,   y, rotation,     ha,     va, color, ratio,            label
0.5, 0.6,        0, center,    top,      ,     1,       OV6
0.3, 0.4,        0, center, middle,      ,     5, Ath1
0.7, 0.4,        0, center, bottom,      ,     5, Bra10
# edges
e, 0, 1
e, 1, 2
```

运行画图：

把自己要画图的物种bed文件连接起来

```
$ cat OV_sca.bed Athaliana.bed  Branap.bed > Ath_OV_Bra.bed   
$ python -m jcvi.graphics.synteny ath.blocks Ath_OV_Bra.bed layout2.txt
```

结果如下
![[Pasted image 20240529201419.png]]