
https://chendianyu.github.io/workflow/2019/09/03/snakemake-tutorial/

## Syntax
```
$ snakemake name-of-your-rule-or-file
```
(Just a few) useful arguments:
- -j - run in parallel when possible
- -f - force to rerun the rule
- -F - force rerun of all the pre-dependencies as well



## Output
It’s not required to list all possible outputs. Just those that you want to monitor or
that are used by a subsequent step as inputs.
```
rule preprocess_bsm:
    ''' Concatenate and normalize input dataset '''
input:
        script = 'preprocess.py',
        file = lambda wildcards: expand(rules.merge_h5_bsm_type.output,
            bsm_type=wildcards.bsm_type),
        qcd = rules.merge_h5_tuples.output.file
output:
        'output/BSM_preprocessed_scaling_{bsm_type}.h5'
shell:
        'python {input.script} --input-file {input.file} \
                               --output-file {output} \
                               --normalize \
                               --expand-dim \
                               --qcd-input {input.qcd}'
```

## Generalizing the read mapping rule
Snakemake能够通过使用通配符实现更加普遍性的规则，如下：
```
rule bwa_map:
    input:
        "data/genome.fa",
        "data/samples/{sample}.fastq"
    output:
        "mapped_reads/{sample}.bam"
    shell:
        "bwa mem {input} | samtools view -Sb - > {output}"
```

Snakemake 会将所有符合*.fastq 的文件处理成*.bam。
执行命令
```
snakemake -np mapped_reads/{A,B,C}.bam
## 等价于 snakemake -np mapped_reads/A.bam mapped_reads/B.bam mapped_reads/C.bam
## 告诉snakemake要生成的target files，就会用A、B、C去替换{sample}
```


```
SAMPLES=["A","B"]
rule bcftools_call:
    input:
        fa="data/genome.fa",
        bam=expand("sorted_reads/{sample}.bam", sample=SAMPLES),
        bai=expand("sorted_reads/{sample}.bam.bai", sample=SAMPLES)
    output:
        "calls/all.vcf"
    shell:
        "samtools mpileup -g -f {input.fa} {input.bam} | "
        "bcftools call -mv - > {output}"
```


## Config files
在上面的代码中，我们是将要跑的样本写在Snakefile，这样的话对于不同的项目，我们每次需要对Snakefile进行修改。对于这种情况，Snakemake提供了配置文件机制。配置文件可以用JSON或YAML语法，例如下面这个config.yaml文件：

```
samples:
    A: data/samples/A.fastq
    B: data/samples/B.fastq
```

而Snakefile文件就可以写成：
```

configfile: "config.yaml"
## Snakemake会加载配置文件并将文件中的内容放入名为config的全局字典中
## 形式为{'samples': {'A': 'data/samples/A.fastq', 'B': 'data/samples/B.fastq'}}

rule bcftools_call:
    input:
        fa="data/genome.fa",
        bam=expand("sorted_reads/{sample}.bam", sample=config["samples"]),
        bai=expand("sorted_reads/{sample}.bam.bai", sample=config["smaples])
    output:
        "calls/all.vcf"
    shell:
        "samtools mpileup -g -f {input.fa} {input.bam} | "
        "bcftools call -mv - > {output}"
## 虽然sample是个字典，但在展开时只会使用key的部分，即A和B
```