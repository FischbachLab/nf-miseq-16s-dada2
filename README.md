README
====================
# 1. Data preparation

## 1) Install Illumina [Basemount](https://help.basespace.illumina.com/cmd-line-interfaces/basespace-cli/introduction-to-basemount)
 Basemount is a tool to mount your BaseSpace Sequence Hub data as a Linux file system

## 2) Mount your BaseSpace account
```{bash}
basemount /path/to/local/storage/Basespace
```
## 3) Copy sequencing files only to a local dir, then syn into s3
```{bash}
# replace project_name with your project name
mkdir -p /path/to/local/storage/myBasespace/project_name/fastqs
cp /path/to/local/storage/Basespace/Projects/project_name/Untitled\ from\ /Samples/*/Files/*  /path/to/local/storage/myBasespace/project_name/fastqs/
```
## 4) Save fastq files into s3
```{bash}
aws s3 sync /path/to/local/storage/myBasespace/project_name/fastqs/  s3://maf-users/MITI/MiSeq/project_name/fastqs/
```


# 2. Run the MiSeq 16S DADA2 pipeline via Nextflow

## Sample commands to run a local job
```{bash}
nextflow run -resume main.nf --project 'MITI-MCB' --input_path 's3://maf-users/MITI/MiSeq' --output_path 's3://genomics-workflow-core/Results/MiSeq-16s-dada2'
```

## Sample commands to submit an aws batch job
```{bash}
aws batch submit-job \
    --job-name nf-miseq-16s-dada2 \
    --job-queue priority-maf-pipelines \
    --job-definition nextflow-production \
    --container-overrides command="FischbachLab/nf-miseq-16s-dada2, \
    "-r", "main", \
    "--project", "MITI-MCB", \
    "--config", "s3://nextflow-pipelines/nf-miseq-16s-dada2/conf/parameters.yaml",\
    "--input_path", "s3://maf-sequencing/Illumina/MiSeq/MITI-MCB", \
    "--output_path", "s3://genomics-workflow-core/Results/MiSeq-16s-dada2" "
```

# Outputs
The full output is saved at
```{bash}
s3://genomics-workflow-core/Results/MiSeq-16s-dada2/MITI-MCB/
```

The summary file can be found at:
```{bash}
s3://genomics-workflow-core/Results/MiSeq-16s-dada2/MITI-MCB/DADA2_summary/Genus_summary.tsv
```
