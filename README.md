Hello world script
====================

A simple script showing the MiSeq 16S DADA2 for the Nextflow framework.
The output directory is at aws s3

```{bash}
nextflow run -resume main.nf --run_id 'MITI-MCB' --input_path 's3://maf-users/MITI/MiSeq' --output_path 's3://genomics-workflow-core/Results/MiSeq-16s-dada'
```

# Sample commands to submit an aws batch job
```{bash}
aws batch submit-job \
    --job-name nf-miseq-16s-dada2 \
    --job-queue priority-maf-pipelines \
    --job-definition nextflow-production \
    --container-overrides command="s3://nextflow-pipelines/nf-miseq-16s-dada2, \
    "--run_id", "MITI-MCB", \
    "--input_path", "s3://maf-users/MITI/MiSeq", \
    "--output_path", "s3://genomics-workflow-core/Results/MiSeq-16s-dada2" "
```
