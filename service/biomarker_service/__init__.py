"""Thin job service around the diagnostic RNA-seq pipeline: upload -> validate -> queue (one job at a time) -> run `nextflow run main.nf` -> results.
Standard library only. The pipeline itself is never modified: the service only calls its documented command and reads its output folder."""
