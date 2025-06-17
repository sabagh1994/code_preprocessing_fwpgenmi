#!/bin/bash

mkdir -p linked_files
ln -s /shared-mounts/sinhas-storage1/mayo/offer_project/recent_preprocessing/pipeline/gencode/protein_coding_genes.list ./linked_files/
ln -s /shared-mounts/sinhas-storage1/mayo/offer_project/pipeline_saba/pipeline/overlap/encode_gencode/TFBS-only-analysis-saba/no-nearest-gene/tfs_protein_coding_200Kb.bed10 ./linked_files/
ln -s /shared-mounts/sinhas-storage1/mayo/offer_project/pipeline_saba/pipeline/overlap/encode_gencode/TFBS-only-analysis-saba/no-nearest-gene/tfs_protein_coding_50Kb.bed10 ./linked_files/
ln -s /shared-mounts/sinhas-storage1/mayo/offer_project/pipeline_saba/pipeline/overlap/encode_gencode/TFBS-only-analysis-saba/no-nearest-gene/tfs_protein_coding_10Kb.bed10 ./linked_files/
ln -s /shared-mounts/sinhas-storage1/mayo/offer_project/pipeline_saba/pipeline/overlap/encode_gencode/TFBS-only-analysis-saba/no-nearest-gene/tfs_protein_coding_1Mb.bed10 ./linked_files
