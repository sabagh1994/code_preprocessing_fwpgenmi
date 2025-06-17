#!/bin/bash

mkdir -p linked_files
ln -s /shared-mounts/sinhas-storage1/mayo/offer_project/recent_preprocessing/pipeline/gencode/protein_coding_genes.list ./linked_files/
ln -s /shared-mounts/sinhas-storage1/mayo/offer_project/OtherCelline_Analysis/overlap/encode-gencode/tfs_protein_coding_200Kb.bed10 ./linked_files
ln -s /shared-mounts/sinhas-storage1/mayo/offer_project/OtherCelline_Analysis/overlap/encode-gencode/tfs_protein_coding_50Kb.bed10 ./linked_files
