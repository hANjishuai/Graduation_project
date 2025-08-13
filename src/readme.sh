# 构建samples.tsv
cd data/Skin/
paste <(ls -1) <(ls -1 | xargs -L1 readlink -f) > ../../db/Pydatabase/sample.tsv
cd ../../
# 运行 get_anndata_from_10xMtx.py
nohup \
python ./src/get_anndata_from_10xMtx.py \
  --samples_table db/Pydatabase/sample.tsv \
  --output output/s00_get_anndata_from_10xMtx/combined.h5ad \
  --delimiter $'\t' \
   > output/s00_get_anndata_from_10xMtx/run.log 2>&1 &

# 运行 get_filteredParam_for_qc.py
nohup \
python ./src/get_filteredParam_for_qc.py \
    -i output/s00_get_anndata_from_10xMtx/combined.h5ad \
    -o output/s01_get_filteredParam_for_qc/qc_adata.h5ad \
    -p db/Pydatabase/s01_get_filteredParam_for_qc/ \
    --plot_dir output/s01_get_filteredParam_for_qc/plots \
    > output/s01_get_filteredParam_for_qc/run.log 2>&1 &

# 更新parameters，再次运行 get_filteredParam_for_qc.py
nohup \
python ./src/get_filteredParam_for_qc.py \
    -i output/s00_get_anndata_from_10xMtx/combined.h5ad \
    -o output/s01_get_filteredParam_for_qc/qc_adata.h5ad \
    -p db/Pydatabase/s01_get_filteredParam_for_qc/ \
    --plot_dir output/s01_get_filteredParam_for_qc/plots \
    --n_genes_max 7500 \
    --mt_pct_max 18 \
    --min_gene_count 3 \
    >> output/s01_get_filteredParam_for_qc/run.log 2>&1 &

#运行 get_PreHarmony_adata.py
nohup \
python ./src/get_PreHarmony_adata.py \
    -i output/s01_get_filteredParam_for_qc/qc_adata.h5ad \
    -o output/s02_get_PreHarmony_adata/pre_harmony.h5ad \
    -p db/Pydatabase/s01_get_filteredParam_for_qc/parameters.json \
    --n_top_genes 10000 \
    --plot_dir output/s02_get_PreHarmony_adata/plots \
    --copy_params \
    > output/s02_get_PreHarmony_adata/run.log 2>&1 &

#运行 get_PostHarmony_adata.py
nohup \
python ./src/get_PostHarmony_adata.py \
    --input output/s02_get_PreHarmony_adata/pre_harmony.h5ad \
    --output output/s03_get_PostHarmony_adata/post_harmony.h5ad \
    --plot_dir output/s03_get_PostHarmony_adata/plots \
    --n_pcs 50 \
    --n_neighbors 50 \
    --res_start 0.0 \
    --res_end 1.0 \
    --res_step 0.2 \
    --leiden_res 0.2 \
    > output/s03_get_PostHarmony_adata/run.log 2>&1 &

#运行 get_wantd_res.py
nohup \
python ./src/get_wantd_res.py \
    -i output/s03_get_PostHarmony_adata/post_harmony.h5ad \
    -o output/s04_get_wantd_res/all_src.h5ad \
    -p output/s04_get_wantd_res/plots \
    > output/s04_get_wantd_res/run.log 2>&1 &

#运行 get_preAnnotate_adata
nohup \
python ./src/get_preAnnotate_adata.py \
    -i output/s04_get_wantd_res/all_src.h5ad \
    -o output/s05_get_preAnnotate_adata/For_annotating.h5ad \
    -s "human" \
    -m "db/Pydatabase/s05_get_preAnnotate_adata/marker_file.csv" \
    -d "db/Pydatabase/s05_get_preAnnotate_adata/diff_genes_file.txt" \
    -p "output/s05_get_preAnnotate_adata/plots" \
    > output/s05_get_preAnnotate_adata/run.log 2>&1 &

##后面可以批量处理
#运行 get_Annotated_adata
nohup \
python ./src/get_Annotated_adata.py \
    -i output/s05_get_preAnnotate_adata/For_annotating.h5ad \
    -o output/s06_get_Annotated_adata/celltype_ker/sub_celltype.h5ad \
    -m db/Pydatabase/s06_get_Annotated_adata/mapping_file.csv \
    --marker-file db/Pydatabase/s06_get_Annotated_adata/marker_file.csv \
    -p output/s06_get_Annotated_adata/celltype_ker/plots \
    -t 'C01_Ker' \
    > output/s06_get_Annotated_adata/celltype_ker/run.log 2>&1 &

nohup \
python ./src/get_Annotated_adata.py \
    -i output/s05_get_preAnnotate_adata/For_annotating.h5ad \
    -o output/s06_get_Annotated_adata/celltype_B/sub_celltype.h5ad \
    -m db/Pydatabase/s06_get_Annotated_adata/mapping_file.csv \
    --marker-file db/Pydatabase/s06_get_Annotated_adata/marker_file.csv \
    -p output/s06_get_Annotated_adata/celltype_B/plots \
    -t 'C09_B' \
    > output/s06_get_Annotated_adata/celltype_B/run.log 2>&1 &

#运行 scanpy2seurat.py
nohup \
python ./src/scanpy2seurat.py \
    --input_h5ad output/s06_get_Annotated_adata/celltype_ker/sub_celltype.h5ad \
    --output_dir output/s07_scanpy2seurat/celltype_ker \
    --reductions "X_pca,X_umap,X_tsne" \
    > output/s07_scanpy2seurat/celltype_ker/run.log 2>&1 &

nohup \
python ./src/scanpy2seurat.py \
    --input_h5ad output/s06_get_Annotated_adata/celltype_B/sub_celltype.h5ad \
    --output_dir output/s07_scanpy2seurat/celltype_B \
    --reductions "X_pca,X_umap,X_tsne" \
    > output/s07_scanpy2seurat/celltype_B/run.log 2>&1 &

#运行 scanpy2seurat2RDS.R
nohup \
Rscript ./src/scanpy2seurat2RDS.R \
    -i output/s07_scanpy2seurat/celltype_ker \
    -o output/s07_scanpy2seurat/celltype_ker/subcell.rds \
    >> output/s07_scanpy2seurat/celltype_ker/run.log 2>&1 &

nohup \
Rscript ./src/scanpy2seurat2RDS.R \
    -i output/s07_scanpy2seurat/celltype_B \
    -o output/s07_scanpy2seurat/celltype_B/subcell.rds \
    >> output/s07_scanpy2seurat/celltype_B/run.log 2>&1 &
