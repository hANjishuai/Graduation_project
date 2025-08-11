# 抗原抗体对接全流程
# step1 读取单细胞组库文件,提取生产性序列并拼接V区（BCR.csv -> pair_BCR.csv）
python ./src/GetBcr_seq_information.py prepare \
  -i ./raw/DLE1/matrix.csv \
  -o ./antibody/DLE/DLE1/prepared_sequences.csv

# step2 配对抗体: 生成重链和轻链配对的FASTA文件 (pair_BCR.csv -> pair_BCR.fasta)
python ./src/GetBcr_seq_information.py pair \
  -i ./antibody/DLE/DLE1/prepared_sequences.csv \
  -o ./antibody/DLE/DLE1/paired_antibodies.fasta

# step3 使用IgFold预测抗体结构 (pair_BCR.fasta -> all_igfold_models/clonetype.pdb clonetype.fasta)
python ./src/GetBcr_seq_information.py fold \
  -f ./antibody/DLE/DLE1/paired_antibodies.fasta \
  -o ./antibody/DLE/DLE1/all_igfold_models \
  --refine \
  --renum \
  -l ./antibody/DLE/DLE1/fold.log

# step4 得到符合haddock3的pdb文件（主要是合并为一条链并重新排序）（clonetype.pdb -> clonetype_renumber.pdb）
mkdir -p ./antibody/DLE/DLE1/all_igfold_models/clonotype1/paratope_meta/
pdb_reres -1 \
  ./antibody/DLE/DLE1/all_igfold_models/clonotype1/clonotype1.pdb | \
  pdb_chain -A | \
  pdb_chainxseg | \
  pdb_tidy \
  > ./antibody/DLE/DLE1/all_igfold_models/clonotype1/paratope_meta/clonetype_renumbered.pdb

# step5 获取免疫组库文件中的cdr位置，作为paratope（clonetype.fasta pair_BCR.csv -> cdr_positions.txt antibody-paratope.act-pass clonetype_renumbered）
python ./src/fasta_tools.py locate-cdr \
  ./antibody/DLE/DLE1/all_igfold_models/clonotype1/clonotype1.fasta \
  -c ./antibody/DLE/DLE1/prepared_sequences.csv \
  -o ./antibody/DLE/DLE1/all_igfold_models/clonotype1/paratope_meta/cdr_positions.txt \
  -a ./antibody/DLE/DLE1/all_igfold_models/clonotype1/paratope_meta/antibody-paratope.act-pass \
  -p ./antibody/DLE/DLE1/all_igfold_models/clonotype1/paratope_meta/clonetype_renumbered

# step6 对抗原蛋白PDB文件进行筛选（0.35）
mkdir -p ./antigen/Skin_celltype/fitted_antigens/ 
# 根据epitoe的Ratio进行筛选
python ./src/Get_epitope.py filter-antigens \
  ./antigen/Skin_celltype \
  --threshold 0.35 

# step7 以抗原蛋白PDB(1EQT_A)为例进行标准化
cat ./antigen/Skin_celltype/fitted-antigens/1EQT_A/1EQT_A_discotope3.pdb | \
  pdb_tidy -strict | \
  pdb_delhetatm | \
  pdb_selaltloc | \
  pdb_keepcoord | \
  pdb_chain -B | \
  pdb_chainxseg | \
  pdb_tidy -strict \
  > ./antigen/Skin_celltype/fitted-antigens/1EQT_A/1EQT_A_clean.pdb

# step8 获取抗原蛋白的潜在的B细胞抗原表位
python ./src/generate_epitope_restraints.py extract-epitope \
  ./antigen/Skin_celltype/fitted_antigens/1EQT_A/1EQT_A_discotope3.csv \
  -p ./antigen/Skin_celltype/fitted_antigens/1EQT_A/1EQT_A_clean.pdb \
  -o ./antigen/Skin_celltype/fitted_antigens/1EQT_A/antigen-epitope.act-pass \
  -g --generate_passive > \
  ./antigen/Skin_celltype/fitted_antigens/1EQT_A/antigen-epitope-positions.txt

#python ./src/generate_epitope_restraints.py generate-passive \
#  ./antigen/Skin_celltype/fitted_antigens/1EQT_A/1EQT_A_clean.pdb \
#  ./antigen/Skin_celltype/fitted_antigens/1EQT_A/antigen-epitope.act-pass

# step9 固定BCR轻重链
haddock3-restraints restrain_bodies \
  ./antibody/DLE/DLE1/all_igfold_models/clonotype1/paratope_meta/clonetype_renumber.pdb > \
  ./antibody/DLE/DLE1/all_igfold_models/clonotype1/paratope_meta/antibody-unambig.tbl


# step10 生成约束文件
haddock3-restraints active_passive_to_ambig \
  ./antibody/DLE/DLE1/all_igfold_models/clonotype1/paratope_meta/antibody-paratope.act-pass \
  ./antigen/Skin_celltype/fitted_antigens/1EQT_A/antigen-epitope.act-pass \
  --segid-one A \
  --segid-two B > \
  ./registry/clonotype1/1EQT_A/ambig-paratope-NMR-epitope.tbl

haddock3-restraints validate_tbl \
  ./registry/clonotype1/1EQT_A/ambig-paratope-NMR-epitope.tbl \
  --silent

# step11 构建cfg文件，准备对接 
mkdir -p ./dock/DLE_DLE1_clonotype94/Skin_celltype_2VXW_B
python ./src/generate_epitope_restraints.py configure-docking \
    ./db/workflows/docking-antibody-antigen.cfg \
    --run_dir ./dock/DLE_DLE1_clonotype94/Skin_celltype_2VXW_B \
    --output ./registry/ambig_restraints/DLE_DLE1_clonotype94/Skin_celltype_2VXW_B/docking-antibody-antigen_updated.cfg \
    --antibody ./antibody/DLE/DLE1/all_igfold_models/clonotype94/paratope_meta/clonetype_renumber.pdb \
    --antigen ./antigen/Skin_celltype/fitted_antigens/2VXW_B/2VXW_B_clean.pdb \
    --ambig_fname ./registry/ambig_restraints/DLE_DLE1_clonotype94/Skin_celltype_2VXW_B/ambig-paratope-NMR-epitope.tbl \
    --unambig_fname ./antibody/DLE/DLE1/all_igfold_models/clonotype94/paratope_meta/antibody-unambig.tbl \
    --sampling 100 \
    --select 20 \
    --top_models 3   

# step12 对接
haddock3 ./registry/ambig_restraints/DLE_DLE1_clonotype94/Skin_celltype_2VXW_B/docking-antibody-antigen_updated.cfg 

# 使用bash命令，批量完成后，继续下面的步骤
python ./src/docking_processor.py process-docking-results \
    ./registry/docking_logs/docking_summary_master.csv \
    --output ./analysis/final_results.csv \
    --verbose #--dry-run

# 热图展示及prism制图文件
Rscript ./src/Plot_docking_results.R \
  -i ./analysis/final_results.csv \
  -o ./analysis
    
    