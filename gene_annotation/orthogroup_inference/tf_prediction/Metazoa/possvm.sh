#!/bin/bash

#####################################################################
# Batch POSSVM Analysis for Transcription Factor Gene Trees
# 
# Input:   iqtree_result/PFAM_ID/cluster_N.treefile
# Output:  possvm_results/
#####################################################################
# 激活环境
source ~/anaconda3/etc/profile.d/conda.sh
conda activate possvm
# ========== 配置参数 ==========
IQTREE_DIR="iq_tree_final"
POSSVM_OUTPUT="possvm_results"
MUS_GENE_NAMES="data/mus_gene_names.csv"
OUTGROUP_FILE="data/outgroup_species.txt"

# POSSVM 参数
ITER_MIDROOT=10
MIN_SUPPORT=50
CUT_GENE_NAMES=100

# 并行处理数（根据你的 CPU 核心数调整）
N_PARALLEL=8

# ========== 创建输出目录 ==========
mkdir -p ${POSSVM_OUTPUT}

# ========== 创建日志文件 ==========
LOG_FILE="${POSSVM_OUTPUT}/possvm_batch.log"
SUMMARY_FILE="${POSSVM_OUTPUT}/possvm_summary.txt"
ERROR_FILE="${POSSVM_OUTPUT}/possvm_errors.txt"

echo "=====================================================" | tee ${LOG_FILE}
echo "Batch POSSVM Analysis Started:  $(date)" | tee -a ${LOG_FILE}
echo "=====================================================" | tee -a ${LOG_FILE}
echo "" | tee -a ${LOG_FILE}

# ========== 检查必需文件 ==========
echo "Checking required files..." | tee -a ${LOG_FILE}

if [[ ! -f ${MUS_GENE_NAMES} ]]; then
    echo "ERROR: Mouse gene names file not found:  ${MUS_GENE_NAMES}" | tee -a ${LOG_FILE}
    exit 1
fi

if [[ !  -f ${OUTGROUP_FILE} ]]; then
    echo "ERROR: Outgroup species file not found: ${OUTGROUP_FILE}" | tee -a ${LOG_FILE}
    exit 1
fi

if [[ ! -d ${IQTREE_DIR} ]]; then
    echo "ERROR: IQ-TREE directory not found: ${IQTREE_DIR}" | tee -a ${LOG_FILE}
    exit 1
fi

echo "✓ All required files found" | tee -a ${LOG_FILE}
echo "" | tee -a ${LOG_FILE}

# ========== 显示外群物种 ==========
echo "Outgroup species:" | tee -a ${LOG_FILE}
cat ${OUTGROUP_FILE} | while read line; do
    echo "  - ${line}" | tee -a ${LOG_FILE}
done
echo "" | tee -a ${LOG_FILE}

# ========== 收集所有树文件 ==========
echo "Collecting gene tree files..." | tee -a ${LOG_FILE}

# 创建临时文件存储任务列表
TASK_LIST=$(mktemp)

total_trees=0
for family_dir in ${IQTREE_DIR}/*/ ; do
    
    if [[ ! -d ${family_dir} ]]; then
        continue
    fi
    
    family=$(basename ${family_dir})
    
    for treefile in ${family_dir}/*.treefile ; do
        
        if [[ ! -f ${treefile} ]]; then
            continue
        fi
        
        # 检查文件是否非空
        if [[ ! -s ${treefile} ]]; then
            echo "  Skipping empty file: ${treefile}" | tee -a ${LOG_FILE}
            continue
        fi
        
        # 提取 cluster 名称
        cluster=$(basename ${treefile} .treefile)
        
        # 生成输出前缀和文件名
        # 格式: PfamID. clusterID. 
        ogprefix="${family}.${cluster}."
        output_prefix="${family}.${cluster}"
        
        # 添加到任务列表
        echo "${treefile}|${family}|${cluster}|${ogprefix}|${output_prefix}" >> ${TASK_LIST}
        
        total_trees=$((total_trees + 1))
    done
done

echo "✓ Found ${total_trees} gene trees to process" | tee -a ${LOG_FILE}
echo "" | tee -a ${LOG_FILE}

if [[ ${total_trees} -eq 0 ]]; then
    echo "ERROR: No tree files found!" | tee -a ${LOG_FILE}
    rm ${TASK_LIST}
    exit 1
fi

# ========== 定义 POSSVM 运行函数 ==========
run_possvm() {
    local treefile=$1
    local family=$2
    local cluster=$3
    local ogprefix=$4
    local output_prefix=$5
    local output_dir=$6
    local mus_names=$7
    local outgroup=$8
    
    # 运行 POSSVM
    python ../../toolkit/possvm-orthology/possvm.py \
        -i ${treefile} \
        -o ${output_dir} \
        -r ${mus_names} \
        -refsps Mmus \
        -itermidroot ${ITER_MIDROOT} \
        -min_support_transfer ${MIN_SUPPORT} \
        -cut_gene_names ${CUT_GENE_NAMES} \
        -outgroup ${outgroup} \
        -ogprefix ${ogprefix} \
        -p ${output_prefix} \
        > ${output_dir}/${output_prefix}.possvm.stdout.log 2> ${output_dir}/${output_prefix}.possvm.stderr.log
    
    local exit_code=$?
    
    # 返回结果
    if [[ ${exit_code} -eq 0 ]]; then
        echo "SUCCESS|${family}|${cluster}"
    else
        echo "FAILED|${family}|${cluster}|exit_code_${exit_code}"
    fi
}

export -f run_possvm
export ITER_MIDROOT MIN_SUPPORT CUT_GENE_NAMES POSSVM_OUTPUT MUS_GENE_NAMES OUTGROUP_FILE

# ========== 并行运行 POSSVM ==========
echo "Running POSSVM on ${total_trees} gene trees (${N_PARALLEL} parallel jobs)..." | tee -a ${LOG_FILE}
echo "This may take a while..." | tee -a ${LOG_FILE}
echo "" | tee -a ${LOG_FILE}

start_time=$(date +%s)

# 使用 GNU parallel 或 xargs 进行并行处理
if command -v parallel &> /dev/null; then
    # 使用 GNU parallel
    cat ${TASK_LIST} | parallel -j ${N_PARALLEL} --colsep '|' \
        run_possvm {1} {2} {3} {4} {5} ${POSSVM_OUTPUT} ${MUS_GENE_NAMES} ${OUTGROUP_FILE}
else
    # 使用简单的后台进程
    processed=0
    running=0
    
    while IFS='|' read -r treefile family cluster ogprefix output_prefix; do
        
        # 等待直到有空闲的进程槽位
        while [[ ${running} -ge ${N_PARALLEL} ]]; do
            sleep 1
            running=$(jobs -r | wc -l)
        done
        
        # 后台运行
        (
            result=$(run_possvm ${treefile} ${family} ${cluster} ${ogprefix} ${output_prefix} ${POSSVM_OUTPUT} ${MUS_GENE_NAMES} ${OUTGROUP_FILE})
            echo ${result}
        ) &
        
        processed=$((processed + 1))
        running=$((running + 1))
        
        # 显示进度
        if [[ $((processed % 10)) -eq 0 ]]; then
            echo "  Processed:  ${processed}/${total_trees}" | tee -a ${LOG_FILE}
        fi
        
    done < ${TASK_LIST}
    
    # 等待所有后台任务完成
    wait
fi

end_time=$(date +%s)
elapsed=$((end_time - start_time))

echo "" | tee -a ${LOG_FILE}
echo "✓ All POSSVM jobs completed" | tee -a ${LOG_FILE}
echo "  Total time: ${elapsed} seconds ($((elapsed / 60)) minutes)" | tee -a ${LOG_FILE}
echo "" | tee -a ${LOG_FILE}

# ========== 统计结果 ==========
echo "=====================================================" | tee -a ${LOG_FILE}
echo "Analysis Summary" | tee -a ${LOG_FILE}
echo "=====================================================" | tee -a ${LOG_FILE}
echo "" | tee -a ${LOG_FILE}

# 统计成功和失败的数量
success_count=$(find ${POSSVM_OUTPUT} -name "*.ortholog_groups.csv" | wc -l)
failed_count=$((total_trees - success_count))

echo "Total gene trees processed: ${total_trees}" | tee -a ${LOG_FILE} | tee ${SUMMARY_FILE}
echo "Successful analyses:         ${success_count}" | tee -a ${LOG_FILE} | tee -a ${SUMMARY_FILE}
echo "Failed analyses:            ${failed_count}" | tee -a ${LOG_FILE} | tee -a ${SUMMARY_FILE}
echo "" | tee -a ${LOG_FILE} | tee -a ${SUMMARY_FILE}

# 按 TF 家族统计
echo "Results by TF family:" | tee -a ${LOG_FILE} | tee -a ${SUMMARY_FILE}
echo "-----------------------------------" | tee -a ${LOG_FILE} | tee -a ${SUMMARY_FILE}
printf "%-20s %10s\n" "Family" "Trees" | tee -a ${LOG_FILE} | tee -a ${SUMMARY_FILE}
echo "-----------------------------------" | tee -a ${LOG_FILE} | tee -a ${SUMMARY_FILE}

for family_dir in ${IQTREE_DIR}/*/ ; do
    family=$(basename ${family_dir})
    count=$(find ${POSSVM_OUTPUT} -name "${family}.*.ortholog_groups.csv" | wc -l)
    if [[ ${count} -gt 0 ]]; then
        printf "%-20s %10d\n" "${family}" "${count}" | tee -a ${LOG_FILE} | tee -a ${SUMMARY_FILE}
    fi
done

echo "-----------------------------------" | tee -a ${LOG_FILE} | tee -a ${SUMMARY_FILE}
echo "" | tee -a ${LOG_FILE} | tee -a ${SUMMARY_FILE}

# 列出输出文件类型
echo "Output file types:" | tee -a ${LOG_FILE} | tee -a ${SUMMARY_FILE}
for ext in ortholog_groups.csv pairs_orthologs.csv ortholog_groups.newick pairs_all.csv; do
    count=$(find ${POSSVM_OUTPUT} -name "*.${ext}" | wc -l)
    echo "  *.${ext}: ${count}" | tee -a ${LOG_FILE} | tee -a ${SUMMARY_FILE}
done
echo "" | tee -a ${LOG_FILE} | tee -a ${SUMMARY_FILE}

# ========== 检查失败的任务 ==========
if [[ ${failed_count} -gt 0 ]]; then
    echo "Checking failed analyses..." | tee -a ${LOG_FILE}
    
    > ${ERROR_FILE}
    
    while IFS='|' read -r treefile family cluster ogprefix output_prefix; do
        
        expected_output="${POSSVM_OUTPUT}/${family}.${cluster}.ortholog_groups.csv"
        
        if [[ ! -f ${expected_output} ]]; then
            echo "FAILED:  ${family} / ${cluster}" | tee -a ${ERROR_FILE}
            
            # 检查错误日志
            stderr_log="${POSSVM_OUTPUT}/${family}.${cluster}.possvm.stderr.log"
            if [[ -f ${stderr_log} && -s ${stderr_log} ]]; then
                echo "  Error log:" | tee -a ${ERROR_FILE}
                head -5 ${stderr_log} | sed 's/^/    /' | tee -a ${ERROR_FILE}
            fi
            echo "" | tee -a ${ERROR_FILE}
        fi
        
    done < ${TASK_LIST}
    
    echo "✓ Failed analyses logged to:  ${ERROR_FILE}" | tee -a ${LOG_FILE}
    echo "" | tee -a ${LOG_FILE}
fi

# ========== 清理 ==========
rm ${TASK_LIST}

# ========== 完成 ==========
echo "=====================================================" | tee -a ${LOG_FILE}
echo "Batch POSSVM Analysis Completed:  $(date)" | tee -a ${LOG_FILE}
echo "=====================================================" | tee -a ${LOG_FILE}
echo "" | tee -a ${LOG_FILE}
echo "Results saved in: ${POSSVM_OUTPUT}/" | tee -a ${LOG_FILE}
echo "Log file:          ${LOG_FILE}" | tee -a ${LOG_FILE}
echo "Summary:          ${SUMMARY_FILE}" | tee -a ${LOG_FILE}

if [[ ${failed_count} -gt 0 ]]; then
    echo "Errors:            ${ERROR_FILE}" | tee -a ${LOG_FILE}
fi

echo "" | tee -a ${LOG_FILE}

# 显示一些示例输出文件
echo "Example output files:" | tee -a ${LOG_FILE}
find ${POSSVM_OUTPUT} -name "*.ortholog_groups.csv" | head -5 | while read file; do
    echo "  ${file}" | tee -a ${LOG_FILE}
done

echo ""
echo "Done!"