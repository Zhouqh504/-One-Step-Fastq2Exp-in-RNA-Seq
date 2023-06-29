#!/bin/bash

# 默认值
thread=4
index_prefix="genome"
genome=""
fq_file=""
annotation=""
length=75
output_dir="output"


# 解析命令行参数
while getopts ":t:i:g:f:a:l:o:h" opt; do
  case $opt in
    t)
      thread=$OPTARG
      ;;
    i)
      index_prefix=$OPTARG
      ;;
    g)
      genome=$OPTARG
      ;;
    f)
      fq_file=$OPTARG
      ;;
    a)
      annotation=$OPTARG
      ;;
	l)
      length=$OPTARG
      ;;
    o)
      output_dir=$OPTARG
      ;;
    h)
      echo "Usage: fastq2Exp.sh [-t thread] [-i index_prefix] [-g genome] [-f fq_file] [-a annotation] [-o output_dir]"
      echo "Options:"
      echo "  -t thread: Number of threads to use (default: 4)"
      echo "  -i index_prefix: Prefix for index files (default: genome)"
      echo "  -g genome: Genome file"
      echo "  -f fq_file: Fastq file list"
      echo "  -a annotation: GFF file"
	  echo "  -l the average read length (default: 75)"
      echo "  -o output_dir: Output directory (default: output)"
      echo "  -h: Display this help message"
      exit 0
      ;;
    :)
      echo "Option -$OPTARG requires an argument." >&2
      exit 1
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      exit 1
      ;;
  esac
done

# 检查必要的参数是否提供
if [ -z "$genome" ] || [ -z "$fq_file" ] || [ -z "$annotation" ]; then
  echo "Missing required arguments. Please provide values for -g, -f, and -a options." >&2
  exit 1
fi

# 创建输出文件夹
if [ ! -d "$output_dir" ]; then
  mkdir "$output_dir"
fi

# 设置输出文件路径
echo -n '' > "$output_dir/sample_info.txt"
sample_info_file="$output_dir/sample_info.txt"
echo -n '' > "$output_dir/input_gtf_list.txt"
input_gtf_list="$output_dir/input_gtf_list.txt"
echo -n '' > "$output_dir/${annotation}_output_merged.gtf"
output_merged_gtf="$output_dir/${annotation}_output_merged.gtf"


# 检查索引文件是否存在
if [ ! -f "${index_prefix}.1.ht2" ]; then
    echo "Index file not found. Building index..."
    hisat2-build -p ${thread} ${genome} ${index_prefix}
else
    echo "Index file found. Skipping index building."
fi

# 读取fq文件中的路径并进行比对
echo >> ${fq_file} #最后一行加换行符

while IFS=$'\t' read -r fq1 fq2; do
    if [ -z "${fq1}" ] && [ -z "${fq2}" ]; then
        echo "Both fq1 and fq2 are empty. Mapping finished."
        break  # 结束循环
    fi

    echo "fq1: ${fq1}"
    echo "fq2: ${fq2}"

    output_prefix=$(basename "${fq1}" .fq.gz)
    bam_file="${output_dir}/${output_prefix}.sorted.bam"

    if [ -f "${bam_file}" ]; then
        echo "Alignment file ${bam_file} already exists. Skipping alignment."
    else
        if [ -z "${fq2}" ]; then
            echo "Running single-end alignment..."
            hisat2 -t --dta -p ${thread} -x ${index_prefix} -U ${fq1} | samtools view -bS - | samtools sort -o "${bam_file}"
        else
            echo "Running paired-end alignment..."
            hisat2 -t --dta -p ${thread} -x ${index_prefix} -1 ${fq1} -2 ${fq2} | samtools view -bS - | samtools sort -o "${bam_file}"
        fi
    fi

    # 进行接下来的操作...
    output_prefix=$(basename "${fq1}" .fq.gz)
    gtf_file="${output_dir}/${output_prefix}.gtf"

    echo "Reconstructing gtf of ${output_prefix}"
    stringtie -p ${thread} -e -G "${annotation}" -o "${gtf_file}" "${bam_file}"
    echo "Reconstructing gtf of ${output_prefix} finish"
    echo  "${gtf_file}" >> "${input_gtf_list}"

    echo "${output_prefix} ${gtf_file}" >> "${sample_info_file}"

done < ${fq_file}

echo "Sample gtf merging"
stringtie --merge -p ${thread} -G "${annotation}" -o "${output_merged_gtf}" "${input_gtf_list}"
echo "Sample gtf merging finished"

echo "Counting"
stringtie_prepcount.py -i "${sample_info_file}" -g "${output_dir}/gene_count_matrix.txt" -t "${output_dir}/transcript_count_matrix.txt" -l ${length}
stringtie_prepfpkm.py -i "${sample_info_file}" -g "${output_dir}/gene_fpkm_matrix.txt" -t "${output_dir}/transcript_fpkm_matrix.txt" -l ${length}
stringtie_preptpm.py -i "${sample_info_file}" -g "${output_dir}/gene_tpm_matrix.txt" -t "${output_dir}/transcript_tpm_matrix.txt" -l ${length}
echo "Counting finished"

# 删除临时文件
if [ -f "${sample_info_file}" ]; then
    rm "${sample_info_file}"
fi
if [ -f "${input_gtf_list}" ]; then
    rm "${input_gtf_list}"
fi
if [ -f "${output_dir}/transcript_fpkm_matrix.csv" ]; then
    rm "${output_dir}/transcript_fpkm_matrix.csv"
fi
if [ -f "${output_dir}/transcript_tpm_matrix.csv" ]; then
    rm "${output_dir}/transcript_tpm_matrix.csv"
fi
if [ -f "${output_dir}/transcript_count_matrix.csv" ]; then
    rm "${output_dir}/transcript_count_matrix.csv"
fi
