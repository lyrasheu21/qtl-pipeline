#!/usr/bin/env bash
# Genotype QC with plink 1.9.
#
#   1. import (bfile or VCF), keep listed donors, restrict to chromosomes
#   2. optionally update sample IDs
#   3. --geno / --maf / --hwe filters
#   4. --test-mishap; drop variants whose haplotype-missingness p < --mishap-p
#   5. write a bgzipped, tabix-indexed VCF (sample IDs = IIDs)
#
# Writes <out>.vcf.gz, <out>.vcf.gz.tbi and <out>_summary.tsv.
set -euo pipefail

in_type="" in_path="" keep="" update_ids="" chroms="" threads=1 out=""
geno="" maf="" hwe="" mishap_p=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    --in-type) in_type="$2" ;;
    --in) in_path="$2" ;;
    --keep) keep="$2" ;;
    --update-ids) update_ids="$2" ;;
    --chr) chroms="$2" ;;
    --geno) geno="$2" ;;
    --maf) maf="$2" ;;
    --hwe) hwe="$2" ;;
    --mishap-p) mishap_p="$2" ;;
    --threads) threads="$2" ;;
    --out) out="$2" ;;
    *) echo "Unknown argument: $1" >&2; exit 1 ;;
  esac
  shift 2
done

work="${out}_work"
mkdir -p "$work"
plink_common=(--keep-allele-order --allow-extra-chr --threads "$threads")

count() { wc -l < "$1" | tr -d ' '; }

# 1. import + keep + chromosomes
if [[ "$in_type" == "bfile" ]]; then input=(--bfile "$in_path"); else input=(--vcf "$in_path" --double-id); fi
keep_arg=()
[[ -n "$keep" ]] && keep_arg=(--keep "$keep")
plink "${input[@]}" "${plink_common[@]}" ${keep_arg[@]+"${keep_arg[@]}"} --chr "$chroms" --make-bed --out "$work/step1_import"
current="$work/step1_import"

# 2. update IDs
if [[ -n "$update_ids" ]]; then
  plink --bfile "$current" "${plink_common[@]}" --update-ids "$update_ids" --make-bed --out "$work/step2_ids"
  current="$work/step2_ids"
fi
n_var_in=$(count "$current.bim")
n_sample_in=$(count "$current.fam")

# 3. missingness / MAF / HWE
plink --bfile "$current" "${plink_common[@]}" --geno "$geno" --maf "$maf" --hwe "$hwe" \
  --make-bed --out "$work/step3_filters"
n_var_filters=$(count "$work/step3_filters.bim")

# 4. haplotype-based test for non-random missingness
plink --bfile "$work/step3_filters" "${plink_common[@]}" --test-mishap --out "$work/step4_mishap"
exclude="$work/mishap_exclude.txt"
if [[ -f "$work/step4_mishap.missing.hap" ]]; then
  awk -v p="$mishap_p" 'NR > 1 && $8 != "NA" && $8 + 0 < p { print $1 }' "$work/step4_mishap.missing.hap" \
    | sort -u > "$exclude"
else
  : > "$exclude"
fi
n_mishap=$(count "$exclude")

# 5. exclude and write VCF
exclude_arg=()
[[ -s "$exclude" ]] && exclude_arg=(--exclude "$exclude")
plink --bfile "$work/step3_filters" "${plink_common[@]}" ${exclude_arg[@]+"${exclude_arg[@]}"} \
  --recode vcf-iid bgz --out "$out"
tabix -f -p vcf "$out.vcf.gz"
n_var_out=$(bcftools view -H "$out.vcf.gz" | wc -l | tr -d ' ')

{
  printf "step\tvariants\tsamples\n"
  printf "input_after_keep_and_chr\t%s\t%s\n" "$n_var_in" "$n_sample_in"
  printf "after_geno_maf_hwe\t%s\t%s\n" "$n_var_filters" "$n_sample_in"
  printf "removed_by_mishap\t%s\tNA\n" "$n_mishap"
  printf "final\t%s\t%s\n" "$n_var_out" "$n_sample_in"
} > "${out}_summary.tsv"

rm -rf "$work"
