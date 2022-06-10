Channel.fromPath("REF_VCF/*.vcf.gz").into{ref_vcf; ref_vcf_2}

raw_vcf = Channel.fromPath("QUERY_VCF/*.vcf.gz")

referenceSequence = file("data/MTBC_ancestor.fa.gz") //pseudotuberculosis

TAXO = file("data/taxonomical_info.txt")
OUTGROUP = file("data/outgroup.txt")

results = file("results")

process samtoolsIndex {
	errorStrategy = "finish"
	label = 'samtools'

	input:
	file referenceSequence

	output:
	tuple file("${referenceSequence.simpleName}.fasta"), file("${referenceSequence.simpleName}.fasta.fai") into refSeq_index

	script:
	"""
	gunzip -c $referenceSequence > ${referenceSequence.simpleName}.fasta
	samtools faidx ${referenceSequence.simpleName}.fasta
	"""
}

process picardIndex {
	label 'GATK'

	input:
	tuple file(fasta), file(fai) from refSeq_index

	output:
	tuple file(fasta), file(fai), file("*.dict") into refSeq_index_dict

	script:
	"""
	gatk CreateSequenceDictionary -R $fasta
	"""
}

process getReferenceSNP {
	label 'bcftools'

	input:
	file(vcf) from ref_vcf

	output:
	file("${vcf.simpleName}.snps.vcf.gz") into snp_vcf, snp_vcf_2

	script:
	"""
	bcftools index $vcf
	bcftools view --types snps -i "QUAL>20 && FORMAT/DP>10" $vcf -Oz -o ${vcf.simpleName}.snps.vcf.gz
	"""
}

process mergeAllVariants {
	label 'bcftools'

	input:
	file(vcf) from snp_vcf.collect()

	output:
	file("ref_snp_pos.bed") into ref_snp_pos

	script:
	"""
	zcat $vcf | grep -v '#' | awk 'OFS="\t" {if(length(\$4) == length(\$5)) print \$1,\$2-1,\$2}' | sort -u > ref_snp_pos.bed
	cp ref_snp_pos.bed $results
	"""
}

process VCF2Fasta_ref {	
	label 'bcftools'

	input:
	file(vcf) from snp_vcf_2
	tuple file(fasta), file(fai), file(dict) from refSeq_index_dict

	output:
	file("${sample}_rebuilt.fasta") into rebuilt_fasta

	script:
	sample = vcf.simpleName
	"""
	bcftools index $vcf
	echo ">"$sample > ${sample}_rebuilt.fasta
	cat $fasta | bcftools consensus $vcf --haplotype 1 | grep -v '>' | tr -d '\n' >> ${sample}_rebuilt.fasta
	echo "" >> ${sample}_rebuilt.fasta
	"""
}

process treeBuild {
	label 'raxml'

	input:
	file(fasta) from rebuilt_fasta.collect()

	output:
	file("MSA.fasta") into reference_alignement
	file("full_tree.raxml.bestTree") into reference_tree
	file("full_tree.raxml.bestModel") into reference_model

	script:
	"""
	cat *.fasta > MSA.fasta
#	goalign subsites -i MSA.fasta --informative > MSA_subsites.fasta

	raxml-ng --all --msa MSA.fasta --model GTR --bs-trees 100 --prefix full_tree --threads ${task.cpus}c

	mkdir -p $results/TREE
	cp MSA.fasta $results
	cp full_tree* $results/TREE
	"""
}


process filterVariants {
	label 'bcftools'

	input:
	file(bed) from ref_snp_pos
	file(vcf) from raw_vcf

	output:
	file("${sample}.filtered.vcf.gz") into target_vcf

	script:
	sample = vcf.simpleName
	"""
	bcftools index $vcf
	bcftools view $vcf --types snps -i "FORMAT/DP>5" -R $bed -Oz -o ${sample}.filtered.vcf.gz
	"""
}

process VCF2Fasta_target {
	label 'bcftools'

	input:
	file(vcf) from target_vcf
	tuple file(fasta), file(fai), file(dict) from refSeq_index_dict

	output:
	file("${sample}_rebuilt.fasta") into sample_to_place

	script:
	sample = vcf.simpleName
	"""
	bcftools index $vcf
	echo ">"$sample > ${sample}_rebuilt.fasta

	cat $fasta | bcftools consensus $vcf --haplotype 1 --missing 'N' --absent 'N' | grep -v '>' >> ${sample}_rebuilt.fasta
	echo "" >> ${sample}_rebuilt.fasta
	"""
}

process placement {
	label 'epa'

	errorStrategy="terminate"

	publishDir "$results", mode: "copy", pattern: '*.jplace'
    publishDir "$results", mode: 'copy', pattern: '*.txt'

	input:
	file(TAXO)
	file(FASTA) from sample_to_place.collect()
	file(MSA) from reference_alignement
	file(TREE) from reference_tree
	file(MODEL) from reference_model

	output:
	file("*.jplace")
	file("output_per_query.txt")

	script:
	"""
	cat $FASTA > QUERY.msa

	epa-ng --tree ${TREE} --ref-msa $MSA --query QUERY.msa -T ${task.cpus} --model $MODEL --no-heur

	gappa examine assign --jplace-path . --file-prefix output_ \\
		--taxon-file $TAXO --per-query-results --best-hit --root-outgroup $OUTGROUP --resolve-missing-paths

	mv output_per_query.tsv output_per_query.txt
	"""
}
