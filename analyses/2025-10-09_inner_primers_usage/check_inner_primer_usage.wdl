version 1.0

workflow CheckInnerPrimerUsage {
    input {
        File bam_file
        File bam_index
        File metadata_file
        String output_prefix = "primer_usage_report"
        Boolean search_3_prime = false
        Boolean reverse_complement = false
        String docker_image = "polumechanos/inner_primer_usage:latest"
        Int cpu = 4
        String memory = "8G"
        Int disk_space = 500
    }

    call AnalyzePrimerUsage {
        input:
            bam_file = bam_file,
            bam_index = bam_index,
            metadata_file = metadata_file,
            output_prefix = output_prefix,
            search_3_prime = search_3_prime,
            reverse_complement = reverse_complement,
            docker_image = docker_image,
            cpu = cpu,
            memory = memory,
            disk_space = disk_space
    }

    output {
        File summary_report = AnalyzePrimerUsage.summary_report
        File detailed_results = AnalyzePrimerUsage.detailed_results
    }
}

task AnalyzePrimerUsage {
    input {
        File bam_file
        File bam_index
        File metadata_file
        String output_prefix
        Boolean search_3_prime
        Boolean reverse_complement
        String docker_image
        Int cpu
        String memory
        Int disk_space
    }

    command <<<
        set -euo pipefail
        
        # Build command with optional flags
        cmd="python /app/check_inner_primer_usage.py ~{bam_file} ~{metadata_file} -o ~{output_prefix}.txt"
        
        if [ "~{search_3_prime}" == "true" ]; then
            cmd="$cmd --3-prime"
        fi
        
        if [ "~{reverse_complement}" == "true" ]; then
            cmd="$cmd --reverse-complement"
        fi
    
        
        echo "Running command: $cmd"
        eval $cmd
        
        # Ensure output files exist
        if [ ! -f "~{output_prefix}.txt" ]; then
            echo "Error: Summary report not generated"
            exit 1
        fi
        
        if [ ! -f "~{output_prefix}_detailed.tsv" ]; then
            echo "Error: Detailed results not generated"
            exit 1
        fi
    >>>

    runtime {
        docker: docker_image
        cpu: cpu
        memory: memory
        disks: "local-disk " + disk_space + " SSD"
        preemptible: 2
    }

    output {
        File summary_report = "~{output_prefix}.txt"
        File detailed_results = "~{output_prefix}_detailed.tsv"
    }
}
