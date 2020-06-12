import "fastqc.wdl" as fastqc
import "multiqc.wdl" as multiqc


workflow FastQCWorkflow{

    Array[File] r1_files
    Array[File] r2_files
    String output_zip_name
    String git_repo_url
    String git_commit_hash

    Array[Pair[File, File]] fastq_pairs = zip(r1_files, r2_files)


    scatter(item in fastq_pairs){

        call trim_reads as read_trimmer{
            input:
                r1 = item.left,
                r2 = item.right
        }

        call fastqc.run_fastqc as fastqc_for_read1 {
            input:
                fastq = read_trimmer.trimmed_r1
        }

        call fastqc.run_fastqc as fastqc_for_read2 {
            input:
                fastq = read_trimmer.trimmed_r2
        }

    }

    call multiqc.create_qc as experimental_qc {
        input:
            r1_fastqc_zips = fastqc_for_read1.fastqc_zip,
            r2_fastqc_zips = fastqc_for_read2.fastqc_zip
    }

    call zip_results {
        input:
            zip_name = output_zip_name,
            multiqc_report = experimental_qc.report,
            r1_fastqc_zips = fastqc_for_read1.fastqc_zip,
            r2_fastqc_zips = fastqc_for_read2.fastqc_zip
    }

    output {
        File zip_out = zip_results.zip_out
    }

    meta {
        workflow_title : "FastQC "
        workflow_short_description : "For running FastQC on a bunch of FastQ"
        workflow_long_description : "Runs FastQC"
    }
}

task trim_reads {

    File r1
    File r2

    # Extract the samplename from the fastq filename
    String sample_name = basename(r1, "_R1.fastq.gz")

    Int disk_size = 200


    command {
        java -jar /opt/software/Trimmomatic-0.39/trimmomatic-0.39.jar PE \
            -trimlog ${sample_name}.trim.log \
            -summary ${sample_name}.trim_summary.log \
            ${r1} ${r2} \
            -baseout ${sample_name}.trimmed.fastq.gz \
            ILLUMINACLIP:/opt/software/Trimmomatic-0.39/adapters/TruSeq3-PE.fa:2:30:10:8:true
    }

    output {
        File trimmed_r1 = "${sample_name}.trimmed_1P.fastq.gz"
        File trimmed_r2 = "${sample_name}.trimmed_2P.fastq.gz"
    }

    runtime {
        docker: "docker.io/blawney/fastqc:v0.0.1"
        cpu: 4
        memory: "12 G"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 0
    }
}


task zip_results {

    String zip_name 
    File multiqc_report
    Array[File] r1_fastqc_zips
    Array[File]? r2_fastqc_zips

    Int disk_size = 100

    command {

        mkdir report
        mkdir report/qc

        mv ${multiqc_report} report/qc/
        mv ${sep=" " r1_fastqc_zips} report/qc/
        mv ${sep=" " r2_fastqc_zips} report/qc/

        zip -r "${zip_name}.zip" report
    }

    output {
        File zip_out = "${zip_name}.zip"
    }

    runtime {
        docker: "docker.io/blawney/fastqc:v0.0.1"
        cpu: 2
        memory: "6 G"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 0
    }
}
