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

        call fastqc.run_fastqc as fastqc_for_read1 {
            input:
                fastq = item.left
        }

        call fastqc.run_fastqc as fastqc_for_read2 {
            input:
                fastq = item.right
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
