version 1.0

workflow extractMitoContigs {

    input {
        File chrM_fa
        String sampleName
        File inputFastaGZ
        File parse_script
        File mitoContig
        Int mat_pat_int
    }
    
    call blastFasta {
        input:
            chrM_fa=chrM_fa,
            sampleName=sampleName,
            inputFastaGZ=inputFastaGZ,
    }

    call parseBlastOutput {
        input:
            sampleName=sampleName,
            blastOutput=blastFasta.blastOutput,
            parse_script=parse_script
    }

    call correctMtAssembly {
        input:
            sampleName=sampleName,
            parsedBlastOutput=parseBlastOutput.parsedBlastOutput,
            unzippedOrigFa=blastFasta.unzippedOrigFa,
            mitoContig=mitoContig,
            mat_pat_int=mat_pat_int
    }

    output {
        File blastOutput       = blastFasta.blastOutput
        File parsedBlastOutput = parseBlastOutput.parsedBlastOutput
        File mitoContigsFn     = correctMtAssembly.mitoContigsFn
        File FinalAssembly     = correctMtAssembly.FinalAssembly
    }
}


task blastFasta {

    input {
        File chrM_fa
        String sampleName
        File inputFastaGZ
        String blastOutputName = "${sampleName}.BlastOutput.txt"

        Int memSizeGB = 4
        Int diskSizeGB = 64
        String dockerImage = "ncbi/blast:latest"
    }

    String inputFasta = basename(inputFastaGZ, ".gz")

    command <<<

        set -o pipefail
        set -e
        set -u
        set -o xtrace

        ## Create blastn database
        makeblastdb \
            -in ~{chrM_fa} \
            -title chrM_ucsc_hg38 \
            -dbtype nucl \
            -out chrM_ref/chrM_ref

        ## gunzip input fasta 
        gunzip -c ~{inputFastaGZ} > ~{inputFasta}
        

        ## Query fasta against MT database
        blastn \
            -query ~{inputFasta} \
            -db chrM_ref/chrM_ref \
            -num_threads 2 \
            -out ~{blastOutputName} \
            -outfmt '6 std qlen slen'
    >>>

    output {
        File blastOutput = blastOutputName
        File unzippedOrigFa = inputFasta
    }

    runtime {
        memory: memSizeGB + " GB"
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: dockerImage
        preemptible: 1
    }
}


task parseBlastOutput {

    input {
        String sampleName
        File blastOutput
        File parse_script
        String parsedBlastOutputName = "${sampleName}.ParsedBlastOutput.txt"

        Int memSizeGB = 4
        Int diskSizeGB = 64
        String dockerImage = "amancevice/pandas:latest"
    }

    command <<<

        set -o pipefail
        set -e
        set -u
        set -o xtrace


        ## gunzip input fasta 
        python ~{parse_script} ~{blastOutput} ~{parsedBlastOutputName}
        
    >>>

    output {
        File parsedBlastOutput = parsedBlastOutputName
    }

    runtime {
        memory: memSizeGB + " GB"
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: dockerImage
        preemptible: 1
    }
}

task correctMtAssembly {

    input {
        String sampleName
        File parsedBlastOutput
        File unzippedOrigFa
        File mitoContig
        Int mat_pat_int

        Int memSizeGB = 4
        Int diskSizeGB = 64
        String dockerImage = "biocontainers/samtools:latest"
    }

    String unzippedOrigFaFai = "${unzippedOrigFa}.fai"
    String mitoContigsFn     = "${sampleName}.ParsedBlastOutput.txt"
    String nonMitoContigs    = "${sampleName}.ParsedBlastOutput.txt"
    String nonMitoAssembly   = "${sampleName}.noMito.fa"
    String renameNonMitoAss  = "${sampleName}.noMito.renamed.fa"
    String FinalAssembly     = "${sampleName}.fa"

    command <<<

        set -o pipefail
        set -e
        set -u
        set -o xtrace

        ## Index original assembly
        samtools faidx ~{unzippedOrigFa}

        ## Pull mito contigs names, then filter out of all contig names to create nonMitoContigs
        cat ~{parsedBlastOutput} | cut -f1 | sed '1d'  > ~{mitoContigsFn}
        cat ~{unzippedOrigFaFai} | cut -f1 | grep -v -f ~{mitoContigsFn} > ~{nonMitoContigs}

        ## Pull nonMitoContigs out of Assembly
        samtools faidx ~{unzippedOrigFa} `cat ~{nonMitoContigs}` > ~{nonMitoAssembly}

        ## Rename contig names to sampleName#1/2#contigName format (1 = paternal, 2 = maternal)
        sed 's/^>/>~{sampleName}#~{mat_pat_int}/' ~{nonMitoAssembly} > ~{renameNonMitoAss}

        ## Now add in the MT assembly from Heng, and zip up the file
        cat ~{renameNonMitoAss} ~{mitoContig} | bgzip -c > ~{FinalAssembly}
        
    >>>

    output {
        File mitoContigsFn = mitoContigsFn
        File FinalAssembly = FinalAssembly
    }

    runtime {
        memory: memSizeGB + " GB"
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: dockerImage
        preemptible: 1
    }
}