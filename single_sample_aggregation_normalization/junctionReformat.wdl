workflow CombineAllJunctions {

    Array[File] inputFiles
    String dockerVersion="4.1"
    String dockerContainer = "vanallenlab/splicejunction_normalization:${dockerVersion}"
    String sampleSetId
    File spliceJunctionReference
    Int? preemptNum
    Int preemptibleCount = select_first([preemptNum, 3])
    Int mem_gb=3

    scatter (spliceJunctionFile in inputFiles) {

        call spliceJunctionReformat { 
        	input: 
        		inputSpliceJunctionFile=spliceJunctionFile, 
                dockerContainer=dockerContainer, 
                preemptibleCount=preemptibleCount,
                mem_gb=mem_gb
                }
    }

    call spliceJunctionCombine { 
    	input: 
        	reformattedFiles=spliceJunctionReformat.reformatJunctionOutput, 
            dockerContainer=dockerContainer, 
            sampleSetId=sampleSetId, 
            preemptibleCount=preemptibleCount,
            mem_gb=mem_gb
            }

    call spliceJunctionNormalization { 
    	input: 
        	inputJunctionsFile=spliceJunctionCombine.allSamplesSpliceJunctions, 
            dockerContainer=dockerContainer, 
            sampleSetId=sampleSetId, 
            spliceReference=spliceJunctionReference, 
            preemptibleCount=preemptibleCount,
            mem_gb=mem_gb
            }

	output {
    	File junctions_all_samples = spliceJunctionNormalization.reformatJunctionOutput
    }
}


task spliceJunctionNormalization {

    File inputJunctionsFile
    String dockerContainer
    String sampleSetId 
    File spliceReference
    Int preemptibleCount
    Int mem_gb

    command <<<
    
    # log resource usage for debugging purposes
       	function runtimeInfo() {
        	echo [$(date)]
        	echo \* CPU usage: $(top -bn 2 -d 0.01 | grep '^%Cpu' | tail -n 1 | awk '{print $2}')%
        	echo \* Memory usage: $(free -m | grep Mem | awk '{ OFMT="%.0f"; print ($3/$2)*100; }')%
        	echo \* Disk usage: $(df | grep cromwell_root | awk '{ print $5 }')
        }
        while true;
        	do runtimeInfo;
           	sleep 15;
       	done &
        
        echo "outputs"
        mkdir outputs

        python /SpliceJunctionNormalization.py -splice_file ${inputJunctionsFile} -transcript_model ${spliceReference} -sample_set_id ${sampleSetId} -output_folder outputs --normalize
    >>>

    output { 
        File reformatJunctionOutput = "outputs/${sampleSetId}_outputs.txt"
    }

    runtime {
        docker: dockerContainer
        preemptible: preemptibleCount
        memory: mem_gb +" GB"
    }

}

task spliceJunctionReformat {

    File inputSpliceJunctionFile
    String dockerContainer
    String base = basename(inputSpliceJunctionFile)
    Int preemptibleCount
    Int mem_gb


    command <<<
    
    # log resource usage for debugging purposes
       	function runtimeInfo() {
        	echo [$(date)]
        	echo \* CPU usage: $(top -bn 2 -d 0.01 | grep '^%Cpu' | tail -n 1 | awk '{print $2}')%
        	echo \* Memory usage: $(free -m | grep Mem | awk '{ OFMT="%.0f"; print ($3/$2)*100; }')%
        	echo \* Disk usage: $(df | grep cromwell_root | awk '{ print $5 }')
        }
        while true;
        	do runtimeInfo;
           	sleep 15;
       	done &
        
        echo "outputs"
        mkdir outputs

        python /spliceJunctionReformat.py -junction_file ${inputSpliceJunctionFile} -output_folder outputs

    >>>

    output {
        File reformatJunctionOutput = "outputs/${base}.splicejunctions.reformatted.tsv"
    }

    runtime {
        docker: dockerContainer
        preemptible: preemptibleCount
        memory: mem_gb +" GB"
    }

}

task spliceJunctionCombine {

    Array[File] reformattedFiles
    String dockerContainer
    String sampleSetId
    Int preemptibleCount
    Int mem_gb

    command <<<
    
    # log resource usage for debugging purposes
       	function runtimeInfo() {
        	echo [$(date)]
        	echo \* CPU usage: $(top -bn 2 -d 0.01 | grep '^%Cpu' | tail -n 1 | awk '{print $2}')%
        	echo \* Memory usage: $(free -m | grep Mem | awk '{ OFMT="%.0f"; print ($3/$2)*100; }')%
        	echo \* Disk usage: $(df | grep cromwell_root | awk '{ print $5 }')
        }
        while true;
        	do runtimeInfo;
           	sleep 15;
       	done &

        echo "mkdir reformatted_junction_dir"
        mkdir reformatted_junction_dir

        echo "mv ${sep = ' ' reformattedFiles} reformatted_junction_dir"
        mv ${sep = ' ' reformattedFiles} reformatted_junction_dir

        echo "outputs"
        mkdir outputs

        python /combineJunctionAcrossSamples.py -input_junction_folder reformatted_junction_dir -output_folder outputs -sample_set_id ${sampleSetId}

        ls -lh outputs

    >>>

    output {
        File allSamplesSpliceJunctions = "outputs/${sampleSetId}_allsplicejunctions.tsv"
    }

    runtime {
        docker: dockerContainer
        preemptible: preemptibleCount
        memory: mem_gb +" GB"
    }

}