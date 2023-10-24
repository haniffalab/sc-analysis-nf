#!/usr/bin/env/ nextflow

nextflow.enable.dsl=2

input = file(params.input)
out = file(params.output)

process MESSAGE {
 
    script:
    """
    echo message.py 'it works!' 
    """
}

process EXPRESSION {
    
    input: path(input)


    output: path(out)


    script:
    """
    echo expression.py 
    """
}

workflow { 
    MESSAGE
    EXPRESSION
}
