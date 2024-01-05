.. _run:

Running
=======

The pipeline offers a workflow that processes files, images, and 
builds a Vitessce config file from the generated files.
Alternatively, the workflow to process files only, or the workflow to process images only  
can be called independently.

Each of these workflows work as entry points that can be specified when running the
pipeline through the command line.

- The ``-c`` Flag directs Nextflow to the configuration file ``nextflow.config``.
- The ``-t`` Flag directs Nextflow to the html template file ``template.html`` to make the log. 
- The ``Process_files`` workflow handles data files and their conversions.
- The ``Process_images`` workflow handles image files and/or label image data and their conversions.

Configurations and data are input through a :ref:`parameters yaml file <parameters_file>`.

To run the pipeline use

.. code-block:: shell

   nextflow run main.nf -c nextflow.config

This will handle all input files, whether they are data files or images, for all datasets
defined.

To view a log of the pipeline use
.. code-block:: shell

   nextflow log modest_mcclintock -t template.html > provenance.html

***************
Leaving this content here for now as a template for future adjustments to the documentation.
***************
You can modify the entry point if you're interested in only getting the converted outputs.
Use ``-entry Process_files`` or ``-entry Process_images`` as you need.

Running using Conda 
-------------------

The default pipeline will run on local executor without any type of environment creation. If you've already setup your conda environment you don't have to do anything else.

However, if you are working on a compute cluster you will need to make sure the conda environment is avaiable and active in your worker nodes. To run the pipeline using a new conda environment use the ``-profile conda`` option:

.. code-block:: shell

   nextflow run main.nf \
            -params-file /path/to/params.yaml \
            -entry Full_pipeline \
            -profile conda

Creating the environment when the pipleine is launched may take a few minutes.

Further reading
---------------

For more information about Docker image pulling/local conda env creation in Nextflow please refer to Nextflow's official docs for `containers <https://www.nextflow.io/docs/latest/container.html>`__ and `conda <https://www.nextflow.io/docs/latest/conda.html>`__.