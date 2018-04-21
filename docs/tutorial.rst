========
Tutorial
========

Upon installation TRuML can be used from the command line by invoking the ``truml`` command with
various options.  All options can be seen with the command:

    :command:`truml -h`

Converting from Kappa to BNGL
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
In order to convert a model written in Kappa syntax to BNGL syntax, run the following command:

    :command:`truml -k model.ka`

where ``model.ka`` is the Kappa model of interest.

Converting from BNGL to Kappa
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Converting from BNGL to Kappa is similarly straightforward:

    :command:`truml -b model.bngl`

Verbose conversion
^^^^^^^^^^^^^^^^^^
The default settings are to output warnings to the terminal.  To enable a slightly more verbose output, describing the
various steps in translation, the ``-v`` flag is available.  Output useful for debugging is available by using the ``-d``
flag.