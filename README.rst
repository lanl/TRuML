TRuML: A translator for rule-based modeling
===========================================

Summary
-------
TRuML is a Python application designed to translate models written in either Kappa or the BioNetGen language (BNGL)
into the other language.  Currently TRuML is only compatible with Python 2 and work is being done to adapt TRuML for
use with Python 3.

Usage
-----
TRuML is invoked on the command line using the command ``truml``.  All options can be seen with the ``-h`` or ``--help``
flags::

    truml -h

Conversion of a model written in BNGL to the Kappa language is done using the command::

    truml -b model.bngl

The corresponding Kappa model will then exist in the same directory as the source BNGL file.  The reverse translation
is done in an equivalent manner::

    truml -k model.ka

Varying levels of verbosity can be achieved.  TRuML uses Python's logging library to track the translation process,
and it's default logging level is WARNING.  The ``-v`` (``--verbose``) and ``-d`` (``--debug``) flags will set the
logging level to INFO and DEBUG, respectively.  The logging output can be redirected to a file with the ``-l`` (``--log_file``)
flag and a corresponding file name.

Note that TRuML only translates the model structure.  It does not consider any simulation or simulation perturbation
directives from either language.

Caveats
-------
The two languages are not interchangeable (there are models in both languages that cannot be exactly translated to the
other).  In many cases, TRuML will raise issues explicitly.  However the languages are both still under development and
unidentified incompatibilities could cause the translation process to crash.  Some cases that cannot be exactly translated
include:

BNGL to Kappa
*************
 - Patterns containing molecules in complex but without explicit binding::

    A().B()

Kappa to BNGL
*************
 - Models containing infinities (``[inf]``)

Finally, the BioNetGen language is interpreted slightly differently for different simulation engines.  The NFsim engine
does not consider molecularity constraints on the right hand side of the rule.  As this is similar to Kappa's purely
local pattern representation (which cannot consider right hand side molecularity), TRuML assumes this convention.
