=================
Advanced concepts
=================

Variables, parameters, and functions
------------------------------------
Both Kappa and BNGL have language features that enable writing mathematical expressions for use in various aspects of
model construction.  These fall into two main categories:
    * Quantities that do not vary with time (static)
    * Quantities that may vary with time (dynamic)

Kappa defines both types as "variables", whereas BNGL distinguishes between the two as "parameters" (static) and
"functions" (dynamic).  In Kappa, variables can be defined in such a way that they are outputted to a file along with
other molecular "observables" of interest to the modeler.  When converting from BNGL to Kappa, the default behavior is
to reassign BNGL functions to be Kappa observables.  If this is not desired (i.e. the user wishes to suppress output of
these dynamic quantities), the flag ``-nf`` can be used to assign BNGL functions to be translated to variables that are
not outputted to file.

Syntax joining molecules and molecularity
-----------------------------------------
In BNGL, there are two operators that may be used to join molecule patterns: the dot (.) and the plus (+).  The dot
operator denotes that the two molecules are in the same complex (regardless of whether they have an explicit bond
between them).  The plus operator enforces separation of the molecules (meaning they are not in the same complex).  The
Kappa language only has the comma operator (,) to join molecule patterns.  Thus molecularity is ambiguous in arbitrary
patterns.  This has consequences for Kappa to BNGL conversion and vice versa

Kappa to BNGL
^^^^^^^^^^^^^
In Kappa rules, molecularity may be specified by a special rate syntax.  If such syntax is specified, conversion to BNGL
is precise.  Otherwise, it is possible that the rule may apply to both intermolecular and intramolecular reactions (if
it is a bond forming rule).  Since it is considered bad form to leave such rules ambiguous in the Kappa language, the
default behavior of TRuML is to assume that bond formation in Kappa rules is intermolecular unless otherwise specified.
To relax this assumption, TRuML provides the ``-dp`` flag, to convert all Kappa bond formation rules to two distinct
BNGL rules: one with the dot operator and one with the plus operator.

BNGL to Kappa
^^^^^^^^^^^^^
Since BNGL rules enforce molecularity by default, a similar assumption is made when converting to Kappa, that
bimolecular bond formation rules (those with the plus operator) can be directly converted to Kappa.  When converting
intramolecular bond formation rules to Kappa, the special rate syntax is used to enforce molecularity.

As mentioned earlier, the dot operator can be used to identify complexes without explicit specification of how they are
bound.  Kappa does not have this capability, and so BNGL observables that use such syntax will not be translated (a
warning will be shown), and rules that involve this syntax will result in a failed conversion