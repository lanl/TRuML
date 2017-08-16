"""Custom exceptions for TRuML"""


class NotAMoleculeException(Exception):
    """Raised when a string is expected to, but does not, conform to molecule syntax"""

    def __init__(self, s):
        self.msg = "%s is not a molecule" % s

    def __str__(self):
        return self.msg


class NotCompatibleException(Exception):
    """Raised when a BNGL string or model cannot be converted to Kappa"""

    def __init__(self, s):
        self.msg = s

    def __str__(self):
        return self.msg


class NotConvertedException(Exception):
    """Raised when a string (if required) has not been converted to Kappa compatible syntax"""

    def __init__(self, s):
        self.msg = "Pattern '%s' not converted" % s

    def __str__(self):
        return self.msg


class NoModelsException(Exception):
    """Raised when no models are detected either in the command line arguments or the filesystem"""

    def __init__(self, s):
        self.msg = s

    def __str__(self):
        return self.msg


class UnknownMoleculeTypeException(Exception):
    """Raised when a Molecule instance has no known corresponding MoleculeDef"""

    def __init__(self, s):
        self.msg = "Cannot find MoleculeDef corresponding to '%s'" % s

    def __str__(self):
        return self.msg


class NoMoleculesException(Exception):
    """Raised when no molecule definitions are found (likely due to incorrect TRuML argument)"""

    def __init__(self, s):
        self.msg = s

    def __str__(self):
        return self.msg
