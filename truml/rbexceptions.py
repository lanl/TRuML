"""truml.rbexceptions: module containing custom exceptions for TRuML"""


class NotAMoleculeException(Exception):
    """Raised when a string is expected to, but does not, conform to molecule syntax"""

    def __init__(self, s):
        self.msg = "%s is not a molecule" % s


class NotCompatibleException(Exception):
    """Raised when a BNGL string or model cannot be converted to Kappa"""

    def __init__(self, s):
        self.msg = s


class NotConvertedException(Exception):
    """Raised when a string (if required) has not been converted to Kappa compatible syntax"""

    def __init__(self):
        self.msg = "Must convert object due to identically named sites"


class NoModelsException(Exception):
    """Raised when no models are detected either in the command line arguments or the filesystem"""

    def __init__(self, s):
        if s is not None:
            self.msg = s
