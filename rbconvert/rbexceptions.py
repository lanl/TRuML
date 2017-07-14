class NotAMoleculeException(Exception):
    """Raised when a string is expected to, but does not, conform to molecule syntax"""

    def __init__(self, s):
        print "%s is not a molecule" % s


class NotCompatibleException(Exception):
    """Raised when a BNGL string or model cannot be converted to Kappa"""

    def __init__(self, s):
        print s


class NotConvertedException(Exception):
    """Raised when a string (if required) has not been converted to Kappa compatible syntax"""

    def __init__(self):
        print "Must convert object due to identically named sites"
