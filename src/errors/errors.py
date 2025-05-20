class CsvLoadError(Exception):
    """Exception raised when loading the CSV file fails."""

    def __init__(self, message: str = "Failed to load CSV file"):
        self.message = message
        super().__init__(self.message)


class GraphLoadError(Exception):
    """Exception raised when loading the Graph file fails."""

    def __init__(self, message: str = "Failed to load Graph JSON file"):
        self.message = message
        super().__init__(self.message)


class ReactionNotFoundError(Exception):
    """Exception raised when a reaction is not found."""

    def __init__(self, rxid: str, message: str = "Reaction with rxid not found"):
        self.rxid = rxid
        self.message = f"{message}: {rxid}"
        super().__init__(self.message)


class InchikeyNotFoundError(Exception):
    """Exception raised when a substance inchikey is not found."""

    def __init__(self, inchikey: str, message: str = "Substance with inchikey not found"):
        self.inchikey = inchikey
        self.message = f"{message}: {inchikey}"
        super().__init__(self.message)


class SubstanceNotFoundInSynthGraphError(Exception):
    """Exception raised when a substance is not found in the synthesis graph."""

    def __init__(self, inchikey: str, message: str = "Substance not found in synthesis graph"):
        self.inchikey = inchikey
        self.message = f"{message}: {inchikey}"
        super().__init__(self.message)


class MultipleInchikeyFoundError(Exception):
    """Exception raised when multiple substance nodes are found when only one is expected."""

    def __init__(
        self,
        inchikey: str,
        message: str = "Multiple substance nodes with inchikey found",
    ):
        self.inchikey = inchikey
        self.message = f"{message}: {inchikey}"
        super().__init__(self.message)


class ReactionYieldPredictionError(Exception):
    """Exception raised when a reaction yield prediction fails."""

    def __init__(self, rxsmiles: str, message: str = "Failed to predict yield"):
        self.rxsmiles = rxsmiles
        self.message = f"{message} for reaction: {rxsmiles}"
        super().__init__(self.message)


class NotCanonicalizableSmilesException(Exception):
    """Exception raised when a SMILES string cannot be canonicalized."""

    def __init__(self, rxsmiles: str):
        self.rxsmiles = rxsmiles
        self.message = f"Failed to canonicalize reaction smiles: {rxsmiles}"
        super().__init__(self.message)


class SynthGraphParsingException(Exception):
    """Exception raised when parsing the synthesis graph fails."""

    def __init__(self, message: str = "Failed to parse synthesis graph"):
        self.message = message
        super().__init__(self.message)


class RxsmilesAtomMappingException(Exception):
    """Exception raised when a reaction smiles has no atom mapping or cannot get one."""

    def __init__(self, rxsmiles: str, message: str = "Failed to get atom mapping"):
        self.rxsmiles = rxsmiles
        self.message = message
        super().__init__(self.message)
