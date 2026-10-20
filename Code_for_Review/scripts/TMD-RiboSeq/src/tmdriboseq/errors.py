"""User-facing exception types."""


class TMDRiboSeqError(Exception):
    """Base exception for expected, user-actionable failures."""


class ConfigurationError(TMDRiboSeqError):
    """Raised when a configuration file cannot be parsed or resolved."""


class ValidationFailed(TMDRiboSeqError):
    """Raised when validation reports one or more fatal errors."""

