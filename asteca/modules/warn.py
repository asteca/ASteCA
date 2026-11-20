"""Warnings issued by ASteCA, shown as ``UserWarning: message`` in lines of at most
``WIDTH`` characters (without the path of the file that issued them and its source
line)"""

import textwrap
import warnings

# Maximum length of the lines of a printed warning
WIDTH = 80


class ASteCAWarning(UserWarning):
    """Warning issued by ASteCA (a ``UserWarning``)"""


def warn(message: str, stacklevel: int = 1) -> None:
    """Issue an :py:class:`ASteCAWarning`, with the line breaks and repeated spaces of
    ``message`` collapsed into single spaces

    :param message: Warning message
    :type message: str
    :param stacklevel: As in ``warnings.warn``, counted from the caller of this
        function: ``1`` (default) is the line that calls :py:func:`warn`
    :type stacklevel: int
    """
    warnings.warn(" ".join(message.split()), ASteCAWarning, stacklevel=stacklevel + 1)


_formatwarning = warnings.formatwarning


def _format(message, category, filename, lineno, line=None) -> str:
    """Short format for :py:class:`ASteCAWarning`, wrapped at ``WIDTH`` characters;
    the original format for any other warning"""
    if issubclass(category, ASteCAWarning):
        return textwrap.fill(f"UserWarning: {message}", WIDTH) + "\n"
    return _formatwarning(message, category, filename, lineno, line)


warnings.formatwarning = _format
