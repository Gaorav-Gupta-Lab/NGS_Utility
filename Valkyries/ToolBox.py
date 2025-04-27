import inspect
import os
from contextlib import suppress

def delete(file_list):
    """
    Delete one or more files.
    :param file_list:
    :return:
    """
    for file in file_list:
        with suppress(FileNotFoundError):
            os.remove(file)

def debug_messenger(reason: str = None):
    """
    Simple call for debugging that will print the filename and line number to the screen.
    :param reason:
    :return:
    """

    if reason is None:
        reason = "Programmer Neglected to Enlighten Us About the Need for Debugging This Section."

    frameinfo = inspect.getframeinfo(inspect.currentframe().f_back)
    print("\033[1;31m***WARNING: Debugging Module {0} at Line {1}.\n\t-->REASON: {2}\033[m"
          .format(frameinfo.filename, frameinfo.lineno, reason))
