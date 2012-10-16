"""
====================
Module **common.py**
====================

Common functions that are useful in many python projects.
"""

###############################################################################
def check_executable(tool_name):
    """Raises an exception if the executable *tool_name* is not found."""
    import os, subprocess
    with open(os.devnull, "w") as fnull:
        try:
            subprocess.Popen([tool_name,'-h'], stdout=fnull, stderr=fnull)
        except OSError:
            raise Exception("The executable '%s' cannot be found in $PATH" % tool_name)

###############################################################################
def wrap_string(s, width):
    """
    Wrap a string to the specified width.

    >>> wrap_string('aaabbb', 3)
    'aaa\\nbbb'
    """
    import textwrap
    wrapper = textwrap.TextWrapper()
    wrapper.width = width
    return '\n'.join(wrapper.wrap(s))

###############################################################################
def andify_strings(list_of_strings):
    """
    Given a list of strings will join them with commas
    and a final "and" word.

    >>> andify_strings(['Apples', 'Oranges', 'Mangos'])
    'Apples, Oranges and Mangos'
    """
    result = ', '.join(list_of_strings)
    comma_index = result.rfind(',')
    if comma_index > -1:
        result = result[:comma_index] + ' and' + result[comma_index+1:]
    return result

################################################################################
def property_cached(f):
    """Decorator for properties evaluated only once.
    It can be used to create a cached property like this::

        class Employee(object):
            @property_cached
            def salary(self):
                return 8000

        bob = Employee()
        print bob.salary
    """
    def get_method(self):
        try:
            return self.__cache__[f]
        except AttributeError:
            self.__cache__ = {}
            x = self.__cache__[f] = f(self)
            return x
        except KeyError:
            x = self.__cache__[f] = f(self)
            return x
    get_method.__doc__ = f.__doc__
    return property(get_method)

################################################################################
def natural_sort(item):
    """
    Sort strings that contain numbers correctly::

        >>> l = ['v1.3.12', 'v1.3.3', 'v1.2.5', 'v1.2.15', 'v1.2.3', 'v1.2.1']
        >>> l.sort(key=natural_sort)
        >>> l.__repr__()
        "['v1.2.1', 'v1.2.3', 'v1.2.5', 'v1.2.15', 'v1.3.3', 'v1.3.12']"
    """
    if item is None: return 0
    import re
    def try_int(s):
        try: return int(s)
        except ValueError: return s
    return map(try_int, re.findall(r'(\d+|\D+)', item))

################################################################################
def shift_list(l, shift, empty=0):
    """
    Shifts the elements of **l** by **shift** positions,
    filling new items with **empty**::

        >>> l = [0, 1, 4, 5, 7, 0]
        >>> shift_list(l, 3)
        [0, 0, 0, 0, 1, 4]
        >>> shift_list(l, -3)
        [5, 7, 0, 0, 0, 0]
        >>> shift_list(l, -8)
        [0, 0, 0, 0, 0, 0]
    """
    result = [empty] * len(l)
    old_index = max(shift, 0)
    new_index = max(-shift, 0)
    offset    = max(len(l) - abs(shift), 0)
    result[old_index:old_index + offset] = l[new_index:new_index + offset]
    return result

################################################################################
class Color(object):
    """Shortcuts for the ANSI escape sequences to control
       formatting, color, etc. on text terminals. Use it like this::

            print Color.red + "Hello world" + Color.end
    """
    # Special #
    end = '\033[0m'
    # Regular #
    blk   = '\033[0;30m' # Black
    red   = '\033[0;31m' # Red
    grn   = '\033[0;32m' # Green
    ylw   = '\033[0;33m' # Yellow
    blu   = '\033[0;34m' # Blue
    pur   = '\033[0;35m' # Purple
    cyn   = '\033[0;36m' # Cyan
    wht   = '\033[0;37m' # White
    # Bold #
    bold  = '\033[1m'
    b_blk = '\033[1;30m' # Black
    b_red = '\033[1;31m' # Red
    b_grn = '\033[1;32m' # Green
    b_ylw = '\033[1;33m' # Yellow
    b_blu = '\033[1;34m' # Blue
    b_pur = '\033[1;35m' # Purple
    b_cyn = '\033[1;36m' # Cyan
    b_wht = '\033[1;37m' # White
    # Light #
    light = '\033[2m'
    l_blk = '\033[2;30m' # Black
    l_red = '\033[2;31m' # Red
    l_grn = '\033[2;32m' # Green
    l_ylw = '\033[2;33m' # Yellow
    l_blu = '\033[2;34m' # Blue
    l_pur = '\033[2;35m' # Purple
    l_cyn = '\033[2;36m' # Cyan
    l_wht = '\033[2;37m' # White
    # Italic #
    italic = '\033[3m'
    i_blk = '\033[3;30m' # Black
    i_red = '\033[3;31m' # Red
    i_grn = '\033[3;32m' # Green
    i_ylw = '\033[3;33m' # Yellow
    i_blu = '\033[3;34m' # Blue
    i_pur = '\033[3;35m' # Purple
    i_cyn = '\033[3;36m' # Cyan
    i_wht = '\033[3;37m' # White
    # Underline #
    underline = '\033[4m'
    u_blk = '\033[4;30m' # Black
    u_red = '\033[4;31m' # Red
    u_grn = '\033[4;32m' # Green
    u_ylw = '\033[4;33m' # Yellow
    u_blu = '\033[4;34m' # Blue
    u_pur = '\033[4;35m' # Purple
    u_cyn = '\033[4;36m' # Cyan
    u_wht = '\033[4;37m' # White
    # Glitter #
    flash = '\033[5m'
    g_blk = '\033[5;30m' # Black
    g_red = '\033[5;31m' # Red
    g_grn = '\033[5;32m' # Green
    g_ylw = '\033[5;33m' # Yellow
    g_blu = '\033[5;34m' # Blue
    g_pur = '\033[5;35m' # Purple
    g_cyn = '\033[5;36m' # Cyan
    g_wht = '\033[5;37m' # White
    # Fill #
    f_blk = '\033[40m'   # Black
    f_red = '\033[41m'   # Red
    f_grn = '\033[42m'   # Green
    f_ylw = '\033[43m'   # Yellow
    f_blu = '\033[44m'   # Blue
    f_pur = '\033[45m'   # Purple
    f_cyn = '\033[46m'   # Cyan
    f_wht = '\033[47m'   # White

################################################################################
if __name__ == "__main__":
    import doctest
    doctest.testmod()