"""mirapy module for astronomical analysis, by Eric G. Suchanek, Ph.D., MIRA"""
# @Author Eric G. Suchanek, Ph.D.
# @License GNU
# @link none
# Microsoft Azure PAT: u3orvr6bwdm7fpmulnhw7bvv63c3njc2s2wqt3etnidb5joxdsga
#
# get rid of the pylint warning about invalid cases for variable names
# pylint: disable=C0103
# opylint: disable=C0111

version_info = ('v', 0, 2019, 10, 'egs')
__VERSION__ = '.'.join(map(str, version_info))
print('mirapy: ', __VERSION__)
