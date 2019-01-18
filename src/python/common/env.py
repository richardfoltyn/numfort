import argparse
from argparse import ArgumentParser

from os.path import expanduser, join
import os.path

projname = 'numfort'
homedir = expanduser('~')
repodir = join(homedir, 'repos', projname)
fortran_root = join(repodir, 'src')


class Namespace(argparse.Namespace):
    """
    Custom Namespace class used to hold parsed command-line arguments
    or default values.
    """

    def __init__(self):
        super(Namespace, self).__init__()

        self.homedir = homedir
        self.fortran_root = fortran_root


def env_setup():
    """
    Perform environment setup, taking into account default values and
    user-provided CLI arguments.
    """

    p = ArgumentParser(description=projname)
    p.add_argument('--repo-dir', action='store', required=False,
                   dest='repodir', default=repodir, help='Git repository root')

    ns = Namespace()
    ns = p.parse_args(namespace=ns)

    ns.fortran_root = join(ns.repodir, 'src')

    return ns

