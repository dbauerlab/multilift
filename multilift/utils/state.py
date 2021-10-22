import logging
from pathlib import Path
from typing import Any, Union

import yaml

from multilift import __prog__
from multilift.utils import open_helper


# Globals #####################################################################


logger = logging.getLogger(__prog__)


# Classes #####################################################################


class State():
    ''' A class to store the current state ('liftover configuration') of
    multilift.  '''

    def __init__(self, file: Path) -> None:
        ''' Init the class and associate it with a file. We always require the
        class to be file-backed so that changes can be immediately dumped. '''
        self.file = file
        self.load()
        self._attr_defaults = {
            'name': None,
            'type': None,
            'filetype': None,
            'md5': None}

    def __repr__(self) -> str:
        return f'State(reference={self.reference}, liftovers={self.liftovers})'

    def __str__(self) -> str:
        ''' Return a YAML representation of the state '''
        return yaml.dump(
            self.data,
            default_flow_style=False, sort_keys=False, explicit_start=True)

    @property
    def reference(self) -> Union[str, None]:
        ''' The name of the reference genome for the current state'''
        return self.data.get('reference', None)

    @property
    def liftovers(self) -> list[str]:
        ''' The name(s) of the liftover genomes for the current state'''
        return self.data.get('liftovers', [])

    @property
    def genomes(self) -> list[str]:
        ''' All genomes loaded for the current state '''
        return ([] if not (ref := self.reference) else [ref]) + self.liftovers

    def load(self) -> None:
        ''' Load the YAML contents of self.file, overwriting current state '''
        with open_helper(self.file) as F:
            self.data = \
                {} if not (data := yaml.load(F, Loader=yaml.Loader)) else data
        self._modified = False

    def dump(self) -> None:
        ''' Dump current state as YAML to self.file, overwriting contents '''
        if self._modified:
            with open_helper(self.file, 'w') as F:
                print(self, file=F, flush=True)
        self._modified = False

    def add_genome(self, genome: str, as_reference: bool=False) -> None:
        ''' Add a genome to the current state and/or set an existing genome
        as the reference '''
        if as_reference:
            if self.reference and self.reference != genome:
                self.data['liftovers'] = \
                    [g for g in self.liftovers + [self.reference] if g != genome]
            self.data['reference'] = genome
        elif genome not in self.liftovers:
            self.data['liftovers'] = self.liftovers + [genome]

    def del_genome(self, genome: str) -> None:
        ''' Delete a genome from the current state and any associated files '''
        if self.reference == genome:
            del(self.data['reference'])
        elif genome in self.liftovers:
            self.data['liftovers'] = [g for g in self.liftovers if g != genome]
        try:
            del(self.data[genome])
        except KeyError:
            pass

    def __getitem__(self, target) -> Any:
        if type(target) is not tuple:
            target = (target, )
        genome, file, attr = target + (None, None, None)[len(target):]
        try:
            if len(target) == 1:
                # return {file: {attr: value, }, } of all files for `genome`
                return self.data[genome]
            if len(target) == 2:
                # return {attr: value, } of all attrs for `genome, file`
                return self.data[genome][file]
            if len(target) == 3:
                # return the value associated with `genome, file, attr` if
                # present else the default value for the attribute
                return self.data[genome][file].get(
                    attr, self._attr_defaults[attr])
        except KeyError:
            logger.error(
                f'The genome "{genome}", file "{file}", or attribute "{attr}"'
                'is not included in this multilift state')

    def __setitem__(self, target, value) -> None:
        genome, file, attr = target
        try:
            self.data[genome][file][attr] = value
        except KeyError:
            if genome not in self.genomes:
                logger.error(
                    f'Genome "{genome}" is not included in this multilift state')
            if genome not in self.data:
                self.data[genome] = {}
            if file not in self.data[genome]:
                self.data[genome][file] = {}
            self.data[genome][file][attr] = value

    def __delitem__(self, target) -> None:
        ''' Delete a genome, file, or attribute, and prune upwards if
        attribute deletion removes all information stored for a file '''
        if type(target) is not tuple:
            target = (target, )
        genome, file, attr = target + (None, None, None)[len(target):]
        try:
            if len(target) == 1:
                # delete the genome and all associated data
                self.del_genome(genome)
            if len(target) == 2:
                # delete a file associated with a genome
                del(self.data[genome][file])
            if len(target) == 3:
                # delete an attribute of a file
                del(self.data[genome][file][attr])
                # prune upwards if no more attrs are stored for file
                if not self.data[genome][file]:
                    del(self[genome][file])
        except KeyError:
            logger.error(
                f'The genome "{genome}", file "{file}", or attribute "{attr}"'
                'is not included in this multilift state')

    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        self.dump()


# Functions ###################################################################





###############################################################################
