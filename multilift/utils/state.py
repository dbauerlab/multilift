import logging
from pathlib import Path, PurePath
from typing import Any, Union

import yaml

from multilift import __prog__
from multilift.utils import basename, file_hash, guess_filetype, open_helper


# Globals #####################################################################


yaml.emitter.Emitter.process_tag = lambda self, *args, **kw: None

Pathish = Union[Path, PurePath, str]

logger = logging.getLogger(__prog__)


# Classes #####################################################################


class MultiliftState():
    ''' A class to store the current state (aka 'liftover configuration') of
    multilift '''

    _attr_defaults = {
        'name': '',
        'type': '',
        'filetype': '',
        'application': '',
        'md5': ''}

    def __init__(self, file: Pathish, dump_on_modify: bool=False,
                overwrite: bool=False) -> None:
        ''' Init the class and associate it with a file. We always require the
        class to be file-backed so that changes can be immediately dumped. '''
        self.file = Path(file)
        self.dump_on_modify = dump_on_modify
        if overwrite:
            self.data = {}
            self.dump()
        else:
            self.load()

    @property
    def genomes(self) -> list[str]:
        ''' The name(s) of the genomes available within the current state '''
        return self.data.get('genomes', [])

    def add_genome(self, genome: str, as_reference: bool=False) -> None:
        ''' Add a genome to the current state '''
        if genome not in self.genomes:
            self.data['genomes'] = self.genomes + [genome]
            self.data[genome] = {}
        if as_reference:
            self.data['reference'] = genome
        if self.dump_on_modify:
            self.dump()

    def del_genome(self, genome: str) -> None:
        ''' Delete a genome from the current state and all associated files '''
        try:
            self.data['genomes'] = [g for g in self.genomes if g != genome]
            if self.reference == genome:
                del(self.data['reference'])
            del(self.data[genome])
        except KeyError:
            pass
        if self.dump_on_modify:
            self.dump()

    @property
    def reference(self) -> Union[str, None]:
        ''' The name of the reference genome for the current state '''
        return self.data.get('reference', None)

    @reference.setter
    def reference(self, genome: str) -> None:
        ''' Set the reference genome for the current state '''
        self.add_genome(genome, True)

    def add_file(self, genome: str, file: Pathish,
            name: str='', filetype: str='',
            application: str='', md5: str='') -> None:
        ''' Add / overwrite a file object for `genome` '''
        file = Path(file)
        filename = str(file)
        if not file.exists():
            logger.error(f'File "{filename}" does not exist')
        meta = {}
        # calculate file md5
        h = file_hash(file)
        if md5 and md5 != h:
            logger.warning(
                f'"{filename}" has a different MD5 digest - has it changed?')
        meta['md5'] = h
        # auto-fill name with file basename
        meta['name'] = name if name else basename(file)
        # auto-fill filetype if guessable
        if not filetype:
            filetype, _ = guess_filetype(file)
        if filetype:
            meta['filetype'] = filetype
        # fill application
        if application:
            meta['application'] = application
        self.data[genome][filename] = meta
        if self.dump_on_modify:
            self.dump()

    def del_file(self, genome: str, file: Pathish) -> None:
        ''' Delete a `genome, file` '''
        file = Path(file)
        filename = str(file)
        try:
            del(self.data[genome][filename])
            if self.dump_on_modify:
                self.dump()
        except KeyError:
            pass

    def __getitem__(self, target) -> Any:
        target = (target, ) if not isinstance(target, tuple) else target
        genome, file, attr = target + (None, None, None)[len(target):]
        try:
            if len(target) == 1:
                # return {file: {attr: value, }, } of all files for `genome`
                return self.data[genome]
            else:
                file = Path(file)
                filename = str(file)
                if len(target) == 2:
                    # return {attr: value, } of all attrs for `genome, file`
                    return self.data[genome][filename]
                if len(target) == 3:
                    # return the value associated with `genome, file, attr` if
                    # present else the default value for the attribute
                    return self.data[genome][filename].get(
                        attr, self._attr_defaults[attr])
        except KeyError:
            logger.error(
                f'The genome "{genome}", file "{file}", or attribute "{attr}" '
                'is not included in this multilift state')

    def __setitem__(self, target, value) -> None:
        genome, file, attr = target
        if attr not in self._attr_defaults:
            logger.error(f'Unsupported attribute "{attr}"')
        file = Path(file)
        filename = str(file)
        try:
            self.data[genome][filename][attr] = value
        except KeyError:
            if genome not in self.genomes:
                logger.error(
                    f'Genome "{genome}" is not included in this multilift state')
            elif filename not in self.data[genome]:
                logger.error(
                    f'File "{filename}" is not included in this multilift state')
        if self.dump_on_modify:
            self.dump()

    def __delitem__(self, target) -> None:
        ''' Delete a genome, file, or attribute '''
        target = (target, ) if not isinstance(target, tuple) else target
        genome, file, attr = target + (None, None, None)[len(target):]
        try:
            if len(target) == 1:
                # delete the genome and all associated data
                self.del_genome(genome)
            else:
                file = Path(file)
                filename = str(file)
                if len(target) == 2:
                    # delete a file associated with a genome
                    del(self.data[genome][filename])
                if len(target) == 3:
                    # delete an attribute of a file
                    del(self.data[genome][filename][attr])
        except KeyError:
            if genome not in self.genomes:
                logger.error(
                    f'Genome "{genome}" is not included in this multilift state')
            elif filename not in self.data[genome]:
                logger.error(
                    f'File "{filename}" is not included in this multilift state')
            else:  # attr not in self.data[genome][filename]
                logger.error(
                    f'The attribute "{attr}" is not recorded for file "{filename}"')
        if self.dump_on_modify:
            self.dump()

    def __repr__(self) -> str:
        return f'MultiliftState({self.genomes})'

    def __str__(self) -> str:
        ''' Return a YAML representation of the state '''
        rep = {}
        if self.genomes:
            rep['genomes'] = self.genomes
        if self.reference:
            rep['reference'] = self.reference
        for genome in self.genomes:
            if self[genome]:
                rep[genome] = self[genome]
        return yaml.dump(
            rep,
            default_flow_style=False, sort_keys=False, explicit_start=True)

    def load(self) -> None:
        ''' Load the YAML contents of self.file, overwriting current state '''
        with open_helper(self.file) as F:
            self.data = {} if not (data := yaml.safe_load(F)) else data
        for genome in self.genomes:
            if genome not in self.data:
                self.data[genome] = {}

    def dump(self) -> None:
        ''' Dump current state as YAML to self.file, overwriting contents '''
        with open_helper(self.file, 'w') as F:
            print(self, file=F, flush=True)

    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        self.dump()


# Functions ###################################################################





###############################################################################
