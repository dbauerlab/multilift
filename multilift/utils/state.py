import logging
from pathlib import Path, PurePath
from typing import Any, Union

import yaml

from multilift import __prog__
from multilift.utils import basename, file_hash, guess_filetype, open_helper


# Globals #####################################################################


Pathish = Union[Path, PurePath, str]

logger = logging.getLogger(__prog__)


# Classes #####################################################################


class MultiliftState():
    ''' A class to store the current state (aka liftover configuration) of
    multilift '''

    _file_attr_defaults = {
        'application': '',
        'genome': '',
        'name': '',
        'type': '',
        'filetype': '',
        'md5': ''}

    def __init__(self, file: Pathish,
                dump_on_modify: bool=False, overwrite: bool=False) -> None:
        ''' Init the class and associate it with a file. We always require the
        class to be file-backed so that changes can be immediately dumped. '''
        self.file = Path(file)
        self.dump_on_modify = dump_on_modify
        if overwrite:
            with open_helper(self.file, 'w'):
                pass
        self.load()

    def add_genome(self, genome: str, as_reference: bool=False) -> None:
        ''' Add a genome to the current state '''
        if genome not in self.genomes:
            self.genomes.append(genome)
        if as_reference:
            self.reference = genome
        if self.dump_on_modify:
            self.dump()

    def del_genome(self, genome: str) -> None:
        ''' Delete a genome from the current state and all associated files '''
        self.genomes = [g for g in self.genomes if g != genome]
        if self.reference == genome:
            self.reference = None
        self._files = {
            k: v for k, v in self._files.items()
            if v.get('genome', _file_attr_defaults['genome']) != genome}
        if self.dump_on_modify:
            self.dump()

    def add_file(self, file: Pathish,
            md5: str='', genome: str='', application: str='',
            name: str='', filetype: str='') -> None:
        ''' Add (/ overwrite) a file to the current state '''
        file = Path(file)
        if not file.exists():
            logger.error(f'File does not exist: {file}')
        if genome and genome not in self.genomes:
            logger.error(f'Genome is not included in this state: {genome}')
        filename = str(file)
        self._files[filename] = {}
        # calculate file md5
        self._files[filename]['md5'] = file_hash(file)
        if md5 and md5 != self._files[filename]['md5']:
            logger.warning(
                f'"File has changed since this state was initiated: {file}')
        # fill name or auto-fill with file basename
        self._files[filename]['name'] = name if name else basename(file)
        # fill genome
        if genome:
            self._files[filename]['genome'] = genome
        # fill filetype and application if guessable
        guessed_filetype, guessed_application = guess_filetype(file)
        filetype = filetype if filetype else guessed_filetype
        if filetype:
            self._files[filename]['filetype'] = filetype
        if application:
            self._files[filename]['application'] = application
        elif filetype:
            self._files[filename]['application'] = guessed_application[0]
        # wrap up
        if self.dump_on_modify:
            self.dump()

    def del_file(self, file: Pathish) -> None:
        ''' Delete a file from the current state '''
        try:
            del(self._files[str(Path(file))])
            if self.dump_on_modify:
                self.dump()
        except KeyError:
            pass

    def files(self, application: str='', genome: str='') -> list:
        if application and genome:
            return [
                k for k, v  in self._files.items()
                if v.get(
                    'application', self._file_attr_defaults['application']
                    ) == application
                and v.get(
                    'genome', self._file_attr_defaults['genome']
                    ) == genome]
        if application:
            return [
                k for k, v  in self._files.items()
                if v.get(
                    'application', self._file_attr_defaults['application']
                    ) == application]
        if genome:
            return [
                k for k, v  in self._files.items()
                if v.get(
                    'genome', self._file_attr_defaults['genome']
                    ) == genome]
        return list(self._files.keys())

    def __getitem__(self, target) -> Any:
        target = (target, ) if not isinstance(target, tuple) else target
        file, attr = target + (None, None)[len(target):]
        file = Path(file)
        filename = str(file)
        try:
            if len(target) == 1:
                # return {attr: value, } dict for `file`
                return self._files[filename]
            # return the value (or default value) for `file, attr`
            return self._files[filename].get(attr, self._file_attr_defaults[attr])
        except KeyError:
            if filename not in self._files:
                logger.error(f'File is not included in this state: {file}')
            logger.error(f'Unsupported attribute: {attr}')

    def __setitem__(self, target, value) -> None:
        ''' Set an attribute for a specified file '''
        file, attr = target
        if attr not in self._file_attr_defaults:
            logger.error(f'Unsupported attribute: {attr}')
        file = Path(file)
        filename = str(file)
        try:
            self._files[filename][attr] = value
        except KeyError:
            logger.error(f'File is not included in this state: {file}')
        if self.dump_on_modify:
            self.dump()

    def __delitem__(self, target) -> None:
        ''' Delete (aka reset to default) an attribute for a specified file or
        delete a file from the state '''
        target = (target, ) if not isinstance(target, tuple) else target
        file, attr = target + (None, None)[len(target):]
        try:
            if len(target) == 1:
                self.del_file(file)
            else:
                del(self._files[str(Path(file))][attr])
            if self.dump_on_modify:
                self.dump()
        except KeyError:
            pass

    def __repr__(self) -> str:
        return f'MultiliftState({self.genomes})'

    def __str__(self) -> str:
        ''' Return a YAML representation of the state '''
        rep = {}
        if self.genomes:
            rep['genomes'] = self.genomes
        if self.reference:
            rep['reference'] = self.reference
        if (alignments := self.files(application='alignment')):
            rep['alignments'] = {}
            for file in alignments:
                rep['alignments'][file] = self[file]
        for genome in self.genomes:
            if (files := self.files(genome=genome)):
                rep[genome] = {}
                for file in files:
                    rep[genome][file] = \
                        {k: v for k, v in self[file].items() if k != 'genome'}
        return yaml.dump(
            rep,
            default_flow_style=False, sort_keys=False, explicit_start=True)

    def load(self) -> None:
        ''' Load the YAML contents of self.file, overwriting current state '''
        with open_helper(self.file) as F:
            yaml_data = {} if not (yaml_data := yaml.safe_load(F)) else yaml_data
        self.genomes = yaml_data.get('genomes', [])
        self.reference = yaml_data.get('reference', None)
        self._files = {}
        try:
            for filename, metadata in yaml_data['alignments'].items():
                self.add_file(filename, **metadata)
        except KeyError:
            pass
        for genome in self.genomes:
            try:
                for filename, metadata in yaml_data[genome].items():
                    self.add_file(filename, genome=genome, **metadata)
            except KeyError:
                pass

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
