from datetime import date
from io import BytesIO
import tarfile
import zipfile

import streamlit as st


def v_space(n: int=1, hline: bool=False) -> None:
    ''' '''
    for i in range(n):
        st.write('')
    if hline:
        st.write('***')


def message(message: str='', level: int=0) -> None:
    ''' Display a message in the `container_message` area on script re-run.
    0: success, 1: info, 2: warning, 3: error '''
    st.session_state.message = message
    st.session_state.message_level = level


def add_to_archive(
    arc: BytesIO, data: BytesIO, name: str, arc_type: str='.tar.gz') -> None:
    ''' '''
    if arc_type == '.tgz':
        info = tarfile.TarInfo(name)
        info.size = len(data.getbuffer())
        arc.addfile(info, data)
    elif arc_type == '.zip':
        info = zipfile.ZipInfo(name, date.timetuple(date.today())[:6])
        info.file_size = len(data.getbuffer())
        arc.writestr(info, data.getbuffer())
