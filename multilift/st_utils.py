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
