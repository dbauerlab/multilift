import streamlit as st

from multilift import __prog_string__


# Globals #####################################################################


multilift_app_path = __file__

__all__ = ['multilift_app_path']


# Functions ###################################################################


def multilift_app():
    st.title(f'This is {__prog_string__}')


###############################################################################


if __name__ == '__main__':
    multilift_app()


###############################################################################
