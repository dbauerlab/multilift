#!/usr/bin/env python3

import logging

from multilift import __prog_string__, parse_args


if __name__ == "__main__":

    args = parse_args()
    logging.info(__prog_string__)
