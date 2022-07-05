# *multilift*: liftover from multiple sequence alignments

`multilift` facilitates the liftover to 'alignment-space' of data recorded for
mutliple related, variant, sequences by calculating or parsing a multiple
sequence alignment of their genomes.

`multilift` runs through a [`streamlit`](https://streamlit.io/) web app GUI and
is designed for rapid deployment on [`streamlit cloud`](https://streamlit.io/cloud),
dedicated servers, or as a lightweight instance running locally.

---

Requirements:

`multilift` makes use of various new Python language features and requires
> Python >= 3.10

One or more of the multiple sequence aligners:
> mafft (recommended)
>
> kalign
>
> clustalo
>
> muscle

See [requirements.txt](requirements.txt) and [packages.txt](packages.txt) for
more information

---

`multilift` is under active development.
