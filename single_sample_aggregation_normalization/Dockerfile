FROM vanallenlab/miniconda:2.7

WORKDIR /

COPY requirements.txt /
RUN conda install --yes -c conda-forge --file requirements.txt

COPY combineJunctionAcrossSamples.py /
COPY spliceJunctionReformat.py /
COPY SpliceJunctionNormalization.py /